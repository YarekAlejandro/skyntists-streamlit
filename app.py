"""Interactive field guide to the Skyntists AstroHack spectral analysis."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

from line_statistics import assess_molecular_evidence, catalog_associations, detect_features

ROOT = Path(__file__).resolve().parent
SPECTRA_DIR = ROOT / "espectros" if (ROOT / "espectros").exists() else ROOT
EVIDENCE_PATH = ROOT / "Entregables" / "Campos_de_evidencia.csv"
if not EVIDENCE_PATH.exists():
    EVIDENCE_PATH = ROOT / "Campos_de_evidencia.csv"
CATALOG_PATH = ROOT / "catalogos" / "procedencia_catalogos.json"
if not CATALOG_PATH.exists():
    CATALOG_PATH = ROOT / "procedencia_catalogos.json"
LINE_CATALOG_PATH = ROOT / "catalogos" / "catalogo_local_prebioticas.csv"
if not LINE_CATALOG_PATH.exists():
    LINE_CATALOG_PATH = ROOT / "catalogo_local_prebioticas.csv"
C_KMS = 299_792.458
V_SYSTEM_KMS = 2.7
POSITIONS = {
    "B_offset": ("B_offset_Jy_beam", "B · posición desplazada"),
    "B_pico": ("B_pico_Jy_beam", "B · pico del continuo"),
    "A": ("A_Jy_beam", "A · extracción alternativa"),
}

st.set_page_config(page_title="Skyntists | AstroHack", page_icon="🔭", layout="wide")
st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=DM+Sans:wght@400;500;600;700&family=Space+Grotesk:wght@500;600;700&display=swap');
.stApp {background: radial-gradient(ellipse at 12% 0%, rgba(20,104,130,.15), transparent 32%), #071522;}
.block-container {max-width: 1500px; padding-top: 1.6rem; padding-bottom: 4rem;}
html, body, [class*="css"] {font-family:'DM Sans',sans-serif;}
h1,h2,h3 {font-family:'Space Grotesk',sans-serif; letter-spacing:-.035em;}
.hero {padding:2rem 2.2rem; border:1px solid rgba(133,222,224,.23); border-radius:22px; background:linear-gradient(118deg,rgba(10,39,61,.98),rgba(13,105,119,.88) 66%,rgba(45,99,117,.78)); box-shadow:0 20px 50px rgba(0,0,0,.22); margin:0 0 1.4rem 0;}
.hero-kicker {font-size:.78rem; font-weight:700; text-transform:uppercase; letter-spacing:.16em; color:#9de5df;}
.hero h1 {color:#fff; font-size:clamp(2.2rem,5vw,4rem); margin:.35rem 0 .55rem;}
.hero p {color:#d3e8ed; font-size:1.12rem; max-width:900px; margin:0;}
.callout {padding:1rem 1.2rem; border:1px solid rgba(116,206,203,.24); border-left:4px solid #55d6c2; background:rgba(18,49,66,.75); border-radius:12px; margin:.8rem 0 1.15rem;}
.small-label {color:#9bb4c0; font-size:.82rem; text-transform:uppercase; letter-spacing:.1em;}
.card {border:1px solid rgba(137,180,195,.18); background:linear-gradient(145deg,rgba(18,41,57,.94),rgba(10,27,41,.9)); border-radius:16px; padding:1.05rem 1.2rem; min-height:128px;}
.card b {color:#f0f7fa; font-size:1.07rem;}
.card p {color:#b6c8d2; margin:.42rem 0 0; font-size:.94rem;}
.muted {color:#b3c6d1;}
div[data-testid="stMetric"] {background:rgba(17,43,59,.75);border:1px solid rgba(130,181,193,.18);padding:1rem 1.1rem;border-radius:15px;}
div[data-testid="stMetricLabel"] {color:#b6c8d2;}
div[data-testid="stMetricValue"] {font-family:'Space Grotesk',sans-serif;}
section[data-testid="stSidebar"] {background:#091a28;}
button[data-baseweb="tab"] {font-weight:600;}
</style>
""", unsafe_allow_html=True)


@st.cache_data(show_spinner=False)
def load_spectrum(window: int) -> pd.DataFrame:
    path = SPECTRA_DIR / f"W{window}_espectros.csv"
    frame = pd.read_csv(path)
    required = {"frecuencia_LSRK_GHz", *(value[0] for value in POSITIONS.values())}
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"Faltan columnas en {path.name}: {', '.join(sorted(missing))}")
    return frame


@st.cache_data(show_spinner=False)
def load_evidence() -> pd.DataFrame:
    if not EVIDENCE_PATH.exists():
        return pd.DataFrame()
    frame = pd.read_csv(EVIDENCE_PATH)
    for col in ("frecuencia_reposo_MHz", "frecuencia_pico_LSRK_MHz", "v_radio_kms", "FWHM_kms", "SNR", "area_Jy_beam_kms", "sigma_Jy_beam", "baseline_Jy_beam"):
        if col in frame:
            frame[col] = pd.to_numeric(frame[col], errors="coerce")
    if {"frecuencia_reposo_MHz", "frecuencia_pico_LSRK_MHz"}.issubset(frame.columns):
        # Positive residual means the fitted center is at higher frequency than
        # the catalog transition shifted by the adopted systemic velocity.
        expected = frame["frecuencia_reposo_MHz"] * (1 - V_SYSTEM_KMS / C_KMS)
        frame["residuo_MHz_vs_vsys"] = frame["frecuencia_pico_LSRK_MHz"] - expected
        frame["v_frecuencia_kms"] = C_KMS * (frame["frecuencia_reposo_MHz"] - frame["frecuencia_pico_LSRK_MHz"]) / frame["frecuencia_reposo_MHz"]
    return frame


@st.cache_data(show_spinner=False)
def load_line_catalog() -> pd.DataFrame:
    if not LINE_CATALOG_PATH.exists():
        return pd.DataFrame()
    frame = pd.read_csv(LINE_CATALOG_PATH, dtype={"entrada": str})
    frame["frecuencia_reposo_MHz"] = pd.to_numeric(frame["frecuencia_reposo_MHz"], errors="coerce")
    frame["incertidumbre_MHz"] = pd.to_numeric(frame["incertidumbre_MHz"], errors="coerce")
    return frame.dropna(subset=["especie", "frecuencia_reposo_MHz"])


@st.cache_data(show_spinner="Comparando el catálogo con desplazamientos nulos…")
def run_interactive_analysis(
    window: int,
    position: str,
    threshold_sigma: float,
    velocity_kms: float,
    velocity_margin_kms: float,
    min_distance: int,
) -> tuple[pd.DataFrame, float, float, pd.DataFrame, pd.DataFrame, tuple[float, float, float]]:
    frame = load_spectrum(window)
    frequency = frame["frecuencia_LSRK_GHz"].to_numpy(dtype=float) * 1000
    intensity = frame[POSITIONS[position][0]].to_numpy(dtype=float)
    peaks, sigma, step = detect_features(
        frequency, intensity, threshold_sigma=threshold_sigma,
        min_distance_channels=min_distance,
    )
    band = {window: (float(frequency.min()), float(frequency.max()), step)}
    peak_map = {window: peaks}
    catalog = load_line_catalog()
    if catalog.empty:
        empty = pd.DataFrame()
        return peaks, sigma, step, empty, empty, band[window]
    summary, _ = assess_molecular_evidence(
        catalog, band, peak_map, velocity_kms=velocity_kms,
        velocity_margin_kms=velocity_margin_kms, min_snr=threshold_sigma,
        simulations=1000, random_state=20261011,
    )
    associations = catalog_associations(
        catalog, band, peak_map, velocity_kms=velocity_kms,
        velocity_margin_kms=velocity_margin_kms, min_snr=threshold_sigma,
    )
    return peaks, sigma, step, summary, associations, band[window]


def plot_layout(fig: go.Figure, height: int = 470) -> go.Figure:
    fig.update_layout(
        height=height, template="plotly_dark", paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(7,21,34,.72)", font=dict(family="DM Sans", color="#d7e4ea"),
        margin=dict(l=22, r=24, t=62, b=28), legend=dict(orientation="h", yanchor="bottom", y=1.02, x=0),
        hoverlabel=dict(bgcolor="#102a3b", bordercolor="#5bd4c4", font_color="#f3fbff"),
    )
    fig.update_xaxes(gridcolor="rgba(155,190,202,.12)", zerolinecolor="rgba(155,190,202,.24)")
    fig.update_yaxes(gridcolor="rgba(155,190,202,.12)", zerolinecolor="rgba(155,190,202,.24)")
    return fig


try:
    evidence = load_evidence()
except Exception as exc:
    st.error(f"No se pudo leer la evidencia: {exc}")
    st.stop()

if evidence.empty:
    st.error("No está disponible Campos_de_evidencia.csv; revisa los archivos del despliegue.")
    st.stop()

confirmed = int(evidence["decision"].astype(str).str.contains("Tentativa", case=False).sum())
species_n = evidence["especie"].nunique() if "especie" in evidence else 0
windows_n = evidence["ventana"].nunique() if "ventana" in evidence else 0
with st.sidebar:
    st.markdown("## Skyntists")
    st.caption("Cuaderno de campo digital · IRAS 16293−2422")
    page = st.radio("Explorar", ["Inicio", "Espectros", "Evidencia molecular", "Demostración Doppler", "Método paso a paso", "Glosario y fuentes"], label_visibility="collapsed", key="wiki_page")
    st.divider()
    st.markdown("**Pregunta guía**")
    st.caption("¿Qué especies moleculares son compatibles con los rasgos espectrales observados?")
    st.caption(f"Datos cargados: {len(evidence)} filas de evidencia · {windows_n} ventanas")

if page == "Inicio":
    st.markdown("""
    <div class="hero">
     <div class="hero-kicker">Guía interactiva · AstroHack SATEM 2026</div>
     <h1>🔭 Skyntists</h1>
     <p>Un recorrido reproducible por los espectros de IRAS 16293−2422: del dato calibrado a una hipótesis molecular, haciendo visibles la evidencia, las decisiones y sus límites.</p>
    </div>
    """, unsafe_allow_html=True)
    m1, m2, m3, m4 = st.columns(4)
    m1.metric("Ventanas ALMA", f"{windows_n}", help="W1 y W2; sus haces y escalas de intensidad son distintos.")
    m2.metric("Canales", "1,918 / ventana")
    m3.metric("Perfiles evaluados", f"{len(evidence)}", help="Incluye candidatos tentativos y el caso de glicina-I insuficiente.")
    m4.metric("Compatibilidad tentativa", f"{confirmed} perfiles", help="No son identificaciones confirmadas ni una probabilidad molecular.")
    st.markdown('<div class="callout"><b>Qué podemos afirmar:</b> seis perfiles muestran compatibilidad espectral tentativa con CH₃CN o CH₃CHO bajo el análisis actual. La coincidencia aislada de Glycine-I es insuficiente. Las líneas pueden mezclarse; el catálogo de apoyo es parcial y no se realizó un modelado LTE conjunto.</div>', unsafe_allow_html=True)

if page == "Inicio":
    st.subheader("Del espectro a la hipótesis")
    st.markdown("La app funciona como wiki del proyecto: combina observación, comparación con catálogos, decisiones de análisis y una explicación explícita de qué no se puede concluir.")
    cols = st.columns(3)
    with cols[0]:
        st.markdown('<div class="card"><span class="small-label">01 · Observar</span><br><b>Espectros calibrados</b><p>Explora W1/W2, elige una extracción y acerca el eje de frecuencia a una transición para inspeccionar el perfil sin suavizado.</p></div>', unsafe_allow_html=True)
    with cols[1]:
        st.markdown('<div class="card"><span class="small-label">02 · Contrastar</span><br><b>Catálogo y ajuste</b><p>Compara el centro ajustado con la frecuencia de catálogo corregida por una velocidad sistémica declarada de 2.7 km/s.</p></div>', unsafe_allow_html=True)
    with cols[2]:
        st.markdown('<div class="card"><span class="small-label">03 · Calificar</span><br><b>Evidencia, no veredicto</b><p>Revisa SNR, ancho, mezcla y alternativas. Una coincidencia de frecuencia por sí sola no identifica una molécula.</p></div>', unsafe_allow_html=True)
    st.markdown("### Lectura rápida del resultado")
    species_summary = evidence.groupby("especie", dropna=False).agg(perfiles=("especie", "size"), decision=("decision", "first"), mediana_SNR=("SNR", "median")).reset_index()
    species_summary["mediana_SNR"] = species_summary["mediana_SNR"].round(1)
    a, b = st.columns([1.15, 1])
    with a:
        fig = px.bar(species_summary, x="especie", y="perfiles", color="especie", text="perfiles", title="Perfiles evaluados por especie", color_discrete_sequence=["#55d6c2", "#6da8ff", "#ffc66d"])
        fig.update_traces(textposition="outside", marker_line_width=0)
        fig.update_layout(showlegend=False, xaxis_title=None, yaxis_title="Número de perfiles")
        st.plotly_chart(plot_layout(fig, 365), width="stretch")
    with b:
        st.markdown("#### Cómo interpretar el SNR")
        st.write("El SNR informado es el contraste local aproximado del perfil respecto al ruido estimado. Un SNR alto señala una estructura destacada en ese intervalo; no mide la probabilidad de que sea una especie concreta.")
        st.info("Identificación robusta requiere varias transiciones independientes, consistencia entre líneas, tratamiento de blends y comparación con alternativas. La app conserva esta distinción en cada visualización.")
    st.markdown("### Ruta de lectura")
    st.write("1. Empieza en **Espectros** y selecciona un rasgo; 2. revisa la comparación de frecuencias en **Evidencia molecular**; 3. sigue las ecuaciones y supuestos en **Método paso a paso**; 4. consulta términos y referencias en **Glosario y fuentes**.")


elif page == "Espectros":
    st.subheader("Espectros de las extracciones calibradas")
    st.write("El modo inicial muestra una sola posición y acerca la vista a una línea candidata. Cambia a banda completa para recuperar el contexto. Usa la leyenda para ocultar o mostrar series.")
    left, right = st.columns([1, 2])
    with left:
        window = st.selectbox("Ventana", [1, 2], format_func=lambda value: f"W{value}")
        position = st.selectbox("Extracción", list(POSITIONS), index=0, format_func=lambda value: POSITIONS[value][1])
        compare = st.checkbox("Comparar tres extracciones", value=False)
        view_mode = st.radio("Escala horizontal", ["Región de línea", "Ventana completa"], horizontal=True)
    df = load_spectrum(window)
    series = list(POSITIONS) if compare else [position]
    freq_key = "frecuencia_LSRK_GHz"
    matching = evidence.loc[evidence["ventana"].eq(f"W{window}")].copy()
    candidate_options = matching.index.tolist()
    if candidate_options:
        candidate_idx = st.selectbox("Línea de referencia", candidate_options, format_func=lambda idx: f"{evidence.loc[idx, 'especie']} · {evidence.loc[idx, 'grupo_rasgo']} · {evidence.loc[idx, 'frecuencia_pico_LSRK_MHz']:.3f} MHz")
        candidate = evidence.loc[candidate_idx]
    else:
        candidate = None
    fig = go.Figure()
    palette = {"B_offset": "#60d7c5", "B_pico": "#76a9fa", "A": "#ffbd69"}
    for key in series:
        value_col, label = POSITIONS[key]
        fig.add_trace(go.Scatter(x=df[freq_key], y=df[value_col], mode="lines", name=label, line=dict(color=palette[key], width=1.35 if len(series)>1 else 1.8), hovertemplate="ν = %{x:.6f} GHz<br>I = %{y:.4f} Jy/beam<extra>%{fullData.name}</extra>"))
    if candidate is not None:
        fit_ghz = float(candidate["frecuencia_pico_LSRK_MHz"])/1000
        rest_shifted_ghz = float(candidate["frecuencia_reposo_MHz"])*(1 - V_SYSTEM_KMS/C_KMS)/1000
        intensity_col = POSITIONS[position][0]
        nearest_i = (df[freq_key] - fit_ghz).abs().idxmin()
        fit_intensity = float(df.loc[nearest_i, intensity_col])
        fig.add_trace(go.Scatter(x=[fit_ghz], y=[fit_intensity], mode="markers", name="Centro ajustado", marker=dict(size=12, color="#fff1a8", line=dict(width=2, color="#102332")), hovertemplate="Centro ajustado<br>ν = %{x:.6f} GHz<extra></extra>"))
        fig.add_vline(x=rest_shifted_ghz, line_dash="dot", line_color="#f2a7ff", line_width=1.7, annotation_text="Catálogo @ 2.7 km/s", annotation_position="top left", annotation_font_color="#f2a7ff")
        if view_mode == "Región de línea":
            fwhm = candidate["FWHM_kms"]
            fwhm = 1.0 if pd.isna(fwhm) else float(fwhm)
            half_span = max(0.006, 5 * fwhm * float(candidate["frecuencia_reposo_MHz"])/(C_KMS*1000))
            fig.update_xaxes(range=[fit_ghz-half_span, fit_ghz+half_span])
        st.caption(f"Selección: {candidate['especie']} · {candidate['grupo_rasgo']}. Línea amarilla: centro del ajuste; línea magenta: frecuencia de catálogo llevada a LSRK con v = {V_SYSTEM_KMS:.1f} km/s. La segunda es una referencia cinemática, no una medida independiente.")
    fig.update_layout(title=f"W{window} · {POSITIONS[position][1] if not compare else 'comparación de posiciones'}", xaxis_title="Frecuencia observada LSRK (GHz)", yaxis_title="Intensidad (Jy/beam)", hovermode="x unified")
    st.plotly_chart(plot_layout(fig, 540), width="stretch")
    st.markdown("#### ¿Qué representa cada extracción?")
    st.write("Las posiciones corresponden a espectros extraídos del producto procesado. B_offset es la posición desplazada recomendada para evaluar emisión menos contaminada por el continuo brillante. B_pico y A sirven como comparación espacial. Los haces difieren entre ventanas; no interpretes diferencias de amplitud entre W1/W2 como cambios de abundancia.")
    st.download_button("Descargar espectro de la ventana", df.to_csv(index=False).encode("utf-8"), f"W{window}_espectros.csv", "text/csv")


elif page == "Evidencia molecular":
    st.subheader("Frecuencia de catálogo y centro del ajuste")
    st.markdown("La comparación directa de frecuencias absolutas es engañosa a esta escala: ambas rondan 349 000 MHz, mientras la separación relevante es de alrededor de 0–2 MHz. Aquí mostramos el **residuo en MHz** después de llevar la frecuencia de reposo a LSRK usando el valor sistémico de 2.7 km/s. Así los puntos quedan en una escala legible y el cero tiene significado físico.")
    species = sorted(evidence["especie"].dropna().unique())
    chosen = st.multiselect("Especies incluidas", species, default=species)
    windows = sorted(evidence["ventana"].dropna().unique())
    chosen_windows = st.multiselect("Ventanas incluidas", windows, default=windows)
    view = evidence[evidence["especie"].isin(chosen) & evidence["ventana"].isin(chosen_windows)].copy()
    if not view.empty:
        fig = px.strip(view, x="especie", y="residuo_MHz_vs_vsys", color="especie", hover_name="grupo_rasgo", hover_data={"frecuencia_reposo_MHz": ":.4f", "frecuencia_pico_LSRK_MHz": ":.4f", "v_frecuencia_kms": ":.3f", "ventana": True, "residuo_MHz_vs_vsys": ":+.4f"}, stripmode="overlay", title="Desplazamiento del centro ajustado frente a la predicción cinemática")
        fig.add_hline(y=0, line_dash="dash", line_color="#f4f7fa", opacity=.85, annotation_text="Coincidencia con v = 2.7 km/s", annotation_position="top left")
        fig.update_traces(jitter=0.22, pointpos=0, marker=dict(size=13, line=dict(width=1, color="#071522")))
        fig.update_layout(xaxis_title="Especie candidata", yaxis_title="ν ajuste − ν catálogo corregida (MHz)", showlegend=False)
        st.plotly_chart(plot_layout(fig, 500), width="stretch")
        st.caption("Residuo positivo: el centro ajustado queda a mayor frecuencia que la predicción. El residuo y v_frecuencia se derivan de estas dos frecuencias; no son mediciones independientes ni incertidumbres. No se dibujan barras de error porque la tabla no contiene incertidumbres del centro del ajuste.")
        chart, explanation = st.columns([1.4, 1])
        with chart:
            plotted = view.dropna(subset=["v_frecuencia_kms", "SNR", "FWHM_kms"])
            if not plotted.empty:
                figv = px.scatter(plotted, x="v_frecuencia_kms", y="SNR", color="especie", symbol="ventana", size="FWHM_kms", hover_name="grupo_rasgo", hover_data={"frecuencia_reposo_MHz": ":.3f", "frecuencia_pico_LSRK_MHz": ":.3f", "FWHM_kms": ":.2f", "mezcla": True}, title="Contraste local y cinemática de los perfiles")
                figv.add_vline(x=V_SYSTEM_KMS, line_dash="dash", line_color="#f5d26c", annotation_text="v sistémica adoptada")
                figv.update_layout(xaxis_title="Velocidad radial derivada (km/s)", yaxis_title="SNR local aproximado")
                st.plotly_chart(plot_layout(figv, 430), width="stretch")
                if len(plotted) < len(view):
                    st.caption("Se excluyen de esta gráfica los perfiles sin SNR o FWHM medidos; esas filas permanecen en el residuo y la tabla.")
        with explanation:
            st.markdown("#### Lectura de la evidencia")
            st.write("En la gráfica de SNR, el tamaño del símbolo representa el FWHM ajustado. Un perfil estrecho y de alto contraste puede ser una buena línea candidata, pero congestión y coincidencias de otras especies reducen la especificidad.")
            st.info("El candidato Glycine-I se conserva para mostrar por qué una coincidencia aislada en una zona congestionada no basta: faltan transiciones independientes que respalden la asignación.")
        table_cols = [c for c in ["especie", "catalogo", "entrada", "frecuencia_reposo_MHz", "frecuencia_pico_LSRK_MHz", "residuo_MHz_vs_vsys", "v_frecuencia_kms", "ventana", "FWHM_kms", "SNR", "mezcla", "alternativas", "decision"] if c in view]
        st.markdown("#### Tabla trazable")
        st.dataframe(view[table_cols], width="stretch", hide_index=True)
        st.download_button("Descargar evidencia filtrada", view.to_csv(index=False).encode("utf-8"), "evidencia_filtrada.csv", "text/csv")
    else:
        st.warning("No hay filas para los filtros seleccionados.")


elif page == "Demostración Doppler":
    st.subheader("Demostración interactiva: ¿se alinean las líneas?")
    st.markdown(
        "Mueve el umbral de detección y la velocidad sistémica. El gráfico vuelve a detectar picos y desplaza las frecuencias de catálogo; "
        "las guías verdes indican transiciones asociadas bajo la tolerancia seleccionada y las ámbar, transiciones sin asociación."
    )
    catalog = load_line_catalog()
    if catalog.empty:
        st.error("No está disponible el catálogo local de frecuencias para esta demostración.")
    else:
        c1, c2, c3 = st.columns(3)
        with c1:
            demo_window = st.selectbox("Ventana espectral", [1, 2], format_func=lambda value: f"W{value}", key="demo_window")
            demo_position = st.selectbox("Extracción", list(POSITIONS), index=0, format_func=lambda value: POSITIONS[value][1], key="demo_position")
        with c2:
            demo_threshold = st.slider("Umbral de prominencia (σ local)", 3.0, 10.0, 5.0, 0.5, key="demo_threshold")
            demo_velocity = st.slider("Velocidad sistémica v (km/s)", 0.0, 6.0, 2.7, 0.05, key="demo_velocity")
        with c3:
            demo_margin = st.slider("Margen cinemático ± (km/s)", 0.1, 2.0, 1.0, 0.1, key="demo_margin")
            demo_distance = st.slider("Separación mínima entre picos (canales)", 1, 12, 3, key="demo_distance")

        species_options = sorted(catalog["especie"].dropna().astype(str).unique())
        demo_species = st.selectbox("Molécula del catálogo que quieres contrastar", species_options, index=species_options.index("CH3CN") if "CH3CN" in species_options else 0)
        peaks, demo_sigma, demo_step, grade_table, all_matches, band = run_interactive_analysis(
            demo_window, demo_position, demo_threshold, demo_velocity, demo_margin, demo_distance
        )
        raw = load_spectrum(demo_window)
        frequency = raw["frecuencia_LSRK_GHz"].to_numpy(dtype=float) * 1000
        intensity = raw[POSITIONS[demo_position][0]].to_numpy(dtype=float)
        low, high, step = band
        molecule_lines = catalog.loc[catalog["especie"].astype(str).eq(demo_species)].copy()
        peak_frequency = peaks["frecuencia_pico_LSRK_MHz"].to_numpy(dtype=float)
        transition_rows = []
        for transition in molecule_lines.itertuples(index=False):
            rest = float(transition.frecuencia_reposo_MHz)
            predicted = rest * (1 - demo_velocity / C_KMS)
            if not low <= predicted <= high:
                continue
            catalog_sigma = float(transition.incertidumbre_MHz) if pd.notna(transition.incertidumbre_MHz) and float(transition.incertidumbre_MHz) > 0 else 0.0
            tolerance = np.sqrt((predicted * demo_margin / C_KMS) ** 2 + catalog_sigma**2 + (step / 2) ** 2)
            tolerance = min(max(float(tolerance), step), 1.5)
            if len(peaks):
                nearest = int(np.argmin(np.abs(peak_frequency - predicted)))
                found_frequency = float(peak_frequency[nearest])
                found_snr = float(peaks.iloc[nearest]["SNR_prominencia_aprox"])
                residual = found_frequency - predicted
                aligned = abs(residual) <= tolerance and found_snr >= demo_threshold
                channel = int(peaks.iloc[nearest]["canal"])
            else:
                found_frequency, found_snr, residual, aligned, channel = np.nan, np.nan, np.nan, False, None
            transition_rows.append({
                "catálogo": str(transition.catalogo),
                "entrada": str(transition.entrada),
                "ν reposo (MHz)": rest,
                "ν predicha LSRK (MHz)": predicted,
                "pico cercano (MHz)": found_frequency,
                "residuo (MHz)": residual,
                "prominencia/σ": found_snr,
                "tolerancia (MHz)": tolerance,
                "grupo de rasgo": f"W{demo_window}-{channel}" if channel is not None else "—",
                "asociación": "alineada" if aligned else "sin asociación",
            })
        transitions = pd.DataFrame(transition_rows)

        summary_row = grade_table.loc[grade_table["especie"].eq(demo_species)] if not grade_table.empty else pd.DataFrame()
        matched_lines = int(transitions["asociación"].eq("alineada").sum()) if not transitions.empty else 0
        distinct_features = int(summary_row["rasgos_distintos_asociados"].iloc[0]) if not summary_row.empty else 0
        q_value = float(summary_row["q_BY_12_especies"].iloc[0]) if not summary_row.empty else np.nan
        grade = str(summary_row["grado_evidencia"].iloc[0]) if not summary_row.empty else "No evaluable"
        metrics = st.columns(4)
        grade_short = {
            "Compatibilidad tentativa": "Tentativa",
            "Sugestiva, ambigua": "Ambigua",
            "Indicio débil; insuficiente para asignar": "Débil",
            "Sin evidencia suficiente en esta prueba": "Insuf.",
        }.get(grade, "No evaluable")
        metrics[0].metric("Ruido σ", f"{demo_sigma:.4f}" if np.isfinite(demo_sigma) else "—", help="Estimación aproximada en Jy/beam a partir de diferencias entre canales.")
        metrics[1].metric("Picos", f"{len(peaks)}", help=f"Máximos con prominencia ≥ {demo_threshold:g}σ local.")
        metrics[2].metric("Líneas", f"{matched_lines}/{len(transitions)}", help="Transiciones del catálogo local asociadas en esta ventana; no es una tasa de identidad molecular.")
        metrics[3].metric("Grado", grade_short, help=f"{grade}. Solo para la ventana elegida; no es probabilidad posterior ni la evaluación conjunta W1+W2.")
        if np.isfinite(q_value):
            st.caption(f"Para {demo_species} en W{demo_window}: {distinct_features} rasgos distintos asociados · q empírico BY = {q_value:.4g} frente a desplazamientos circulares, corregido por las 12 especies del catálogo local.")

        fig = go.Figure()
        fig.add_trace(go.Scatter(
            x=frequency, y=intensity, mode="lines", name="Espectro observado",
            line=dict(color="#60d7c5", width=1.4),
            hovertemplate="ν=%{x:.5f} MHz<br>I=%{y:.4f} Jy/beam<extra></extra>",
        ))
        if len(peaks):
            fig.add_trace(go.Scatter(
                x=peak_frequency, y=peaks["pico_Jy_beam"], mode="markers", name="Pico sobre umbral",
                marker=dict(color="#ffc66d", size=7, line=dict(color="#071522", width=1)),
                customdata=np.column_stack([peaks["SNR_prominencia_aprox"], peaks["canal"]]),
                hovertemplate="Pico detectado · canal %{customdata[1]:.0f}<br>ν=%{x:.5f} MHz<br>prominencia/σ=%{customdata[0]:.1f}<extra></extra>",
            ))
        if not transitions.empty:
            y_span = float(np.nanmax(intensity) - np.nanmin(intensity)) or 1.0
            rug_y = float(np.nanmax(intensity)) + y_span * 0.08
            marker_y = np.linspace(rug_y, rug_y + y_span * 0.10, len(transitions))
            colors = ["#55d6c2" if value == "alineada" else "#ffbd69" for value in transitions["asociación"]]
            marker_hover = [[
                str(row["entrada"]),
                float(row["residuo (MHz)"]) if pd.notna(row["residuo (MHz)"]) else np.nan,
                float(row["prominencia/σ"]) if pd.notna(row["prominencia/σ"]) else np.nan,
                str(row["asociación"]),
            ] for _, row in transitions.iterrows()]
            fig.add_trace(go.Scatter(
                x=transitions["ν predicha LSRK (MHz)"], y=marker_y, mode="markers", name="Transición de catálogo desplazada",
                marker=dict(color=colors, size=9, symbol="line-ns-open", line=dict(width=2)),
                customdata=marker_hover,
                hovertemplate="Entrada %{customdata[0]}<br>ν predicha=%{x:.5f} MHz<br>residuo=%{customdata[1]} MHz<br>prominencia/σ=%{customdata[2]}<br>%{customdata[3]}<extra></extra>",
            ))
            for _, line in transitions.iterrows():
                color = "#55d6c2" if line["asociación"] == "alineada" else "#ffbd69"
                fig.add_vline(x=float(line["ν predicha LSRK (MHz)"]), line_color=color, line_width=1, line_dash="dot", opacity=0.42)
            if len(transitions) > 1:
                option_labels = [f"{row['entrada']} · {row['ν reposo (MHz)']:.4f} MHz reposo" for _, row in transitions.iterrows()]
                selected_line = st.selectbox("Centrar la vista en una transición", range(len(transitions)), format_func=lambda idx: option_labels[idx], key="demo_focus_line")
                center = float(transitions.iloc[selected_line]["ν predicha LSRK (MHz)"])
                focus_width = st.slider("Semiancho de la vista local (MHz)", 0.5, 15.0, 3.0, 0.5, key="demo_focus_width")
                fig.update_xaxes(range=[center - focus_width, center + focus_width])
        fig.update_layout(title=f"W{demo_window} · {demo_species} · v = {demo_velocity:.2f} km/s", xaxis_title="Frecuencia LSRK (MHz)", yaxis_title="Intensidad (Jy/beam)", hovermode="x unified")
        st.plotly_chart(plot_layout(fig, 565), width="stretch")
        if not transitions.empty:
            st.markdown("#### Transiciones cubiertas del catálogo")
            st.dataframe(transitions, width="stretch", hide_index=True)
        st.info(
            "El p/q empírico compara el número de rasgos asociados con patrones de catálogo desplazados circularmente dentro de la misma ventana y corrige por las 12 especies evaluadas. "
            "Es condicional a los deslizadores seleccionados y no corrige por explorar muchas velocidades; úsalo para sensibilidad visual, no para optimizar v y reportar significancia. **No es el porcentaje de probabilidad de que la molécula esté presente**. "
            "Los canales pueden estar correlacionados, el catálogo es parcial, las intensidades no se modelan y los blends requieren inspección."
        )


elif page == "Método paso a paso":
    st.subheader("De los datos calibrados al estado de evidencia")
    st.write("Esta sección documenta qué transforma el análisis y qué queda pendiente. Abre los pasos para ver las ecuaciones y los supuestos.")
    with st.expander("1 · Datos de entrada y marco de referencia", expanded=True):
        st.write("Se exploran dos ventanas espectrales ALMA del proyecto ADS/JAO.ALMA#2013.1.00278.S, productos procesados ARI-L. Las frecuencias están en LSRK y las intensidades en Jy/beam. Los CSV de la app son extracciones 1D; no contienen por sí solos la información espacial completa de los cubos FITS.")
    with st.expander("2 · Frecuencia de reposo a frecuencia observada"):
        st.latex(r"\nu_{\rm pred,LSRK}=\nu_{\rm rest}\left(1-\frac{v_{\rm sys}}{c}\right),\qquad v_{\rm sys}=2.7\;\mathrm{km\,s^{-1}}")
        st.write("Se usa la aproximación no relativista para velocidades pequeñas frente a c. La corrección se aplica una sola vez porque el espectro ya está etiquetado en LSRK. El valor sistémico sirve como referencia común; la velocidad inferida de cada rasgo puede diferir.")
    with st.expander("3 · Residuo que se representa"):
        st.latex(r"\Delta\nu=\nu_{\rm fit,LSRK}-\nu_{\rm pred,LSRK}")
        st.write("El gráfico principal centra cada candidata en una escala local en MHz. El residuo indica concordancia con la velocidad sistémica adoptada, pero no es independiente de los datos de ajuste ni una incertidumbre estadística.")
    with st.expander("4 · Ruido robusto, ajuste y SNR"):
        st.write("La sigma aproximada usa la dispersión robusta (MAD) de primeras diferencias dividida entre √2. El perfil se ajusta localmente con componente gaussiana y baseline. El SNR es contraste local aproximado. Correlación entre canales, baseline imperfecto, blends y selección de ventanas pueden alterar su valor; no debe interpretarse como significancia global.")
    with st.expander("5 · Área integrada y agrupación"):
        st.write("El área se aproxima integrando el perfil sobre el ancho seleccionado; su error aproximado no incluye correlación, calibración ni incertidumbre de modelo. Varias transiciones pueden asociarse al mismo grupo de rasgo: 107 asociaciones del catálogo se agruparon en 64 rasgos, por lo que asociación no significa detección independiente.")
    with st.expander("6 · Qué exige una identificación más fuerte"):
        st.write("Buscar múltiples transiciones no confundidas de la misma especie, modelar blends y alternativas, reportar incertidumbres de frecuencia y parámetros de ajuste, probar robustez a selección de baseline/ventana y comparar intensidades con un modelo de excitación/transferencia radiativa. El inventario de contaminantes también debe ser más completo.")
    with st.expander("7 · Grado de compatibilidad de la demostración"):
        st.latex(r"p_{\rm circ}=\frac{1+\#\{T_b\geq T_{\rm obs}\}}{1+B},\qquad q_{\rm BY}=\mathrm{ajuste\ por\ especies}")
        st.write("El estadístico T cuenta rasgos espectrales distintos asociados a una especie. El nulo desplaza circularmente el patrón de frecuencias dentro de cada ventana y vuelve a contar coincidencias; así aproxima cuántos alineamientos pueden surgir por densidad de picos/catálogo. Los p empíricos se ajustan con Benjamini–Yekutieli entre las especies del catálogo para controlar la tasa de falsos descubrimientos bajo dependencia general.")
        st.write("La regla operativa marca compatibilidad tentativa si hay al menos tres rasgos distintos, q≤0.05, mediana de velocidad a no más del margen elegido, dispersión robusta de velocidad ≤0.5 km/s y menos de la mitad de los rasgos compartidos con otras especies. Es una regla explícita de este proyecto, no un umbral universal; la etiqueta se mantiene tentativa porque aquí faltan intensidades modeladas, catálogo completo y separación de blends.")
        st.warning("El valor q no es la probabilidad de que una especie sea falsa ni P(molécula | datos). En la demostración interactiva el q es condicional a los parámetros seleccionados y no corrige por recorrer el deslizador de velocidad; usa el cuaderno, con la velocidad de referencia fijada, para la tabla de grado. El nulo tampoco simula intensidades, líneas ausentes esperadas, catálogo incompleto, covarianza espectral ni transferencia radiativa.")
    st.markdown("### ¿Por qué importan los coeficientes de Einstein?")
    st.write("El coeficiente Aᵤₗ (s⁻¹) describe la probabilidad de emisión espontánea por unidad de tiempo para una transición entre niveles u y l. Junto a la degeneración estadística, energía del nivel superior, temperatura y función de partición, ayuda a predecir intensidades relativas en equilibrio local (LTE) y a modelar transferencia radiativa fuera de LTE. Por sí solo no identifica una línea: hacen falta parámetros consistentes y un modelo de la emisión.")
    st.warning("El catálogo local no entrega Aᵤₗ para estas filas. La columna log10_I_300K_nm2_MHz representa intensidad de catálogo bajo condiciones definidas; no es log10(Aᵤₗ). No rellenamos Aᵤₗ con valores inventados.")


else:
    st.subheader("Glosario, datos y procedencia")
    terms = {
        "LSRK": "Marco cinemático local estándar de reposo. La frecuencia observada ya fue corregida al marco de referencia indicado por los productos.",
        "Frecuencia de reposo": "Frecuencia de una transición molecular tabulada en el catálogo, antes de aplicar el desplazamiento Doppler del gas.",
        "Centro del ajuste": "Frecuencia central del perfil usado para describir una estructura observada; no equivale automáticamente a una identificación molecular.",
        "FWHM": "Ancho a mitad de altura de un perfil gaussiano ajustado, aquí expresado en km/s.",
        "SNR": "Razón señal/ruido local aproximada para el rasgo. No representa una probabilidad de identidad o significancia global.",
        "q empírico BY": "Valor p empírico de alineamiento contra desplazamientos circulares, ajustado por pruebas múltiples con Benjamini–Yekutieli. No es una probabilidad posterior de identidad molecular.",
        "Mezcla (blend)": "Solapamiento espectral de varias transiciones dentro del ancho instrumental o del perfil físico.",
        "Aᵤₗ (Einstein)": "Probabilidad de emisión espontánea por unidad de tiempo de la transición u→l, en s⁻¹.",
    }
    for term, definition in terms.items():
        with st.expander(term):
            st.write(definition)
    st.markdown("### Procedencia y referencias")
    st.markdown("""
- **Observaciones:** ALMA, proyecto `ADS/JAO.ALMA#2013.1.00278.S`, productos procesados ARI-L; espectros LSRK en Jy/beam.
- **Contexto científico:** Jørgensen et al. (2016), PILS, *Astronomy & Astrophysics*, 595, A117 ([arXiv:1607.08733](https://arxiv.org/abs/1607.08733)).
- **Catálogos moleculares:** [CDMS](https://cdms.astro.uni-koeln.de/classic/entries/) y [JPL](https://spec.jpl.nasa.gov/).
- **Control de falsos descubrimientos:** [Miller et al. (2001), aplicación en análisis astrofísico](https://doi.org/10.1086/324109); [Benjamini & Yekutieli (2001), dependencia entre pruebas](https://doi.org/10.1214/aos/1013699998).
- **Contexto de identificación molecular:** [Jørgensen et al. (2016), PILS](https://arxiv.org/abs/1607.08733).
- **Archivos:** `W1_espectros.csv`, `W2_espectros.csv`, `Campos_de_evidencia.csv`, `procedencia_catalogos.json`.
- **Reproducibilidad:** código de análisis y ajustes detallados en el cuaderno `herramientas/Inicio_AstroHack.ipynb` del repositorio de proyecto.
""")
    st.warning("La lista de catálogo de apoyo es parcial y no constituye un inventario exhaustivo de contaminantes. Los datos públicos de esta app son espectros extraídos y una tabla de evidencia: no incluyen cubos FITS ni la presentación PowerPoint.")

st.divider()
st.caption("Skyntists · AstroHack SATEM 2026 · Análisis exploratorio. La compatibilidad espectral descrita aquí es tentativa y depende de los datos, catálogos y supuestos documentados.")

