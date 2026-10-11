"""Interactive field guide to the Skyntists AstroHack spectral analysis."""
from __future__ import annotations

from pathlib import Path

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

ROOT = Path(__file__).resolve().parent
SPECTRA_DIR = ROOT / "espectros" if (ROOT / "espectros").exists() else ROOT
EVIDENCE_PATH = ROOT / "Entregables" / "Campos_de_evidencia.csv"
if not EVIDENCE_PATH.exists():
    EVIDENCE_PATH = ROOT / "Campos_de_evidencia.csv"
CATALOG_PATH = ROOT / "catalogos" / "procedencia_catalogos.json"
if not CATALOG_PATH.exists():
    CATALOG_PATH = ROOT / "procedencia_catalogos.json"
C_KMS = 299_792.458
V_SYSTEM_KMS = 2.7
POSITIONS = {
    "B_offset": ("B_offset_Jy_beam", "B · posición desplazada"),
    "B_pico": ("B_pico_Jy_beam", "B · pico del continuo"),
    "A": ("A_Jy_beam", "A · extracción alternativa"),
}

st.set_page_config(page_title="Skyntists | AstroHack", layout="wide")
st.markdown("""
<style>
.stApp {background:#f3f5f7; color:#1c2c3a;}
.block-container {max-width:1440px; padding-top:1.4rem; padding-bottom:3.5rem;}
html, body, [class*="css"] {font-family:Inter, "Segoe UI", Arial, sans-serif; color:#1c2c3a;}
h1,h2,h3 {font-family:Inter, "Segoe UI", Arial, sans-serif; color:#17324d; letter-spacing:-.025em;}
h2, h3 {font-weight:650;}
.hero {padding:1.65rem 1.9rem; border:1px solid #d9e1e8; border-left:5px solid #176b87; border-radius:8px; background:#fff; box-shadow:0 3px 14px rgba(20,43,62,.045); margin:0 0 1.15rem;}
.hero-kicker {font-size:.75rem; font-weight:700; text-transform:uppercase; letter-spacing:.13em; color:#176b87;}
.hero h1 {color:#17324d; font-size:clamp(2rem,4vw,3.35rem); margin:.3rem 0 .45rem;}
.hero p {color:#526474; font-size:1.04rem; max-width:940px; margin:0; line-height:1.65;}
.callout {padding:.9rem 1.1rem; border:1px solid #d7e3e9; border-left:4px solid #176b87; background:#fff; border-radius:5px; margin:.75rem 0 1.1rem; color:#33495b; line-height:1.6;}
.small-label {color:#607487; font-size:.75rem; text-transform:uppercase; letter-spacing:.09em; font-weight:700;}
.card {border:1px solid #dce3e9; background:#fff; border-radius:7px; padding:1rem 1.1rem; min-height:128px; box-shadow:0 2px 8px rgba(20,43,62,.035);}
.card b {color:#17324d; font-size:1.04rem;}
.card p {color:#526474; margin:.42rem 0 0; font-size:.92rem; line-height:1.55;}
.muted {color:#607487;}
div[data-testid="stMetric"] {background:#fff; border:1px solid #dce3e9; padding:.9rem 1rem; border-radius:7px; box-shadow:0 2px 8px rgba(20,43,62,.035);}
div[data-testid="stMetricLabel"] {color:#607487;}
div[data-testid="stMetricValue"] {font-family:Inter, "Segoe UI", Arial, sans-serif; color:#17324d; font-weight:650;}
section[data-testid="stSidebar"] {background:#e9eef2; border-right:1px solid #d8e0e6;}
section[data-testid="stSidebar"] [data-testid="stMarkdownContainer"] p {color:#3e5263;}
button[data-baseweb="tab"] {font-weight:600;}
div[data-testid="stDataFrame"] {border:1px solid #dce3e9; border-radius:6px; overflow:hidden;}
.stDownloadButton button {border-radius:5px; font-weight:600;}
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


def plot_layout(fig: go.Figure, height: int = 470) -> go.Figure:
    fig.update_layout(
        height=height, template="plotly_white", paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="#ffffff", font=dict(family="Inter, Segoe UI, Arial, sans-serif", color="#263b4d", size=12),
        margin=dict(l=28, r=26, t=64, b=34), legend=dict(orientation="h", yanchor="bottom", y=1.02, x=0),
        hoverlabel=dict(bgcolor="#ffffff", bordercolor="#176b87", font_color="#1c2c3a"),
    )
    fig.update_xaxes(gridcolor="#e7edf1", zerolinecolor="#aab8c3", linecolor="#b8c4cd", ticks="outside")
    fig.update_yaxes(gridcolor="#e7edf1", zerolinecolor="#aab8c3", linecolor="#b8c4cd", ticks="outside")
    return fig


try:
    evidence = load_evidence()
except Exception as exc:
    st.error(f"No se pudo leer la evidencia: {exc}")
    st.stop()

st.markdown("""
<div class="hero">
 <div class="hero-kicker">Guía interactiva · AstroHack SATEM 2026</div>
 <h1>Skyntists</h1>
 <p>Un recorrido reproducible por los espectros de IRAS 16293−2422: del dato calibrado a una hipótesis molecular, haciendo visibles la evidencia, las decisiones y sus límites.</p>
</div>
""", unsafe_allow_html=True)

if evidence.empty:
    st.error("No está disponible Campos_de_evidencia.csv; revisa los archivos del despliegue.")
    st.stop()

confirmed = int(evidence["decision"].astype(str).str.contains("Tentativa", case=False).sum())
species_n = evidence["especie"].nunique() if "especie" in evidence else 0
windows_n = evidence["ventana"].nunique() if "ventana" in evidence else 0
m1, m2, m3, m4 = st.columns(4)
m1.metric("Ventanas ALMA", f"{windows_n}", help="W1 y W2; sus haces y escalas de intensidad son distintos.")
m2.metric("Canales", "1,918 / ventana")
m3.metric("Perfiles evaluados", f"{len(evidence)}", help="Incluye candidatos tentativos y el caso de glicina-I insuficiente.")
m4.metric("Compatibilidad tentativa", f"{confirmed} perfiles", help="No son identificaciones confirmadas ni una probabilidad molecular.")

st.markdown('<div class="callout"><b>Qué podemos afirmar:</b> seis perfiles muestran compatibilidad espectral tentativa con CH₃CN o CH₃CHO bajo el análisis actual. La coincidencia aislada de Glycine-I es insuficiente. Las líneas pueden mezclarse; el catálogo de apoyo es parcial y no se realizó un modelado LTE conjunto.</div>', unsafe_allow_html=True)

with st.sidebar:
    st.markdown("## Skyntists")
    st.caption("Cuaderno de campo digital · IRAS 16293−2422")
    page = st.radio("Explorar", ["Inicio", "Espectros", "Evidencia molecular", "Método paso a paso", "Glosario y fuentes"], label_visibility="collapsed", key="wiki_page")
    st.divider()
    st.markdown("**Pregunta guía**")
    st.caption("¿Qué especies moleculares son compatibles con los rasgos espectrales observados?")
    st.caption(f"Datos cargados: {len(evidence)} filas de evidencia · {windows_n} ventanas")


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
        fig = px.bar(species_summary, x="especie", y="perfiles", color="especie", text="perfiles", title="Perfiles evaluados por especie", color_discrete_sequence=["#176b87", "#7656a5", "#bd6b2d"])
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
    palette = {"B_offset": "#087f8c", "B_pico": "#3568a8", "A": "#bd6b2d"}
    for key in series:
        value_col, label = POSITIONS[key]
        fig.add_trace(go.Scatter(x=df[freq_key], y=df[value_col], mode="lines", name=label, line=dict(color=palette[key], width=1.35 if len(series)>1 else 1.8), hovertemplate="ν = %{x:.6f} GHz<br>I = %{y:.4f} Jy/beam<extra>%{fullData.name}</extra>"))
    if candidate is not None:
        fit_ghz = float(candidate["frecuencia_pico_LSRK_MHz"])/1000
        rest_shifted_ghz = float(candidate["frecuencia_reposo_MHz"])*(1 - V_SYSTEM_KMS/C_KMS)/1000
        intensity_col = POSITIONS[position][0]
        nearest_i = (df[freq_key] - fit_ghz).abs().idxmin()
        fit_intensity = float(df.loc[nearest_i, intensity_col])
        fig.add_trace(go.Scatter(x=[fit_ghz], y=[fit_intensity], mode="markers", name="Centro ajustado", marker=dict(size=11, color="#d49300", line=dict(width=2, color="#ffffff")), hovertemplate="Centro ajustado<br>ν = %{x:.6f} GHz<extra></extra>"))
        fig.add_vline(x=rest_shifted_ghz, line_dash="dot", line_color="#7656a5", line_width=1.7, annotation_text="Catálogo @ 2.7 km/s", annotation_position="top left", annotation_font_color="#684894")
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
        fig.add_hline(y=0, line_dash="dash", line_color="#65798a", opacity=.9, annotation_text="Coincidencia con v = 2.7 km/s", annotation_position="top left")
        fig.update_traces(jitter=0.22, pointpos=0, marker=dict(size=13, line=dict(width=1, color="#071522")))
        fig.update_layout(xaxis_title="Especie candidata", yaxis_title="ν ajuste − ν catálogo corregida (MHz)", showlegend=False)
        st.plotly_chart(plot_layout(fig, 500), width="stretch")
        st.caption("Residuo positivo: el centro ajustado queda a mayor frecuencia que la predicción. El residuo y v_frecuencia se derivan de estas dos frecuencias; no son mediciones independientes ni incertidumbres. No se dibujan barras de error porque la tabla no contiene incertidumbres del centro del ajuste.")
        chart, explanation = st.columns([1.4, 1])
        with chart:
            plotted = view.dropna(subset=["v_frecuencia_kms", "SNR", "FWHM_kms"])
            if not plotted.empty:
                figv = px.scatter(plotted, x="v_frecuencia_kms", y="SNR", color="especie", symbol="ventana", size="FWHM_kms", hover_name="grupo_rasgo", hover_data={"frecuencia_reposo_MHz": ":.3f", "frecuencia_pico_LSRK_MHz": ":.3f", "FWHM_kms": ":.2f", "mezcla": True}, title="Contraste local y cinemática de los perfiles")
                figv.add_vline(x=V_SYSTEM_KMS, line_dash="dash", line_color="#bd6b2d", annotation_text="v sistémica adoptada")
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
- **Archivos:** `W1_espectros.csv`, `W2_espectros.csv`, `Campos_de_evidencia.csv`, `procedencia_catalogos.json`.
- **Reproducibilidad:** código de análisis y ajustes detallados en el cuaderno `herramientas/Inicio_AstroHack.ipynb` del repositorio de proyecto.
""")
    st.warning("La lista de catálogo de apoyo es parcial y no constituye un inventario exhaustivo de contaminantes. Los datos públicos de esta app son espectros extraídos y una tabla de evidencia: no incluyen cubos FITS ni la presentación PowerPoint.")

st.divider()
st.caption("Skyntists · AstroHack SATEM 2026 · Análisis exploratorio. La compatibilidad espectral descrita aquí es tentativa y depende de los datos, catálogos y supuestos documentados.")
