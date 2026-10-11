from __future__ import annotations

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

from config import C_KMS, POSITIONS, V_SYSTEM_KMS
from data import load_evidence, load_spectrum
from ui import apply_theme, page_header, plot_layout, render_summary_metrics

st.set_page_config(page_title="Skyntists | AstroHack", layout="wide")
apply_theme()

try:
    evidence = load_evidence()
except Exception as exc:
    st.error(f"No se pudo leer la evidencia: {exc}")
    st.stop()

if evidence.empty:
    st.error("No está disponible Campos_de_evidencia.csv; revisa los archivos del despliegue.")
    st.stop()

render_summary_metrics(evidence)

with st.sidebar:
    st.markdown("## Skyntists")
    st.caption("Cuaderno de campo digital · IRAS 16293−2422")
    page = st.radio(
        "Explorar",
        ["Inicio", "Espectros", "Evidencia molecular", "Método paso a paso", "Glosario y fuentes"],
        label_visibility="collapsed",
        key="wiki_page",
    )
    st.divider()
    st.markdown("**Pregunta guía**")
    st.caption("¿Qué especies moleculares son compatibles con los rasgos espectrales observados?")
    st.caption(f"Datos cargados: {len(evidence)} filas de evidencia · {evidence['ventana'].nunique()} ventanas")

def species_summary_table(evidence: pd.DataFrame) -> pd.DataFrame:
    if "especie" not in evidence.columns:
        return pd.DataFrame()
    summary = evidence.groupby("especie", dropna=False).agg(
        perfiles=("especie", "size"),
        decision=("decision", "first"),
    ).reset_index()
    if "SNR" in evidence.columns:
        med = (
            evidence.groupby("especie", dropna=False)["SNR"]
            .median()
            .round(1)
            .rename("mediana_SNR")
            .reset_index()
        )
        summary = summary.merge(med, on="especie", how="left")
    return summary

def render_home():
    page_header("Inicio")
    st.markdown(
        '<div class="callout"><b>Resultado actual:</b> seis perfiles muestran compatibilidad espectral tentativa con CH₃CN o CH₃CHO. La coincidencia aislada de Glycine-I es insuficiente. La base de datos queda conservada como evidencia, no como identificación concluida.</div>',
        unsafe_allow_html=True,
    )
    st.markdown("La app funciona como wiki del proyecto: combina observación, comparación con catálogos, decisiones de análisis y una explicación explícita de qué no se puede concluir.")
    cols = st.columns(3)
    with cols[0]:
        st.markdown('<div class="card"><span class="small-label">01 · Observar</span><br><b>Espectros calibrados</b><p>Explora W1/W2, elige una extracción y acerca el eje de frecuencia a una línea candidata.</p></div>', unsafe_allow_html=True)
    with cols[1]:
        st.markdown('<div class="card"><span class="small-label">02 · Contrastar</span><br><b>Catálogo y ajuste</b><p>Compara el centro ajustado con la frecuencia de catálogo corregida por una velocidad sistémica.</p></div>', unsafe_allow_html=True)
    with cols[2]:
        st.markdown('<div class="card"><span class="small-label">03 · Calificar</span><br><b>Evidencia, no veredicto</b><p>Revisa SNR, ancho, mezcla y alternativas. Una coincidencia de frecuencia por sí sola no basta.</p></div>', unsafe_allow_html=True)

    st.markdown("### Inventario de rasgos por especie")
    species_summary = species_summary_table(evidence)
    if not species_summary.empty:
        a, b = st.columns([1.15, 1])
        with a:
            fig = px.bar(
                species_summary,
                x="especie",
                y="perfiles",
                color="especie",
                text="perfiles",
                title="Perfiles evaluados por especie",
                color_discrete_sequence=["#176b87", "#7656a5", "#bd6b2d"],
            )
            fig.update_traces(textposition="outside", marker_line_width=0)
            fig.update_layout(showlegend=False, xaxis_title=None, yaxis_title="Número de perfiles")
            st.plotly_chart(plot_layout(fig, 365), width="stretch")
        with b:
            st.markdown("#### Cómo interpretar el SNR")
            st.write("El SNR informado es el contraste local aproximado del perfil respecto al ruido estimado. Un SNR alto señala una estructura destacada en ese intervalo; no mide la probabilidad de que la asignación sea correcta.")
            st.info("Identificación robusta requiere varias transiciones independientes, consistencia entre líneas, tratamiento de blends y comparación con alternativas. La app conserva esta distinción para no sobreinterpretar la evidencia.")
    if {"especie", "ventana", "decision"}.issubset(evidence.columns):
        by_window = evidence.groupby(["ventana", "especie", "decision"], dropna=False).size().reset_index(name="perfiles")
        fig_status = px.bar(
            by_window,
            x="especie",
            y="perfiles",
            color="decision",
            facet_col="ventana",
            barmode="stack",
            hover_data={"decision": True, "ventana": True},
            title="Estado de evidencia · separado por ventana",
        )
        fig_status.update_layout(xaxis_title=None, yaxis_title="Número de perfiles", legend_title="Evaluación")
        st.plotly_chart(plot_layout(fig_status, 390), width="stretch")
        st.caption("La separación por ventana evita sumar observaciones con haces y escalas instrumentales diferentes.")
    st.markdown("### Ruta de lectura")
    st.write("1. Empieza en **Espectros** y selecciona un rasgo; 2. revisa la comparación de frecuencias en **Evidencia molecular**; 3. sigue las ecuaciones y supuestos en **Método paso a paso**; 4. consulta **Glosario y fuentes** cuando necesites definir términos o ubicaciones.")

def render_spectra():
    page_header("Espectros")
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
        candidate_idx = st.selectbox(
            "Línea de referencia",
            candidate_options,
            format_func=lambda idx: (
                f"{evidence.loc[idx, 'especie']} · {evidence.loc[idx, 'grupo_rasgo']} · "
                f"{evidence.loc[idx, 'entrada']}"
            ),
        )
        candidate = evidence.loc[candidate_idx]
    else:
        candidate = None

    fig = go.Figure()
    palette = {"B_offset": "#55d6c2", "B_pico": "#80b7ff", "A": "#ffbd69"}

    for key in series:
        value_col, label = POSITIONS[key]
        fig.add_trace(
            go.Scatter(
                x=df[freq_key],
                y=df[value_col],
                mode="lines",
                name=label,
                line=dict(color=palette[key], width=1.35 if len(series) > 1 else 1.8),
                hovertemplate="ν = %{x:.6f} GHz<br>I = %{y:.4f} Jy/beam<extra>%{fullData.name}</extra>",
            )
        )

    if candidate is not None:
        fit_ghz = float(candidate["frecuencia_pico_LSRK_MHz"]) / 1000
        rest_shifted_ghz = float(candidate["frecuencia_reposo_MHz"]) * (1 - V_SYSTEM_KMS / C_KMS) / 1000
        intensity_col = POSITIONS[position][0]
        nearest_i = (df[freq_key] - fit_ghz).abs().idxmin()
        fit_intensity = float(df.loc[nearest_i, intensity_col])
        fig.add_trace(
            go.Scatter(
                x=[fit_ghz],
                y=[fit_intensity],
                mode="markers",
                name="Centro ajustado",
                marker=dict(size=11, color="#d49300", line=dict(width=2, color="#ffffff")),
                hovertemplate="Centro ajustado: %{x:.6f} GHz<extra></extra>",
            )
        )
        fig.add_vline(
            x=rest_shifted_ghz,
            line_dash="dot",
            line_color="#c7a5ff",
            line_width=1.7,
            annotation_text="Catálogo @ 2.7 km/s",
            annotation_position="top left",
            annotation_font_color="#d9c5ff",
        )

        if view_mode == "Región de línea":
            fwhm = candidate["FWHM_kms"]
            fwhm = 1.0 if pd.isna(fwhm) else float(fwhm)
            half_span = max(0.006, 5 * fwhm * float(candidate["frecuencia_reposo_MHz"]) / (C_KMS * 1000))
            fig.update_xaxes(range=[fit_ghz - half_span, fit_ghz + half_span])

        st.caption(
            f"Selección: {candidate['especie']} · {candidate['grupo_rasgo']}. Marcador dorado: centro del ajuste; línea violeta: frecuencia de catálogo trasladada a LSRK con v = {V_SYSTEM_KMS} km/s."
        )

    fig.update_layout(
        title=f"W{window} · {POSITIONS[position][1] if not compare else 'comparación de posiciones'}",
        xaxis_title="Frecuencia observada LSRK (GHz)",
        yaxis_title="Intensidad (Jy/beam)",
    )

    st.plotly_chart(plot_layout(fig, 540), width="stretch")
    st.markdown("#### ¿Qué representa cada extracción?")
    st.write("Las posiciones corresponden a espectros extraídos del producto procesado. B_offset es la posición desplazada recomendada para evaluar emisión menos contaminada por el continuo brillante; B_pico sigue el máximo del continuo; A es una extracción alternativa.")
    st.download_button("Descargar espectro de la ventana", df.to_csv(index=False).encode("utf-8"), f"W{window}_espectros.csv", "text/csv")

def render_evidence():
    page_header("Evidencia molecular")
    st.markdown("La comparación directa de frecuencias absolutas es engañosa a esta escala: ambas rondan 349 000 MHz, mientras la separación relevante es de alrededor de 0–2 MHz. Aquí mostramos la diferencia respecto a la predicción cinemática y la relación con SNR.")
    species = sorted(evidence["especie"].dropna().unique())
    chosen = st.multiselect("Especies incluidas", species, default=species)
    windows = sorted(evidence["ventana"].dropna().unique())
    chosen_windows = st.multiselect("Ventanas incluidas", windows, default=windows)

    view = evidence[evidence["especie"].isin(chosen) & evidence["ventana"].isin(chosen_windows)].copy()
    if view.empty:
        st.warning("No hay filas para los filtros seleccionados.")
        return

    fig = px.strip(
        view,
        x="especie",
        y="residuo_MHz_vs_vsys",
        color="especie",
        hover_name="grupo_rasgo",
        hover_data={
            "frecuencia_reposo_MHz": ":.4f",
            "frecuencia_pico_LSRK_MHz": ":.4f",
            "v_frecuencia_kms": ":.2f",
            "ventana": True,
        },
    )
    fig.add_hline(y=0, line_dash="dash", line_color="#c7a5ff", opacity=.9, annotation_text="Coincidencia con v = 2.7 km/s", annotation_position="top left", annotation_font_color="#d9c5ff")
    fig.update_traces(jitter=0.22, pointpos=0, marker=dict(size=13, line=dict(width=1, color="#08111c")))
    fig.update_layout(xaxis_title="Especie candidata", yaxis_title="ν ajuste − ν catálogo corregida (MHz)", showlegend=False)
    st.plotly_chart(plot_layout(fig, 500), width="stretch")
    st.caption("Residuo positivo: el centro ajustado queda a mayor frecuencia que la predicción. El residuo y v_frecuencia se derivan de estas dos frecuencias; no son mediciones independientes ni una probabilidad molecular.")

    chart, explanation = st.columns([1.4, 1])
    with chart:
        plotted = view.dropna(subset=["v_frecuencia_kms", "SNR", "FWHM_kms"])
        if not plotted.empty:
            figv = px.scatter(
                plotted,
                x="v_frecuencia_kms",
                y="SNR",
                color="especie",
                symbol="ventana",
                size="FWHM_kms",
                hover_name="grupo_rasgo",
                hover_data={
                    "frecuencia_reposo_MHz": ":.3f",
                    "frecuencia_pico_LSRK_MHz": ":.3f",
                    "ventana": True,
                },
            )
            figv.add_vline(x=V_SYSTEM_KMS, line_dash="dash", line_color="#bd6b2d", annotation_text="v sistémica adoptada")
            figv.update_layout(xaxis_title="Velocidad radial derivada (km/s)", yaxis_title="SNR local aproximado")
            st.plotly_chart(plot_layout(figv, 430), width="stretch")
            if len(plotted) < len(view):
                st.caption("Se excluyen de esta gráfica los perfiles sin SNR o FWHM medidos; esas filas permanecen en el residuo y la tabla.")
    with explanation:
        st.markdown("#### Lectura de la evidencia")
        st.write("En la gráfica de SNR, el tamaño del símbolo representa el FWHM ajustado. Un perfil estrecho y de alto contraste puede ser una buena línea candidata, pero congestión y coincidencias con líneas del entorno pueden falsear la lectura.")
        st.info("El candidato Glycine-I se conserva para mostrar por qué una coincidencia aislada en una zona congestionada no basta: faltan transiciones independientes que respalden la asignación.")

    table_cols = [
        c for c in [
            "especie",
            "catalogo",
            "entrada",
            "frecuencia_reposo_MHz",
            "frecuencia_pico_LSRK_MHz",
            "residuo_MHz_vs_vsys",
            "v_frecuencia_kms",
            "ventana",
            "FWHM_kms",
            "SNR",
            "mezcla",
            "decision",
            "grupo_rasgo",
        ] if c in view.columns
    ]
    st.markdown("#### Tabla trazable")
    st.dataframe(view[table_cols], width="stretch", hide_index=True)
    st.download_button("Descargar evidencia filtrada", view.to_csv(index=False).encode("utf-8"), "evidencia_filtrada.csv", "text/csv")

def render_method():
    page_header("Método paso a paso")
    st.write("Esta sección documenta qué transforma el análisis y qué queda pendiente. Abre los pasos para ver las ecuaciones y los supuestos.")
    with st.expander("1 · Datos de entrada y marco de referencia", expanded=True):
        st.write("Se exploran dos ventanas espectrales ALMA del proyecto ADS/JAO.ALMA#2013.1.00278.S, productos procesados ARI-L. Las frecuencias están en LSRK y las intensidades en Jy/beam. Los CSV contienen la señal observada y la evidencia de análisis.")
    with st.expander("2 · Frecuencia de reposo a frecuencia observada"):
        st.latex(r"\nu_{\rm pred,LSRK}=\nu_{\rm rest}\left(1-\frac{v_{\rm sys}}{c}\right),\qquad v_{\rm sys}=2.7\;\mathrm{km\,s^{-1}}")
        st.write("Se usa la aproximación no relativista para velocidades pequeñas frente a c. La corrección se aplica una sola vez porque el espectro ya está etiquetado en LSRK. El valor sistémico se toma como valor de trabajo, no como identificación.")
    with st.expander("3 · Residuo que se representa"):
        st.latex(r"\Delta\nu=\nu_{\rm fit,LSRK}-\nu_{\rm pred,LSRK}")
        st.write("El gráfico principal centra cada candidata en una escala local en MHz. El residuo indica concordancia con la velocidad sistémica adoptada, pero no es independiente de los datos de entrada ni una medida absoluta de certeza.")
    with st.expander("4 · Ruido robusto, ajuste y SNR"):
        st.write("La sigma aproximada usa la dispersión robusta (MAD) de primeras diferencias dividida entre √2. El perfil se ajusta localmente con componente gaussiana y baseline. El SNR es contra la señal local estimada, no una probabilidad global.")
    with st.expander("5 · Área integrada y agrupación"):
        st.write("El área se aproxima integrando el perfil sobre el ancho seleccionado; su error aproximado no incluye correlación, calibración ni incertidumbre de modelo. Varias transiciones pueden agruparse por especie y ventana para comparar patrones.")
    with st.expander("6 · Qué exige una identificación más fuerte"):
        st.write("Buscar múltiples transiciones no confundidas de la misma especie, modelar blends y alternativas, reportar incertidumbres de frecuencia y parámetros de ajuste, probar robustez a selección de baseline y validar con comparaciones independientes.")
    st.markdown("### ¿Por qué importan los coeficientes de Einstein?")
    st.write("El coeficiente Aᵤₗ (s⁻¹) describe la probabilidad de emisión espontánea por unidad de tiempo para una transición entre niveles u y l. Junto a la degeneración estadística, energia y población, define la intensidad relativa de una línea.")
    st.warning("El catálogo local no entrega Aᵤₗ para estas filas. La columna log10_I_300K_nm2_MHz representa intensidad de catálogo bajo condiciones definidas; no es log10(Aᵤₗ). No rellenamos esa interpretación como si fuera una medición física directa.")

def render_glossary():
    page_header("Glosario y fuentes")
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
    st.warning("La lista de catálogo de apoyo es parcial y no constituye un inventario exhaustivo de contaminantes. Los datos públicos de esta app son espectros extraídos y una tabla de evidencia: no es un catálogo completo ni una clasificación definitiva.")

def main():
    pages = {
        "Inicio": render_home,
        "Espectros": render_spectra,
        "Evidencia molecular": render_evidence,
        "Método paso a paso": render_method,
        "Glosario y fuentes": render_glossary,
    }
    pages[page]()

if __name__ == "__main__":
    main()
