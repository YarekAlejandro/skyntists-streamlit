"""Streamlit explorer for the Skyntists AstroHack spectral evidence."""
from __future__ import annotations

from pathlib import Path

import pandas as pd
import streamlit as st

ROOT = Path(__file__).resolve().parent
SPECTRA_DIR = ROOT / "espectros" if (ROOT / "espectros").exists() else ROOT
EVIDENCE_PATH = ROOT / "Entregables" / "Campos_de_evidencia.csv"
if not EVIDENCE_PATH.exists():
    EVIDENCE_PATH = ROOT / "Campos_de_evidencia.csv"

st.set_page_config(
    page_title="Skyntists · AstroHack",
    page_icon="🔭",
    layout="wide",
    initial_sidebar_state="expanded",
)

# Restrained scientific styling: high contrast, explicit units, no chart smoothing.
st.markdown("""
<style>
.block-container {max-width: 1440px; padding-top: 2rem; padding-bottom: 3rem;}
.hero {padding: 1.4rem 1.6rem; border-radius: 16px; background: linear-gradient(115deg,#102a43,#137c8b); color: white; margin-bottom: 1.25rem;}
.hero h1 {color:white; margin-bottom:.25rem; font-size:2.25rem;}
.hero p {color:#d8edf0; margin:0; font-size:1.05rem;}
.note {padding: .85rem 1rem; border-left: 4px solid #d97745; background:#fff7ed; border-radius:6px;}
</style>
<div class="hero"><h1>🔭 Skyntists · AstroHack SATEM 2026</h1><p>Explorador de espectros y evidencia molecular de IRAS 16293−2422</p></div>
""", unsafe_allow_html=True)

@st.cache_data(show_spinner=False)
def load_spectrum(window: int) -> pd.DataFrame:
    path = SPECTRA_DIR / f"W{window}_espectros.csv"
    frame = pd.read_csv(path)
    required = {"frecuencia_LSRK_GHz", "B_offset_Jy_beam", "B_pico_Jy_beam", "A_Jy_beam"}
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"Faltan columnas en {path.name}: {', '.join(sorted(missing))}")
    return frame

@st.cache_data(show_spinner=False)
def load_evidence() -> pd.DataFrame:
    if not EVIDENCE_PATH.exists():
        return pd.DataFrame()
    frame = pd.read_csv(EVIDENCE_PATH)
    for col in ("frecuencia_reposo_MHz", "frecuencia_pico_LSRK_MHz", "FWHM_kms", "SNR", "area_Jy_beam_kms"):
        if col in frame:
            frame[col] = pd.to_numeric(frame[col], errors="coerce")
    return frame

try:
    evidence = load_evidence()
except Exception as exc:
    st.error(f"No se pudo cargar la tabla de evidencia: {exc}")
    st.stop()

c1, c2, c3, c4 = st.columns(4)
c1.metric("Ventanas", "W1 + W2")
c2.metric("Canales", "1,918 por ventana")
c3.metric("Rasgos agrupados", "64", help="34 en W1 y 30 en W2; 107 asociaciones de catálogo no equivalen a detecciones independientes.")
c4.metric("Estado", "Tentativo", help="Las coincidencias son candidatas, no identificaciones moleculares confirmadas.")

st.markdown('<div class="note"><b>Lectura científica:</b> seis perfiles se consideran compatibles tentativamente con CH₃CN y CH₃CHO. Glicina-I se clasifica como insuficiente. El SNR es contraste local aproximado, no probabilidad de identidad molecular.</div>', unsafe_allow_html=True)

tabs = st.tabs(["Espectros", "Candidatos", "Método y límites", "Fuentes"])
with tabs[0]:
    st.subheader("Explora las extracciones calibradas")
    left, right = st.columns([1, 2])
    with left:
        window = st.selectbox("Ventana espectral", [1, 2], format_func=lambda x: f"W{x}")
        position = st.selectbox("Posición", ["B_offset", "B_pico", "A"], index=0)
        show_all = st.checkbox("Comparar las tres posiciones", value=True)
    df = load_spectrum(window)
    spectrum_columns = [
        "B_offset_Jy_beam", "B_pico_Jy_beam", "A_Jy_beam"
    ] if show_all else [f"{position}_Jy_beam"]
    st.line_chart(
        df.set_index("frecuencia_LSRK_GHz")[spectrum_columns],
        x_label="Frecuencia observada LSRK (GHz)",
        y_label="Intensidad (Jy/beam)",
    )
    st.caption("Las ventanas tienen haces distintos: sus amplitudes no se comparan directamente como abundancias. No se suavizan ni desplazan los datos.")

with tabs[1]:
    st.subheader("Tabla de evidencia")
    if evidence.empty:
        st.warning("No se encontró Entregables/Campos_de_evidencia.csv.")
    else:
        species = sorted(evidence["especie"].dropna().unique()) if "especie" in evidence else []
        chosen = st.multiselect("Especies", species, default=species)
        view = evidence[evidence["especie"].isin(chosen)] if chosen else evidence.iloc[0:0]
        if {"frecuencia_reposo_MHz", "frecuencia_pico_LSRK_MHz"}.issubset(view.columns) and not view.empty:
            st.subheader("Frecuencia de catálogo y centro del ajuste")
            st.scatter_chart(
                view,
                x="frecuencia_reposo_MHz",
                y="frecuencia_pico_LSRK_MHz",
                color="especie" if "especie" in view else None,
            )
        cols = [c for c in ["especie","catalogo","entrada","frecuencia_reposo_MHz","ventana","frecuencia_pico_LSRK_MHz","v_radio_kms","FWHM_kms","SNR","mezcla","decision"] if c in view.columns]
        st.dataframe(view[cols], use_container_width=True, hide_index=True)
        st.download_button("Descargar CSV de evidencia", view.to_csv(index=False).encode("utf-8"), "Campos_de_evidencia_filtrado.csv", "text/csv")

with tabs[2]:
    st.subheader("Cómo leer los resultados")
    st.markdown("""
- **Datos:** productos procesados ARI-L del proyecto `ADS/JAO.ALMA#2013.1.00278.S`, marco LSRK, Jy/beam; B_offset es la extracción recomendada por la guía.
- **Comparación Doppler:** para B se parte de 2.7 km/s y se aplica una sola vez (`ν_obs = ν_rest (1 − v_radio/c)`).
- **Ruido y ajuste:** MAD de diferencias entre canales y perfil gaussiano local con baseline son aproximaciones descriptivas. La correlación entre canales, blends y regiones con líneas limitan el SNR.
- **Conteo:** varias transiciones pueden coincidir con un mismo rasgo. 107 asociaciones se agrupan en 64 rasgos; no son 107 detecciones.
- **Conclusión:** las señales apoyan hipótesis tentativas. No se infieren abundancias, temperatura ni presencia de vida.

**Einstein Aij.** El coeficiente (s⁻¹) caracteriza la probabilidad de emisión espontánea y es útil, junto con degeneraciones, energías y función de partición, para modelar excitación y transferencia radiativa. El catálogo local no proporciona Aij; `log10_I_300K_nm2_MHz` no es `log10(Aij)`. La app no infiere ni simula valores ausentes.
""")

with tabs[3]:
    st.subheader("Procedencia y referencias")
    st.markdown("""
- Datos ALMA: `ADS/JAO.ALMA#2013.1.00278.S`, productos procesados ARI-L.
- Jørgensen et al. (2016), PILS, *A&A*, 595, A117. [arXiv:1607.08733](https://arxiv.org/abs/1607.08733).
- [Catálogo CDMS](https://cdms.astro.uni-koeln.de/classic/entries/) · [Catálogo JPL](https://spec.jpl.nasa.gov/).
- Identificadores, hashes y procedencia local: `procedencia_catalogos.json`.

El catálogo de apoyo es parcial y no constituye un inventario exhaustivo de contaminantes espectrales. La app visualiza los archivos entregados; la ejecución reproducible y los ajustes detallados siguen en `herramientas/Inicio_AstroHack.ipynb`.
""")

with st.sidebar:
    st.markdown("### Skyntists")
    st.markdown("**Pregunta:** ¿qué especies son compatibles con los rasgos observados?")
    st.caption("Análisis científico exploratorio · AstroHack SATEM 2026")
    if not evidence.empty and "decision" in evidence:
        st.caption(f"Filas de evidencia: {len(evidence)}")
