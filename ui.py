from __future__ import annotations

import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

from config import PAGE_COPY


CSS = """
<style>
:root {color-scheme:dark; --ink:#e7eef7; --muted:#a9b9c9; --panel:#111e2c; --line:#26394c; --teal:#55d6c2; --blue:#80b7ff; --gold:#ffd166; --violet:#c7a5ff;}
.stApp {background:radial-gradient(ellipse at 78% -10%,rgba(39,112,143,.20),transparent 35%),#08111c;color:var(--ink);}
.block-container {max-width:1440px;padding-top:1.5rem;padding-bottom:3.5rem;}
html,body,[class*="css"] {font-family:Inter,"Segoe UI",Arial,sans-serif;color:var(--ink);}
h1,h2,h3 {font-family:Inter,"Segoe UI",Arial,sans-serif;color:#f0f5fb;letter-spacing:-.025em;}
h2,h3 {font-weight:650;}
.hero {position:relative;overflow:hidden;padding:1.8rem 2rem;border:1px solid #294157;border-left:4px solid var(--teal);border-radius:12px;background:linear-gradient(115deg,rgba(15,32,48,.98),rgba(18,36,55,.96));}
.callout {padding:.9rem 1.1rem;border:1px solid #294157;border-left:4px solid var(--teal);background:#101d2a;border-radius:7px;margin:.75rem 0 1.1rem;color:#d3dfeb;line-height:1.65;}
.small-label {color:#91a7ba;font-size:.75rem;text-transform:uppercase;letter-spacing:.09em;font-weight:700;}
.card {border:1px solid #2a4054;background:linear-gradient(145deg,#142333,#101b28);border-radius:9px;padding:1rem 1.1rem;min-height:128px;box-shadow:0 5px 18px rgba(0,0,0,.12);}
.card b {color:#ecf4fb;font-size:1.04rem;}
.card p {color:#b2c2d0;margin:.42rem 0 0;font-size:.92rem;line-height:1.6;}
.muted {color:var(--muted);}
div[data-testid="stMetric"] {background:#111e2c;border:1px solid #293e52;padding:.9rem 1rem;border-radius:9px;box-shadow:0 5px 18px rgba(0,0,0,.12);}
div[data-testid="stMetricLabel"] {color:#a9b9c9;}
div[data-testid="stMetricValue"] {font-family:Inter,"Segoe UI",Arial,sans-serif;color:#f0f5fb;font-weight:650;}
section[data-testid="stSidebar"] {background:#0d1926;border-right:1px solid #25384b;}
section[data-testid="stSidebar"] [data-testid="stMarkdownContainer"] p {color:#b7c7d5;}
button[data-baseweb="tab"] {font-weight:600;}
div[data-testid="stDataFrame"] {border:1px solid #293e52;border-radius:7px;overflow:hidden;}
.stDownloadButton button {border-radius:6px;font-weight:600;}
div[data-testid="stAlert"] {border:1px solid #34495c;}
</style>
"""


def apply_theme() -> None:
    st.markdown(CSS, unsafe_allow_html=True)


def page_header(page: str) -> None:
    kicker, title, intro = PAGE_COPY[page]
    st.markdown(
        f'<div class="section-head"><div class="eyebrow">{kicker}</div><h1>{title}</h1><p class="page-intro">{intro}</p></div>',
        unsafe_allow_html=True,
    )


def plot_layout(fig: go.Figure, height: int = 470) -> go.Figure:
    fig.update_layout(
        height=height,
        template="plotly_dark",
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="#0c1723",
        font=dict(family="Inter, Segoe UI, Arial, sans-serif", color="#dce7f1", size=12),
        margin=dict(l=28, r=26, t=64, b=34),
        legend=dict(orientation="h", yanchor="bottom", y=1.02, x=0),
        hoverlabel=dict(bgcolor="#152536", bordercolor="#55d6c2", font_color="#f0f5fb"),
    )
    fig.update_xaxes(gridcolor="#203244", zerolinecolor="#63788c", linecolor="#52677a", ticks="outside")
    fig.update_yaxes(gridcolor="#203244", zerolinecolor="#63788c", linecolor="#52677a", ticks="outside")
    return fig


def render_summary_metrics(evidence) -> None:
    confirmed = 0
    if "decision" in evidence.columns:
        confirmed = int(evidence["decision"].astype(str).str.contains("Tentativa", case=False).sum())

    species_n = evidence["especie"].nunique() if "especie" in evidence.columns else 0
    windows_n = evidence["ventana"].nunique() if "ventana" in evidence.columns else 0

    m1, m2, m3, m4 = st.columns(4)
    m1.metric("Ventanas ALMA", f"{windows_n}", help="W1 y W2; sus haces y escalas de intensidad son distintos.")
    m2.metric("Canales", "1,918 / ventana")
    m3.metric("Perfiles evaluados", f"{len(evidence)}", help="Incluye candidatos tentativos y el caso de glicina-I insuficiente.")
    m4.metric("Compatibilidad tentativa", f"{confirmed} perfiles", help="No son identificaciones confirmadas ni una probabilidad molecular.")


def build_species_summary(evidence):
    if "especie" not in evidence.columns:
        return None

    summary = (
        evidence.groupby("especie", dropna=False)
        .agg(perfiles=("especie", "size"), decision=("decision", "first"))
        .reset_index()
    )

    if "SNR" in evidence.columns:
        medians = (
            evidence.groupby("especie", dropna=False)["SNR"]
            .median()
            .round(1)
            .rename("mediana_SNR")
            .reset_index()
        )
        summary = summary.merge(medians, on="especie", how="left")

    return summary
