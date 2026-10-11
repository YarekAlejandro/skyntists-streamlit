# Skyntists · AstroHack SATEM 2026

Wiki interactiva y explorador científico para inspeccionar los espectros de IRAS 16293−2422 y la evidencia usada al evaluar candidatas moleculares.

**Aplicación en vivo:** [skyntists-astrohack.streamlit.app](https://skyntists-astrohack.streamlit.app/) · **Código y datos de despliegue:** [GitHub](https://github.com/YarekAlejandro/skyntists-streamlit).

## Cómo recorrer la wiki

- **Inicio:** pregunta científica, síntesis de resultados y ruta de lectura.
- **Espectros:** W1/W2 con selector de extracción, zoom alrededor de un centro ajustado, modo de banda completa y referencia de catálogo corregida a LSRK.
- **Evidencia molecular:** filtros por especie/ventana, residuo de frecuencia contra la velocidad sistémica adoptada, exploración de SNR/ancho y tabla descargable.
- **Método paso a paso:** sistema de referencia, corrección Doppler, residuo, estimación de ruido, ajustes, agrupación, límites y contexto de los coeficientes de Einstein.
- **Glosario y fuentes:** términos, procedencia de datos, catálogos y referencias.

Los tooltips conservan las frecuencias originales. Las gráficas no ocultan que algunos resultados son tentativos ni inventan incertidumbres ausentes en los archivos.

## Ejecutar localmente

```bash
python -m venv .venv
# Windows: .venv\\Scripts\\activate
# macOS/Linux: source .venv/bin/activate
pip install -r requirements.txt
streamlit run app.py
```

Para despliegue, conectar este repositorio público en [Streamlit Community Cloud](https://share.streamlit.io/) y seleccionar `app.py`; el servicio instala `requirements.txt`.

## Archivos de entrada

- `W1_espectros.csv`, `W2_espectros.csv`: extracciones procesadas, frecuencias LSRK (GHz), intensidades (Jy/beam).
- `Campos_de_evidencia.csv`: centros, velocidades derivadas, anchos, SNR local, mezcla y decisiones exploratorias.
- `procedencia_catalogos.json`: procedencia y trazabilidad de catálogos.
- `app.py`: aplicación wiki en Streamlit.

Este repositorio es el paquete público mínimo de visualización. No contiene cubos FITS ni la presentación PowerPoint. El análisis reproducible detallado permanece en el repositorio del proyecto, en `herramientas/Inicio_AstroHack.ipynb`.

## Lectura científica y limitaciones

Seis perfiles son compatibles tentativamente con CH₃CN o CH₃CHO en el análisis actual; Glycine-I queda como insuficiente por la congestión y la falta de confirmación mediante otras transiciones independientes. 107 asociaciones de catálogo se agruparon en 64 rasgos; no son 107 detecciones independientes.

La comparación usa una velocidad sistémica de referencia de 2.7 km/s para llevar las frecuencias de reposo a LSRK mediante la aproximación no relativista. El residuo y la velocidad por frecuencia son magnitudes derivadas de los mismos valores de catálogo/ajuste. El SNR es un contraste local aproximado, no una probabilidad molecular ni una significancia global. El inventario espectral es parcial y existen blends; no se infieren abundancias, condiciones físicas ni presencia de vida.

El catálogo local no proporciona coeficientes Einstein Aᵤₗ para estas filas. `log10_I_300K_nm2_MHz` no equivale a `log10(Aᵤₗ)` y la app no sustituye datos faltantes por estimaciones.

## Procedencia

Datos ALMA del proyecto `ADS/JAO.ALMA#2013.1.00278.S`, productos procesados ARI-L. Contexto: Jørgensen et al. (2016), PILS, *A&A*, 595, A117 ([arXiv:1607.08733](https://arxiv.org/abs/1607.08733)); [CDMS](https://cdms.astro.uni-koeln.de/classic/entries/) y [JPL](https://spec.jpl.nasa.gov/).
