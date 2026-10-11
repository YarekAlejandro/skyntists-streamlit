# Skyntists · AstroHack SATEM 2026

App pública para explorar espectros W1/W2 de IRAS 16293−2422 y evidencia molecular.

## Despliegue

Esta raíz se publica en GitHub y se conecta en [Streamlit Community Cloud](https://share.streamlit.io/) con `app.py` como entrypoint. Cloud instala `requirements.txt`.

## Archivos

- `app.py`: aplicación Streamlit.
- `W1_espectros.csv`, `W2_espectros.csv`: espectros LSRK en Jy/beam.
- `Campos_de_evidencia.csv`: mediciones del notebook y decisiones tentativas.
- `procedencia_catalogos.json`: trazabilidad del catálogo local.

Datos de espectro y tabla de evidencia se hacen públicos con el código para que la visualización global sea reproducible. Este paquete no contiene cubos FITS ni presentación.

## Resultado y límites

Seis perfiles son compatibles tentativamente con CH₃CN/CH₃CHO; glicina-I permanece insuficiente. El SNR y los ajustes son descriptivos; blends y catálogo parcial limitan la identificación. No se infiere presencia de vida. El catálogo no incluye Aij; `log10_I_300K_nm2_MHz` no equivale a `log10(Aij)`.

Datos ALMA `ADS/JAO.ALMA#2013.1.00278.S`, productos ARI-L. Jørgensen et al. (2016), PILS, *A&A* 595, A117 ([arXiv:1607.08733](https://arxiv.org/abs/1607.08733)); [CDMS](https://cdms.astro.uni-koeln.de/classic/entries/) y [JPL](https://spec.jpl.nasa.gov/).
