# Gait_Analisis_SODvsWT

Análisis de la marcha lateral en cinta de ratones WT y SOD1 entre P30 y P85. El punto de partida son las coordenadas de seis puntos del miembro posterior (cresta ilíaca, cadera, rodilla, tobillo, articulación metatarsofalángica y punta del dedo) estimadas con DeepLabCut. A partir de ellas se miden los ángulos articulares cuadro a cuadro, se detectan los ciclos de marcha, se calculan las variables de cada ciclo y se comparan los dos genotipos.

Todo el análisis se hace a partir de los CSV de DeepLabCut; no se necesitan los videos.

## Contenido

```
datos/
  csv/<día>/                 CSV de DeepLabCut, uno por animal y día (71 sesiones)
  INDICE_videos.csv          cuadros por segundo de cada sesión
scripts/
  rutas.py                   rutas de entrada y salida
  marcha_v2.py               coordenadas -> ángulos -> ciclos -> variables por ciclo
  validar_60fps.py           validación del análisis a 30 fps con las sesiones grabadas a 60 fps
  qc_plot.py                 figura de control de la detección de ciclos
  estadistica_v2.py          resumen por animal, pruebas por día y modelo mixto
  exportar_v2.py             tablas por variable, figuras y libro de Excel
  fig_rodilla.py             figura de la goniometría de rodilla
  R/                         modelo mixto en R (lme4, lmerTest, emmeans)
ejecutar_analisis.py         encadena todos los pasos
requirements.txt             librerías de Python
```

## Requisitos

Python 3.13 con las librerías de `requirements.txt` (las versiones indicadas son con las que se trabajó):

| Librería | Versión |
|---|---|
| numpy | 2.5.3 |
| pandas | 3.0.5 |
| scipy | 1.18.1 |
| statsmodels | 0.15.0 |
| matplotlib | 3.11.2 |
| openpyxl | 3.1.5 |

R 4.3 o posterior, solo para el análisis estadístico final, con los paquetes `lme4`, `lmerTest`, `emmeans`, `ggplot2` y `pbkrtest`. El script `scripts/R/instalar_paquetes_R.R` los instala si faltan.

Probado en Linux, en local y en una máquina limpia de GitHub Actions. Las rutas no dependen del sistema operativo. No hace falta GPU.

## Cómo ejecutarlo

```
pip install -r requirements.txt
python ejecutar_analisis.py
```

Tarda alrededor de un minuto y medio. Si `Rscript` no está disponible, el script ejecuta solo la parte de Python y lo avisa; también se puede omitir la parte de R con `python ejecutar_analisis.py --sin-r`. En ese caso, el paso de R se corre aparte:

```
Rscript scripts/R/instalar_paquetes_R.R
Rscript scripts/R/analisis_rodilla.R resultados/Analisis_marcha_v2
Rscript scripts/R/todas_las_variables_LRT.R resultados/Analisis_marcha_v2
```

Los pasos de Python también se pueden correr uno por uno desde la raíz del repositorio, en este orden:

```
python scripts/marcha_v2.py --crudo                 # versión principal: sin suavizar, 30 fps
python scripts/marcha_v2.py                         # sensibilidad: suavizado Savitzky-Golay de 5 cuadros
python scripts/marcha_v2.py --crudo --sin-decimar   # sin reducir las sesiones de 60 fps
python scripts/marcha_v2.py --sin-decimar
python scripts/validar_60fps.py
python scripts/qc_plot.py
python scripts/estadistica_v2.py > resultados/intermedios/estadistica_salida.txt
python scripts/exportar_v2.py
python scripts/fig_rodilla.py
```

## Datos de entrada

`datos/csv/<día>/` contiene un CSV por sesión, con el formato habitual de DeepLabCut: tres filas de encabezado (`scorer`, `bodyparts`, `coords`) y una fila por cuadro, con `x`, `y` y `likelihood` de los marcadores `crest`, `hip`, `knee`, `ankle`, `foot` y `toe`. Las coordenadas están en píxeles de la imagen original. Los tres primeros caracteres del nombre del archivo son el número del animal.

Hay 12 animales y 6 días (P30, P37, P44, P51, P65 y P85), con un total de 71 sesiones: el animal 866 no tiene sesión en P30. Los genotipos están definidos en `marcha_v2.py` (`GENO`):

- WT: 821, 860, 861, 863, 866, 867
- SOD1: 856, 857, 859, 862, 864, 865

Las sesiones de P85 y la 821_P51 provienen de una iteración posterior del modelo de DeepLabCut; por eso la estadística incluye también un análisis sin P85.

`datos/INDICE_videos.csv` (separador `;`) indica los cuadros por segundo reales de cada sesión. Dos sesiones (856_P30 y 859_P30) se grabaron a 60 fps y se reducen a 30 fps para que todas se midan igual; el resto está a 30 fps o a 600/19 fps (unos 31,6). P51 está en 1280x720 y se escala a 1920x1080.

## Método

El encabezado de `scripts/marcha_v2.py` describe el procedimiento completo, paso por paso, con los criterios numéricos. En resumen:

1. Se leen las coordenadas y se descartan los puntos con `likelihood` menor que 0,70. Solo se interpolan huecos de hasta 2 cuadros. No se suaviza.
2. Las coordenadas se pasan a un plano cartesiano común (y hacia arriba y el animal mirando siempre hacia +x).
3. En cada cuadro se calculan los ángulos de rodilla, tobillo, pie y cadera con `atan2` a partir de los vectores que salen de cada articulación. En las cuatro articulaciones, un ángulo menor significa más flexión. Además se reconstruye una rodilla triangulada a partir de los largos medianos de fémur y tibia.
4. El período de la zancada de cada sesión se obtiene por autocorrelación de la altura del dedo.
5. Un paso se acepta solo si coinciden la subida del dedo y la flexión de rodilla y de cadera, cada una sobre un umbral del 30 % de la mediana de la sesión.
6. Un ciclo va de un paso aceptado al siguiente y es válido si dura entre 0,6 y 1,6 veces el período y el dedo se eleva más del 30 % de la mediana de la sesión.
7. Para cada ciclo válido se calculan las variables de la tabla siguiente. La elevación del dedo se pasa a milímetros con la velocidad de la cinta (20 cm/s), que se usa como referencia de escala en el apoyo.
8. La unidad experimental es el animal: se promedian los ciclos válidos de cada animal y día.

### Variables

| Código | Variable |
|---|---|
| V1 | Excursión de tobillo |
| V2 | Excursión de pie (MTP) |
| V3a, V3b | Excursión de rodilla (marcador de DeepLabCut y rodilla triangulada) |
| V4a, V4b, V4c | Rodilla: extensión máxima, flexión máxima y ángulo medio |
| V4d, V4e, V4f | Tobillo: extensión máxima, flexión máxima y ángulo medio |
| V4g, V4h, V4i | Pie: extensión máxima, flexión máxima y ángulo medio |
| V5a | Frecuencia de paso |
| V5b, V5c | Tiempo de apoyo y de oscilación |
| X1 | Proporción de apoyo (duty) |
| V6 | Elevación del dedo (mm) |

### Estadística

Para cada variable se ajusta un modelo lineal mixto, `y ~ genotipo * día + (1 | animal)`, y se obtienen pruebas de razón de verosimilitud (ajuste por máxima verosimilitud) para la interacción genotipo x día y para el efecto del genotipo, con corrección de Benjamini-Hochberg entre variables. Por día se comparan los genotipos con la prueba de Welch, Mann-Whitney exacto y g de Hedges, con corrección de Holm sobre los 6 días. Hay análisis de sensibilidad sin P85, con los 10 ciclos más extremos de cada sesión (los mayores en excursiones y extensión máxima, los menores en flexión máxima), con el suavizado y con la calibración de la cinta por sesión.

En R (`scripts/R`) se repite el modelo de la goniometría de rodilla con `lme4` y `lmerTest`: prueba F tipo III con grados de libertad de Kenward-Roger, contrastes por día con `emmeans`, y los supuestos (Shapiro-Wilk de los residuos y Brown-Forsythe). El segundo script compara las pruebas de razón de verosimilitud de R con las de Python para todas las variables.

## Resultados generados

Todo se escribe en `resultados/`, que no se versiona:

- `resultados/intermedios/`: tablas de sesiones, ciclos y pasos candidatos de cada versión del análisis (`sesiones_*.csv`, `ciclos_*.csv`, `eventos_*.csv`), `validacion_30_vs_60fps.csv` y la salida de la estadística.
- `resultados/Analisis_marcha_v2/`: `Resultados_marcha_v2.xlsx`, `RESUMEN_todas_las_variables.csv` y una carpeta por grupo de variables con la tabla por animal y día, la comparación WT vs SOD1 por día y la figura. La carpeta `QC/` guarda la validación y los controles, y `R/` los resultados del análisis en R.

## Reproducibilidad

El análisis no tiene componentes aleatorios, salvo el desplazamiento horizontal de los puntos en algunas figuras, que usa una semilla fija. Con las versiones indicadas arriba, dos ejecuciones sobre los mismos CSV producen tablas idénticas.

Cada cambio en el repositorio dispara en GitHub Actions (`.github/workflows/reproducir.yml`) una ejecución completa, con Python y R, sobre una máquina limpia. La ejecución falla si algún paso termina con error o si no se generan las salidas esperadas.

## Autoría

Martín Volpi Canales.
