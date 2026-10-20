# -*- coding: utf-8 -*-
"""Rutas del análisis. Todos los scripts las toman de aquí.

Todo es relativo a la carpeta raíz del repositorio (un nivel por encima de scripts/):
    datos/csv/<día>/*.csv           CSV de DeepLabCut (entrada)
    datos/INDICE_videos.csv         cuadros por segundo de cada sesión (entrada)
    resultados/intermedios/         tablas de ciclos y sesiones, estadistica.pkl (se generan)
    resultados/Analisis_marcha_v2/  tablas finales, figuras y Excel (se generan)

Para usar otra carpeta raíz sin modificar el código, definir la variable de entorno MARCHA_RAIZ.
"""
import os

AQUI = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.environ.get('MARCHA_RAIZ', os.path.abspath(os.path.join(AQUI, '..')))
DATOS = os.path.join(RAIZ, 'datos')
CSV = os.path.join(DATOS, 'csv')
INDICE = os.path.join(DATOS, 'INDICE_videos.csv')
INTER = os.path.join(RAIZ, 'resultados', 'intermedios')
OUT = os.path.join(RAIZ, 'resultados', 'Analisis_marcha_v2')

for _d in (INTER, OUT):
    os.makedirs(_d, exist_ok=True)
