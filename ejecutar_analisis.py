# -*- coding: utf-8 -*-
"""Ejecuta el análisis completo a partir de los CSV de DeepLabCut (datos/csv).

    python ejecutar_analisis.py            Python y, si está instalado, R
    python ejecutar_analisis.py --sin-r    solo la parte de Python

Los pasos son los mismos que se pueden correr a mano desde la carpeta scripts/; aquí solo se
encadenan en el orden correcto. Salida en resultados/.
"""
import argparse
import os
import shutil
import subprocess
import sys
import time

RAIZ = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(RAIZ, 'scripts')
INTER = os.path.join(RAIZ, 'resultados', 'intermedios')
OUT = os.path.join(RAIZ, 'resultados', 'Analisis_marcha_v2')
ENV = dict(os.environ, PYTHONUTF8='1', MPLBACKEND='Agg', PYTHONDONTWRITEBYTECODE='1')


def correr(titulo, cmd, salida=None):
    print(titulo, flush=True)
    t0 = time.time()
    os.makedirs(os.path.dirname(salida), exist_ok=True) if salida else None
    destino = open(salida, 'w', encoding='utf-8') if salida else subprocess.DEVNULL
    try:
        r = subprocess.run(cmd, cwd=RAIZ, env=ENV, stdout=destino, stderr=subprocess.PIPE,
                           text=True, encoding='utf-8')
    finally:
        if salida:
            destino.close()
    if r.returncode != 0:
        print(r.stderr, file=sys.stderr)
        sys.exit(f'El paso terminó con error: {" ".join(cmd)}')
    print(f'    listo ({time.time() - t0:.0f} s)', flush=True)


def py(nombre, *args):
    return [sys.executable, os.path.join(SCRIPTS, nombre), *args]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--sin-r', action='store_true', help='omite el análisis en R')
    a = ap.parse_args()

    if not os.path.isdir(os.path.join(RAIZ, 'datos', 'csv')):
        sys.exit('No encuentro datos/csv')

    correr('[1/9] Ciclos y variables: sin suavizar, 30 fps (versión principal)', py('marcha_v2.py', '--crudo'))
    correr('[2/9] Ciclos y variables: suavizado SG5, 30 fps', py('marcha_v2.py'))
    correr('[3/9] Ciclos y variables: sin suavizar, 60 fps sin reducir', py('marcha_v2.py', '--crudo', '--sin-decimar'))
    correr('[4/9] Ciclos y variables: suavizado SG5, 60 fps sin reducir', py('marcha_v2.py', '--sin-decimar'))
    correr('[5/9] Validación 30 vs 60 fps', py('validar_60fps.py'))
    correr('[6/9] Figura de control de la detección de ciclos', py('qc_plot.py'))
    correr('[7/9] Estadística (resumen por animal, pruebas por día, modelo mixto)', py('estadistica_v2.py'),
           os.path.join(INTER, 'estadistica_salida.txt'))
    correr('[8/9] Tablas, carpetas por variable y Excel', py('exportar_v2.py'))
    correr('[9/9] Figura de rodilla', py('fig_rodilla.py'))

    if a.sin_r:
        print('Se omite el análisis en R (--sin-r).')
        return
    rscript = shutil.which('Rscript')
    if rscript is None:
        print('No se encontró Rscript: se omite el análisis en R. Ver README (sección Requisitos).')
        return
    correr('[R 1/3] Paquetes de R', [rscript, os.path.join(SCRIPTS, 'R', 'instalar_paquetes_R.R')])
    correr('[R 2/3] Rodilla: modelo mixto, Kenward-Roger, contrastes por día',
           [rscript, os.path.join(SCRIPTS, 'R', 'analisis_rodilla.R'), OUT])
    correr('[R 3/3] Las 18 variables: LRT en R',
           [rscript, os.path.join(SCRIPTS, 'R', 'todas_las_variables_LRT.R'), OUT],
           os.path.join(OUT, 'R', 'comparacion_18_variables_R_vs_Python.txt'))
    print('Terminado. Resultados en', os.path.relpath(OUT, RAIZ))


if __name__ == '__main__':
    main()
