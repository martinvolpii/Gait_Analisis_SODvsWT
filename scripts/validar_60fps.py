# -*- coding: utf-8 -*-
"""Validación de la resolución temporal y del suavizado con las dos sesiones grabadas a 60 fps.

856_P30 y 859_P30 son las únicas sesiones con video a 60 fps. marcha_v2.py las analiza normalmente
a 30 fps (toma un cuadro de cada dos) para que todas las sesiones tengan la misma resolución.
Aquí se comparan, para esas dos sesiones, cuatro versiones del mismo cálculo:
    sg5_30    30 fps con suavizado Savitzky-Golay 5      (ciclos_sg5.csv)
    crudo_30  30 fps sin suavizado  <- versión principal (ciclos_crudo.csv)
    sg5_60    60 fps con suavizado                       (ciclos_sg5_60fps.csv)
    crudo_60  60 fps sin suavizado  <- referencia        (ciclos_crudo_60fps.csv)
y se expresa cada una como % de diferencia frente a crudo_60.

Requiere haber corrido antes marcha_v2.py cuatro veces:
    python marcha_v2.py --crudo
    python marcha_v2.py
    python marcha_v2.py --crudo --sin-decimar
    python marcha_v2.py --sin-decimar
Salida: resultados/intermedios/validacion_30_vs_60fps.csv
"""
import os
import pandas as pd
import rutas

SESIONES = ['856_P30', '859_P30']
VARS = ['rodilla_exc', 'rodilla_tri_exc', 'tobillo_exc', 'pie_exc', 'rodilla_max', 'rodilla_min', 'tobillo_max',
        'tobillo_min', 'pie_max', 'pie_min', 'rodilla_med', 'tobillo_med', 'pie_med', 'dur_s', 'apoyo_s', 'oscil_s',
        'duty', 'elev_px']
VERSIONES = [('sg5_30', 'ciclos_sg5.csv'), ('crudo_30', 'ciclos_crudo.csv'),
             ('sg5_60', 'ciclos_sg5_60fps.csv'), ('crudo_60', 'ciclos_crudo_60fps.csv')]


def main():
    C = {k: pd.read_csv(os.path.join(rutas.INTER, f)) for k, f in VERSIONES}
    filas = []
    for s in SESIONES:
        for v in VARS:
            r = dict(sesion=s, var=v)
            for k, _ in VERSIONES:
                c = C[k]
                x = c[(c.sesion == s) & c.valido][v].dropna()
                r[k] = x.mean()
                r['n_' + k] = len(x)
            for k in ['sg5_30', 'crudo_30', 'sg5_60']:
                r[f'%{k}_vs_crudo60'] = 100 * (r[k] - r['crudo_60']) / r['crudo_60']
            filas.append(r)
    D = pd.DataFrame(filas)
    D.to_csv(os.path.join(rutas.INTER, 'validacion_30_vs_60fps.csv'), index=False)
    print(D[['sesion', 'var', '%sg5_30_vs_crudo60', '%crudo_30_vs_crudo60', '%sg5_60_vs_crudo60']].round(1).to_string(index=False))


if __name__ == '__main__':
    main()
