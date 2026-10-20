# -*- coding: utf-8 -*-
"""Figuras (una por variable), carpetas por variable, Excel consolidado y LEEME."""
import numpy as np, pandas as pd, os, shutil
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

import rutas
S = rutas.INTER
OUT = rutas.OUT
DIAS = ['P30', 'P37', 'P44', 'P51', 'P65', 'P85']
COL = {'WT': '#2a78d6', 'SOD1': '#eb6834'}; MK = {'WT': 'o', 'SOD1': 's'}
INK, MUTED, GRID = '#1f1f1e', '#6b6b68', '#e6e6e3'
CARPETA = {'1': '01_Excursion_tobillo', '2': '02_Excursion_pie_MTP', '3': '03_Excursion_rodilla',
           '4': '04_Goniometria_rodilla_tobillo_pie', '5': '05_Frecuencia_apoyo_oscilacion', '6': '06_Elevacion_pata'}
DEF = {
    'V1': 'máx − mín del ángulo rodilla–tobillo–pie(MTP) dentro de cada ciclo',
    'V2': 'máx − mín del ángulo tobillo–pie(MTP)–dedo dentro de cada ciclo (ángulo con signo, sin plegado en 180°)',
    'V3a': 'máx − mín del ángulo cadera–rodilla–tobillo, rodilla de DeepLabCut',
    'V3b': 'máx − mín del ángulo de rodilla con la rodilla reconstruida desde cadera y tobillo (largos fijos de fémur y tibia = mediana de la sesión)',
    'V4a': 'ángulo máximo de rodilla en el ciclo (mayor extensión)', 'V4b': 'ángulo mínimo de rodilla en el ciclo (mayor flexión)',
    'V4c': 'ángulo medio de rodilla en el ciclo', 'V4d': 'ángulo máximo de tobillo (mayor extensión)',
    'V4e': 'ángulo mínimo de tobillo (mayor flexión)', 'V4f': 'ángulo medio de tobillo',
    'V4g': 'ángulo máximo de pie/MTP', 'V4h': 'ángulo mínimo de pie/MTP', 'V4i': 'ángulo medio de pie/MTP',
    'V5a': '1 / duración media del ciclo (ciclo = entre dos pasos confirmados consecutivos: dedo arriba + flexión de cadera y rodilla)',
    'V5b': 'tiempo del ciclo con el dedo moviéndose hacia atrás respecto de la cadera (apoyo)',
    'V5c': 'tiempo del ciclo con el dedo moviéndose hacia adelante respecto de la cadera (oscilación)',
    'V6': 'máx − mín de la altura del dedo en el ciclo, en mm (px a escala 1080p / px por cm calibrado con la cinta)',
    'X1': 'apoyo / duración del ciclo (variable extra; comprimida a 30 fps)',
}


def fig_variable(A, D, code, name, col, unit, path):
    fig, ax = plt.subplots(figsize=(7.2, 4.4), dpi=150)
    fig.patch.set_facecolor('#fcfcfb'); ax.set_facecolor('#fcfcfb')
    x = np.arange(len(DIAS)); off = {'WT': -0.12, 'SOD1': 0.12}
    rng = np.random.default_rng(1)
    for g in ['WT', 'SOD1']:
        m, se = [], []
        for d in DIAS:
            v = A[(A.dia == d) & (A.genotipo == g)][col].dropna().values
            m.append(v.mean() if len(v) else np.nan); se.append(v.std(ddof=1) / np.sqrt(len(v)) if len(v) > 1 else np.nan)
            xx = x[DIAS.index(d)] + off[g] + rng.uniform(-0.04, 0.04, len(v))
            ax.scatter(xx, v, s=16, marker=MK[g], facecolor='none', edgecolor=COL[g], linewidth=0.9, alpha=0.75, zorder=2)
        m, se = np.array(m), np.array(se)
        ax.errorbar(x + off[g], m, yerr=se, color=COL[g], lw=2, marker=MK[g], ms=7, mfc=COL[g], mec='#fcfcfb',
                    mew=1.5, capsize=3, zorder=3, label=f'{g} (media ± EE)')
        ax.annotate(g, (x[-1] + off[g] + 0.12, m[-1]), color=INK, fontsize=9, va='center')
    ymax = np.nanmax(A[col]); ymin = np.nanmin(A[col]); span = ymax - ymin
    for i, r in D.iterrows():
        if r.p_welch < 0.05:
            lab = f'p={r.p_welch:.3f}' + ('\nHolm<0,05' if r.p_welch_holm_6dias < 0.05 else '')
            ax.text(x[DIAS.index(r.dia)], ymax + 0.04 * span, lab, ha='center', va='bottom', fontsize=8, color=INK)
    ax.set_ylim(ymin - 0.08 * span, ymax + 0.22 * span)
    ax.set_xticks(x); ax.set_xticklabels(DIAS); ax.set_xlim(-0.5, len(DIAS) - 0.2)
    ax.set_ylabel(f'{name} ({unit})' if unit else name, color=INK)
    ax.set_title(f'{code}  {name}', loc='left', fontsize=11, color=INK)
    ax.grid(axis='y', color=GRID, lw=0.8); ax.set_axisbelow(True)
    for s in ['top', 'right']: ax.spines[s].set_visible(False)
    for s in ['left', 'bottom']: ax.spines[s].set_color('#bdbdb8')
    ax.tick_params(colors=MUTED)
    ax.legend(frameon=False, fontsize=8, loc='lower left')
    fig.text(0.01, 0.01, 'Puntos = animales (media de sus ciclos). P30–P65: snapshot 1000; P85 (y 821_P51): iteración 5. Pruebas por día: Welch.',
             fontsize=6.5, color=MUTED)
    fig.tight_layout(rect=(0, 0.03, 1, 1)); fig.savefig(path, facecolor=fig.get_facecolor()); plt.close(fig)


def main():
    o = pd.read_pickle(os.path.join(S, 'estadistica.pkl'))
    A, R, pdias, sens = o['A'], o['R'], o['pdias'], o['sens']
    # No se borra la carpeta: se sobrescriben los archivos (asi no se pierde nada que el usuario haya agregado).
    os.makedirs(OUT, exist_ok=True)
    for code, row in R.set_index('codigo').iterrows():
        sub = os.path.join(OUT, CARPETA[row.grupo[0]]); os.makedirs(sub, exist_ok=True)
        col, name, unit = row.columna, row.variable, row.unidad
        wide = A.pivot_table(index=['genotipo', 'animal'], columns='dia', values=col, observed=False).reindex(columns=DIAS)
        wide = wide.sort_index(level=0, ascending=False).round(4)
        stem = f'{code}_' + col
        wide.to_csv(os.path.join(sub, f'{stem}_por_animal.csv'), sep=';', decimal=',', encoding='utf-8-sig')
        pdias[code].round(5).to_csv(os.path.join(sub, f'{stem}_WT_vs_SOD1_por_dia.csv'), sep=';', decimal=',', index=False, encoding='utf-8-sig')
        fig_variable(A, pdias[code], code, name, col, unit, os.path.join(sub, f'{stem}.png'))
    # resumen ordenado
    keep = ['codigo', 'grupo', 'variable', 'unidad', 'media_WT', 'media_SOD1', 'efecto_SOD1_vs_WT', 'IC95_inf', 'IC95_sup',
            'p_genotipo_principal', 'q_FDR_genotipo_principal', 'p_interaccion_GxD', 'q_FDR_interaccion_GxD',
            'p_genotipo_global', 'q_FDR_genotipo_global', 'p_interaccion_sin_P85', 'p_genotipo_principal_sin_P85',
            'dias_p<0.05', 'dia_p_min', 'p_min_por_dia', 'p_min_por_dia_holm', 'ICC_animal', 'n_obs', 'n_animales']
    Rs = R[keep].sort_values('p_interaccion_GxD').copy()
    Rs.insert(3, 'definicion', Rs.codigo.map(DEF))
    Rs.round(5).to_csv(os.path.join(OUT, 'RESUMEN_todas_las_variables.csv'), sep=';', decimal=',', index=False, encoding='utf-8-sig')
    # Excel
    xl = os.path.join(OUT, 'Resultados_marcha_v2.xlsx')
    with pd.ExcelWriter(xl, engine='openpyxl') as w:
        Rs.round(5).to_excel(w, sheet_name='Resumen', index=False)
        for code, row in R.set_index('codigo').iterrows():
            col = row.columna
            wide = A.pivot_table(index=['genotipo', 'animal'], columns='dia', values=col, observed=False).reindex(columns=DIAS).sort_index(level=0, ascending=False)
            sh = f'{code}'
            wide.round(4).to_excel(w, sheet_name=sh, startrow=2)
            pdias[code].round(5).to_excel(w, sheet_name=sh, startrow=18, index=False)
            ws = w.sheets[sh]
            ws.cell(1, 1, f'{code} {row.variable} ({row.unidad}) - {DEF[code]}')
            ws.cell(17, 1, 'WT vs SOD1 por día (Welch, Mann-Whitney exacto, Hedges g; Holm sobre 6 días; q FDR sobre 17 variables x 6 días)')
            ws.cell(27, 1, f'Modelo mixto y ~ genotipo*día + (1|animal): p interacción = {row.p_interaccion_GxD:.4g} (q {row.get("q_FDR_interaccion_GxD", np.nan):.3g}); '
                           f'p genotipo (principal) = {row.p_genotipo_principal:.4g}; efecto SOD1−WT = {row.efecto_SOD1_vs_WT:.3g} {row.unidad} '
                           f'[{row.IC95_inf:.3g}, {row.IC95_sup:.3g}]; sin P85: p interacción = {row.p_interaccion_sin_P85:.4g}')
        sens.round(5).to_excel(w, sheet_name='Sensibilidad', index=False)
        pd.read_csv(os.path.join(S, 'validacion_30_vs_60fps.csv')).round(4).to_excel(w, sheet_name='Validacion_30_vs_60fps', index=False)
        o['ses'].round(4).to_excel(w, sheet_name='Sesiones_QC', index=False)
        A.round(5).to_excel(w, sheet_name='Por_animal_y_dia', index=False)
        pd.read_csv(os.path.join(S, 'ciclos_crudo.csv')).round(5).to_excel(w, sheet_name='Ciclos', index=False)
    # QC
    q = os.path.join(OUT, 'QC'); os.makedirs(q, exist_ok=True)
    for f in ['qc_ciclos.png', 'validacion_30_vs_60fps.csv', 'sesiones_crudo.csv', 'ciclos_crudo.csv', 'eventos_crudo.csv', 'ciclos_sg5.csv']:
        shutil.copy(os.path.join(S, f), q)
    print('px/cm', o['px_cm']); print(sorted(os.listdir(OUT)))


if __name__ == '__main__':
    main()
