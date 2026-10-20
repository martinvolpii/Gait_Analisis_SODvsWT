# -*- coding: utf-8 -*-
"""Estadistica por variable (cada una por separado): resumen por animal y dia, pruebas por dia,
modelo mixto genotipo x dia con intercepto aleatorio por animal, correccion de multiplicidad."""
import numpy as np, pandas as pd, os, warnings
from scipy import stats
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
import rutas
warnings.filterwarnings('ignore')

S = rutas.INTER
DIAS = ['P30', 'P37', 'P44', 'P51', 'P65', 'P85']

# (codigo, grupo, nombre legible, columna en ciclos, unidad, tipo de resumen top10: 'max'|'min'|None)
VARS = [
    ('V1', '1 Excursion de tobillo', 'Excursión de tobillo', 'tobillo_exc', '°', 'max'),
    ('V2', '2 Excursion de pie (MTP)', 'Excursión de pie (MTP)', 'pie_exc', '°', 'max'),
    ('V3a', '3 Excursion de rodilla', 'Excursión de rodilla (DLC)', 'rodilla_exc', '°', 'max'),
    ('V3b', '3 Excursion de rodilla', 'Excursión de rodilla (triangulada)', 'rodilla_tri_exc', '°', 'max'),
    ('V4a', '4 Goniometria', 'Rodilla: extensión máxima', 'rodilla_max', '°', 'max'),
    ('V4b', '4 Goniometria', 'Rodilla: flexión máxima (ángulo mínimo)', 'rodilla_min', '°', 'min'),
    ('V4c', '4 Goniometria', 'Rodilla: ángulo medio', 'rodilla_med', '°', None),
    ('V4d', '4 Goniometria', 'Tobillo: extensión máxima', 'tobillo_max', '°', 'max'),
    ('V4e', '4 Goniometria', 'Tobillo: flexión máxima (ángulo mínimo)', 'tobillo_min', '°', 'min'),
    ('V4f', '4 Goniometria', 'Tobillo: ángulo medio', 'tobillo_med', '°', None),
    ('V4g', '4 Goniometria', 'Pie: extensión máxima', 'pie_max', '°', 'max'),
    ('V4h', '4 Goniometria', 'Pie: flexión máxima (ángulo mínimo)', 'pie_min', '°', 'min'),
    ('V4i', '4 Goniometria', 'Pie: ángulo medio', 'pie_med', '°', None),
    ('V5a', '5 Frecuencia, apoyo y oscilacion', 'Frecuencia de paso', 'frecuencia_hz', 'Hz', None),
    ('V5b', '5 Frecuencia, apoyo y oscilacion', 'Tiempo de apoyo', 'apoyo_s', 's', None),
    ('V5c', '5 Frecuencia, apoyo y oscilacion', 'Tiempo de oscilación', 'oscil_s', 's', None),
    ('V6', '6 Elevacion de la pata', 'Elevación del dedo', 'elev_mm', 'mm', 'max'),
]
EXTRA = [('X1', '5 Frecuencia, apoyo y oscilacion', 'Duty factor (apoyo/ciclo)', 'duty', '', None)]


def resumen_por_animal(cyc, ses, px_cm_global, top10=False):
    c = cyc[cyc.valido].copy()
    c['elev_mm'] = c.elev_px / px_cm_global * 10.0
    filas = []
    for (sesion, dia, animal, geno), g in c.groupby(['sesion', 'dia', 'animal', 'genotipo']):
        r = dict(sesion=sesion, dia=dia, animal=str(animal), genotipo=geno, n_ciclos=len(g))
        for code, grp, name, col, unit, tk in VARS + EXTRA:
            if col == 'frecuencia_hz':
                d = g.dur_s.dropna(); r[col] = 1.0 / d.mean() if len(d) else np.nan; r['n_' + col] = len(d)
                continue
            x = g[col].dropna()
            r['n_' + col] = len(x)
            if top10 and tk is not None:
                x = x.sort_values(ascending=(tk == 'min')).iloc[:10]
            r[col] = x.mean() if len(x) else np.nan
        filas.append(r)
    A = pd.DataFrame(filas)
    A['dia'] = pd.Categorical(A.dia, DIAS, ordered=True)
    return A.sort_values(['dia', 'genotipo', 'animal']).reset_index(drop=True)


def hedges_g(a, b):
    na, nb = len(a), len(b)
    if na < 2 or nb < 2: return np.nan
    sp = np.sqrt(((na - 1) * np.var(a, ddof=1) + (nb - 1) * np.var(b, ddof=1)) / (na + nb - 2))
    if sp == 0: return np.nan
    d = (np.mean(b) - np.mean(a)) / sp
    return d * (1 - 3 / (4 * (na + nb) - 9))


def por_dia(A, col):
    rows = []
    for dia in DIAS:
        w = A[(A.dia == dia) & (A.genotipo == 'WT')][col].dropna().values
        s = A[(A.dia == dia) & (A.genotipo == 'SOD1')][col].dropna().values
        r = dict(dia=dia, n_WT=len(w), n_SOD1=len(s), media_WT=np.mean(w) if len(w) else np.nan,
                 de_WT=np.std(w, ddof=1) if len(w) > 1 else np.nan, media_SOD1=np.mean(s) if len(s) else np.nan,
                 de_SOD1=np.std(s, ddof=1) if len(s) > 1 else np.nan)
        r['dif_SOD1_menos_WT'] = r['media_SOD1'] - r['media_WT']
        r['dif_%'] = 100 * r['dif_SOD1_menos_WT'] / r['media_WT'] if r['media_WT'] else np.nan
        if len(w) >= 2 and len(s) >= 2:
            r['p_welch'] = stats.ttest_ind(s, w, equal_var=False).pvalue
            r['p_mann_whitney'] = stats.mannwhitneyu(s, w, alternative='two-sided', method='exact').pvalue
            r['hedges_g'] = hedges_g(w, s)
        else:
            r['p_welch'] = r['p_mann_whitney'] = r['hedges_g'] = np.nan
        rows.append(r)
    D = pd.DataFrame(rows)
    ok = D.p_welch.notna()
    D['p_welch_holm_6dias'] = np.nan
    if ok.any():
        D.loc[ok, 'p_welch_holm_6dias'] = multipletests(D.loc[ok, 'p_welch'], method='holm')[1]
    return D


def lrt(m_full, m_red, df):
    stat = 2 * (m_full.llf - m_red.llf)
    return max(stat, 0.0), stats.chi2.sf(max(stat, 0.0), df)


def ajustar(formula, d, reml):
    # Prueba los 4 optimizadores y se queda con el de mayor verosimilitud (v2.1, 19-10-2026).
    # Antes se quedaba con el primero que convergia; en V4c y V4d ese no era el optimo
    # (comprobado contra lme4 en R) y el p del efecto principal de genotipo salia distinto.
    best, last = None, None
    for meth in ['bfgs', 'powell', 'nm', 'lbfgs']:
        try:
            r = smf.mixedlm(formula, d, groups=d['animal']).fit(reml=reml, method=meth)
            if np.isfinite(r.llf) and (best is None or r.llf > best.llf + 1e-9):
                best = r
        except Exception as e:
            last = e
    if best is None:
        raise RuntimeError(f'no converge: {last}')
    return best


def modelo_mixto(A, col, dias=DIAS):
    d = A[A.dia.isin(dias)][['animal', 'genotipo', 'dia', col]].dropna().copy()
    d['dia'] = d.dia.astype(str); d = d.rename(columns={col: 'y'})
    d['genotipo'] = pd.Categorical(d.genotipo, ['WT', 'SOD1'])
    nd = d.dia.nunique()
    out = dict(n_obs=len(d), n_animales=d.animal.nunique())
    try:
        full = ajustar('y ~ C(genotipo) * C(dia)', d, False)
        add = ajustar('y ~ C(genotipo) + C(dia)', d, False)
        nul = ajustar('y ~ C(dia)', d, False)
        out['chi2_interaccion'], out['p_interaccion_GxD'] = lrt(full, add, nd - 1)
        out['chi2_genotipo_principal'], out['p_genotipo_principal'] = lrt(add, nul, 1)
        out['chi2_genotipo_global'], out['p_genotipo_global'] = lrt(full, nul, nd)
        # efecto principal estimado (modelo aditivo, REML)
        addr = ajustar('y ~ C(genotipo) + C(dia)', d, True)
        k = [p for p in addr.params.index if 'genotipo' in p][0]
        out['efecto_SOD1_vs_WT'] = addr.params[k]
        ci = addr.conf_int().loc[k]; out['IC95_inf'], out['IC95_sup'] = ci[0], ci[1]
        vc = float(addr.cov_re.iloc[0, 0]); out['var_animal'] = vc; out['var_residual'] = addr.scale
        out['ICC_animal'] = vc / (vc + addr.scale) if (vc + addr.scale) > 0 else np.nan
    except Exception as e:
        out['error'] = str(e)[:120]
    return out


def main():
    ses = pd.read_csv(os.path.join(S, 'sesiones_crudo.csv'))
    px_cm = float(np.nanmedian(ses.px_por_cm_cinta))
    cyc = pd.read_csv(os.path.join(S, 'ciclos_crudo.csv'))
    cyc_sg = pd.read_csv(os.path.join(S, 'ciclos_sg5.csv'))
    A = resumen_por_animal(cyc, ses, px_cm)
    A10 = resumen_por_animal(cyc, ses, px_cm, top10=True)
    Asg = resumen_por_animal(cyc_sg, ses, px_cm)
    # elevacion con calibracion por sesion (sensibilidad)
    cal = ses.set_index('sesion').px_por_cm_cinta
    Acal = A.copy(); Acal['elev_mm'] = A.elev_mm * px_cm / A.sesion.map(cal).values
    res, pdias, sens = [], {}, []
    for code, grp, name, col, unit, tk in VARS + EXTRA:
        D = por_dia(A, col); pdias[code] = D
        m = modelo_mixto(A, col); m_sin85 = modelo_mixto(A, col, DIAS[:-1])
        r = dict(codigo=code, grupo=grp, variable=name, columna=col, unidad=unit,
                 media_WT=A[A.genotipo == 'WT'][col].mean(), media_SOD1=A[A.genotipo == 'SOD1'][col].mean(), **m)
        r['p_interaccion_sin_P85'] = m_sin85.get('p_interaccion_GxD'); r['p_genotipo_principal_sin_P85'] = m_sin85.get('p_genotipo_principal')
        j = D.p_welch.idxmin() if D.p_welch.notna().any() else None
        r['dia_p_min'] = D.loc[j, 'dia'] if j is not None else None
        r['p_min_por_dia'] = D.loc[j, 'p_welch'] if j is not None else np.nan
        r['p_min_por_dia_holm'] = D.loc[j, 'p_welch_holm_6dias'] if j is not None else np.nan
        r['dias_p<0.05'] = ','.join(D[D.p_welch < 0.05].dia.astype(str))
        res.append(r)
        # sensibilidades
        s = dict(codigo=code, variable=name)
        for lab, AA in [('top10', A10), ('suavizado_SG5', Asg)] + ([('elev_calib_por_sesion', Acal)] if col == 'elev_mm' else []):
            if lab == 'top10' and tk is None: continue
            mm = modelo_mixto(AA, col); DD = por_dia(AA, col)
            s[f'{lab}_p_interaccion'] = mm.get('p_interaccion_GxD'); s[f'{lab}_p_genotipo'] = mm.get('p_genotipo_principal')
            s[f'{lab}_efecto'] = mm.get('efecto_SOD1_vs_WT')
            s[f'{lab}_dias_p<0.05'] = ','.join(DD[DD.p_welch < 0.05].dia.astype(str))
        sens.append(s)
    R = pd.DataFrame(res)
    prim = R.codigo.isin([v[0] for v in VARS])
    for c in ['p_interaccion_GxD', 'p_genotipo_principal', 'p_genotipo_global']:
        R[c.replace('p_', 'q_FDR_')] = np.nan
        R.loc[prim, c.replace('p_', 'q_FDR_')] = multipletests(R.loc[prim, c].fillna(1), method='fdr_bh')[1]
    # FDR sobre todas las pruebas por dia (17 variables x 6 dias)
    allp = pd.concat([pdias[c].assign(codigo=c) for c in R[prim].codigo])
    allp['q_FDR_102'] = multipletests(allp.p_welch.fillna(1), method='fdr_bh')[1]
    for c in R[prim].codigo:
        pdias[c] = pdias[c].merge(allp[allp.codigo == c][['dia', 'q_FDR_102']], on='dia', how='left')
    R['p_mejor'] = R[['p_interaccion_GxD', 'p_genotipo_principal']].min(axis=1)
    pd.to_pickle(dict(A=A, A10=A10, Asg=Asg, R=R, pdias=pdias, sens=pd.DataFrame(sens), px_cm=px_cm, ses=ses),
                 os.path.join(S, 'estadistica.pkl'))
    pd.set_option('display.width', 250)
    print('px/cm global', round(px_cm, 2))
    print(R[['codigo', 'variable', 'media_WT', 'media_SOD1', 'p_interaccion_GxD', 'q_FDR_interaccion_GxD', 'p_genotipo_principal',
             'q_FDR_genotipo_principal', 'efecto_SOD1_vs_WT', 'dias_p<0.05', 'p_min_por_dia', 'p_min_por_dia_holm',
             'p_interaccion_sin_P85', 'ICC_animal']].round(4).to_string())
    print(pd.DataFrame(sens).round(4).to_string())


if __name__ == '__main__':
    main()
