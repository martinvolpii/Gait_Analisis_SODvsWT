# -*- coding: utf-8 -*-
"""
marcha_v2.py — De las coordenadas de DeepLabCut a las variables de marcha (miembro posterior, vista lateral, cinta)
Versión 2.1. Estudio longitudinal WT vs SOD1, P30–P85. Autor: Martín Volpi Canales.

Este encabezado describe qué hace el script y por qué, en el mismo orden en que lo hace el código.
Cada paso está marcado en el código con un comentario "# PASO n".

ENTRADA Y SALIDA
  Entrada: un CSV de DeepLabCut por sesión (datos/csv/<día>/<animal>_<día>DLC_....csv). Cada fila es un
           cuadro del video; para cada uno de 6 marcadores hay x, y (píxeles) y likelihood (0–1):
             crest = cresta ilíaca · hip = cadera · knee = rodilla · ankle = tobillo
             foot  = articulación metatarsofalángica (MTP) · toe = punta del dedo
           Además, datos/INDICE_videos.csv da los cuadros por segundo reales de cada sesión.
  Salida (resultados/intermedios/):
             sesiones_<v>.csv  una fila por sesión (período, umbrales, calibración, nº de pasos y ciclos)
             ciclos_<v>.csv    una fila por ciclo, con todas las variables
             eventos_<v>.csv   una fila por paso candidato, con sus 3 señales y si se confirmó
           <v> = crudo (versión principal) | sg5 (suavizado, sensibilidad) | *_60fps (validación)
  Uso:  python marcha_v2.py --crudo        (versión principal; ejecutar_analisis.py corre las 4 versiones)
        python marcha_v2.py                (con suavizado Savitzky-Golay de 5 cuadros)
        opción --sin-decimar               (no reduce a 30 fps las sesiones grabadas a 60 fps)

PASO 1 — LEER LOS PUNTOS (función cargar)
  Se leen x, y y likelihood de los 6 marcadores, cuadro por cuadro.

PASO 2 — LIMPIEZA (comienzo de procesar_sesion)
  - Un punto con likelihood < 0,70 (PCUT) es dato faltante en ese cuadro; nunca se inventa.
  - Exclusiones fijas revisadas a mano: 866_P85 cuadros 0 y 2–11 (sin animal); 862_P85 cuadros 0–1 y
    foot/toe del cuadro 21 (error de tracking de 414 px).
  - Videos de 60 fps (856_P30, 859_P30): se toma un cuadro de cada dos para quedar en 30 fps, como el resto.
  - Huecos de 1 o 2 cuadros (MAXGAP) se rellenan con una recta; huecos más largos quedan vacíos (interp_cortos).
  - Sin suavizado en la versión principal (--crudo). A 30 fps una zancada dura ~7 cuadros y un filtro la
    aplasta: con los videos de 60 fps se midió que Savitzky-Golay 5 reduce las excursiones 33–45 %.
    El suavizado (suavizar) solo se usa como análisis de sensibilidad.

PASO 3 — PLANO CARTESIANO (procesar_sesion)
  En la imagen, el origen está arriba a la izquierda y la y crece hacia ABAJO. Se pasa a un plano común:
  - y -> -y: la altura crece hacia arriba.
  - Si la cresta queda a la izquierda de la cadera (mediana de la sesión), el animal mira a la izquierda:
    x -> -x, para que el animal mire siempre hacia +x. No cambia ningún ángulo ni distancia, solo el
    sentido de giro, que importa para el signo del ángulo.
  - P51 se grabó en 1280x720: sus coordenadas se multiplican por 1,5 para quedar en la escala 1080p.
  Las coordenadas absolutas no importan: solo se usan diferencias entre puntos y cambios en el tiempo.

PASO 4 — ÁNGULOS EN CADA CUADRO (angulo_signado y el diccionario ang)
  En el vértice V (la articulación) se arman dos vectores hacia los marcadores vecinos P y D:
      u = P - V,  w = D - V
      theta = atan2(u_x*w_y - u_y*w_x,  u_x*w_x + u_y*w_y)      en grados, llevado a [0, 360)
  atan2 usa el producto cruz (sentido de giro) y el producto punto (alineación), así que, a diferencia
  del arcocoseno (0–180°), el ángulo puede pasar de 180° sin "rebotar" (importa en el pie).
      rodilla  = theta(cadera, RODILLA, tobillo)
      tobillo  = 360 - theta(rodilla, TOBILLO, pie)
      pie      = 360 - theta(tobillo, PIE, dedo)
      cadera   = 360 - theta(cresta, CADERA, rodilla)
  El 360 - theta se usa porque tobillo, pie y cadera se doblan hacia el lado contrario que la rodilla.
  Con eso, en las cuatro articulaciones un ángulo menor significa más flexión.
  Ejemplo de la convención: con V = (0, 0), P = (0, 1) y D = (1, 0), u = (0, 1) y w = (1, 0); el producto
  cruz vale -1 y el punto 0, así que atan2 da -90° y, llevado a [0, 360), theta = 270°. Si se intercambian
  P y D, theta = 90°.
  La altura del dedo de cada cuadro es su y en el plano (h).

PASO 5 — RODILLA TRIANGULADA (triangular_rodilla)
  El marcador de rodilla se desliza sobre la piel. Se reconstruye la rodilla como el cruce de dos círculos:
  radio Lf (mediana cadera-rodilla) centrado en la cadera y radio Lt (mediana rodilla-tobillo) centrado en
  el tobillo; de los dos cruces se toma el del lado de la rodilla marcada. Si la pata "no alcanza"
  (cadera-tobillo > Lf + Lt) se pone sobre la línea cadera-tobillo. Solo se usa para V3b.

PASO 6 — CALIBRACIÓN CON LA CINTA (bloque px_cm)
  La banda corre a 20 cm/s y en el apoyo el pie viaja con ella. En los cuadros de apoyo (dedo yendo hacia
  atrás respecto de la cadera y en el 40 % más bajo de su altura) se mide la velocidad horizontal del dedo
  en px/s; su mediana / 20 = píxeles por cm de la sesión. La elevación en mm usa la mediana de las 71
  sesiones (se calcula en estadistica_v2.py). La ubicación de la banda en la imagen no hace falta aquí:
  todas las variables son diferencias dentro de un mismo ciclo.

PASO 7 — RITMO DE LA ZANCADA, PERÍODO T (periodo_autocorr)
  Autocorrelación de la altura del dedo: se compara la señal consigo misma desplazada L cuadros; T es el
  primer máximo entre 0,12 y 0,60 s (correlación > 0,1), refinado con una parábola. T es la referencia de
  cada sesión para todo lo que sigue.

PASO 8 — DETECTAR PASOS: 3 SEÑALES QUE DEBEN COINCIDIR (bloque "deteccion multisenal")
  Un paso cuenta solo si, en la oscilación, el dedo sube y la rodilla y la cadera se flexionan (su ángulo
  baja). En los datos, el mínimo de cadera cae en el cuadro del pico del dedo y el de rodilla un cuadro antes.
  a) Candidatos: picos de la altura del dedo separados al menos 0,7*T (find_peaks).
  b) En cada candidato p:
       subida del dedo     = h[p] - mínimo de h en los 0,8*T cuadros anteriores               (px)
       flexión de rodilla  = máximo del ángulo en los 0,8*T anteriores
                             - mínimo del ángulo entre p-0,3*T y p+0,15*T                     (grados)
       flexión de cadera   = lo mismo con el ángulo de cadera                                 (grados)
  c) Umbral de cada señal = 30 % de su mediana en la sesión (pisos: rodilla 5°, cadera 3°).
     Confirmado = las tres señales sobre su umbral. Señal sin tracking = no confirmado.

PASO 9 — CICLOS (bucle de ciclos)
  Ciclo = de un paso confirmado al siguiente (de un máximo del dedo al siguiente). Su duración se mide
  entre los picos refinados con una parábola (pico_parabolico). Es válido si:
     dura entre 0,6*T y 1,6*T (si dura más, falta un paso intermedio), y
     la elevación del dedo supera 0,3 x la mediana de la sesión (si no, no hubo paso real).

PASO 10 — VARIABLES DE CADA CICLO (bucle de ciclos)
  Para cada ángulo, sobre los cuadros del ciclo (ambos picos incluidos; si falta un dato, ese ángulo no se
  calcula en ese ciclo): máximo (extensión máxima), mínimo (flexión máxima), media (ángulo medio) y
  excursión (máximo - mínimo).
     V1  excursión de tobillo        V2  excursión de pie (MTP)
     V3a excursión de rodilla        V3b excursión de rodilla triangulada
     V4a/b/c rodilla máx/mín/media   V4d/e/f tobillo máx/mín/media   V4g/h/i pie máx/mín/media
     V6  elevación del dedo = máx - mín de h en el ciclo (px; en mm al dividir por px/cm x 10)
  Apoyo y oscilación (V5b, V5c): r = x_dedo - x_cadera, v = dr/dt. Apoyo = v < 0 (la banda lleva el pie
  hacia atrás respecto del cuerpo); oscilación = v > 0. v se interpola en 400 puntos entre los picos:
     duty = fracción con v < 0,  apoyo = duty x duración,  oscilación = (1 - duty) x duración.
  Frecuencia (V5a) = 1 / duración media de los ciclos válidos (se calcula en estadistica_v2.py).

DESPUÉS (estadistica_v2.py)
  Unidad experimental = el animal: se promedian los ciclos válidos de cada animal y día. Luego, Welch +
  Mann-Whitney + Hedges g por día (Holm sobre 6 días) y modelo mixto y ~ genotipo x día + (1 | animal).
"""
import numpy as np, pandas as pd, glob, os, re, json, sys
import rutas
from scipy.signal import savgol_filter, find_peaks

BP = ['crest', 'hip', 'knee', 'ankle', 'foot', 'toe']
PCUT = 0.70
MAXGAP = 2
GENO = {**{a: 'WT' for a in ['821', '860', '861', '863', '866', '867']},
        **{a: 'SOD1' for a in ['856', '857', '859', '862', '864', '865']}}
DIAS = ['P30', 'P37', 'P44', 'P51', 'P65', 'P85']


def leer_fps():
    if not os.path.exists(rutas.INDICE):
        sys.exit(f'Falta {rutas.INDICE} (cuadros por segundo de cada sesion).')
    d = pd.read_csv(rutas.INDICE, sep=';', encoding='utf-8-sig')
    d = d[d.modelo == 'snapshot_1000']
    out = {}
    for _, r in d.iterrows():
        f = r.fps
        if isinstance(f, str) and '/' in f:
            a, b = f.split('/'); f = float(a) / float(b)
        out[r.sesion] = float(f)
    return out


def cargar(path):
    # PASO 1: leer x, y y likelihood de los 6 marcadores (CSV de DLC con 3 filas de encabezado)
    d = pd.read_csv(path, header=[1, 2], index_col=0)
    X = {}
    for b in BP:
        X[b] = np.c_[d[(b, 'x')].values, d[(b, 'y')].values].astype(float)
        X[b + '_l'] = d[(b, 'likelihood')].values.astype(float)
    return X, len(d)


def interp_cortos(v, maxgap=MAXGAP):
    """interpola linealmente solo corridas internas de NaN de largo <= maxgap"""
    v = v.copy(); n = len(v); isn = np.isnan(v)
    if isn.all() or not isn.any():
        return v
    i = 0
    while i < n:
        if isn[i]:
            j = i
            while j < n and isn[j]:
                j += 1
            if i > 0 and j < n and (j - i) <= maxgap:
                v[i:j] = np.interp(np.arange(i, j), [i - 1, j], [v[i - 1], v[j]])
            i = j
        else:
            i += 1
    return v


def suavizar(v, w=5, o=2):
    """SG sobre cada tramo continuo sin NaN de largo >= w"""
    out = v.copy(); ok = ~np.isnan(v); n = len(v); i = 0
    while i < n:
        if ok[i]:
            j = i
            while j < n and ok[j]:
                j += 1
            if j - i >= w:
                out[i:j] = savgol_filter(v[i:j], w, o)
            i = j
        else:
            i += 1
    return out


def angulo_signado(P, V, D):
    """angulo en V entre V->P y V->D, en (0,360), con signo segun la orientacion (x ya espejado)"""
    u = P - V; w = D - V
    cr = u[:, 0] * w[:, 1] - u[:, 1] * w[:, 0]
    dt = (u * w).sum(1)
    a = np.degrees(np.arctan2(cr, dt))
    return np.mod(a, 360.0)


def periodo_autocorr(h, fps):
    """periodo dominante (cuadros) de la senal h con NaN, autocorrelacion normalizada por pares validos"""
    x = h - np.nanmean(h); n = len(x)
    lo, hi = int(np.floor(0.12 * fps)), int(np.ceil(0.60 * fps))
    ac = []
    for L in range(0, min(hi + 2, n - 5)):
        a, b = x[:n - L], x[L:]
        m = ~np.isnan(a) & ~np.isnan(b)
        ac.append(np.sum(a[m] * b[m]) / max(m.sum(), 1) if m.sum() > 10 else np.nan)
    ac = np.array(ac)
    if len(ac) <= lo + 1 or np.isnan(ac[0]) or ac[0] == 0:
        return np.nan
    ac = ac / ac[0]
    cand = [L for L in range(max(lo, 1), len(ac) - 1)
            if ac[L] >= ac[L - 1] and ac[L] >= ac[L + 1] and ac[L] > 0.1]
    if not cand:
        return np.nan
    L = cand[0]
    # refinamiento parabolico
    y0, y1, y2 = ac[L - 1], ac[L], ac[L + 1]
    den = y0 - 2 * y1 + y2
    return L + (0.5 * (y0 - y2) / den if den != 0 else 0.0)


def pico_parabolico(h, p):
    if 0 < p < len(h) - 1 and not np.isnan(h[p - 1]) and not np.isnan(h[p + 1]):
        y0, y1, y2 = h[p - 1], h[p], h[p + 1]; den = y0 - 2 * y1 + y2
        if den != 0:
            return p + 0.5 * (y0 - y2) / den
    return float(p)


def triangular_rodilla(H, K, A, Lf, Lt):
    d = np.linalg.norm(A - H, axis=1)
    u = (A - H) / d[:, None]
    perp = np.c_[-u[:, 1], u[:, 0]]
    a = (Lf ** 2 - Lt ** 2 + d ** 2) / (2 * d)
    h2 = Lf ** 2 - a ** 2
    estirada = h2 < 0
    h = np.sqrt(np.clip(h2, 0, None))
    P = H + a[:, None] * u
    K1 = P + h[:, None] * perp; K2 = P - h[:, None] * perp
    # lado de la rodilla trackeada; si falta, el lado mayoritario de la sesion
    s = np.sign((A[:, 0] - H[:, 0]) * (K[:, 1] - H[:, 1]) - (A[:, 1] - H[:, 1]) * (K[:, 0] - H[:, 0]))
    smaj = np.sign(np.nansum(s)) or 1.0
    s = np.where(np.isnan(s) | (s == 0), smaj, s)
    s1 = np.sign((A[:, 0] - H[:, 0]) * (K1[:, 1] - H[:, 1]) - (A[:, 1] - H[:, 1]) * (K1[:, 0] - H[:, 0]))
    Kt = np.where((s1 == s)[:, None], K1, K2)
    return Kt, estirada


def procesar_sesion(path, dia, animal, fps, suave=True, decimar60=True):
    X, n0 = cargar(path)
    ses = f'{animal}_{dia}'
    escala = 1.5 if dia == 'P51' else 1.0   # PASO 3: P51 esta en 720p -> escala 1080p
    # PASO 2: limpieza. Exclusiones fijas revisadas a mano
    nanrows = []
    if ses == '866_P85': nanrows = [0] + list(range(2, 12))
    if ses == '862_P85': nanrows = [0, 1]
    for b in BP:
        X[b] = X[b] * escala
        X[b][X[b + '_l'] < PCUT] = np.nan   # PASO 2: likelihood < 0,70 = dato faltante
        X[b][nanrows] = np.nan
    if ses == '862_P85':
        X['foot'][21] = np.nan; X['toe'][21] = np.nan
    if decimar60 and fps > 50:   # PASO 2: videos de 60 fps -> 30 fps (un cuadro de cada dos)
        for b in BP:
            X[b] = X[b][::2]
        fps = fps / 2
    n = len(X['toe'])
    # PASO 3: plano cartesiano. x espejado para que el animal mire a +x; y hacia arriba
    dirx = np.sign(np.nanmedian(X['crest'][:, 0] - X['hip'][:, 0]))
    if dirx == 0 or np.isnan(dirx):
        dirx = np.sign(np.nanmedian(X['knee'][:, 0] - X['hip'][:, 0])) or 1.0
    C = {}
    for b in BP:
        xy = X[b].copy(); xy[:, 0] = xy[:, 0] * dirx
        xy[:, 1] = -xy[:, 1]  # y hacia arriba
        for k in range(2):
            v = interp_cortos(xy[:, k])   # PASO 2: rellena huecos de <= 2 cuadros
            xy[:, k] = suavizar(v) if suave else v
        C[b] = xy
    # PASO 5: largos de segmento (mediana) y rodilla triangulada (cruce de dos circulos)
    Lf = np.nanmedian(np.linalg.norm(C['knee'] - C['hip'], axis=1))
    Lt = np.nanmedian(np.linalg.norm(C['ankle'] - C['knee'], axis=1))
    Ktri, estirada = triangular_rodilla(C['hip'], C['knee'], C['ankle'], Lf, Lt)
    # PASO 4: angulos en cada cuadro. Convencion: angulo interno (lado de la flexion); tobillo, pie y cadera se miden en el sentido
    # opuesto a la rodilla (con el animal mirando a +x e y hacia arriba) -> 360 - theta
    ang = {
        'rodilla': angulo_signado(C['hip'], C['knee'], C['ankle']),
        'rodilla_tri': angulo_signado(C['hip'], Ktri, C['ankle']),
        'tobillo': 360.0 - angulo_signado(C['knee'], C['ankle'], C['foot']),
        'pie': 360.0 - angulo_signado(C['ankle'], C['foot'], C['toe']),
        'cadera': 360.0 - angulo_signado(C['crest'], C['hip'], C['knee']),
    }
    h = C['toe'][:, 1].copy()              # altura del dedo (y en el plano)
    T = periodo_autocorr(h, fps)           # PASO 7: periodo de la zancada
    rel = C['toe'][:, 0] - C['hip'][:, 0]
    vrel = np.gradient(rel)  # px/cuadro (NaN se propaga)
    # PASO 6: calibracion con la cinta: en apoyo el dedo viaja con la cinta (20 cm/s)
    # (vrel < 0 = el dedo retrocede respecto de la cadera = apoyo)
    vabs = np.gradient(C['toe'][:, 0]) * fps  # px/s
    if np.isfinite(h).sum() > 10:
        hlo = np.nanpercentile(h, 40)
        m = (vrel < 0) & (h <= hlo) & np.isfinite(vabs)
        px_cm = np.nanmedian(np.abs(vabs[m])) / 20.0 if m.sum() >= 5 else np.nan
    else:
        px_cm = np.nan
    info = dict(sesion=ses, dia=dia, animal=animal, genotipo=GENO[animal], fps=fps, cuadros=n,
                dirx=dirx, T_cuadros=T, T_s=T / fps if T == T else np.nan, Lf_px=Lf, Lt_px=Lt,
                px_por_cm_cinta=px_cm,
                pct_rodilla_estirada_tri=100 * np.nanmean(estirada[~np.isnan(Ktri[:, 0])]) if (~np.isnan(Ktri[:, 0])).any() else np.nan)
    if not np.isfinite(T):
        return info, pd.DataFrame(), dict(ang=ang, h=h, peaks=[], C=C, fps=fps, eventos=pd.DataFrame())
    # PASO 8a: candidatos = picos de la altura del dedo separados >= 0,7 T
    hf = h.copy()
    pk0, _ = find_peaks(np.nan_to_num(hf, nan=np.nanmin(hf)), distance=max(1, int(np.floor(0.7 * T))))
    pk0 = [p for p in pk0 if not np.isnan(hf[p])]
    # ---- PASO 8b: deteccion multisenal: elevacion del dedo + descenso (flexion) de cadera y rodilla ----
    # En los datos el minimo de cadera cae en el mismo cuadro del pico del dedo y el de rodilla 1 cuadro antes.
    pre = int(np.ceil(0.8 * T)); wl = int(np.ceil(0.3 * T)); wr = int(np.ceil(0.15 * T))
    ev = []
    for p in pk0:
        r = dict(cuadro_pico=p)
        a0 = max(0, p - pre)
        for nom, sig in [('dedo', h), ('rodilla', ang['rodilla']), ('cadera', ang['cadera'])]:
            antes = sig[a0:p]; cerca = sig[max(0, p - wl):p + wr + 1]
            if np.isfinite(antes).sum() < 2 or np.isfinite(cerca).sum() < 1:
                r['cambio_' + nom] = np.nan; continue
            if nom == 'dedo':
                r['cambio_' + nom] = sig[p] - np.nanmin(antes)          # subida del dedo (px)
            else:
                r['cambio_' + nom] = np.nanmax(antes) - np.nanmin(cerca)  # descenso del angulo (grados)
        ev.append(r)
    E = pd.DataFrame(ev)
    umb = {}   # PASO 8c: umbral de cada senal = 30 % de su mediana en la sesion (con pisos)
    for nom, piso in [('dedo', 0.0), ('rodilla', 5.0), ('cadera', 3.0)]:
        col = 'cambio_' + nom
        med = np.nanmedian(E[col]) if len(E) and E[col].notna().any() else np.nan
        umb[nom] = max(piso, 0.3 * med) if np.isfinite(med) else np.nan
        E['ok_' + nom] = E[col] >= umb[nom]
        E['eval_' + nom] = E[col].notna()
    E['confirmado'] = E.ok_dedo & E.ok_rodilla & E.ok_cadera if len(E) else []
    pk = [int(p) for p in E.loc[E.confirmado, 'cuadro_pico']] if len(E) else []   # pasos confirmados (las 3 senales)
    tp = [pico_parabolico(hf, p) for p in pk]
    filas = []   # PASO 9: ciclos = de un paso confirmado al siguiente
    elev_all = []
    for i in range(len(pk) - 1):
        a, b = pk[i], pk[i + 1]
        seg = slice(a, b + 1)
        hh = h[seg]
        elev_all.append(np.nanmax(hh) - np.nanmin(hh) if np.isfinite(hh).sum() >= 2 else np.nan)
    med_elev = np.nanmedian(elev_all) if len(elev_all) else np.nan
    for i in range(len(pk) - 1):
        a, b = pk[i], pk[i + 1]
        dur_fr = tp[i + 1] - tp[i]
        r = dict(sesion=ses, dia=dia, animal=animal, genotipo=GENO[animal], ciclo=i, ini=a, fin=b,
                 dur_s=dur_fr / fps)
        r['valido_dur'] = (0.6 * T <= dur_fr <= 1.6 * T)   # criterio 1: duracion entre 0,6 T y 1,6 T
        seg = slice(a, b + 1)  # incluye ambos picos (ciclo cerrado)
        hh = h[seg]
        r['elev_px'] = (np.max(hh) - np.min(hh)) if not np.isnan(hh).any() else np.nan
        # criterio 2: elevacion del dedo > 0,3 x mediana de la sesion
        r['valido_elev'] = bool(np.isfinite(r['elev_px']) and r['elev_px'] > 0.3 * med_elev)
        for k, v in ang.items():   # PASO 10: maximo, minimo, media y excursion de cada angulo en el ciclo
            av = v[seg]
            if np.isnan(av).any():
                r[k + '_max'] = r[k + '_min'] = r[k + '_med'] = r[k + '_exc'] = np.nan
            else:
                r[k + '_max'] = av.max(); r[k + '_min'] = av.min(); r[k + '_med'] = av.mean()
                r[k + '_exc'] = av.max() - av.min()
        # PASO 10: apoyo/oscilacion con interpolacion lineal de vrel entre cuadros (apoyo = vrel < 0)
        vv = vrel[a:b + 1]
        if np.isnan(vv).any():
            r['apoyo_s'] = r['oscil_s'] = r['duty'] = np.nan
        else:
            t = np.arange(a, b + 1, dtype=float)
            # recorte a los tiempos de pico refinados
            t0, t1 = tp[i], tp[i + 1]
            tt = np.linspace(t0, t1, 400)
            vi = np.interp(tt, t, vv)
            neg = (vi < 0).mean()
            r['apoyo_s'] = neg * dur_fr / fps
            r['oscil_s'] = (1 - neg) * dur_fr / fps
            r['duty'] = neg
        filas.append(r)
    cyc = pd.DataFrame(filas)
    if len(cyc):
        cyc['valido'] = cyc.valido_dur & cyc.valido_elev
    info.update(n_picos_dedo=len(E), n_eventos_sin_rodilla=int((E.eval_rodilla & ~E.ok_rodilla).sum()) if len(E) else 0,
                n_eventos_sin_cadera=int((E.eval_cadera & ~E.ok_cadera).sum()) if len(E) else 0,
                n_eventos_cadera_no_evaluable=int((~E.eval_cadera).sum()) if len(E) else 0,
                n_eventos_confirmados=int(E.confirmado.sum()) if len(E) else 0,
                umbral_dedo_px=umb.get('dedo'), umbral_rodilla_grados=umb.get('rodilla'), umbral_cadera_grados=umb.get('cadera'),
                n_ciclos_detectados=len(cyc), n_ciclos_validos=int(cyc['valido'].sum()) if len(cyc) else 0,
                elev_mediana_px=med_elev)
    if len(E):
        E.insert(0, 'sesion', ses)
    return info, cyc, dict(ang=ang, h=h, peaks=tp, C=C, fps=fps, Ktri=Ktri, eventos=E)


def sesiones():
    fpsmap = leer_fps()
    out = []
    n = sum(len(glob.glob(os.path.join(rutas.CSV, d, '*.csv'))) for d in DIAS)
    if n == 0:
        sys.exit(f'No hay CSV de DeepLabCut en {os.path.join(rutas.CSV, "<dia>")} (P30 ... P85).')
    for dia in DIAS:
        for f in sorted(glob.glob(os.path.join(rutas.CSV, dia, '*.csv'))):
            animal = os.path.basename(f)[:3]
            ses = f'{animal}_{dia}'
            fps = fpsmap.get(ses, 30.0)
            out.append((f, dia, animal, fps))
    return out


if __name__ == '__main__':
    suave = '--crudo' not in sys.argv
    dec = '--sin-decimar' not in sys.argv
    tag = ('sg5' if suave else 'crudo') + ('' if dec else '_60fps')
    infos, cycs, evs = [], [], []
    for f, dia, animal, fps in sesiones():
        info, cyc, aux = procesar_sesion(f, dia, animal, fps, suave=suave, decimar60=dec)
        infos.append(info); cycs.append(cyc); evs.append(aux['eventos'])
    I = pd.DataFrame(infos); Cy = pd.concat(cycs, ignore_index=True)
    I.to_csv(os.path.join(rutas.INTER, f'sesiones_{tag}.csv'), index=False)
    Cy.to_csv(os.path.join(rutas.INTER, f'ciclos_{tag}.csv'), index=False)
    pd.concat(evs, ignore_index=True).to_csv(os.path.join(rutas.INTER, f'eventos_{tag}.csv'), index=False)
    print(I[['sesion', 'T_s', 'n_picos_dedo', 'n_eventos_sin_rodilla', 'n_eventos_sin_cadera', 'n_eventos_cadera_no_evaluable', 'n_eventos_confirmados', 'n_ciclos_detectados', 'n_ciclos_validos', 'umbral_rodilla_grados', 'umbral_cadera_grados']].round(3).to_string())
