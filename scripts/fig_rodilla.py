import numpy as np, pandas as pd, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import os, rutas
o=pd.read_pickle(os.path.join(rutas.INTER,'estadistica.pkl')); A=o['A']; R=o['R'].set_index('codigo'); P=o['pdias']
DIAS=['P30','P37','P44','P51','P65','P85']; COL={'WT':'#2a78d6','SOD1':'#eb6834'}; MK={'WT':'o','SOD1':'s'}
INK,MUTED,GRID,BG='#1f1f1e','#6b6b68','#e6e6e3','#fcfcfb'
pan=[('V4a','rodilla_max','Extensión máxima'),('V4c','rodilla_med','Ángulo medio'),('V4b','rodilla_min','Flexión máxima (ángulo mínimo)')]
fig,axs=plt.subplots(1,3,figsize=(15,5.2),dpi=150,sharex=True); fig.patch.set_facecolor(BG)
x=np.arange(6); off={'WT':-0.12,'SOD1':0.12}; rng=np.random.default_rng(1)
for ax,(code,col,tit) in zip(axs,pan):
    ax.set_facecolor(BG)
    for g in ['WT','SOD1']:
        m,se=[],[]
        for i,d in enumerate(DIAS):
            v=A[(A.dia==d)&(A.genotipo==g)][col].dropna().values
            m.append(v.mean()); se.append(v.std(ddof=1)/np.sqrt(len(v)))
            ax.scatter(i+off[g]+rng.uniform(-.04,.04,len(v)),v,s=14,marker=MK[g],facecolor='none',edgecolor=COL[g],lw=.9,alpha=.7,zorder=2)
        ax.errorbar(x+off[g],m,yerr=se,color=COL[g],lw=2,marker=MK[g],ms=7,mfc=COL[g],mec=BG,mew=1.5,capsize=3,zorder=3,label=('WT (n = 6; en P30, 5)' if g=='WT' else 'SOD1 (n = 6)'))
    D=P[code]; lo,hi=np.nanmin(A[col]),np.nanmax(A[col]); sp=hi-lo
    for _,r in D.iterrows():
        if r.p_welch<0.05:
            ax.text(DIAS.index(r.dia),hi+0.04*sp,f'p = {r.p_welch:.3f}'+('\nHolm 6 días = %.3f'%r.p_welch_holm_6dias if r.p_welch_holm_6dias<0.05 else ''),ha='center',va='bottom',fontsize=8,color=INK)
    ax.set_ylim(lo-0.06*sp,hi+0.22*sp); ax.set_xticks(x); ax.set_xticklabels(DIAS)
    rr=R.loc[code]
    ax.set_title(f'{tit}\ninteracción genotipo × día p = {rr.p_interaccion_GxD:.4f}  (q FDR = {rr.q_FDR_interaccion_GxD:.3f})',loc='left',fontsize=10,color=INK)
    ax.grid(axis='y',color=GRID,lw=.8); ax.set_axisbelow(True)
    for s in ['top','right']: ax.spines[s].set_visible(False)
    for s in ['left','bottom']: ax.spines[s].set_color('#bdbdb8')
    ax.tick_params(colors=MUTED)
axs[0].set_ylabel('Ángulo de rodilla, cadera–rodilla–tobillo (°)\nmenor = más flexionada',color=INK)
axs[0].legend(frameon=False,fontsize=8,loc='lower left')
fig.suptitle('Ángulo de rodilla longitudinal, WT vs SOD1 (media ± EE; puntos = animales)',x=0.01,ha='left',fontsize=12,color=INK)
fig.text(0.01,0.01,'Media de todos los ciclos válidos (detección con 3 señales), sin suavizado. P30–P65 snapshot 1000; P85 y 821_P51 iteración 5. Pruebas por día: Welch; modelo mixto genotipo × día + (1|animal).',fontsize=7,color=MUTED)
fig.tight_layout(rect=(0,0.03,1,0.95)); fig.savefig(os.path.join(rutas.OUT,'04_Goniometria_rodilla_tobillo_pie','V4_rodilla_longitudinal_3_paneles.png'),facecolor=BG)
