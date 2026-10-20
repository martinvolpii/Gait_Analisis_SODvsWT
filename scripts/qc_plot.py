import sys, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import os, rutas; import marcha_v2 as M
ses=dict((f'{a}_{d}',(f,d,a,fps)) for f,d,a,fps in M.sesiones())
pick=['821_P30','862_P51','866_P85']
fig,axs=plt.subplots(len(pick),2,figsize=(14,3.2*len(pick)),gridspec_kw=dict(width_ratios=[1,1]))
for i,s in enumerate(pick):
    f,d,a,fps=ses[s]; info,cyc,aux=M.procesar_sesion(f,d,a,fps,suave=False)
    h=aux['h']; t=np.arange(len(h))/aux['fps']; C=aux['C']
    vrel=np.gradient(C['toe'][:,0]-C['hip'][:,0])
    ax=axs[i,0]; ax.plot(t,h,color='#2a78d6',lw=1.5,label='altura del dedo (px)')
    for p in aux['peaks']: ax.axvline(p/aux['fps'],color='#888',lw=0.6)
    st=(vrel<0)
    ax.fill_between(t,np.nanmin(h),np.nanmax(h),where=st,color='#eb6834',alpha=0.15,step='mid',label='apoyo (dedo hacia atrás)')
    ax.set_title(f'{s}: picos = límites de ciclo; T={info["T_s"]:.3f} s; válidos {info["n_ciclos_validos"]}/{info["n_ciclos_detectados"]}',fontsize=9)
    ax.legend(fontsize=7,loc='upper right'); ax.set_xlim(0,min(t[-1],3))
    ax=axs[i,1]
    for k,c in [('rodilla','#2a78d6'),('rodilla_tri','#1baf7a'),('tobillo','#eb6834'),('pie','#e87ba4')]:
        ax.plot(t,aux['ang'][k],lw=1.2,color=c,label=k)
    for p in aux['peaks']: ax.axvline(p/aux['fps'],color='#888',lw=0.6)
    ax.set_xlim(0,min(t[-1],3)); ax.legend(fontsize=7,ncol=4); ax.set_ylabel('grados')
plt.tight_layout(); plt.savefig(os.path.join(rutas.INTER,'qc_ciclos.png'),dpi=90)
