"""Source-bound publication figures; all curves/points trace to proofs/results."""
import os,json
from pathlib import Path
from fractions import Fraction
ROOT=Path(__file__).resolve().parents[1]
os.environ.setdefault('MPL_IGNORE_SYSTEM_FONTS','1')
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'.plot-cache'))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,
                     'axes.spines.right':False,'svg.fonttype':'path','savefig.bbox':'tight'})
from scipy.special import i0e
out=ROOT/'paper'
k=np.linspace(.001,4,500)
fig,ax=plt.subplots(figsize=(6.4,2.6))
L=np.sqrt(k/(k+2));T=np.sqrt(k/(4*(k+1)))
ax.fill_between(k,0,L,color='#dce7f5',alpha=.65)
ax.plot(k,L,color='#204e82',lw=2,label='Longitudinal envelope')
ax.plot(k,T,color='#a35b22',lw=2,ls='--',label='Transverse envelope')
ax.set(xlabel='Reduced wave number k / N',ylabel='Attainable imaginary-axis growth γ / N',xlim=(0,4),ylim=(0,.9))
ax.grid(alpha=.15);ax.legend(frameon=False,loc='lower right')
fig.tight_layout();fig.savefig(out/'sharp-bands.svg');fig.savefig(out/'sharp-bands.png',dpi=180);plt.close(fig)
smooth=json.loads((ROOT/'results/smooth.json').read_text())
fig,(a,b)=plt.subplots(1,2,figsize=(6.4,2.6),gridspec_kw={'width_ratios':[1.1,1]})
u=np.linspace(-1,1,800)
colors=['#577ca0','#2d6183','#163e59']
for kap,col in zip((8,16,32),colors):
    v=kap*np.sqrt((1-u*u)*8/9)
    f=kap/(1-np.exp(-2*kap))*(.75*np.exp(kap*u/3-kap+v)*i0e(v)+.25*np.exp(-kap*u-kap))
    a.plot(u,f,label=f'κ = {kap}',color=col,lw=1.4)
a.axhline(.5,color='#777777',ls=':',lw=1,label='Stable uniform pair')
a.set(xlabel='Polar cosine u',ylabel='Electron species density per du',xlim=(-1,1));a.legend(frameon=False,fontsize=7)
ks=[];lo=[];hi=[]
for r in smooth['cases']:
    if r['status']=='certified_growing_mode':
        ks.append(r['kappa']);lo.append(float(Fraction(r['gamma_interval'][0])));hi.append(float(Fraction(r['gamma_interval'][1])))
mid=(np.array(lo)+np.array(hi))/2
b.errorbar(ks,mid,yerr=[mid-np.array(lo),np.array(hi)-mid],fmt='o-',color='#204e82',capsize=3,lw=1)
b.axhline(1/np.sqrt(3),color='#a35b22',ls='--',lw=1,label='Sharp measure envelope')
b.set(xscale='log',xlabel='Smooth-kernel concentration κ',ylabel='Certified growing mode γ / N',ylim=(.2,.63))
b.set_xticks(ks,[str(x) for x in ks]);b.minorticks_off();b.legend(frameon=False,fontsize=7,loc='lower right')
for ax in (a,b):ax.grid(alpha=.15)
fig.tight_layout();fig.savefig(out/'smooth-completions.svg');fig.savefig(out/'smooth-completions.png',dpi=180);plt.close(fig)
