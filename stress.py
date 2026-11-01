"""Predetermined exact boundary/scaling cases and checked grid Farkas rays."""
from fractions import Fraction as F
from pathlib import Path
import json,time,hashlib
import numpy as np
from scipy.optimize import linprog
from moments import QC,classify,classify_numeric,root_target
from independent_verifier import verify
from baselines import grid_lp
from screen import encode
ROOT=Path(__file__).resolve().parent

def grid_farkas(N1,F1,N2,F2,w,k,sector,n):
    target=root_target(N1,F1,N2,F2,w,k,sector);z=w/k
    if target is None:return None
    cols=[]
    for species in (0,1):
        for j in range(n):
            u=F(2*j,n-1)-1;s=1/(z-u);sign=1 if species==0 else -1
            cols.append([F(species==0),u if species==0 else F(0),F(species==1),u if species==1 else F(0),sign*s.re,sign*s.im])
    A=np.array([[float(v) for v in c] for c in cols]);b=list(map(F,(N1,F1,N2,F2)))+[target.re,target.im]
    candidate=linprog(np.array(list(map(float,b))),A_ub=-A,b_ub=np.zeros(2*n),bounds=[(-1,1)]*6,method='highs')
    if not candidate.success:return None
    for d in (10**6,10**9,10**12):
        y=[F(float(v)).limit_denominator(d) for v in candidate.x]
        floor=min(sum(v*t for v,t in zip(c,y)) for c in cols)
        correction=max(F(0),-floor);y[0]+=correction;y[2]+=correction
        products=[sum(v*t for v,t in zip(c,y)) for c in cols]
        objective=sum(v*t for v,t in zip(b,y))
        if min(products)>=0 and objective<0:
            return {'dual':y,'minimum_column_product':min(products),'negative_rhs_product':objective}
    return None

def main():
    rows=[];t=time.perf_counter()
    offsets=(F(0),)+tuple(s*d for d in (F(1,100),F(1,10000),F(1,10**8),F(1,10**12)) for s in (-1,1))
    for sector,k0,gamma0 in (('longitudinal',F(2,3),F(1,2)),('transverse',F(1,3),F(1,4))):
        for scale in (F(1,10**50),F(1,10**6),F(1),F(10**6),F(10**50)):
            for delta in offsets:
                p={'N1':scale,'F1':F(0),'N2':scale,'F2':F(0),'w':QC(0,scale*(gamma0+delta)),
                   'k':scale*k0,'sector':sector}
                r=classify(p['N1'],p['F1'],p['N2'],p['F2'],p['w'],p['k'],sector);verify(p,r)
                expected='possible' if delta<=0 else 'impossible'
                if r['status']!=expected:raise RuntimeError('sharp threshold contradiction')
                row={'problem':p,'delta':delta,'expected':expected,'result':r}
                if scale==1:
                    grids=[]
                    for n in (32,128,512):
                        g=grid_lp(scale,0,scale,0,p['w'],p['k'],sector,n)
                        if g['status']!='exact_grid_witness':g['farkas']=grid_farkas(scale,0,scale,0,p['w'],p['k'],sector,n)
                        grids.append({'nodes':n,**g})
                    row['grids']=grids
                    try:row['numeric_result']=classify_numeric(scale,0,scale,0,p['w'],p['k'],sector)
                    except Exception as exc:row['numeric_result']={'status':'unresolved','error_class':type(exc).__name__}
                    if row['numeric_result']['status']!='unresolved':verify(p,row['numeric_result'])
                rows.append(row)
        print(sector,'completed',len(rows),flush=True)
    report={'scope':'designed exact-boundary/scaling adversarial cases, not random physical generalization',
            'cases':rows,'passed':True,'elapsed_seconds':time.perf_counter()-t,
            'source_hashes':{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in ('stress.py','ADVERSARIAL_PLAN.md','moments.py','geometry.py','independent_verifier.py')}}
    (ROOT/'results'/'stress.json').write_text(json.dumps(encode(report),indent=2)+'\n')
    print('stress passed',len(rows),flush=True)
if __name__=='__main__':main()
