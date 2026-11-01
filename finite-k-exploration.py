"""Bounded exploratory LP attack on a fixed-k equal-Q conjecture.

The endpoint/central atom family supplies an attainable pure-imaginary
rate g2=2(1-a)ka-k²a² below k*=(1-a)/a and (1-a)² above it.
The full all-Re claim is NOT assumed. This file searches specified points
1.001 times above that rate on a finite angular grid; feasibility gives
a candidate counterexample, infeasibility does not prove continuous
infeasibility. No production geometry is imported.
"""
from pathlib import Path
import hashlib,json,time
import numpy as np
from scipy.optimize import linprog

pairs=[(1/3,.5),(.1,.5),(.7,.1),(.1,2)]
real_parts=[0,.1,.25,.5,1]
u=np.linspace(-1,1,257)
rows=[]
for a,k in pairs:
 cap=(1-a)**2 if k >= (1-a)/a else 2*(1-a)*k*a-k*k*a*a
 gamma=np.sqrt(cap)*1.001
 for x in real_parts:
  h=(1-u*u)/(x+1j*gamma-k*u)
  zero=np.zeros(len(u));one=np.ones(len(u))
  A=np.array([np.r_[one,zero],np.r_[u,zero],np.r_[u*u,zero],
              np.r_[zero,one],np.r_[zero,u],np.r_[zero,u*u],
              np.r_[h.real,-h.real],np.r_[h.imag,-h.imag]])
  b=np.array([1,0,a,1,0,a,-1,0])
  t=time.perf_counter();r=linprog(np.zeros(2*len(u)),A_eq=A,b_eq=b,bounds=(0,None),method='highs')
  row=dict(a=a,k=k,real=x,gamma=gamma,attainer_gamma_squared=cap,grid_nodes=len(u),
           lp_status=int(r.status),seconds=time.perf_counter()-t)
  if r.success:
   row['residual']=float(np.max(np.abs(A@r.x-b)))
   row['atoms']=[dict(species=s,u=float(u[j]),weight=float(r.x[s*len(u)+j]))
     for s in range(2) for j in range(len(u)) if r.x[s*len(u)+j]>1e-10]
  rows.append(row)
p=Path('results/finite-k-exploration.json');p.write_text(json.dumps(dict(scope=__doc__,rows=rows,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),indent=2)+'\n')
print(json.dumps({'cases':len(rows),'grid_feasible_above_family':sum(r['lp_status']==0 for r in rows),'artifact':str(p)}))
