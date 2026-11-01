"""Frozen continuous comparison, independent of any angular grid.

Generic degree-three Hausdorff SOCP via Clarabel, exact certificate extraction,
and ordinary floating Weyl geometry. No exact-classifier fallback in baselines.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json
import math
import platform
import random
import statistics
import sys
import time

import clarabel
import numpy as np
import scipy
from scipy import sparse

from moments import (QC, classify, root_target, hausdorff_affine, shifted_scaled,
                     attaining, dual_from_directions)
from geometry import disks, Quad
from independent_verifier import verify

ROOT = Path(__file__).resolve().parent
REPEATS = 3
MARGIN_TOL = 1e-8
FLOAT_CONTAIN_TOL = 1e-9

def encode(x):
    if isinstance(x, F): return str(x)
    if isinstance(x, QC): return {'re': encode(x.re), 'im': encode(x.im)}
    if isinstance(x, Quad): return {'a': str(x.a), 'b': str(x.b), 'd': str(x.d)}
    if isinstance(x, dict): return {str(k): encode(v) for k, v in x.items()}
    if isinstance(x, (tuple, list)): return [encode(v) for v in x]
    return x

def decode_problem(p):
    return dict(N1=F(p['N1']), F1=F(p['F1']), N2=F(p['N2']), F2=F(p['F2']),
                w=QC(F(p['w']['re']), F(p['w']['im'])), k=F(p['k']), sector=p['sector'])

def suite():
    cases=[]
    for filename, group in (('screen.json','development_screen'), ('stress.json','designed_stress')):
        data=json.loads((ROOT/'results'/filename).read_text())
        rows=data['cases']
        for i,row in enumerate(rows):
            cases.append((f'{group}:{i}', group, decode_problem(row['problem'])))
    rng=random.Random(271828)
    masses=[F(1,2),F(1),F(3,2),F(2)]
    flux_ratios=[F(-1),F(-3,4),F(-1,3),F(0),F(1,3),F(3,4),F(1)]
    real=[F(i,3) for i in range(-9,10)]
    imag=[F(1,100),F(1,20),F(1,5),F(1,2),F(1),F(2)]
    wave=[F(1,4),F(1,2),F(1),F(2)]
    for i in range(128):
        n1,n2=rng.choice(masses),rng.choice(masses)
        p=dict(N1=n1,F1=n1*rng.choice(flux_ratios),N2=n2,F2=n2*rng.choice(flux_ratios),
               w=QC(rng.choice(real),rng.choice(imag)),k=rng.choice(wave),
               sector=rng.choice(['longitudinal','transverse']))
        cases.append((f'fresh_synthetic:{i}','fresh_synthetic',p))
    return cases

def normalized_lmis(p):
    target=root_target(p['N1'],p['F1'],p['N2'],p['F2'],p['w'],p['k'],p['sector'])
    if target is None: return None,None,None,None
    z=p['w']/p['k']; scale=p['N1']/z.im
    if scale<=0: raise ValueError('planned comparator requires positive first-species mass')
    hs=[shifted_scaled(h,QC(),scale) for h in hausdorff_affine(p['N1'],p['F1'],z)]
    hs += [shifted_scaled(h,target,scale) for h in hausdorff_affine(p['N2'],p['F2'],z)]
    hs=[tuple(tuple(x/max(F(1),max(abs(t) for part in h for t in part)) for x in part)
              for part in h) for h in hs]
    return z,target,scale,hs

def cone_vector(q): return (q[0]+q[2],q[0]-q[2],2*q[1])

def conic_proposal(p):
    z,target,scale,hs=normalized_lmis(p)
    if target is None:
        return {'status':'impossible','proof_kind':'zero_A_nonzero_B'},None
    matrices=[];rhs=[]
    for h in hs:
        magnitude=max(abs(t) for part in h for t in part)
        # v already scales S by N1/Im(z); up/down matrix normalization now
        # removes the remaining arbitrary common physical unit amplitude.
        C,A,B=tuple(tuple(t/magnitude for t in part) for part in h) if magnitude else h
        cv,av,bv=map(cone_vector,(C,A,B))
        matrices.extend([[-float(av[j]),-float(bv[j]),2. if j==0 else 0.] for j in range(3)])
        rhs.extend(map(float,cv))
    matrices.extend([[1.,0.,0.],[-1.,0.,0.],[0.,1.,0.],[0.,-1.,0.]])
    rhs.extend([1.,1.,0.,1.])
    settings=clarabel.DefaultSettings()
    settings.verbose=False;settings.max_iter=200
    for name in ('tol_gap_abs','tol_gap_rel','tol_feas','tol_infeas_abs','tol_infeas_rel'):
        setattr(settings,name,1e-10)
    if hasattr(settings,'max_threads'): settings.max_threads=1
    solver=clarabel.DefaultSolver(sparse.csc_matrix((3,3)),np.array([0.,0.,-1.]),
              sparse.csc_matrix(matrices),np.array(rhs),
              [clarabel.SecondOrderConeT(3) for _ in range(4)]+[clarabel.NonnegativeConeT(4)],settings)
    sol=solver.solve()
    v=np.asarray(sol.x,dtype=float)
    if len(v)!=3 or not np.all(np.isfinite(v)):
        return {'status':'unresolved','solver_status':str(sol.status),'reason':'nonfinite or absent solution'},None
    margin=float(v[2])
    suggestion='possible' if margin>MARGIN_TOL else ('impossible' if margin < -MARGIN_TOL else 'unresolved')
    diagnostics={'solver_status':str(sol.status),'solver_margin':margin,
                 'floating_suggestion':suggestion,'solver_iterations':int(sol.iterations),
                 'solver_internal_seconds':float(sol.solve_time)}
    xy=tuple(F(float(q)) for q in v[:2])
    for denom in (10**6,10**9,10**12,None):
        q=tuple(a.limit_denominator(denom) if denom else a for a in xy)
        s1=QC(scale*q[0],scale*q[1]);s2=s1-target
        a=attaining(p['N1'],p['F1'],z,s1);b=attaining(p['N2'],p['F2'],z,s2)
        if a is not None and b is not None:
            return {'status':'possible','s1':s1,'s2':s2,'witnesses':[a,b],**diagnostics},sol
    dual=dual_from_directions(hs,xy)
    if dual is not None:
        return {'status':'impossible','proof_kind':'rational_dual','dual':dual,**diagnostics},sol
    return {'status':'unresolved',**diagnostics},sol

def conic(p): return conic_proposal(p)[0]

def floating_geometry(p):
    target=root_target(p['N1'],p['F1'],p['N2'],p['F2'],p['w'],p['k'],p['sector'])
    if target is None:return {'numerical_status':'impossible','margin':None}
    z=p['w']/p['k']
    ds=disks(p['N1'],p['F1'],z)+[((c[0]+target.re,c[1]+target.im),r2)
                              for c,r2 in disks(p['N2'],p['F2'],z)]
    # Exact common rescaling makes common unit changes irrelevant before float conversion.
    scale=max(F(1),*(abs(t) for c,r2 in ds for t in c),*(abs(r2) for c,r2 in ds))
    # sqrt(scale^2) is deliberately unnecessary: any positive common coordinate scale works.
    cc=np.array([[float(a/scale) for a in c] for c,r2 in ds])
    rr=np.array([math.sqrt(float(r2/(scale*scale))) for c,r2 in ds])
    spread=max(float(np.max(np.abs(cc))),float(np.max(rr)),1e-300)
    cc/=spread;rr/=spread
    candidates=[c.copy() for c in cc]
    for i in range(4):
        for j in range(i+1,4):
            d=cc[j]-cc[i];dist=float(np.dot(d,d))
            if dist==0: continue
            aa=(rr[i]*rr[i]-rr[j]*rr[j]+dist)/(2*dist)
            disc=rr[i]*rr[i]/dist-aa*aa
            if disc < -FLOAT_CONTAIN_TOL: continue
            rad=math.sqrt(max(0.,disc));perp=np.array([-d[1],d[0]])
            candidates.extend([cc[i]+aa*d+rad*perp,cc[i]+aa*d-rad*perp])
    best=max(float(np.min(rr-np.linalg.norm(cc-c,axis=1))) for c in candidates)
    return {'numerical_status':'possible' if best >= -FLOAT_CONTAIN_TOL else 'impossible',
            'margin':best,'containment_tolerance':FLOAT_CONTAIN_TOL}

def summarize(rows):
    out={}
    for group in ['all','development_screen','designed_stress','fresh_synthetic']:
        rs=[r for r in rows if group=='all' or r['group']==group]
        summary={'count':len(rs)}
        for method in ['exact','conic','floating_geometry']:
            key='numerical_status' if method=='floating_geometry' else 'status'
            counts={s:sum(r[method][key]==s for r in rs) for s in ['possible','impossible','unresolved']}
            counts['median_generation_seconds']=statistics.median(r['timing'][method] for r in rs)
            summary[method]=counts
        summary['conic_exact_disagreements']=sum(r['conic']['status']!='unresolved' and r['conic']['status']!=r['exact']['status'] for r in rs)
        summary['float_exact_disagreements']=sum(r['floating_geometry']['numerical_status']!=r['exact']['status'] for r in rs)
        common=[r for r in rs if r['conic']['status']!='unresolved' and r['exact']['status']!='unresolved']
        summary['jointly_certified']=len(common)
        if common:
            summary['median_conic_over_exact_generation_ratio_joint']=statistics.median(r['timing']['conic']/r['timing']['exact'] for r in common)
        out[group]=summary
    return out

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',default='results/continuous-baseline.json');args=parser.parse_args()
    cases=suite();tasks=list(range(len(cases)));rng=random.Random(1729)
    methods={'exact':lambda p:classify(**p),'conic':conic,'floating_geometry':floating_geometry}
    # One untimed smoke per backend, fixed in advance to the first case; no tuning from outcomes.
    for name,fn in methods.items(): fn(cases[0][2])
    samples={i:{name:[] for name in methods} for i in tasks};answers={}
    started=time.perf_counter()
    for rep in range(REPEATS):
        rng.shuffle(tasks)
        for i in tasks:
            order=list(methods);rng.shuffle(order)
            for name in order:
                t=time.perf_counter();res=methods[name](cases[i][2]);dt=time.perf_counter()-t
                samples[i][name].append(dt)
                if rep==0:answers[i,name]=res
                else:
                    key='numerical_status' if name=='floating_geometry' else 'status'
                    if res[key]!=answers[i,name][key]:raise RuntimeError('nondeterministic classification')
        print(f'repetition {rep+1}/{REPEATS} complete ({len(cases)} cases)',flush=True)
    rows=[]
    for i,(name,group,p) in enumerate(cases):
        checks={}
        for method in ['exact','conic']:
            answer=answers[i,method]
            if answer['status']!='unresolved':
                t=time.perf_counter();verify(p,answer);checks[method]=time.perf_counter()-t
        if answers[i,'conic']['status']!='unresolved' and answers[i,'conic']['status']!=answers[i,'exact']['status']:
            raise RuntimeError('independently certified disagreement')
        rows.append({'id':name,'group':group,'problem':p,**{m:answers[i,m] for m in methods},
                     'generation_seconds':samples[i], 'timing':{m:statistics.median(v) for m,v in samples[i].items()},
                     'independent_check_seconds':checks})
    result={'scope':'controlled continuous mathematical comparison; reused development/stress plus deterministic synthetic cases',
       'plan_sha256':hashlib.sha256((ROOT/'independent-continuous-review.md').read_bytes()).hexdigest(),
       'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
       'input_sha256':{f:hashlib.sha256((ROOT/'results'/f).read_bytes()).hexdigest() for f in ['screen.json','stress.json']},
       'runtime':{'python':sys.version,'platform':platform.platform(),'numpy':np.__version__,'scipy':scipy.__version__,'clarabel':clarabel.__version__},
       'settings':{'repeats':REPEATS,'threads':1,'conic_tolerances':1e-10,'margin_tolerance':MARGIN_TOL,
                   'floating_containment_tolerance':FLOAT_CONTAIN_TOL,'case_seed':271828,'order_seed':1729},
       'summary':summarize(rows),'cases':rows,'elapsed_seconds':time.perf_counter()-started}
    path=ROOT/args.output;path.parent.mkdir(exist_ok=True);path.write_text(json.dumps(encode(result),indent=2)+'\n')
    print(json.dumps(result['summary'],indent=2),flush=True)
if __name__=='__main__':main()
