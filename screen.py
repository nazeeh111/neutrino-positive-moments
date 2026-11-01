"""Predeclared504-point screen with exact independent certificate checking."""
from fractions import Fraction as F
from pathlib import Path
import json,hashlib,time
import importlib.util
from moments import QC,classify,classify_numeric
from geometry import Quad
from independent_verifier import verify
from baselines import grid_lp,zero_mode
ROOT=Path(__file__).resolve().parent

def encode(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,Quad):return {'a':str(x.a),'b':str(x.b),'d':str(x.d)}
    if isinstance(x,QC):return {'re':encode(x.re),'im':encode(x.im)}
    if isinstance(x,dict):return {str(k):encode(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)):return [encode(v) for v in x]
    return x

def main():
    cases=(('equal_zero_flux',(F(1),F(0),F(1),F(0))),
           ('fiducial_axial',(F(1),F(1,3),F(1),F(-1,3))),
           ('single_species',(F(1),F(1,3),F(0),F(0))))
    real=(F(-3,2),F(-1),F(-1,2),F(0),F(1,2),F(1),F(3,2))
    imag=(F(1,20),F(1,10),F(1,5),F(2,5),F(4,5),F(6,5))
    rows=[];started=time.perf_counter();(ROOT/'results').mkdir(exist_ok=True)
    for name,m in cases:
        for k in (F(1),F(1,2)):
            for x in real:
                for y in imag:
                    for sector in ('longitudinal','transverse'):
                        w=QC(x,y);data=dict(zip(('N1','F1','N2','F2'),m));data.update(w=w,k=k,sector=sector)
                        t=time.perf_counter();answer=classify(*m,w,k,sector);algorithm_time=time.perf_counter()-t
                        t=time.perf_counter()
                        if answer['status']!='unresolved':verify(data,answer)
                        check_time=time.perf_counter()-t
                        t=time.perf_counter()
                        try:numeric=classify_numeric(*m,w,k,sector)
                        except Exception as exc:numeric={'status':'unresolved','error_class':type(exc).__name__}
                        numeric_time=time.perf_counter()-t
                        if numeric['status']!='unresolved':verify(data,numeric)
                        if numeric['status']!='unresolved' and answer['status']!='unresolved' and numeric['status']!=answer['status']:
                            raise RuntimeError('exact backends disagree')
                        grids=[]
                        for n in (32,128,512):
                            t=time.perf_counter();r=grid_lp(*m,w,k,sector,n);dt=time.perf_counter()-t
                            if r['status']=='exact_grid_witness' and answer['status']=='impossible':
                                raise RuntimeError('continuous impossibility contradicts exact grid witness')
                            grids.append({'nodes_per_species':n,'seconds':dt,**r})
                        rows.append({'case':name,'problem':data,'result':answer,'algorithm_seconds':algorithm_time,
                                     'independent_check_seconds':check_time,'numeric_result':numeric,
                                     'numeric_seconds':numeric_time,'grids':grids})
            print(name,str(k),'completed',len(rows),flush=True)
            (ROOT/'results'/'screen_partial.json').write_text(json.dumps(encode({'cases':rows}),indent=2)+'\n')
    certified=sum(r['result']['status']!='unresolved' for r in rows)
    zero_attained=len({(str(r['problem']['k']),str(r['problem']['w'].re),str(r['problem']['w'].im))
        for r in rows if r['case']=='equal_zero_flux' and r['result']['status']=='possible'})
    single_bad=sum(r['case']=='single_species' and r['result']['status']=='possible' for r in rows)
    result={'scope':'504 predeclared axial fixed-frequency existential classifications',
            'cases':rows,'certified_decisions':certified,'equal_zero_flux_attained_frequencies':zero_attained,
            'single_species_false_positives':single_bad,
            'gate_passed':certified>=.9*len(rows) and zero_attained>=4 and single_bad==0,
            'elapsed_seconds':time.perf_counter()-started,
            'zero_mode_baselines':{name:zero_mode(*m) for name,m in cases},
            'conic_baseline_available':importlib.util.find_spec('clarabel') is not None,
            'source_hashes':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
                for p in ('moments.py','geometry.py','baselines.py','screen.py','independent_verifier.py','PREREGISTRATION.md')}}
    (ROOT/'results'/'screen.json').write_text(json.dumps(encode(result),indent=2)+'\n')
    print('gate',result['gate_passed'],'certified',certified,'zero_flux_attained',zero_attained,flush=True)
if __name__=='__main__':main()
