"""Compare the complete envelope with ALL three previously proved upper caps.
Exact rational arithmetic; this is a bound comparison, not an impact benchmark.
"""
from fractions import Fraction as F
from pathlib import Path
import sys,json,hashlib
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from finite_band_envelope import envelope
ROOT=Path(__file__).resolve().parent

def cubic(a,k):
    lo,hi=F(0),F(1)
    for _ in range(100):
        m=(lo+hi)/2;p=m**3-a*m*m+a*k*k*m-a*k*k
        if p<0:lo=m
        else:hi=m
    return lo,hi

def run():
    rows=[];a=F(1,3)
    ks=sorted(set([F(j,50) for j in range(1,151)]+[F(1,10**j) for j in range(1,19)]))
    for k in ks:
        result=envelope(1,a,k,bits=100)
        full2=F(result['exact_growth_squared'])
        m1sq=k/(k+2);allksq=(1-a)**2
        clo,chi=cubic(a,k)
        baseline_lo=min(m1sq,allksq,clo*clo);baseline_hi=min(m1sq,allksq,chi*chi)
        # The two-cap comparison is an exact analytic upper bound on all3.
        two_cap=min(m1sq,allksq)
        if not full2<=two_cap or not 8*two_cap<=9*full2:
            raise RuntimeError('Isotropic analytic inequality failed')
        if baseline_hi<full2:raise RuntimeError('Old upper cap below sharp value')
        rows.append({'k':str(k),'exact_sharp_growth_squared':str(full2),'m1_upper_squared':str(m1sq),'all_k_upper_squared':str(allksq),'combined_upper_squared_interval':[str(baseline_lo),str(baseline_hi)],'old_over_new_ratio_squared_interval':[str(baseline_lo/full2),str(baseline_hi/full2)]})
    out={'scope':'IsotropicQ=N/3, normalizedN1longitudinal; transverse followsk->2k andgrowth/2; no physicaldata or solverruntime claim','cases':len(rows),'analytic_maximum_old_over_new_ratio_squared':'9/8','maximum_relative_tightening_expression':'1-sqrt(8/9)','maximum_relative_tightening_approx':1-(8/9)**.5,'records':rows,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (ROOT/'isotropic-bound-comparison.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k!='records'},indent=2))
if __name__=='__main__':run()
