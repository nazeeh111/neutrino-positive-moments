"""Validated positive smooth angular completions with identical M0,M1,M2."""
from fractions import Fraction as F
from pathlib import Path
import json,time,hashlib
from flint import arb,acb,ctx,fmpq
ROOT=Path(__file__).resolve().parent

def rb(q):
    q=F(q);return arb(fmpq(q.numerator,q.denominator))

def density(u,kappa):
    # Density with respect to du after integrating the azimuth, mass one.
    c=arb(kappa)/(2*arb(kappa).sinh())
    ring=(u*kappa/3).exp()*(kappa*kappa*(1-u*u)*arb(2)/9).hypgeom_0f1(1)
    pole=(-u*kappa).exp()
    return c*(ring*arb(3)/4+pole/4)

def integrate(f,a,b):
    out=acb.integral(lambda u,analytic:f(u),a,b,rel_tol=arb('1e-25'),abs_tol=arb('1e-25'),
                     eval_limit=200000,depth_limit=30)
    if not out.is_finite():raise ArithmeticError('integration enclosure failed')
    return out

def dispersion(gamma,kappa):
    g=rb(gamma)
    val=integrate(lambda u:u*(density(u,kappa)-density(-u,kappa))/(g*g+u*u),0,1)
    return 2*(g*g+1)*val.real-1

def main():
    ctx.prec=160;rows=[];started=time.perf_counter()
    samples=tuple(map(F,('1/20','1/10','1/5','3/10','2/5','1/2','3/5','7/10')))
    for kappa in (2,4,8,16,32,64):
        row={'kappa':kappa};t=time.perf_counter()
        try:
            moment=[integrate(lambda u,j=j:u**j*density(u,kappa),-1,1) for j in range(3)]
            row['moment_enclosures']=[str(x) for x in moment]
            row['moment_checks']=[moment[0].real.contains(1),moment[1].real.contains(0),
                                   (moment[2].real-rb(F(1,3))).contains(0)]
            if not all(row['moment_checks']):raise ArithmeticError('smooth moment crosscheck failed')
            values=[dispersion(g,kappa) for g in samples]
            row['sample_signs']=[{'gamma':str(g),'value':str(v)} for g,v in zip(samples,values)]
            bracket=next(((samples[j],samples[j+1]) for j in range(len(samples)-1)
                          if values[j]>0 and values[j+1]<0),None)
            if bracket is None:row['status']='no_certified_bracket'
            else:
                lo,hi=bracket
                while hi-lo>F(1,10000):
                    mid=(lo+hi)/2;v=dispersion(mid,kappa)
                    if v>0:lo=mid
                    elif v<0:hi=mid
                    else:raise ArithmeticError('unresolved dispersion sign')
                row.update(status='certified_growing_mode',gamma_interval=[str(lo),str(hi)],
                           endpoint_values=[str(dispersion(lo,kappa)),str(dispersion(hi,kappa))])
        except (ArithmeticError,ValueError) as exc:row.update(status='refused',reason=str(exc))
        row['seconds']=time.perf_counter()-t;rows.append(row)
        print(kappa,row['status'],row.get('gamma_interval'),flush=True)
    qualifying=sum(r['kappa']<=32 and r['status']=='certified_growing_mode' and F(r['gamma_interval'][0])>F(1,4) for r in rows)
    report={'scope':'smooth positive axial completions, same first/second moments, k=N=1',
            'precision_bits':ctx.prec,'cases':rows,'qualifying_cases':qualifying,'gate_passed':qualifying>=2,
            'elapsed_seconds':time.perf_counter()-started,
            'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'plan_sha256':hashlib.sha256((ROOT/'SMOOTH_VALIDATION_PLAN.md').read_bytes()).hexdigest()}
    (ROOT/'results'/'smooth.json').write_text(json.dumps(report,indent=2)+'\n')
    print('smooth gate',report['gate_passed'],flush=True)
if __name__=='__main__':main()
