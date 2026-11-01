"""Independent counter-flux axial mode reference, fixed physical K=j/8.

No production moment, disk classifier or numerical eigenvalue solver is used
for the continuous reference. Sign-changing imaginary-frequency roots are
bracketed on a fixed logarithmic gamma list, then solved by SciPy brentq.
This is numerical corroboration, NOT a root census or interval certificate.
The published Minerbo closure mode also has an exact rational gamma^2.
"""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,platform,sys
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import roots_legendre

ROOT=Path(__file__).resolve().parent
FLUX=F(1,3)
CURRENT=2*FLUX
CHI=F(17,45)  # Froustey et al. Eq12 at f=1/3, not exact entropy chi.
BETA=brentq(lambda b: 1/np.tanh(b)-1/b-1/3,.01,10,xtol=1e-14)
SOURCES={
 'linear_input':'https://arxiv.org/html/2101.02745#S3.SS1; Eq8',
 'shift':'https://arxiv.org/html/2311.11968v3#S3.SS1; Eq23',
 'closure':'https://arxiv.org/html/2311.11968v3#S2.SS2; Eqs12,19,24-27',
 'reported_peak':'https://arxiv.org/html/2311.11968v3#S4; Eq39'}
GAMMA_GRID=np.geomspace(1e-7,3.,160)


def densities(u,kind):
    if kind=='linear':return (1+u)/2,(1-u)/2
    norm=BETA/(2*np.sinh(BETA))
    return norm*np.exp(BETA*u),norm*np.exp(-BETA*u)


def signed_density(u,kind):
    if kind=='linear':return u
    return BETA/np.sinh(BETA)*np.sinh(BETA*u)


def integrals_quad(gamma,k,kind):
    # Pair +/-u before integrating, so odd cancellation is analytical.
    a,ae=quad(lambda u:u*signed_density(u,kind)/(gamma*gamma+k*k*u*u),0.,1.,
              epsabs=2e-12,epsrel=2e-12,limit=250)
    b,be=quad(lambda u:u**3*signed_density(u,kind)/(gamma*gamma+k*k*u*u),0.,1.,
              epsabs=2e-12,epsrel=2e-12,limit=250)
    return np.array([-2*k*a,-2j*gamma*a,-2*k*b]),max(ae,be)


def factors(v):
    i0,i1,i2=v
    return ((1+i0)*(-1+i2)-i1*i1,-1+(i0-i2)/2)


def value(gamma,k,kind,sector):
    v,_=integrals_quad(gamma,k,kind)
    return float(factors(v)[sector].real)


def direct_legendre(gamma,k,kind,n):
    u,w=roots_legendre(n)
    g=signed_density(u,kind)
    return np.array([np.dot(w,u**j*g/(1j*gamma-k*u)) for j in range(3)])


def solve_case(K,kind):
    k=float(K-CURRENT)
    result=[]
    for sector,name in ((0,'longitudinal'),(1,'transverse')):
        vals=[value(g,k,kind,sector) for g in GAMMA_GRID]
        roots=[]
        for i in range(len(vals)-1):
            if vals[i]*vals[i+1]>=0:continue
            lo,hi=GAMMA_GRID[i:i+2]
            gamma=brentq(lambda g:value(g,k,kind,sector),lo,hi,xtol=5e-15,rtol=1e-14)
            if roots and abs(gamma-roots[-1]['gamma'])<1e-10:continue
            original,_=integrals_quad(gamma,k,kind)
            check512=direct_legendre(gamma,k,kind,512)
            check256=direct_legendre(gamma,k,kind,256)
            roots.append({'gamma':gamma,'bracket':[float(lo),float(hi)],
                          'continuous_residual':float(abs(factors(original)[sector])),
                          'independent_512_residual':float(abs(factors(check512)[sector])),
                          'integrals_quad_vs_512_maxerr':float(np.max(np.abs(original-check512))),
                          'integrals_256_vs_512_maxerr':float(np.max(np.abs(check256-check512)))})
        result.append({'sector':name,'roots':roots,
                       'status':'sign_changing_imaginary_roots_found' if roots else
                                'no_sign_changing_root_found_in_fixed_gamma_search'})
    return {'physical_K':str(K),'reduced_k':str(K-CURRENT),'completion':kind,'branches':result}


def closure_case(K):
    k=K-CURRENT
    g2=(k+2*FLUX)*(2*FLUX-CHI*k)
    # Independent four-component matrix assembled from Eq24--27.
    f,c,kf=float(FLUX),float(CHI),float(k)
    matrix=np.array([[-1,kf+f,1,-f],[c*kf-f,c,f,-c],
                     [-1,f,1,kf-f],[-f,c,c*kf+f,-c]])
    ev=np.linalg.eigvals(matrix)
    positive=[float(e.imag) for e in ev if e.imag>1e-9]
    pred=np.sqrt(float(g2)) if g2>0 else 0.
    residual=abs(max(positive,default=0.)-pred)
    return {'physical_K':str(K),'reduced_k':str(k),'sector':'longitudinal',
            'gamma_squared':str(g2),'growing':g2>0,'gamma':pred,
            'independent_matrix_peak_error':residual,
            'method':'published polynomial Minerbo/FT quantum closure, not full angular ME dispersion'}


def main():
    ks=[F(j,8) for j in range(-8,9)]
    moment_reference={}
    for kind in ('linear','maximum_entropy'):
        moment_reference[kind]={}
        for species in (0,1):
            vals=[quad(lambda u:u**j*densities(u,kind)[species],-1,1,
                       epsabs=2e-13,epsrel=2e-13)[0] for j in range(3)]
            moment_reference[kind]['plus' if species==0 else 'minus']=vals
    rows=[solve_case(k,kind) for kind in ('linear','maximum_entropy') for k in ks]
    closures=[closure_case(k) for k in ks]
    cm=max(c['independent_matrix_peak_error'] for c in closures)
    if cm>1e-12:raise RuntimeError('closure formula/matrix mismatch')
    d={'scope':'Fixed 17 physical K points, only sign-changing imaginary roots in predeclared gamma range; no all-root census and no interval-certified roots.',
       'physical_K_grid':[str(k) for k in ks], 'shift':'k_reduced=K_physical−2/3; omega unchanged because signed mass is zero',
       'normalization':'N_plus=N_minus=1, mu0=sqrt(2) G_F N; c=1, hbar=1',
       'gamma_search':{'minimum':1e-7,'maximum':3.,'log_nodes':160},
       'beta_maximum_entropy':BETA,'chi_polynomial_Minerbo':str(CHI),
       'moments_numerical':moment_reference,'continuous_modes':rows,'closure_modes':closures,
       'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
       'sources':SOURCES,'python':sys.version,'platform':platform.platform(),
       'method':'SciPy quad + brentq; independently cross-checked direct 256/512-point Gauss-Legendre original integrals; closure exact gamma² independently compared with 4x4 eigenvalues.'}
    out=ROOT/'results/counterflux-modes.json';out.write_text(json.dumps(d,indent=2)+'\n')
    print(json.dumps({'output':str(out),'beta':BETA,'moments':moment_reference,
          'root_counts':[{'completion':kind,'sector':name,'count':sum(len(b['roots']) for r in rows if r['completion']==kind for b in r['branches'] if b['sector']==name)} for kind in ('linear','maximum_entropy') for name in ('longitudinal','transverse')],
          'max_closure_matrix_error':cm,
          'max_direct512_residual':max(q['independent_512_residual'] for r in rows for b in r['branches'] for q in b['roots'])},indent=2))

if __name__=='__main__':main()
