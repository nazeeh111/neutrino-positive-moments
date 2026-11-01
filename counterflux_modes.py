"""Independent counter-flux axial mode reference, fixed physical K=j/8.

No production moment, disk classifier or numerical eigenvalue solver is used
for the continuous reference. Sign-changing imaginary-frequency roots are
bracketed on a fixed logarithmic gamma list, then solved by SciPy brentq.
The root search is not a census. Saved LINEAR-completion roots additionally
receive Arb sign/IVT brackets; entropy roots remain numerical references.
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
    # Subtract the regular G(u)/u value at zero. The remaining adaptive
    # integrand has no narrow Lorentzian layer when gamma is very small.
    h0=1. if kind=='linear' else BETA*BETA/np.sinh(BETA)
    def h(u):
        if kind=='linear':return 1.
        if abs(u)<1e-8:return h0*(1+(BETA*u)**2/6)
        return signed_density(u,kind)/u
    if k==0:
        a,ae=quad(lambda u:u*u*h(u)/(gamma*gamma),0.,1.,
                  epsabs=2e-12,epsrel=2e-12)
    else:
        a0=h0/(k*k)*(1-gamma/abs(k)*np.arctan(abs(k)/gamma))
        rem,ae=quad(lambda u:u*u*(h(u)-h0)/(gamma*gamma+k*k*u*u),0.,1.,
                    epsabs=2e-12,epsrel=2e-12,limit=250)
        a=a0+rem
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
                     [-1,-f,1,kf+f],[f,c,c*kf-f,-c]])
    ev=np.linalg.eigvals(matrix)
    positive=[float(e.imag) for e in ev if e.imag>1e-9]
    pred=np.sqrt(float(g2)) if g2>0 else 0.
    # At K=0 the zero eigenvalue is defective. Floating eigvals can split it
    # into a fake conjugate pair ~sqrt(eps); use coefficient identity there.
    poly=np.poly(matrix)
    expected=np.polymul([1.,0.,-c*kf*kf],[1.,0.,float(g2)])
    poly_error=float(np.max(np.abs(poly-expected)))
    residual=abs(max(positive,default=0.)-pred) if g2>0 else 0.
    return {'physical_K':str(K),'reduced_k':str(k),'sector':'longitudinal',
            'gamma_squared':str(g2),'growing':g2>0,'gamma':pred,
            'independent_matrix_peak_error':residual,'independent_matrix_characteristic_maxerr':poly_error,
            'floating_eigenvalue_max_imaginary':max(positive,default=0.),
            'degenerate_zero_mode':g2==0,
            'method':'published polynomial Minerbo/FT quantum closure, not full angular ME dispersion'}


def validate_linear_brackets(rows):
    """Arb sign enclosures plus IVT, using the analytical linear-G integral.

    This certifies existence of each saved positive root, not uniqueness or
    absence of off-axis roots. No angular discretization enters the proof.
    """
    from flint import arb,ctx,fmpq
    ctx.prec=160
    def ball(q):
        q=F(q)
        return arb(fmpq(q.numerator,q.denominator))
    def original_factor(gamma,k,sector):
        g,kb=ball(gamma),ball(k)
        if k==0:
            A=1/(3*g*g);B=1/(5*g*g)
        else:
            kk=abs(kb)
            A=(1-g/kk*(kk/g).atan())/(kb*kb)
            B=(ball(F(1,3))-g*g*A)/(kb*kb)
        if sector=='longitudinal':
            return (1-2*kb*A)*(-1-2*kb*B)+4*g*g*A*A
        return -1-kb*(A-B)
    certified=[]
    for row in rows:
        if row['completion']!='linear':continue
        k=F(row['reduced_k'])
        for branch in row['branches']:
            for root in branch['roots']:
                base=int(np.floor(root['gamma']*10**8))
                lo,hi=F(base,10**8),F(base+1,10**8)
                vl,vh=(original_factor(g,k,branch['sector']) for g in (lo,hi))
                if not ((vl>0 and vh<0) or (vl<0 and vh>0)):
                    raise ArithmeticError('Arb did not validate the proposed linear-root bracket')
                cert={'physical_K':row['physical_K'],'reduced_k':str(k),
                      'sector':branch['sector'],'gamma_interval':[str(lo),str(hi)],
                      'endpoint_values':[str(vl),str(vh)],'precision_bits':160,
                      'method':'closed-form original dispersion, Arb strict signs and IVT'}
                root['validated_gamma_interval']=cert['gamma_interval']
                certified.append(cert)
    return certified


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
    qt=(1-CHI)/2
    tm=((qt,-qt),(qt,-qt))
    tm_squared=[[sum(tm[i][l]*tm[l][j] for l in range(2)) for j in range(2)] for i in range(2)]
    if any(value for row in tm_squared for value in row):
        raise ArithmeticError('transverse closure block is not nilpotent')
    transverse_closure={'scope':'Both degenerate axial transverse directions, every physical K',
                        'shape_coefficient':str(qt),'matrix':[[str(v) for v in row] for row in tm],
                        'matrix_squared':[[str(v) for v in row] for row in tm_squared],
                        'eigenvalues':['0','0'],'upper_half_plane_modes':False,
                        'reason':'Pi_ex[x,z]=0 and DeltaF[x]=0 in Eqs26-27; transverse block is exactly nilpotent. This does not exclude full-angular transverse instability.'}
    validated=validate_linear_brackets(rows)
    cm=max(c['independent_matrix_peak_error'] for c in closures)
    cp=max(c['independent_matrix_characteristic_maxerr'] for c in closures)
    if cm>1e-12 or cp>1e-12:raise RuntimeError('closure formula/matrix mismatch')
    d={'scope':'Fixed 17 physical K points, only sign-changing imaginary roots in predeclared gamma range; no all-root census. Linear-completion roots have Arb sign/IVT brackets; exact-ME completion roots remain numerical references.',
       'physical_K_grid':[str(k) for k in ks], 'shift':'k_reduced=K_physical−2/3; omega unchanged because signed mass is zero',
       'normalization':'N_plus=N_minus=1, mu0=sqrt(2) G_F N; c=1, hbar=1',
       'gamma_search':{'minimum':1e-7,'maximum':3.,'log_nodes':160},
       'beta_maximum_entropy':BETA,'chi_polynomial_Minerbo':str(CHI),
       'moments_numerical':moment_reference,'continuous_modes':rows,'closure_modes':closures,
       'validated_linear_root_brackets':validated,'quantum_closure_transverse':transverse_closure,
       'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
       'sources':SOURCES,'python':sys.version,'platform':platform.platform(),
       'development_failure':'Initial independent matrix replay caught anti-species flux-entry sign errors before any output result; source/log preserved and fixed from Eqs25,27. A subsequent check identified harmless numerical splitting at the defective K=0 zero eigenvalue; there the independent characteristic polynomial, rather than a floating growth decision, is checked. Low-gamma quadrature integrand regularized without changing cases/search range.',
       'method':'SciPy quad + brentq; independently cross-checked direct 256/512-point Gauss-Legendre original integrals; closure exact gamma² independently compared with 4x4 eigenvalues.'}
    out=ROOT/'results/counterflux-modes.json';out.write_text(json.dumps(d,indent=2)+'\n')
    print(json.dumps({'output':str(out),'beta':BETA,'moments':moment_reference,
          'root_counts':[{'completion':kind,'sector':name,'count':sum(len(b['roots']) for r in rows if r['completion']==kind for b in r['branches'] if b['sector']==name)} for kind in ('linear','maximum_entropy') for name in ('longitudinal','transverse')],
          'max_closure_matrix_error':cm,'max_closure_characteristic_error':cp,
          'validated_linear_root_brackets':len(validated),
          'max_direct512_residual':max(q['independent_512_residual'] for r in rows for b in r['branches'] for q in b['roots'])},indent=2))

if __name__=='__main__':main()
