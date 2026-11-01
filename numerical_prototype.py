"""Exact rational Hausdorff/Stieltjes feasibility for axial neutrino roots.

Classical moment theory is not claimed new. Numerical optimization proposes;
Fraction arithmetic checks attainability or separating certificates.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import combinations
import math
import numpy as np
from scipy.optimize import minimize

@dataclass(frozen=True)
class QC:
    re:F=F(0)
    im:F=F(0)
    def __post_init__(self):
        object.__setattr__(self,'re',F(self.re));object.__setattr__(self,'im',F(self.im))
    @staticmethod
    def cast(x): return x if isinstance(x,QC) else QC(F(x))
    def __add__(self,o):
        o=self.cast(o);return QC(self.re+o.re,self.im+o.im)
    __radd__=__add__
    def __neg__(self): return QC(-self.re,-self.im)
    def __sub__(self,o): return self+-self.cast(o)
    def __rsub__(self,o): return self.cast(o)+-self
    def __mul__(self,o):
        o=self.cast(o);return QC(self.re*o.re-self.im*o.im,self.re*o.im+self.im*o.re)
    __rmul__=__mul__
    def __truediv__(self,o):
        o=self.cast(o);d=o.re**2+o.im**2
        if d==0: raise ZeroDivisionError
        return self*QC(o.re/d,-o.im/d)
    def __rtruediv__(self,o): return self.cast(o)/self
    def __complex__(self): return complex(float(self.re),float(self.im))

def la(a,b): return tuple(x+y for x,y in zip(a,b))
def ls(a,c): return tuple(c*x for x in a)

def hausdorff_affine(N,flux,z):
    """Return H±=(C,A,B); s=x+iy, matrices stored (00,01,11)."""
    N,flux=F(N),F(flux)
    if N<0 or abs(flux)>N or z.im<=0: raise ValueError('invalid positive-moment data')
    t0=(F(0),F(0),-1/z.im)
    t1=la(ls(t0,z.re),(F(0),F(-1),F(0)))
    t2=la((N,F(0),F(0)),la(ls(t1,2*z.re),ls(t0,-(z.re*z.re+z.im*z.im))))
    t3=la((flux,F(0),F(0)),la(ls(t2,2*z.re),ls(t1,-(z.re*z.re+z.im*z.im))))
    out=[]
    for sign in (1,-1):
        entries=(la(t0,ls(t1,sign)),la(t1,ls(t2,sign)),la(t2,ls(t3,sign)))
        out.append(tuple(tuple(e[j] for e in entries) for j in range(3)))
    return out

def eval_aff(h,x,y): return tuple(c+a*x+b*y for c,a,b in zip(*h))
def psd(q): return q[0]>=0 and q[2]>=0 and q[0]*q[2]-q[1]*q[1]>=0

def reweighted(N,flux,z,s):
    t0=-s.im/z.im;t1=z.re*t0-s.re
    t2=F(N)+2*z.re*t1-(z.re*z.re+z.im*z.im)*t0
    t3=F(flux)+2*z.re*t2-(z.re*z.re+z.im*z.im)*t1
    return t0,t1,t2,t3

def attaining(N,flux,z,s):
    """Exact algebraic quadrature specification, or None if not feasible."""
    hs=hausdorff_affine(N,flux,z)
    if not all(psd(eval_aff(h,s.re,s.im)) for h in hs): return None
    t0,t1,t2,t3=reweighted(N,flux,z,s)
    if N==0:
        return {'kind':'zero'} if s==QC() else None
    if t0<=0: return None
    det=t0*t2-t1*t1
    if det==0:
        u=t1/t0
        if -1<=u<=1 and t2==t0*u*u and t3==t0*u*u*u:
            return {'kind':'one','u':u,'nu_mass':t0}
        return None
    if det<0: return None
    # monic p2(u)=u²+a*u+b; exact Gaussian quadrature of nu.
    b=(-t2*t2+t1*t3)/det
    a=(-t0*t3+t1*t2)/det
    disc=a*a-4*b;center=-a/2;mean=t1/t0
    if not (disc>0 and -1<=center<=1 and disc/4<=(1-abs(center))**2
            and (mean-center)**2<=disc/4):
        return None
    return {'kind':'two','p2_a':a,'p2_b':b,'nu_t0':t0,'nu_t1':t1}

def root_target(N1,F1,N2,F2,w,k,sector):
    g0=F(N1)-F(N2);g1=F(F1)-F(F2);k=F(k)
    if k<=0 or w.im<=0: raise ValueError('requires k>0 and Im omega>0')
    if sector=='longitudinal':
        A=w*w-k*k+w*g0-k*g1
        B=k*k+w*g0+k*g1+g0*g0
    elif sector=='transverse':
        A=k*k-w*w;B=2*k*k-w*g0-k*g1
    else: raise ValueError('unknown sector')
    if A==QC():
        if B==QC(): raise ArithmeticError('unexpected upper-half-plane identity')
        return None
    return k*B/A # target for s1-s2, not J=(s1-s2)/k

def shifted_scaled(h,shift,scale):
    C,A,B=h
    return (tuple(c-a*shift.re-b*shift.im for c,a,b in zip(C,A,B)),
            tuple(a*scale for a in A),tuple(b*scale for b in B))

def determinant3(m):
    return (m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1])
           -m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0])
           +m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]))

def dual_from_directions(hs,xy):
    """Rational PSD rank-one directions plus exact nonnegative weights."""
    pool=[]
    for j,h in enumerate(hs):
        m=eval_aff(h,*xy)
        _,vs=np.linalg.eigh(np.array([[float(m[0]),float(m[1])],[float(m[1]),float(m[2])]]))
        for col in range(2):
            v=tuple(F(float(x)).limit_denominator(10**6) for x in vs[:,col])
            q=(v[0]**2,2*v[0]*v[1],v[1]**2)
            vals=tuple(sum(qi*hi for qi,hi in zip(q,hpart)) for hpart in h)
            pool.append((j,v,vals))
    # Valid normalized Stieltjes bounds: |Re(s)/scale|<=1, -1<=Im(s)/scale<=0.
    bounds=((F(1),F(1),F(0)),(F(1),F(-1),F(0)),
            (F(1),F(0),F(1)),(F(0),F(0),F(-1)))
    for j,h in enumerate(bounds): pool.append((len(hs)+j,None,h))
    pool.sort(key=lambda r:float(r[2][0]+r[2][1]*xy[0]+r[2][2]*xy[1]))
    for n in (1,2,3):
        for rows in combinations(pool,n):
            if n==1:
                lambdas=(F(1),) if rows[0][2][1:]==(0,0) else None
            elif n==2:
                va,vb=rows[0][2],rows[1][2]
                dim=1 if va[1]!=vb[1] else 2
                if va[dim]==vb[dim]: continue
                t=-vb[dim]/(va[dim]-vb[dim]);lambdas=(t,1-t)
            else:
                m=[[F(1)]*3,[r[2][1] for r in rows],[r[2][2] for r in rows]]
                d=determinant3(m)
                if d==0: continue
                lambdas=[]
                for j in range(3):
                    a=[row[:] for row in m]
                    for i in range(3): a[i][j]=F(i==0)
                    lambdas.append(determinant3(a)/d)
            if lambdas is None or min(lambdas)<0: continue
            total=tuple(sum(l*r[2][j] for l,r in zip(lambdas,rows)) for j in range(3))
            if total[1:]==(0,0) and total[0]<0:
                return {'terms':[{'constraint':r[0],'vector':r[1],'weight':l}
                                 for r,l in zip(rows,lambdas)],'negative_constant':total[0]}
    return None

def classify(N1,F1,N2,F2,w,k,sector='longitudinal'):
    N1,F1,N2,F2,k=map(F,(N1,F1,N2,F2,k))
    for N,f in ((N1,F1),(N2,F2)):
        if N<0 or abs(f)>N: raise ValueError('nonrealizable species moments')
    target=root_target(N1,F1,N2,F2,w,k,sector)
    if target is None: return {'status':'impossible','proof_kind':'zero_A_nonzero_B'}
    z=w/k
    fixed1=QC() if N1==0 else (N1/(z-F1/N1) if abs(F1)==N1 else None)
    fixed2=QC() if N2==0 else (N2/(z-F2/N2) if abs(F2)==N2 else None)
    if fixed1 is not None or fixed2 is not None:
        s1=fixed1 if fixed1 is not None else fixed2+target
        s2=s1-target
        a=attaining(N1,F1,z,s1);b=attaining(N2,F2,z,s2)
        if a is not None and b is not None:
            return {'status':'possible','s1':s1,'s2':s2,'witnesses':[a,b]}
        violations=[]
        for species,(N,f,st) in enumerate(((N1,F1,s1),(N2,F2,s2))):
            for j,h in enumerate(hausdorff_affine(N,f,z)):
                if not psd(eval_aff(h,st.re,st.im)): violations.append((species,j))
        if violations:
            return {'status':'impossible','proof_kind':'fixed_species','s1':s1,'s2':s2,'violations':violations}
        return {'status':'unresolved','reason':'degenerate witness reconstruction'}
    scale=N1/z.im
    hs=[shifted_scaled(h,QC(),scale) for h in hausdorff_affine(N1,F1,z)]
    hs += [shifted_scaled(h,target,scale) for h in hausdorff_affine(N2,F2,z)]
    # Positive scaling preserves exact LMI feasibility and improves conditioning.
    hs=[tuple(tuple(x/max(F(1),max(abs(t) for part in h for t in part)) for x in part)
              for part in h) for h in hs]
    arrays=np.array([[[[float(a[0]),float(a[1])],[float(a[1]),float(a[2])]] for a in h] for h in hs])
    def cons(v,gradient=False):
        vals=[];grads=[]
        for C,A,B in arrays:
            ev,evec=np.linalg.eigh(C+v[0]*A+v[1]*B);q=evec[:,0]
            vals.append(ev[0]-v[2]);grads.append([q@A@q,q@B@q,-1.])
        return np.asarray(grads if gradient else vals)
    cold=complex(N1/(z-F1/N1))/float(scale)
    start=np.array([np.clip(cold.real,-1,1),np.clip(cold.imag,-1,0),-1.])
    result=minimize(lambda v:-v[2],start,jac=lambda v:np.array([0.,0.,-1.]),
                    constraints=[{'type':'ineq','fun':cons,'jac':lambda v:cons(v,True)}],
                    bounds=[(-1,1),(-1,0),(-10,10)],method='SLSQP',
                    options={'ftol':1e-11,'maxiter':250})
    xy=tuple(F(float(v)) for v in result.x[:2])
    diagnostics={'solver_success':bool(result.success),'solver_margin':float(result.x[2]),
                 'solver_iterations':int(result.nit)}
    for denom in (10**6,10**9,10**12,None):
        point=tuple(v.limit_denominator(denom) if denom else v for v in xy)
        s1=QC(scale*point[0],scale*point[1]);s2=s1-target
        a=attaining(N1,F1,z,s1);b=attaining(N2,F2,z,s2)
        if a is not None and b is not None:
            return {'status':'possible','s1':s1,'s2':s2,'witnesses':[a,b],**diagnostics}
    dual=dual_from_directions(hs,xy)
    if dual is not None:
        return {'status':'impossible','proof_kind':'rational_dual','dual':dual,**diagnostics}
    return {'status':'unresolved',**diagnostics}
