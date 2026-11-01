"""Exact four-disk geometry over real quadratic fields, no numerical solver."""
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import combinations
from math import isqrt,sqrt

def _sgn(x): return (x>0)-(x<0)

@dataclass(frozen=True,eq=False)
class Quad:
    a:F=F(0)
    b:F=F(0)
    d:F=F(0)
    def __post_init__(self):
        a,b,d=F(self.a),F(self.b),F(self.d)
        if d<0: raise ValueError('negative real radical')
        n,q=isqrt(d.numerator),isqrt(d.denominator)
        if b==0 or d==0: b,d=F(0),F(0)
        elif n*n==d.numerator and q*q==d.denominator: a,b,d=a+b*F(n,q),F(0),F(0)
        object.__setattr__(self,'a',a);object.__setattr__(self,'b',b);object.__setattr__(self,'d',d)
    def pair(self,o):
        if not isinstance(o,Quad): o=Quad(F(o))
        if self.d!=o.d and self.b and o.b: raise ValueError('mixed quadratic fields')
        return o,self.d if self.b else o.d
    def __add__(self,o):
        try:o,d=self.pair(o)
        except TypeError:return NotImplemented
        return Quad(self.a+o.a,self.b+o.b,d)
    __radd__=__add__
    def __neg__(self):return Quad(-self.a,-self.b,self.d)
    def __sub__(self,o):return self+-o if isinstance(o,Quad) else self+-F(o)
    def __rsub__(self,o):return -self+o
    def __mul__(self,o):
        try:o,d=self.pair(o)
        except TypeError:return NotImplemented
        return Quad(self.a*o.a+self.b*o.b*d,self.a*o.b+self.b*o.a,d)
    __rmul__=__mul__
    def __truediv__(self,o):
        o,d=self.pair(o);den=o.a*o.a-o.b*o.b*d
        if den==0:raise ZeroDivisionError
        return self*Quad(o.a/den,-o.b/den,d)
    def __rtruediv__(self,o):return Quad(F(o))/self
    def __pow__(self,n):
        if not isinstance(n,int) or n<0: raise ValueError('requires nonnegative integer power')
        out=Quad(1)
        for _ in range(n):out=out*self
        return out
    def sign(self):
        sa,sb=_sgn(self.a),_sgn(self.b)
        if not sb:return sa
        if not sa:return sb
        if sa==sb:return sa
        c=_sgn(self.a*self.a-self.b*self.b*self.d)
        return c if sa>0 else -c
    def __eq__(self,o):
        try:return (self-o).sign()==0
        except (TypeError,ValueError):return False
    def __lt__(self,o):return (self-o).sign()<0
    def __le__(self,o):return (self-o).sign()<=0
    def __gt__(self,o):return (self-o).sign()>0
    def __ge__(self,o):return (self-o).sign()>=0
    def __abs__(self):return self if self.sign()>=0 else -self
    def __float__(self):return float(self.a)+float(self.b)*sqrt(float(self.d))

def disks(N,M,z):
    N,M=F(N),F(M);x,y=z.re,z.im
    if N<0 or abs(M)>N or y<=0:raise ValueError('invalid positive moments')
    out=[]
    for sign in (1,-1):
        a=(x+sign)**2+y*y
        lx=N*(2*x+sign)-M
        ly=-((1+sign*x)*(N+sign*M)+2*N*y*y)/y
        out.append(((lx/(2*a),ly/(2*a)),(N+sign*M)**2/(4*y*y*a)))
    return out

def inside(point,disk):
    (x,y),r2=disk
    return (point[0]-x)**2+(point[1]-y)**2<=r2

def candidates(ds):
    for i,(c,r2) in enumerate(ds):yield ('center',i),tuple(Quad(t) for t in c)
    for i,j in combinations(range(len(ds)),2):
        (c,r2),(e,s2)=ds[i],ds[j];dx=e[0]-c[0];dy=e[1]-c[1];D=dx*dx+dy*dy
        if D==0:continue
        p=(r2-s2+D)/(2*D);d=r2/D-p*p
        if d<0:continue
        for sign in ((1,) if d==0 else (1,-1)):
            yield ('pair',i,j,sign),(Quad(c[0]+p*dx,-sign*dy,d),Quad(c[1]+p*dy,sign*dx,d))

def feasible_candidates(ds):
    return [(name,p) for name,p in candidates(ds) if all(inside(p,d) for d in ds)]
