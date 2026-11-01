"""Exact axial root compatibility with supplied number moments N, M and Q.

Classical degree-four Hausdorff theory gives two ellipses per species.
Coordinates S=p+i*sqrt(gamma2)*q use a rational metric, so rational squared
rates need no rounding. Positive results mean existence of SOME completion.
General completions can require three atoms; no two-atom claim is made.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import combinations
import argparse
import json
from geometry import Quad


@dataclass(frozen=True)
class ComplexBasis:
    p: object
    q: object
    g: F

    def __post_init__(self):
        for key in ('p', 'q'):
            v = getattr(self, key)
            object.__setattr__(self, key, v if isinstance(v, Quad) else F(v))
        object.__setattr__(self, 'g', F(self.g))
        if self.g <= 0:
            raise ValueError('gamma2 must be positive')

    def cast(self, value):
        if isinstance(value, ComplexBasis):
            if value.g != self.g:
                raise ValueError('different imaginary bases')
            return value
        return ComplexBasis(value, 0, self.g)

    def __add__(self, other):
        o = self.cast(other)
        return ComplexBasis(self.p+o.p, self.q+o.q, self.g)
    __radd__ = __add__
    def __neg__(self): return ComplexBasis(-self.p, -self.q, self.g)
    def __sub__(self, other): return self+-self.cast(other)
    def __rsub__(self, other): return self.cast(other)+-self
    def __mul__(self, other):
        o = self.cast(other)
        return ComplexBasis(self.p*o.p-self.g*self.q*o.q,
                            self.p*o.q+self.q*o.p, self.g)
    __rmul__ = __mul__
    def __truediv__(self, other):
        o = self.cast(other)
        norm = o.p*o.p+self.g*o.q*o.q
        if norm == 0:
            raise ZeroDivisionError('zero complex denominator')
        return self*ComplexBasis(o.p/norm, -o.q/norm, self.g)
    def __rtruediv__(self, other): return self.cast(other)/self
    def zero(self): return self.p == 0 and self.q == 0


def parse(data):
    names = ('N1','F1','Q1','N2','F2','Q2','wr','gamma2','k')
    n1,m1,q1,n2,m2,q2,wr,g,k = (F(data[s]) for s in names)
    if g <= 0 or k <= 0:
        raise ValueError('requires gamma2>0 and reduced k>0')
    for n,m,q in ((n1,m1,q1), (n2,m2,q2)):
        if n < 0 or q < 0 or q > n or m*m > n*q:
            raise ValueError('nonrealizable positive number moments')
    sector = data.get('sector','longitudinal')
    if sector not in ('longitudinal','transverse'):
        raise ValueError('unknown sector')
    return (n1,m1,q1),(n2,m2,q2),wr,g,k,sector


def sector_target(sp1, sp2, wr, g, k, sector):
    n1,m1,_ = sp1
    n2,m2,_ = sp2
    w = ComplexBasis(wr,1,g)
    g0,g1 = n1-n2,m1-m2
    if sector == 'longitudinal':
        a = w*w-k*k+w*g0-k*g1
        b = k*k+w*g0+k*g1+g0*g0
    else:
        a = k*k-w*w
        b = 2*k*k-w*g0-k*g1
    if a.zero():
        if b.zero():
            raise ArithmeticError('unexpected UHP identically-zero factor')
        return None
    return k*b/a


def species_ellipses(species, wr, g, k):
    n,m,q = species
    z = ComplexBasis(wr/k,1/k,g)
    delta = n*q-m*m
    if n == 0:
        return [((F(0),F(0)),F(0))]*2
    if delta == 0:
        s = n/(z-m/n)
        return [((s.p,s.q),F(0))]*2
    if q == n:
        s = ((n+m)/2)/(z-1)+((n-m)/2)/(z+1)
        return [((s.p,s.q),F(0))]*2
    x,rho = wr/k,(wr*wr+g)/(k*k)
    d = q-2*m*x+n*rho
    if d <= 0:
        raise ArithmeticError('positive mass and UHP frequency require D>0')
    gram = ((n*(n*x-m)/d, -(2*n*n/k+delta*k/g)/(2*d)),
            delta*delta*k*k/(4*g*d*d))
    remainder = n-q
    den = z*z-1
    c = (n*z+m+ComplexBasis(0,remainder*k/(2*g),g))/den
    support = ((c.p,c.q),remainder*remainder*k*k /
               (4*g*(den.p*den.p+g*den.q*den.q)))
    return [gram,support]


def ellipse_candidates(ellipses, g):
    for index,(center,_) in enumerate(ellipses):
        yield ('center',index),tuple(Quad(t) for t in center)
    for i,j in combinations(range(len(ellipses)),2):
        (c,r2),(d,s2) = ellipses[i],ellipses[j]
        dx,dy = d[0]-c[0],d[1]-c[1]
        distance = dx*dx+g*dy*dy
        if distance == 0:
            continue
        along = (r2-s2+distance)/(2*distance)
        radical = (r2/distance-along*along)/g
        if radical < 0:
            continue
        for sign in ((1,) if radical == 0 else (1,-1)):
            yield ('pair',i,j,sign),(
                Quad(c[0]+along*dx,-sign*g*dy,radical),
                Quad(c[1]+along*dy,sign*dx,radical))


def contains(point, ellipse, g):
    center,radius2 = ellipse
    return (point[0]-center[0])**2+g*(point[1]-center[1])**2 <= radius2


def classify(data):
    sp1,sp2,wr,g,k,sector = parse(data)
    target = sector_target(sp1,sp2,wr,g,k,sector)
    if target is None:
        return dict(status='impossible',proof_kind='m2_four_disks')
    e1 = species_ellipses(sp1,wr,g,k)
    e2 = [((c[0]+target.p,c[1]+target.q),r)
          for c,r in species_ellipses(sp2,wr,g,k)]
    ellipses = e1+e2
    # Prefer algebraic boundary candidates; centers cover disk-containment cases.
    points = list(ellipse_candidates(ellipses,g))
    points.sort(key=lambda item:item[0][0] == 'center')
    for tag,point in points:
        if all(contains(point,e,g) for e in ellipses):
            s1 = ComplexBasis(*point,g)
            s2 = s1-target
            return dict(status='possible',s1=[s1.p,s1.q],s2=[s2.p,s2.q],
                        candidate=list(tag),proof_kind='degree4_hausdorff')
    return dict(status='impossible',proof_kind='m2_four_disks')


def jsonable(value):
    if isinstance(value,F): return str(value)
    if isinstance(value,Quad):
        return {k:jsonable(getattr(value,k)) for k in ('a','b','d')}
    if isinstance(value,dict): return {k:jsonable(v) for k,v in value.items()}
    if isinstance(value,(tuple,list)): return [jsonable(v) for v in value]
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input',help='JSON file with N1,F1,Q1,N2,F2,Q2,wr,gamma2,k,sector')
    args = parser.parse_args()
    data = json.loads(__import__('pathlib').Path(args.input).read_text())
    print(json.dumps(jsonable(classify(data)),indent=2))


if __name__ == '__main__': main()
