"""Stdlib-only independent supplied-second-moment frequency proof checker.

No production geometry or frequency implementation is imported. A transform
is encoded as [p,q] meaning p+i*sqrt(gamma2)*q at z=(wr+i*sqrt(gamma2))/k.
Each coordinate may be rational or an independent_verifier.R serialization.
See independent-m2-locus.md for the constructive moment-completion proof.
"""
from dataclasses import dataclass
from fractions import Fraction as Q
from itertools import combinations

from independent_verifier import InvalidCertificate, R, as_real, rational, require


@dataclass(frozen=True)
class CQ:
    p: object
    q: object
    g: Q

    def __post_init__(self):
        object.__setattr__(self, 'p', as_real(self.p))
        object.__setattr__(self, 'q', as_real(self.q))
        object.__setattr__(self, 'g', rational(self.g))
        require(self.g > 0, 'positive imaginary square required')

    def coerce(self, other):
        if isinstance(other, CQ):
            require(other.g == self.g, 'inconsistent imaginary basis')
            return other
        return CQ(other, 0, self.g)

    def __add__(self, other):
        o = self.coerce(other)
        return CQ(self.p+o.p, self.q+o.q, self.g)
    __radd__ = __add__

    def __neg__(self): return CQ(-self.p, -self.q, self.g)
    def __sub__(self, other): return self+-self.coerce(other)
    def __rsub__(self, other): return self.coerce(other)+-self

    def __mul__(self, other):
        o = self.coerce(other)
        return CQ(self.p*o.p-self.g*self.q*o.q,
                  self.p*o.q+self.q*o.p, self.g)
    __rmul__ = __mul__

    def __truediv__(self, other):
        o = self.coerce(other)
        norm = o.p*o.p+self.g*o.q*o.q
        require(norm > 0, 'zero complex denominator')
        return self*CQ(o.p/norm, -o.q/norm, self.g)

    def __rtruediv__(self, other): return self.coerce(other)/self
    def zero(self): return self.p == 0 and self.q == 0


def parse(data):
    keys = ('N1', 'F1', 'Q1', 'N2', 'F2', 'Q2', 'wr', 'gamma2', 'k')
    n1, f1, q1, n2, f2, q2, wr, g, k = (rational(data[s]) for s in keys)
    require(k > 0 and g > 0, 'requires k>0 and gamma2>0')
    for n, f, q in ((n1,f1,q1), (n2,f2,q2)):
        require(n >= 0, 'negative species mass')
        if n == 0:
            require(f == 0 and q == 0, 'zero mass requires zero moments')
        else:
            require(q <= n and n*q >= f*f, 'nonrealizable species moments')
    sector = data.get('sector', 'longitudinal')
    require(sector in ('longitudinal', 'transverse'), 'unknown sector')
    return (n1,f1,q1), (n2,f2,q2), wr, g, k, sector


def read_s(value, g):
    require(isinstance(value, (list,tuple)) and len(value) == 2,
            'transform requires [p,q]')
    return CQ(value[0], value[1], g)


def moments(species, wr, g, k, s):
    n, f, q = species
    x, rho = wr/k, (wr*wr+g)/(k*k)
    t0 = -k*s.q
    t1 = x*t0-s.p
    t2 = n+2*x*t1-rho*t0
    t3 = f+2*x*t2-rho*t1
    t4 = q+2*x*t3-rho*t2
    return t0,t1,t2,t3,t4


def support_matrices(species, wr, g, k, s):
    t = moments(species,wr,g,k,s)
    return ([[t[i+j] for j in range(3)] for i in range(3)],
            [[t[i+j]-t[i+j+2] for j in range(2)] for i in range(2)])


def determinant(a):
    if len(a) == 1: return a[0][0]
    if len(a) == 2: return a[0][0]*a[1][1]-a[0][1]*a[1][0]
    require(len(a) == 3, 'unsupported determinant size')
    return sum(((-1)**j)*a[0][j]*determinant(
        [[a[i][h] for h in range(3) if h != j] for i in (1,2)])
        for j in range(3))


def psd(matrix):
    return all(determinant([[matrix[i][j] for j in indices] for i in indices]) >= 0
               for size in range(1,len(matrix)+1)
               for indices in combinations(range(len(matrix)),size))


def feasible_s(species, wr, g, k, s):
    return all(psd(a) for a in support_matrices(species,wr,g,k,s))


def original_factor(sp1, sp2, wr, g, k, sector, s1, s2):
    """Evaluate the original, unmultiplied factor through direct recurrences."""
    w = CQ(wr,1,g)
    i0 = (s1-s2)/k
    i1 = (w*i0-(sp1[0]-sp2[0]))/k
    i2 = (w*i1-(sp1[1]-sp2[1]))/k
    return ((1+i0)*(-1+i2)-i1*i1 if sector == 'longitudinal'
            else -1+(i0-i2)/2)


def affine_target(sp1, sp2, wr, g, k, sector):
    """Infer the root's affine transform equation by direct factor evaluations.

    The factor is affine in s1-s2 by polynomial cancellation. This evaluation
    does not import or copy the production sector-coefficient formulas.
    """
    zero, one = CQ(0,0,g), CQ(1,0,g)
    constant = original_factor(sp1,sp2,wr,g,k,sector,zero,zero)
    slope = original_factor(sp1,sp2,wr,g,k,sector,one,zero)-constant
    if slope.zero():
        return None, constant.zero()
    return -constant/slope, False


def unique_transform(species, wr, g, k):
    n,f,q = species
    z = CQ(wr/k,1/k,g)
    if n == 0: return CQ(0,0,g)
    if n*q == f*f: return n/(z-f/n)
    if q == n:
        return ((n-f)/2)/(z+1)+((n+f)/2)/(z-1)
    return None


def determinant_ellipses(sp1, sp2, wr, g, k, target):
    """Recover all four necessary ellipses via six determinant evaluations."""
    ellipses = []
    for species, shift in ((sp1,CQ(0,0,g)), (sp2,target)):
        for which in (0,1):
            def f(p,q):
                s = CQ(p,q,g)-shift
                return determinant(support_matrices(species,wr,g,k,s)[which])
            f0 = f(0,0)
            lp, lq = (f(1,0)-f(-1,0))/2, (f(0,1)-f(0,-1))/2
            ap, aq = (f(1,0)+f(-1,0))/2-f0, (f(0,1)+f(0,-1))/2-f0
            mixed = f(1,1)-f0-lp-lq-ap-aq
            require(mixed == 0 and ap < 0 and aq == g*ap,
                    'determinant is not the expected ellipse')
            center = (-lp/(2*ap), -lq/(2*aq))
            radius2 = center[0]**2+g*center[1]**2-f0/ap
            require(radius2 >= 0, 'negative necessary ellipse radius')
            ellipses.append((center,radius2))
    return ellipses


def candidates(ellipses, g):
    for center,_ in ellipses: yield center
    for (c,r), (other,s) in combinations(ellipses,2):
        dx,dy = other[0]-c[0],other[1]-c[1]
        d2 = dx*dx+g*dy*dy
        if d2 == 0: continue
        a = (r-s+d2)/(2*d2)
        t2 = (r-a*a*d2)/(g*d2)
        if t2 < 0: continue
        t = R(0,1,t2)
        p,q = c[0]+a*dx,c[1]+a*dy
        yield p-g*dy*t,q+dx*t
        yield p+g*dy*t,q-dx*t


def in_ellipse(point, ellipse, g):
    c,r = ellipse
    return (point[0]-c[0])**2+g*(point[1]-c[1])**2 <= r


def _verify(data, result):
    sp1,sp2,wr,g,k,sector = parse(data)
    status = result['status']
    require(status in ('possible','impossible'), 'unresolved is not a certificate')
    if status == 'possible':
        s1,s2 = read_s(result['s1'],g),read_s(result['s2'],g)
        require(original_factor(sp1,sp2,wr,g,k,sector,s1,s2).zero(),
                'original dispersion factor does not vanish')
        require(feasible_s(sp1,wr,g,k,s1) and feasible_s(sp2,wr,g,k,s2),
                'a species fails full degree-four Hausdorff positivity')
        return True
    require(result.get('proof_kind') == 'm2_four_disks', 'unknown impossibility proof')
    target, identically_zero = affine_target(sp1,sp2,wr,g,k,sector)
    if target is None:
        require(not identically_zero, 'factor vanishes for all admissible measures')
        return True
    u1,u2 = unique_transform(sp1,wr,g,k),unique_transform(sp2,wr,g,k)
    if u1 is not None or u2 is not None:
        s1 = u1 if u1 is not None else u2+target
        s2 = s1-target
        require(not (feasible_s(sp1,wr,g,k,s1) and feasible_s(sp2,wr,g,k,s2)),
                'unique-moment completion supplies a feasible root')
        return True
    ellipses = determinant_ellipses(sp1,sp2,wr,g,k,target)
    for point in candidates(ellipses,g):
        require(not all(in_ellipse(point,e,g) for e in ellipses),
                'claimed empty intersection contains an exact candidate')
    return True


def verify(data, result):
    """Return True or raise InvalidCertificate; explicit checks survive -O."""
    try:
        return _verify(data,result)
    except InvalidCertificate:
        raise
    except (KeyError,TypeError,ValueError,ZeroDivisionError,IndexError) as exc:
        raise InvalidCertificate('malformed M2 certificate: '+str(exc)) from exc


def verify_m1_possible(data, result):
    """Check an M0/M1 positive witness at rational gamma2, with no Q claim.

    The input may omit Q1/Q2. Supplied Q values are deliberately ignored:
    this checks membership of the larger information class only. Transforms
    retain the same [p,q] imaginary-basis representation as the M2 checker.
    """
    try:
        # Endpoint second moments are always realizable for admissible N,F;
        # they only let the common parser validate the remaining input.
        parsed = dict(data,Q1=data['N1'],Q2=data['N2'])
        sp1,sp2,wr,g,k,sector = parse(parsed)
        require(result['status'] == 'possible', 'M1 API accepts positive proofs only')
        s1,s2 = read_s(result['s1'],g),read_s(result['s2'],g)
        require(original_factor(sp1,sp2,wr,g,k,sector,s1,s2).zero(),
                'original dispersion factor does not vanish')
        for sp,s in ((sp1,s1),(sp2,s2)):
            t = moments(sp,wr,g,k,s)
            for sign in (-1,1):
                h = [[t[i+j]+sign*t[i+j+1] for j in range(2)] for i in range(2)]
                require(psd(h), 'M1 species fails degree-three Hausdorff positivity')
        return True
    except InvalidCertificate:
        raise
    except (KeyError,TypeError,ValueError,ZeroDivisionError,IndexError) as exc:
        raise InvalidCertificate('malformed M1 positive certificate: '+str(exc)) from exc
