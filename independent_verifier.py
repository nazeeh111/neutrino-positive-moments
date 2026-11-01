"""Independent exact proof checker. Standard library only; no solver/production import.

Input: data={N1,F1,N2,F2,w,k,sector}; w/s values may be objects with
Fraction .re/.im, {'re':..., 'im':...}, or length-two pairs. Fractions may
also be serialized as integer strings or 'numerator/denominator'.
verify(data, result) raises InvalidCertificate on invalid input/proof, and
returns True only for certified possible/impossible. An unresolved proposal
is deliberately not a certificate. Solver diagnostics have no authority.
"""
from dataclasses import dataclass
from fractions import Fraction as Q
from math import isqrt
from itertools import combinations


class InvalidCertificate(ValueError):
    pass


def require(condition, message):
    if not condition:
        raise InvalidCertificate(message)


def rational(value):
    require(isinstance(value, (int, str, Q)) and not isinstance(value, bool),
            'expected exact rational')
    return Q(value)


def rational_sqrt(value):
    if value < 0:
        return None
    p,q = isqrt(value.numerator),isqrt(value.denominator)
    return Q(p,q) if p*p == value.numerator and q*q == value.denominator else None


@dataclass(frozen=True, eq=False)
class R:
    """Independent ordered real quadratic arithmetic; no floating decisions.

    Equivalent radical presentations are aligned by an exact rational square
    ratio, avoiding reliance on square-free factorization or production code.
    """
    a: Q = Q(0)
    b: Q = Q(0)
    d: Q = Q(0)

    def __post_init__(self):
        a,b,d = rational(self.a),rational(self.b),rational(self.d)
        require(d >= 0, 'negative real radical')
        root = rational_sqrt(d)
        if b == 0 or d == 0:
            b,d = Q(0),Q(0)
        elif root is not None:
            a,b,d = a+b*root,Q(0),Q(0)
        object.__setattr__(self,'a',a)
        object.__setattr__(self,'b',b)
        object.__setattr__(self,'d',d)

    def aligned(self, other):
        o = as_real(other)
        if not isinstance(o,R):
            o = R(o)
        if self.b and o.b and self.d != o.d:
            ratio = rational_sqrt(o.d/self.d)
            require(ratio is not None, 'independent quadratic fields in one expression')
            o = R(o.a,o.b*ratio,self.d)
        return o,self.d if self.b else o.d

    def __add__(self,other):
        o,d = self.aligned(other)
        return R(self.a+o.a,self.b+o.b,d)
    __radd__ = __add__

    def __neg__(self):
        return R(-self.a,-self.b,self.d)

    def __sub__(self,other):
        return self + -as_real(other)

    def __rsub__(self,other):
        return -self + other

    def __mul__(self,other):
        o,d = self.aligned(other)
        return R(self.a*o.a+self.b*o.b*d,self.a*o.b+self.b*o.a,d)
    __rmul__ = __mul__

    def __truediv__(self,other):
        o,d = self.aligned(other)
        norm = o.a*o.a-o.b*o.b*d
        require(norm != 0, 'zero quadratic denominator')
        return self * R(o.a/norm,-o.b/norm,d)

    def __rtruediv__(self,other):
        return R(rational(other))/self

    def __pow__(self,power):
        require(type(power) is int and power >= 0, 'invalid quadratic exponent')
        out = R(1)
        for _ in range(power):
            out = out*self
        return out

    def sign(self):
        if self.b == 0:
            return (self.a > 0)-(self.a < 0)
        if self.a == 0:
            return (self.b > 0)-(self.b < 0)
        if self.a > 0 and self.b > 0:
            return 1
        if self.a < 0 and self.b < 0:
            return -1
        difference = self.a*self.a-self.b*self.b*self.d
        return ((difference > 0)-(difference < 0))*(1 if self.a > 0 else -1)

    def __eq__(self,other):
        try:
            return (self-other).sign() == 0
        except (InvalidCertificate,TypeError,ValueError):
            return False

    def __lt__(self,other): return (self-other).sign() < 0
    def __le__(self,other): return (self-other).sign() <= 0
    def __gt__(self,other): return (self-other).sign() > 0
    def __ge__(self,other): return (self-other).sign() >= 0
    def __abs__(self): return self if self.sign() >= 0 else -self


def as_real(value):
    if isinstance(value,R):
        return value
    if isinstance(value,dict) and set(value) == {'a','b','d'}:
        return R(value['a'],value['b'],value['d'])
    if all(hasattr(value,key) for key in ('a','b','d')):
        return R(value.a,value.b,value.d)
    return rational(value)


@dataclass(frozen=True)
class C:
    re: Q = Q(0)
    im: Q = Q(0)

    def __post_init__(self):
        object.__setattr__(self, 're', as_real(self.re))
        object.__setattr__(self, 'im', as_real(self.im))

    def __add__(self, other):
        o = complex_rational(other)
        return C(self.re + o.re, self.im + o.im)
    __radd__ = __add__

    def __neg__(self):
        return C(-self.re, -self.im)

    def __sub__(self, other):
        return self + -complex_rational(other)

    def __rsub__(self, other):
        return complex_rational(other) + -self

    def __mul__(self, other):
        o = complex_rational(other)
        return C(self.re * o.re - self.im * o.im,
                 self.re * o.im + self.im * o.re)
    __rmul__ = __mul__

    def __truediv__(self, other):
        o = complex_rational(other)
        d = o.re * o.re + o.im * o.im
        require(d > 0, 'zero complex denominator')
        return self * C(o.re / d, -o.im / d)

    def __rtruediv__(self, other):
        return complex_rational(other) / self


def complex_rational(value):
    if isinstance(value, C):
        return value
    if isinstance(value, dict):
        return C(value['re'], value['im'])
    if isinstance(value, (tuple, list)) and len(value) == 2:
        return C(*value)
    if hasattr(value, 're') and hasattr(value, 'im'):
        return C(value.re, value.im)
    return C(as_real(value))


def parse(data):
    n1, f1, n2, f2, k = (rational(data[x]) for x in ('N1', 'F1', 'N2', 'F2', 'k'))
    w = complex_rational(data['w'])
    require(isinstance(w.re,Q) and isinstance(w.im,Q), 'frequency must be rational')
    sector = data.get('sector', 'longitudinal')
    require(n1 >= 0 and n2 >= 0 and abs(f1) <= n1 and abs(f2) <= n2,
            'nonrealizable positive species moments')
    require(k > 0 and w.im > 0, 'requires positive k and upper-half-plane frequency')
    require(sector in ('longitudinal', 'transverse'), 'unknown sector')
    return n1, f1, n2, f2, w, k, sector


def sector_coefficients(n1, f1, n2, f2, w, k, sector):
    """Direct recurrence elimination: I1=(w I0-g0)/k, I2=(w I1-g1)/k.

    D_L=(1+I0)(-1+I2)-I1²; T=-1+(I0-I2)/2.
    k² D_L = A I0-B, 2k² T = A I0-B, I0=(s1-s2)/k.
    """
    g0, g1 = n1-n2, f1-f2
    if sector == 'longitudinal':
        return w*w + w*g0 - k*(k+g1), k*(k+g1) + w*g0 + g0*g0
    return k*k-w*w, 2*k*k-w*g0-k*g1


def transformed_moments(n, f, z, s):
    """nu=mu/((u-Re z)^2+(Im z)^2), t_j=integral u^j dnu."""
    t0 = -s.im / z.im
    t1 = z.re*t0-s.re
    t2 = n + 2*z.re*t1-(z.re*z.re+z.im*z.im)*t0
    t3 = f + 2*z.re*t2-(z.re*z.re+z.im*z.im)*t1
    return (t0, t1, t2, t3)


def support_matrices(n, f, z, s):
    t = transformed_moments(n, f, z, s)
    # Gram matrices of (1,u), weighted by (1+u) and (1-u).
    return tuple((t[0]+e*t[1], t[1]+e*t[2], t[2]+e*t[3]) for e in (1, -1))


def positive_semidefinite(m):
    return m[0] >= 0 and m[2] >= 0 and m[0]*m[2] >= m[1]*m[1]


def check_witness(n, f, z, s, witness):
    """Validate a positive one/two-atom measure without irrational evaluation.

    Two nodes c+-sqrt(d)/2 have weights determined by total and mean. The
    checked support/mean inequalities prove both nodes and weights valid.
    Polynomial reduction recovers moments through degree three exactly.
    Multiplication of the weights by |z-u|² produces the positive mu.
    """
    kind = witness['kind']
    if kind == 'zero':
        require(n == f == 0 and s == C(), 'invalid zero witness')
        return True
    if kind == 'one':
        u, mass = (as_real(witness[x]) for x in ('u', 'nu_mass'))
        require(-1 <= u <= 1 and mass > 0, 'invalid one-atom support/weight')
        t = tuple(mass*u**j for j in range(4))
    elif kind == 'two':
        a, b, t0, t1 = (as_real(witness[x]) for x in ('p2_a','p2_b','nu_t0','nu_t1'))
        d, center = a*a-4*b, -a/2
        require(t0 > 0 and d > 0, 'invalid two-atom mass/discriminant')
        require(abs(center) <= 1 and d/4 <= (1-abs(center))**2,
                'quadrature node outside support')
        require((t1/t0-center)**2 <= d/4, 'negative quadrature weight')
        t2 = -a*t1-b*t0
        t = (t0, t1, t2, -a*t2-b*t1)
    else:
        raise InvalidCertificate('unknown witness kind')
    require(t == transformed_moments(n, f, z, s), 'witness moment/transform mismatch')
    require(all(positive_semidefinite(m) for m in support_matrices(n, f, z, s)),
            'witness violates necessary Hausdorff positivity')
    return True


def normalized_constraints(n1, f1, n2, f2, z, target):
    """Independent affine interpolation in normalized coordinates s1=scale*(x+iy).

    This exactly reconstructs the certificate's positive row normalization.
    |s1|<=N1/Im(z) implies |x|<=1, -1<=y<=0; these box inequalities
    are necessary for any positive mu1, even on singular moment boundaries.
    """
    scale = n1/z.im
    require(scale > 0, 'invalid dual scaling')
    constraints = []
    for n, f, shift in ((n1, f1, C()), (n2, f2, target)):
        at = [support_matrices(n, f, z, scale*v-shift)
              for v in (C(), C(1), C(0,1))]
        for j in range(2):
            base = at[0][j]
            dx = tuple(v-u for u,v in zip(base, at[1][j]))
            dy = tuple(v-u for u,v in zip(base, at[2][j]))
            norm = max(Q(1), *(abs(v) for row in (base,dx,dy) for v in row))
            constraints.append(tuple(tuple(v/norm for v in row) for row in (base,dx,dy)))
    return constraints


def check_dual(n1, f1, n2, f2, z, target, proof):
    hs = normalized_constraints(n1, f1, n2, f2, z, target)
    box = ((Q(1),Q(1),Q(0)), (Q(1),Q(-1),Q(0)),
           (Q(1),Q(0),Q(1)), (Q(0),Q(0),Q(-1)))
    total = [Q(0),Q(0),Q(0)]
    require(len(proof['terms']) > 0, 'empty dual')
    for term in proof['terms']:
        j, weight = term['constraint'], rational(term['weight'])
        require(type(j) is int and 0 <= j < 8 and weight >= 0, 'invalid dual index/weight')
        if j < 4:
            v = tuple(rational(x) for x in term['vector'])
            require(len(v) == 2, 'invalid PSD direction')
            a, b = v
            row = tuple(a*a*m[0]+2*a*b*m[1]+b*b*m[2] for m in hs[j])
        else:
            require(term.get('vector') is None, 'box constraint has PSD vector')
            row = box[j-4]
        for r in range(3):
            total[r] += weight*row[r]
    require(total[1:] == [0,0] and total[0] < 0, 'dual does not separate')
    require(total[0] == rational(proof['negative_constant']), 'dual constant mismatch')
    return True


def independent_disks(n1,f1,n2,f2,z,target):
    """Recover circle polynomials from independently evaluated determinants.

    Each determinant is quadratic in (Re s,Im s). Interpolation recovers
    its coefficients, rather than trusting the production center formulas.
    Necessity of these disks for PSD suffices for the emptiness proof.
    """
    result=[]
    for n,f,shift in ((n1,f1,C()),(n2,f2,target)):
        for j in range(2):
            def determinant_at(x,y):
                p,q,r=support_matrices(n,f,z,C(x,y)-shift)[j]
                return p*r-q*q
            constant=determinant_at(0,0)
            xp,xm,yp,ym=(determinant_at(*p) for p in ((1,0),(-1,0),(0,1),(0,-1)))
            lx,ly=(xp-xm)/2,(yp-ym)/2
            xx,yy=(xp+xm-2*constant)/2,(yp+ym-2*constant)/2
            require(xx == yy and xx < 0, 'determinant is not a negative circle quadratic')
            require(determinant_at(1,1) == constant+lx+ly+xx+yy,
                    'unexpected mixed determinant coefficient')
            center=(-lx/(2*xx),-ly/(2*xx))
            radius2=sum((v*v for v in center),Q(0))-constant/xx
            require(radius2 >= 0, 'negative disk radius')
            result.append((center,radius2))
    return result


def independent_candidates(disks):
    for center,_ in disks:
        yield tuple(R(v) for v in center)
    for (center,radius2),(other,other_radius2) in combinations(disks,2):
        dx,dy=other[0]-center[0],other[1]-center[1]
        distance2=dx*dx+dy*dy
        if distance2 == 0:
            continue
        # Radical axis foot; h² is the perpendicular coordinate squared.
        alpha=(radius2-other_radius2+distance2)/(2*distance2)
        height2=radius2/distance2-alpha*alpha
        if height2 < 0:
            continue
        footx,footy=center[0]+alpha*dx,center[1]+alpha*dy
        for orientation in ((1,) if height2 == 0 else (-1,1)):
            yield (R(footx,-orientation*dy,height2),R(footy,orientation*dx,height2))


def in_disk(point,disk):
    center,radius2=disk
    return (point[0]-center[0])**2+(point[1]-center[1])**2 <= radius2


def check_four_disk_emptiness(n1,f1,n2,f2,z,target):
    disks=independent_disks(n1,f1,n2,f2,z,target)
    # Closed disks intersect iff a center or a pair-boundary crossing is
    # feasible. Zero radius centers, nesting, coincidence and tangency are
    # included; proof of this geometric theorem is separate from sampling.
    for point in independent_candidates(disks):
        require(not all(in_disk(point,disk) for disk in disks),
                'claimed empty disk intersection contains an exact candidate')
    return True


def verify(data, result):
    """Check one exact classification, raising InvalidCertificate on refusal."""
    try:
        return _verify(data, result)
    except InvalidCertificate:
        raise
    except (ValueError, TypeError, KeyError, IndexError, ZeroDivisionError, AttributeError) as e:
        raise InvalidCertificate(str(e)) from e


def _verify(data, result):
    n1,f1,n2,f2,w,k,sector = parse(data)
    a,b = sector_coefficients(n1,f1,n2,f2,w,k,sector)
    status = result['status']
    if status == 'impossible' and result.get('proof_kind') == 'zero_A_nonzero_B':
        require(a == C() and b != C(), 'invalid zero-A proof')
        return True
    require(a != C(), 'missing zero-A proof')
    target, z = k*b/a, w/k
    if status == 'possible':
        s1,s2 = (complex_rational(result[x]) for x in ('s1','s2'))
        require(s1-s2 == target, 'sector root identity fails')
        require(len(result['witnesses']) == 2, 'need both species witnesses')
        for n,f,s,h in zip((n1,n2),(f1,f2),(s1,s2),result['witnesses']):
            check_witness(n,f,z,s,h)
        return True
    require(status == 'impossible', 'unresolved is not a certificate')
    kind = result['proof_kind']
    if kind == 'four_disks':
        return check_four_disk_emptiness(n1,f1,n2,f2,z,target)
    if kind == 'rational_dual':
        # The submitted scaling is meaningful only when N1 is nonzero.
        return check_dual(n1,f1,n2,f2,z,target,result['dual'])
    require(kind == 'fixed_species', 'unknown impossibility proof')
    fixed = []
    for n,f in ((n1,f1),(n2,f2)):
        fixed.append(C() if n == 0 else n/(z-f/n) if abs(f) == n else None)
    require(any(s is not None for s in fixed), 'no fixed species')
    expected1 = fixed[0] if fixed[0] is not None else fixed[1]+target
    expected2 = expected1-target
    require(complex_rational(result['s1']) == expected1 and
            complex_rational(result['s2']) == expected2, 'wrong forced transform')
    actual = []
    for i,(n,f,s) in enumerate(((n1,f1,expected1),(n2,f2,expected2))):
        for j,m in enumerate(support_matrices(n,f,z,s)):
            if not positive_semidefinite(m):
                actual.append((i,j))
    require(actual and all(tuple(v) in actual for v in result['violations']) and
            len(result['violations']) > 0, 'invalid fixed-species violation')
    return True


def direct_dispersion(atoms1, atoms2, w, k):
    """Exact azimuthally averaged direct beam tensor, independent of targets."""
    w, k = complex_rational(w), rational(k)
    ints = [sum((sign*rational(m)*rational(u)**j/(w-k*rational(u))
                 for atoms,sign in ((atoms1,1),(atoms2,-1)) for u,m in atoms), C())
            for j in range(3)]
    i0,i1,i2 = ints
    return (1+i0)*(-1+i2)-i1*i1, -1+(i0-i2)/2
