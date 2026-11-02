"""Exact algebra for the nonaxial purely imaginary cap; no off-axis claim."""
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import hashlib
import json
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from independent_verifier import C
from exact_tangent_certificate import run as sharp_check


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def add(*polys):
    out = {}
    for p in polys:
        for ij, a in p.items():
            out[ij] = out.get(ij, F(0))+a
    return {ij: a for ij, a in out.items() if a}


def scale(p, a):
    return {ij: F(a)*v for ij, v in p.items() if a*v}


def mul(p, q):
    out = {}
    for (i, j), a in p.items():
        for (k, l), b in q.items():
            ij = i+k, j+l
            out[ij] = out.get(ij, F(0))+a*b
    return {ij: a for ij, a in out.items() if a}


def polynomial_checks():
    one, y, z = {(0, 0): F(1)}, {(1, 0): F(1)}, {(0, 1): F(1)}
    y2, z2, yz = mul(y, y), mul(z, z), mul(y, z)
    s, d, e = add(one, scale(z2, -1)), add(z2, y2), add(one, y2)
    B = add(y, mul(e, z), scale(mul(y, z2), -1))
    zm = add(z, scale(y, -1))
    # P00 and mixed-entry identities, with positive denominators cleared.
    require(add(mul(add(s, scale(yz, 2)), d), scale(mul(yz, e), -2))
            == mul(mul(zm, zm), s), 'Upper diagonal identity failed')
    require(add(mul(add(one, scale(y2, 2), scale(yz, -1)), d), scale(mul(y2, e), -2))
            == mul(zm, B), 'Mixed entry identity failed')
    positive_L = add(d, scale(mul(yz, s), 2))
    square = add(z, mul(y, s))
    require(positive_L == add(mul(square, square), mul(mul(y2, z2), add(scale(one, 2), scale(z2, -1)))),
            'Parallel diagonal nonnegative decomposition failed')
    require(add(mul(e, positive_L), scale(mul(B, B), -1)) == mul(mul(y2, s), d),
            'Two-by-two determinant polynomial identity failed')


def determinant(a):
    if len(a) == 1:
        return a[0][0]
    return sum((-1)**j*a[0][j]*determinant([row[:j]+row[j+1:] for row in a[1:]])
               for j in range(len(a)))


def positive(a, strict=False):
    for size in range(1, len(a)+1):
        for ids in combinations(range(len(a)), size):
            value = determinant([[a[i][j] for j in ids] for i in ids])
            if value < 0 or strict and value == 0:
                return False
    return True


def rotation(q):
    w, x, y, z = map(F, q)
    n = w*w+x*x+y*y+z*z
    return [[1-2*(y*y+z*z)/n, 2*(x*y-z*w)/n, 2*(x*z+y*w)/n],
            [2*(x*y+z*w)/n, 1-2*(x*x+z*z)/n, 2*(y*z-x*w)/n],
            [2*(x*z-y*w)/n, 2*(y*z+x*w)/n, 1-2*(x*x+y*y)/n]]


def atoms(q, mass):
    base = [((F(x, 3), F(y, 3), F(1, 3)), F(3, 16))
            for x in (-2, 2) for y in (-2, 2)] + [((F(0), F(0), F(-1)), F(1, 4))]
    rot = rotation(q)
    out = [(tuple(sum(row[j]*v[j] for j in range(3)) for row in rot), mass*w)
           for v, w in base]
    require(sum(w for v, w in out) == mass, 'Wrong mass')
    for v, w in out:
        require(sum(a*a for a in v) == 1 and w > 0, 'Invalid atom')
    for i in range(3):
        require(sum(w*v[i] for v, w in out) == 0, 'Nonzero species flux')
        for j in range(3):
            require(sum(w*v[i]*v[j] for v, w in out) == (mass/3 if i == j else 0),
                    'Wrong species tensor')
    return out


def point_matrix(v, y):
    x, u, z = v
    d, c = z*z+y*y, F(3, 2)/(1+y*y)
    a, b = c*(1+2*y*y), -c*y
    P = [[F(0) for j in range(3)] for i in range(3)]
    P[0][0] = c*(1-z*z+2*y*z)-3*y*z/d
    for i, vi in enumerate((x, u), 1):
        P[0][i] = P[i][0] = (-3*y*y/d+a+b*z)*vi
        for j, vj in enumerate((x, u), 1):
            P[i][j] = (F(3, 2)*z*z if i == j else 0)+(F(3, 2)+3*y*z/d)*vi*vj
    require(positive(P), 'Point matrix not positive semidefinite')
    return P


def matrix_check(plus, minus, mass, gamma, y):
    k, omega, r = gamma/y, C(0, gamma), C(0, y)
    Pi = [[C(1 if i == j == 0 else -1 if i == j else 0) for j in range(4)] for i in range(4)]
    for distribution, sign in ((plus, 1), (minus, -1)):
        for v, weight in distribution:
            four = (F(1), v[2], v[0], v[1])
            for i in range(4):
                for j in range(4):
                    Pi[i][j] += sign*weight*four[i]*four[j]/(omega-k*v[2])
    S, T = Pi[0][0]-1, Pi[0][2:]
    require(Pi[0][1] == r*S and Pi[1][1] == -1+r*r*S, 'Longitudinal reduction failed')
    require(Pi[1][2:] == [r*a for a in T], 'Mixed reduction failed')
    M = [[Pi[0][0]-1+1/(1-r*r), *T],
         [T[0], Pi[2][2], Pi[2][3]], [T[1], Pi[3][2], Pi[3][3]]]
    J = [C(0, 1), C(1), C(1)]
    H = [[-J[i]*M[i][j]*J[j] for j in range(3)] for i in range(3)]
    herm = [[(H[i][j]+C(H[j][i].re, -H[j][i].im))/2 for j in range(3)] for i in range(3)]
    require(all(a.im == 0 for row in herm for a in row), 'Unexpected Hermitian imaginary entries')
    herm = [[a.re for a in row] for row in herm]
    effective = [(v, w/(2*mass)) for v, w in plus]
    effective += [(tuple(-a for a in v), w/(2*mass)) for v, w in minus]
    t, diagonal = 2*mass/(3*gamma), [1/(1+y*y), F(1), F(1)]
    expected = [[(1-t)*diagonal[i] if i == j else F(0) for j in range(3)] for i in range(3)]
    for v, weight in effective:
        P = point_matrix(v, y)
        for i in range(3):
            for j in range(3):
                expected[i][j] += t*weight*P[i][j]
    require(herm == expected, 'Hermitian part differs from positive certificate')
    if gamma > 2*mass/3:
        require(positive(herm, strict=True), 'Hermitian part not positive definite above cap')


def run():
    polynomial_checks()
    pairs = [((1, 0, 0, 0), (0, 1, 0, 0)),
             ((1, 2, 3, 4), (4, -1, 2, 3)),
             ((2, 1, -3, 1), (1, 3, 2, -2))]
    count = above = 0
    for mass in (F(1), F(7, 5)):
        for qa, qb in pairs:
            plus, minus = atoms(qa, mass), atoms(qb, mass)
            for ratio in (F(1, 3), F(2, 3), F(1), F(7, 3)):
                gamma = mass*ratio
                for y in (F(1, 10), F(1, 3), F(1), F(3)):
                    matrix_check(plus, minus, mass, gamma, y)
                    count += 1
                    above += int(gamma > 2*mass/3)
    require(sharp_check()['passed'], 'Existing full-polarization attainer failed')
    return {'passed': True, 'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'bivariate_polynomial_identities': 4, 'full_matrix_identity_cases': count,
            'strict_positive_cases_above_cap': above, 'point_psd_checks': 10*count,
            'sharp_attainer_replayed': True,
            'statement': 'gamma<=2N/3 for purely imaginary roots with common isotropic full moments, without axial or inversion symmetry',
            'limits': 'Not nonzero-real-frequency roots; primary checked, no independent review or priority claim'}


if __name__ == '__main__':
    print(json.dumps(run(), indent=2))
