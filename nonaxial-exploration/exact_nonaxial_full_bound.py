"""Exact full-frequency nonaxial cap certificate; primary checks, not review."""
from fractions import Fraction as F
from math import comb
from pathlib import Path
import hashlib
import json
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from independent_verifier import C
from exact_nonaxial_imaginary_bound import atoms, positive


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def add(*polys):
    out = {}
    for p in polys:
        for key, value in p.items():
            out[key] = out.get(key, F(0))+value
    return {key: value for key, value in out.items() if value}


def scale(p, factor):
    return {key: value*factor for key, value in p.items() if value*factor}


def mul(p, q):
    out = {}
    for i, a in p.items():
        for j, b in q.items():
            key = tuple(x+y for x, y in zip(i, j))
            out[key] = out.get(key, F(0))+a*b
    return {key: value for key, value in out.items() if value}


def reduce_constraints(p):
    # Variables (a,b,h,t,z); substitute b²=1-a² and h²=1+t².
    out = {}
    for (a, b, h, t, z), value in p.items():
        for i in range(b//2+1):
            for j in range(h//2+1):
                key = (a+2*i, b % 2, h % 2, t+2*j, z)
                out[key] = out.get(key, F(0))+value*comb(b//2, i)*(-1)**i*comb(h//2, j)
    return {key: value for key, value in out.items() if value}


def identity(left, right, label):
    require(not reduce_constraints(add(left, scale(right, -1))), label)


def polynomial_checks():
    one = {(0, 0, 0, 0, 0): F(1)}
    variables = []
    for j in range(5):
        powers = [0]*5
        powers[j] = 1
        variables.append({tuple(powers): F(1)})
    a, b, h, t, z = variables
    square = lambda p: mul(p, p)
    a2, b2, h2, t2, z2 = map(square, variables)
    x, y = mul(a, h), mul(b, t)
    v = add(z, scale(x, -1))
    s = add(one, scale(z2, -1))
    D = add(square(v), square(y))
    A = add(v, scale(y, -1))
    K = add(mul(b, h), mul(a, t))
    Kminus = add(mul(b, h), scale(mul(a, t), -1))
    ell = add(scale(mul(t, D), -1), mul(h, add(mul(Kminus, v), mul(K, y))))
    identity(add(D, scale(mul(y, v), -2)), square(A), 'Upper-diagonal square failed')
    mixed = add(mul(add(mul(b, add(one, scale(t2, 2))), scale(mul(t, z), -1)), D),
                scale(mul(mul(mul(b, t2), h), add(h, scale(mul(a, z), -1))), -2))
    identity(mixed, mul(A, ell), 'Mixed-entry factorization failed')
    rhs = add(mul(D, add(mul(square(K), z2), mul(add(square(K), scale(t2, -1)), s))),
              scale(mul(mul(mul(h2, y), v), s), 2))
    identity(square(ell), rhs, 'Rank-one block identity failed')
    identity(add(square(K), square(Kminus)), scale(add(b2, t2), 2), 'Averaged transverse constant failed')
    identity(add(h2, scale(b2, -1), scale(t2, F(-1, 2))), add(a2, scale(t2, F(1, 2))),
             'Nonnegative transverse remainder failed')


def point_matrix(v, a, b, h, t):
    vx, vy, z = v
    y = b*t
    D = (z-a*h)**2+y*y
    s = 1-z*z
    K = b*h+a*t
    e = (K*K-t*t)/(h*h)
    alpha, beta = b*(1+2*t*t)/h, -t/h
    F00 = -(z-a*h)*s/D
    P = [[F(0) for _ in range(3)] for _ in range(3)]
    P[0][0] = F(3, 2)*s+3*y*F00
    for i, vi in enumerate((vx, vy), 1):
        P[0][i] = P[i][0] = F(3, 2)*(alpha+beta*z)*vi-3*y*t*(h-a*z)*vi/D
        for j, vj in enumerate((vx, vy), 1):
            P[i][j] = (F(3, 2)*K*K*z*z/(h*h) if i == j else 0)
            P[i][j] += (F(3, 2)*e+3*y*(z-a*h)/D)*vi*vj
    require(positive(P), 'Point matrix not positive semidefinite')
    return P


def polarization(plus, minus, omega, k):
    Pi = [[C(1 if i == j == 0 else -1 if i == j else 0) for j in range(4)] for i in range(4)]
    for distribution, sign in ((plus, 1), (minus, -1)):
        for v, weight in distribution:
            four = (F(1), v[2], v[0], v[1])
            for i in range(4):
                for j in range(4):
                    Pi[i][j] += sign*weight*four[i]*four[j]/(omega-k*v[2])
    return Pi


def matrix_check(plus, minus, mass, gamma, a, b, h, t):
    x, y, p, q = a*h, b*t, b*h, -a*t
    r, root = C(x, y), C(p, q)
    require(root*root == 1-r*r and p > 0, 'Wrong square-root branch')
    k = gamma/y
    Pi = polarization(plus, minus, r*k, k)
    S, T = Pi[0][0]-1, Pi[0][2:]
    require(Pi[0][1] == r*S and Pi[1][1] == -1+r*r*S, 'Full reduction failed')
    require(Pi[1][2:] == [r*v for v in T], 'Mixed Ward identity failed')
    M = [[1/(1-r*r)+S, *T], [T[0], Pi[2][2], Pi[2][3]], [T[1], Pi[3][2], Pi[3][3]]]
    J = [C(0, 1)*root, C(1), C(1)]
    H = [[-J[i]*M[i][j]*J[j] for j in range(3)] for i in range(3)]
    herm = [[(H[i][j]+C(H[j][i].re, -H[j][i].im))/2 for j in range(3)] for i in range(3)]
    require(all(value.im == 0 for row in herm for value in row), 'Nonreal Hermitian part')
    herm = [[value.re for value in row] for row in herm]
    tau = 2*mass/(3*gamma)
    remainder = (a*a+t*t/2)/(h*h)
    expected = [[(1-tau)+(tau*remainder if i else 0) if i == j else F(0)
                 for j in range(3)] for i in range(3)]
    for distribution, sign in ((plus, 1), (minus, -1)):
        for v, weight in distribution:
            reflected = v if sign == 1 else tuple(-u for u in v)
            P = point_matrix(reflected, sign*a, b, h, t)
            for i in range(3):
                for j in range(3):
                    expected[i][j] += tau*weight*P[i][j]/(2*mass)
    require(herm == expected, 'Full-frequency Hermitian decomposition failed')
    if gamma > 2*mass/3:
        require(positive(herm, strict=True), 'Strict positivity above cap failed')


def oscillating_attainer():
    plus = atoms((1, 0, 0, 0), F(1))
    minus = []
    for axis in range(3):
        for sign in (-1, 1):
            v = [F(0)]*3
            v[axis] = F(sign)
            minus.append((tuple(v), F(1, 6)))
    omega, k = C(F(2, 3), F(2, 3)), F(4)
    Pi = polarization(plus, minus, omega, k)
    vector = [C(1), -omega/k, C(), C()]
    require(all(sum((row[j]*vector[j] for j in range(4)), C()) == C() for row in Pi),
            'Original full polarization rejects oscillating attainer')
    require(Pi[2][2] != C() and Pi[3][3] != C(), 'Unexpected transverse root')
    # Directly verify the second species' full moments as well.
    require(sum(w for v, w in minus) == 1, 'Attainer mass wrong')
    for i in range(3):
        require(sum(w*v[i] for v, w in minus) == 0, 'Attainer flux wrong')
        for j in range(3):
            require(sum(w*v[i]*v[j] for v, w in minus) == (F(1, 3) if i == j else 0),
                    'Attainer tensor wrong')
    return {'k': '4', 'real_frequency': '2/3', 'growth': '2/3', 'full_polarization_null_vector_exact': True}


def run():
    polynomial_checks()
    pairs = [((1, 0, 0, 0), (0, 1, 0, 0)), ((1, 2, 3, 4), (4, -1, 2, 3))]
    count = above = off_axis = 0
    for circle in (F(1, 2), F(3, 4), F(1), F(4, 3), F(2)):
        a, b = (1-circle*circle)/(1+circle*circle), 2*circle/(1+circle*circle)
        for hyper in (F(1, 10), F(1, 3), F(1, 2)):
            h, t = (1+hyper*hyper)/(1-hyper*hyper), 2*hyper/(1-hyper*hyper)
            require(a*a+b*b == 1 and h*h-t*t == 1, 'Parameter constraints failed')
            for mass in (F(1), F(7, 5)):
                for qa, qb in pairs:
                    plus, minus = atoms(qa, mass), atoms(qb, mass)
                    for ratio in (F(1, 3), F(2, 3), F(1)):
                        gamma = mass*ratio
                        matrix_check(plus, minus, mass, gamma, a, b, h, t)
                        count += 1
                        above += int(gamma > 2*mass/3)
                        off_axis += int(a != 0)
    return {'passed': True, 'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'polynomial_identities_modulo_parameter_constraints': 5,
            'full_matrix_identity_cases': count, 'nonzero_real_frequency_cases': off_axis,
            'strict_positive_cases_above_cap': above, 'point_psd_checks': 10*count,
            'oscillating_attainer': oscillating_attainer(),
            'statement': 'Sharp gamma<=2N/3 for all real frequency parts and wave vectors, independent positive species with common isotropic full moments',
            'limits': 'Primary checked, no independent review or priority claim; not the nonaxial fixed-k envelope or a utility pass'}


if __name__ == '__main__':
    print(json.dumps(run(), indent=2))
