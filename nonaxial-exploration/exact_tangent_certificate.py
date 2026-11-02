"""Exact first-variation certificate at k=2, omega=2i/3; not a global bound."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from independent_verifier import C


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def mul(a, b):
    out = [F(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] += x*y
    return out


def run():
    omega, k = C(0, F(2, 3)), F(2)
    plus = [((F(x, 3), F(y, 3), F(1, 3)), F(3, 16))
            for x in (-2, 2) for y in (-2, 2)] + [((F(0), F(0), F(-1)), F(1, 4))]
    minus = [(tuple(-a for a in v), m) for v, m in plus]
    for atoms in (plus, minus):
        require(sum(m for v, m in atoms) == 1, 'Wrong mass')
        for v, m in atoms:
            require(m > 0 and sum(a*a for a in v) == 1, 'Invalid sphere atom')
        for i in range(3):
            require(sum(m*v[i] for v, m in atoms) == 0, 'Nonzero flux')
            for j in range(3):
                require(sum(m*v[i]*v[j] for v, m in atoms) == (F(1, 3) if i == j else 0),
                        'Nonisotropic tensor')
    signed = [(F(1, 3), F(3, 4)), (F(-1), F(1, 4)),
              (F(-1, 3), F(-3, 4)), (F(1), F(-1, 4))]
    I = [sum((m*z**j/(omega-k*z) for z, m in signed), C()) for j in range(3)]
    Ip = [sum((-m*z**j/((omega-k*z)*(omega-k*z))
               for z, m in signed), C()) for j in range(3)]
    require(I == [C(F(-9, 10)), C(0, F(-3, 10)), C(F(1, 10))], 'Wrong transform')
    D = (1+I[0])*(-1+I[2])-I[1]*I[1]
    Dp = Ip[0]*(-1+I[2])+(1+I[0])*Ip[2]-2*I[1]*Ip[1]
    transverse = -1+(I[0]-I[2])/2
    Pi = [[C(1 if i == j == 0 else -1 if i == j else 0)
           for j in range(4)] for i in range(4)]
    for atoms, sign in ((plus, 1), (minus, -1)):
        for v, m in atoms:
            four = (F(1), *v)
            for i in range(4):
                for j in range(4):
                    Pi[i][j] += sign*m*four[i]*four[j]/(omega-k*v[2])
    expected = [[C() for j in range(4)] for i in range(4)]
    expected[0][0], expected[3][3] = 1+I[0], -1+I[2]
    expected[0][3] = expected[3][0] = I[1]
    expected[1][1] = expected[2][2] = transverse
    require(Pi == expected, 'Full nonaxial base polarization differs')
    require(D == C() and Dp == C(0, F(3, 2)), 'Root or derivative incorrect')
    require(transverse == C(F(-3, 2)), 'Other polarization factors singular')
    require(transverse*transverse*Dp == C(0, F(27, 8)), 'Full root not simple')
    # Coefficients establish the derivative as a rational function, not samples.
    require([Dp*omega, -k*Dp] == [C(-1), C(0, -3)], 'Wrong denominator polynomial')
    require([-1+I[2], -2*I[1], 1+I[0]] ==
            [C(F(-9, 10)), C(0, F(3, 5)), C(F(1, 10))], 'Wrong numerator polynomial')
    # 5(1+9z²)(1+3z/5-z²)-3z(11-z²)
    # = 5(1-3z)²(1-z²), nonnegative on [-1,1].
    lhs = mul([5, 3, -5], [1, 0, 9])
    lhs[1] -= 33
    lhs[3] += 3
    rhs = [5*x for x in mul([1, -6, 9], [1, 0, -1])]
    require(lhs == rhs == [5, -30, 40, 30, -45], 'Dual polynomial identity failed')
    h = lambda z: 3*z*(11-z*z)/(10*(1+9*z*z))
    for j in range(-20, 21):
        z = F(j, 20)
        derivative = -((-1+I[2])-2*I[1]*z+(1+I[0])*z*z)/((omega-k*z)*Dp)
        require(derivative.im == h(z), 'Derivative implementation differs from formula')
        require((1+F(3, 5)*z-z*z)/2-h(z) >= 0, 'Negative dual slack')
    base_plus = [(F(1, 3), F(3, 4)), (F(-1), F(1, 4))]
    require([sum(m*z**j for z, m in base_plus) for j in range(3)] ==
            [F(1), F(0), F(1, 3)], 'Wrong base moments')
    require(sum(m*h(z) for z, m in base_plus) == F(1, 3), 'No dual contact')
    # Equality supports {-1,1/3,1}; moment equations uniquely fix the weights.
    moment_matrix = [[F(z)**j for z in (-1, F(1, 3), 1)] for j in range(3)]
    weights = [F(1, 4), F(3, 4), F(0)]
    require([sum(a*b for a, b in zip(row, weights)) for row in moment_matrix]
            == [F(1), F(0), F(1, 3)], 'Equality weights wrong')
    determinant = (moment_matrix[1][1]-moment_matrix[1][0])*(moment_matrix[1][2]-moment_matrix[1][0])*(moment_matrix[1][2]-moment_matrix[1][1])
    require(determinant == F(16, 9), 'Equality support not uniquely determined')
    return {'passed': True, 'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'base_frequency': '2i/3', 'k': '2', 'full_determinant_derivative': '27i/8',
            'dual_slack': '(1-3z)^2*(1-z^2)/(2*(1+9z^2))',
            'conclusion': 'Nonpositive first growth derivative for each moment-preserving mixture path from this base, including nonaxial and independently varied species',
            'limits': 'Not a global nonaxial cap, finite-change or higher-order result; no independent review completed'}


if __name__ == '__main__':
    print(json.dumps(run(), indent=2))
