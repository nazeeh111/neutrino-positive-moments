"""Primary exact feasibility checks against original polarization equations."""
from fractions import Fraction as F
from itertools import permutations, product
from pathlib import Path
import hashlib
import json
import random
from general_moment_certificate import C, moments, evaluate, mean_Q, constant_rayleigh_squared
from exact_nonaxial_imaginary_bound import determinant, positive
from exact_nonaxial_full_bound import point_matrix


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def encode(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, C):
        return {'real': encode(value.re), 'imag': encode(value.im)}
    if isinstance(value, list):
        return [encode(x) for x in value]
    if isinstance(value, dict):
        return {key: encode(x) for key, x in value.items()}
    return value


def original_polarization(plus, minus, omega, k):
    Pi = [[C(1 if i == j == 0 else -1 if i == j else 0) for j in range(4)] for i in range(4)]
    for distribution, sign in ((plus, 1), (minus, -1)):
        for v, weight in distribution:
            four = [F(1), v[2], v[0], v[1]]
            for i in range(4):
                for j in range(4):
                    Pi[i][j] += sign*weight*four[i]*four[j]/(omega-k*v[2])
    return Pi


def verify_identity(plus, minus, k, a, b, h, t):
    mp, mm = moments(plus), moments(minus)
    result = evaluate(mp, mm, k, a, b, h, t)
    require(result['Delta'] != C(), 'Unexpected singular fixture')
    omega, Delta = result['omega'], result['Delta']
    Pi = original_polarization(plus, minus, omega, k)
    g0 = mp['N']-mm['N']
    g = [x-y for x, y in zip(mp['F'], mm['F'])]
    Omega, Kz = omega+g0, k+g[2]
    # Map reduced coordinates (q0,bx,by) to the original covector.
    transform = [[k*Kz/Delta, omega*g[0]/Delta, omega*g[1]/Delta],
                 [-k*Omega/Delta, -k*g[0]/Delta, -k*g[1]/Delta],
                 [C(), C(1), C()], [C(), C(), C(1)]]
    projection = [[sum((Pi[i][s]*transform[s][j] for s in range(4)), C())
                   for j in range(3)] for i in (0, 2, 3)]
    M = [[value for value in row] for row in result['M0']]
    Psum = [[F(0) for _ in range(3)] for _ in range(3)]
    for distribution, sign in ((plus, 1), (minus, -1)):
        for v, weight in distribution:
            three = [F(1), v[0], v[1]]
            P = point_matrix(v if sign == 1 else tuple(-x for x in v), sign*a, b, h, t)
            for i in range(3):
                for j in range(3):
                    M[i][j] += sign*weight*three[i]*three[j]/(omega-k*v[2])
                    Psum[i][j] += weight*P[i][j]
    require(M == projection, 'Reduced matrix differs from original row equations')
    require(determinant(Pi) == -Delta*determinant(M)/(k*k), 'Original determinant relation failed')
    J = result['J']
    H = [[-J[i]*M[i][j]*J[j] for j in range(3)] for i in range(3)]
    for i in range(3):
        for j in range(3):
            actual = (H[i][j]+C(H[j][i].re, -H[j][i].im))/2
            expected = result['lower_matrix'][i][j]+Psum[i][j]/(3*omega.im)
            require(actual == C(expected), 'Full Hermitian decomposition failed')
    if result['status'] == 'excluded':
        require(determinant(Pi) != C(), 'Excluded a root of the original matrix')
    return result


def run():
    nodes = []
    for axis in range(3):
        for sign in (1, -1):
            v = [F(0)]*3
            v[axis] = F(sign)
            nodes.append(tuple(v))
    plus = list(zip(nodes, map(F, ['1/5', '1/10', '1/5', '1/10', '1/5', '1/5'])))
    minus = list(zip(nodes, map(F, ['1/10', '1/5', '1/10', '1/10', '1/5', '1/10'])))
    fixture = verify_identity(plus, minus, F(99, 25), F(0), F(1), F(101, 99), F(20, 99))
    require(fixture['status'] == 'excluded', 'Frozen strict-improvement fixture failed')
    norm_lower = constant_rayleigh_squared(moments(plus), moments(minus))
    require(norm_lower == F(663, 800) and fixture['omega'].im**2 < norm_lower,
            'Norm baseline comparison failed')
    south, north = (F(0), F(0), F(-1)), (F(0), F(0), F(1))
    root = verify_identity([(south, F(1))], [(north, F(1))], F(6, 5), F(0), F(1), F(5, 3), F(4, 3))
    require(root['omega'] == C(0, F(8, 5)), 'Wrong exact root frequency')
    require(determinant(original_polarization([(south, F(1))], [(north, F(1))], root['omega'], root['k'])) == C(),
            'Two-beam fixture is not an original root')
    require(root['status'] == 'unresolved', 'False exclusion of exact growing root')
    # Adversarial shortcut: equal-current base matrix would falsely exclude it.
    qp = mean_Q(moments([(south, F(1))]), F(0), F(1), F(5, 3), F(4, 3))
    qm = mean_Q(moments([(north, F(1))]), F(0), F(1), F(5, 3), F(4, 3), reflected=True)
    wrong = [[(F(1) if i == j else 0)-(qp[i][j]+qm[i][j])/(3*root['omega'].im)
              for j in range(3)] for i in range(3)]
    require(positive(wrong, strict=True), 'Negative control no longer tests the known invalid shortcut')
    singular = evaluate(moments([(south, F(25, 24))]), moments([(north, F(25, 24))]),
                        F(3, 4), F(0), F(1), F(5, 3), F(4, 3))
    require(singular['Delta'] == C() and singular['status'] == 'unresolved', 'Singular branch mishandled')
    zero = evaluate(moments([]), moments([]), F(1), F(0), F(1), F(5, 3), F(4, 3))
    require(zero['status'] == 'excluded', 'Zero-interaction control failed')
    grid = set(nodes)
    for raw in set(permutations((2, 2, 1))):
        for signs in product((-1, 1), repeat=3):
            grid.add(tuple(F(x*s, 3) for x, s in zip(raw, signs)))
    grid = sorted(grid)
    rng = random.Random(2026110208)
    statuses = {'excluded': 0, 'unresolved': 0}
    unequal = nonparallel = off_axis = 0
    for index in range(40):
        pp = [(v, F(rng.randint(1, 9), 23)) for v in rng.sample(grid, 5)]
        pm = [(v, F(rng.randint(1, 9), 23)) for v in rng.sample(grid, 7)]
        circle = (F(1, 2), F(1), F(2))[index % 3]
        hyper = (F(1, 10), F(1, 3), F(1, 2))[(index//3) % 3]
        a, b = (1-circle*circle)/(1+circle*circle), 2*circle/(1+circle*circle)
        h, t = (1+hyper*hyper)/(1-hyper*hyper), 2*hyper/(1-hyper*hyper)
        answer = verify_identity(pp, pm, F(rng.randint(1, 12), 5), a, b, h, t)
        statuses[answer['status']] += 1
        mp, mm = moments(pp), moments(pm)
        unequal += int(mp['N'] != mm['N'])
        nonparallel += int(any(mp['F'][i]*mm['F'][j] != mp['F'][j]*mm['F'][i] for i in range(3) for j in range(3)))
        off_axis += int(a != 0)
    root_dir = Path(__file__).resolve().parent
    return encode({'passed': True, 'random_original_matrix_identity_cases': 40,
                   'unequal_mass_cases': unequal, 'nonparallel_flux_cases': nonparallel,
                   'nonzero_real_frequency_cases': off_axis, 'statuses_not_an_accuracy_benchmark': statuses,
                   'fixed_fixture': {'status': fixture['status'], 'k': fixture['k'], 'omega': fixture['omega'],
                                     'lower_matrix': fixture['lower_matrix'], 'norm_squared_rayleigh_lower': norm_lower,
                                     'moments_plus': moments(plus), 'moments_minus': moments(minus)},
                   'exact_growing_root_not_excluded': True, 'invalid_equal_current_shortcut_detected': True,
                   'singular_delta_returns_unresolved': True, 'zero_species_control': True,
                   'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                   'certificate_source_sha256': hashlib.sha256((root_dir/'general_moment_certificate.py').read_bytes()).hexdigest(),
                   'plan_sha256': hashlib.sha256((root_dir/'GENERAL_CERTIFICATE_PLAN.md').read_bytes()).hexdigest(),
                   'scope': 'Primary exact feasibility checks; not independent review, complete classification, physical-use pass or priority evidence'})


if __name__ == '__main__':
    print(json.dumps(run(), indent=2))
