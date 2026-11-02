"""Exact mixture/scaling witness for the full indistinguishable growth range.

This corroborates a corollary of the full nonaxial cap. It is not an
independent review or a benchmark on astrophysical data.
"""
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


def merge(atoms):
    out = {}
    for v, weight in atoms:
        out[v] = out.get(v, F(0))+weight
    return {v: weight for v, weight in out.items() if weight}


def moments(measure, mass):
    require(all(w > 0 and sum(a*a for a in v) == 1 for v, w in measure.items()), 'Invalid atom')
    require(sum(measure.values()) == mass, 'Wrong species mass')
    for i in range(3):
        require(sum(w*v[i] for v, w in measure.items()) == 0, 'Nonzero species flux')
        for j in range(3):
            require(sum(w*v[i]*v[j] for v, w in measure.items()) == (mass/3 if i == j else 0),
                    'Nonisotropic species tensor')


def run():
    base = [((F(x, 3), F(y, 3), F(1, 3)), F(3, 16))
            for x in (-2, 2) for y in (-2, 2)] + [((F(0), F(0), F(-1)), F(1, 4))]
    opposite = [(tuple(-a for a in v), w) for v, w in base]
    common = []
    for axis in range(3):
        for sign in (-1, 1):
            v = [F(0)]*3
            v[axis] = F(sign)
            common.append((tuple(v), F(1, 6)))
    require(set(v for v, w in base).isdisjoint(v for v, w in opposite), 'Base species overlap')
    records = []
    for mass in (F(1), F(7, 5)):
        original = merge([(v, mass*w) for v, w in base]+[(v, -mass*w) for v, w in opposite])
        for epsilon in (F(0), F(1, 100), F(1, 10), F(1, 3), F(1, 2), F(9, 10), F(1)):
            plus = merge([(v, mass*(1-epsilon)*w) for v, w in common]
                         +[(v, mass*epsilon*w) for v, w in base])
            minus = merge([(v, mass*(1-epsilon)*w) for v, w in common]
                          +[(v, mass*epsilon*w) for v, w in opposite])
            moments(plus, mass)
            moments(minus, mass)
            nodes = sorted(set(plus)|set(minus))
            g = [plus.get(v, F(0))-minus.get(v, F(0)) for v in nodes]
            require(all(weight == epsilon*original.get(v, F(0)) for v, weight in zip(nodes, g)),
                    'Common component did not cancel exactly')
            jp, jm = sum(max(w, F(0)) for w in g), sum(max(-w, F(0)) for w in g)
            require(jp == jm == epsilon*mass, 'Wrong positive/negative crossing masses')
            growth, k = 2*epsilon*mass/3, 2*epsilon*mass
            if epsilon:
                omega, ratio = C(0, growth), C(0, F(1, 3))
                vector = [(1-ratio*v[2])/(omega-k*v[2]) for v in nodes]
                require(any(value != C() for value in vector), 'Zero witness vector')
                for i, v in enumerate(nodes):
                    action = k*v[2]*vector[i]
                    action -= sum(((1-sum(a*b for a, b in zip(v, u)))*weight*value
                                   for u, weight, value in zip(nodes, g, vector)), C())
                    require(action == omega*vector[i], 'Original beam operator residual nonzero')
                require(jp/growth == F(3, 2), 'Wrong proxy-to-actual ratio')
            else:
                require(all(w == 0 for w in g), 'Stable endpoint still has signed interaction')
            records.append({'mass': str(mass), 'epsilon': str(epsilon), 'same_supplied_moments': True,
                            'maximum_growth': str(growth), 'attaining_k': str(k),
                            'positive_crossing_mass': str(jp), 'negative_crossing_magnitude': str(jm)})
    return {'passed': True, 'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'exact_cases': len(records), 'records': records,
            'upper_bound_basis': 'Apply the full nonaxial theorem to effective species epsilon*nu_s; the common component cancels from the actual operator',
            'analytical_range': '[0,2N/3] at identical supplied (N,0,NI/3) moments',
            'optimal_worst_case_absolute_error_for_a_moments_only_point_estimate': 'N/3',
            'scope': 'Arbitrary positive measures, all real wave vectors; primary-checked corollary, not a new practical-use pass or independent review'}


if __name__ == '__main__':
    print(json.dumps(run(), indent=2))
