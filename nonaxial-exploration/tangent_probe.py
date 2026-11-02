"""Discovery-only 30-node tangent LP. The exact certificate supersedes it."""
from itertools import permutations, product
from pathlib import Path
import hashlib
import json
import numpy as np
import scipy
from scipy.linalg import eig
from scipy.optimize import linprog


def run():
    labels = {tuple(a*b for a, b in zip(v, signs))
              for v in set(permutations((2, 2, 1)))
              for signs in product((-1, 1), repeat=3)}
    for axis in range(3):
        for sign in (-1, 1):
            v = [0, 0, 0]
            v[axis] = 3*sign
            labels.add(tuple(v))
    labels = sorted(labels)
    v = np.array(labels)/3
    inverse = np.array([labels.index(tuple(-a for a in t)) for t in labels])
    w = np.array([3/16 if t[2] == 1 and abs(t[0]) == abs(t[1]) == 2
                  else 1/4 if t == (0, 0, -3) else 0 for t in labels])
    K = 1-v@v.T
    A = np.diag(2*v[:, 2])-K*(w-w[inverse])[None, :]
    values, left, right = eig(A, left=True, right=True)
    index = np.argmin(abs(values-2j/3))
    q, l = right[:, index], left[:, index]
    denominator = np.vdot(l, q)
    gradient = np.array([(np.vdot(l, -K[:, j]*q[j]+K[:, inverse[j]]*q[inverse[j]])
                          / denominator).imag for j in range(len(v))])
    M = np.array([np.ones(len(v)), *v.T, v[:, 0]**2, v[:, 1]**2,
                  v[:, 0]*v[:, 1], v[:, 0]*v[:, 2], v[:, 1]*v[:, 2]])
    b = np.array([1, 0, 0, 0, 1/3, 1/3, 0, 0, 0])
    result = linprog(-gradient, A_eq=M, b_eq=b, bounds=(0, None), method='highs')
    if not result.success:
        raise RuntimeError(result.message)
    report = {'scope': 'Numerical discovery only, fixed 30 nodes and inversion-paired species',
              'nodes': len(v), 'base_root': [float(values[index].real), float(values[index].imag)],
              'spectral_separation': float(min(abs(np.delete(values, index)-values[index]))),
              'lp_success': bool(result.success), 'moment_residual': float(max(abs(M@result.x-b))),
              'directional_gain': float(gradient@(result.x-w)),
              'candidate_support': [{'integer_direction': labels[j], 'weight': float(a)}
                                    for j, a in enumerate(result.x) if a > 1e-10],
              'numpy': np.__version__, 'scipy': scipy.__version__,
              'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'plan_sha256': hashlib.sha256((Path(__file__).parent/'TANGENT_PLAN.md').read_bytes()).hexdigest()}
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    run()
