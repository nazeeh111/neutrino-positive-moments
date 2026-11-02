"""Exact audit of a restricted published reconstruction branch, not FLASH.

Grohs et al., arXiv:2501.07540v2, Eqs. 25-31. Off-diagonal inputs
are strictly positive real and increasing; the circular branch is not used.
Run with --output PATH to retain the deterministic report. Standard library.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def mm(a, b):
    if a * b <= 0:
        return F(0)
    return min(a, b) if a > 0 else max(a, b)


def matrix(x, z):
    return [[(1 + z) / 2, x / 2], [x / 2, (1 - z) / 2]]


def multiply(a, b):
    return [[sum(a[i][k] * b[k][j] for k in range(2))
             for j in range(2)] for i in range(2)]


def rotate(a):
    # U = B/sqrt(5); no floating square root is needed in U A U^dagger.
    b = [[2, -1], [1, 2]]
    bt = [[2, 1], [-1, 2]]
    require(multiply(b, bt) == [[5, 0], [0, 5]], 'Nonunitary rotation')
    return [[v / 5 for v in row] for row in multiply(multiply(b, a), bt)]


def physical(a):
    return (a[0][1] == a[1][0] and a[0][0] >= 0 and a[1][1] >= 0
            and a[0][0] * a[1][1] >= a[0][1] ** 2
            and a[0][0] + a[1][1] == 1)


def reconstruct(cells, side=1):
    left, center, right = cells
    # This guard is essential: do not silently implement the general algorithm
    # by applying real minmod to arbitrary complex or decreasing off-diagonals.
    vals = [a[0][1] for a in cells]
    require(0 < vals[0] < vals[1] < vals[2], 'Outside audited line branch')
    require(all(a[0][1] == a[1][0] for a in cells), 'Non-Hermitian input')
    require(max(vals[1] - vals[0], vals[2] - vals[1]) < vals[2] - vals[0],
            'Published ordering check failed')
    out = [[center[i][j] + side * mm(right[i][j] - center[i][j],
                                    center[i][j] - left[i][j]) / 2
            for j in range(2)] for i in range(2)]
    require(out[0][1] > 0, 'Published modulus sign correction would be needed')
    return out


def subtract(a, b):
    return [[a[i][j] - b[i][j] for j in range(2)] for i in range(2)]


def frobenius2(a):
    return sum(v * v for row in a for v in row)


def purity(a):
    square = multiply(a, a)
    return square[0][0] + square[1][1]


def profile(t):
    return matrix(F(3, 10) + F(3, 20) * t + F(1, 20) * t * t,
                  F(3, 10) + F(3, 20) * t - F(1, 20) * t * t)


def check_case(h):
    cells = [profile(-h), profile(F(0)), profile(h)]
    rotated_cells = [rotate(a) for a in cells]
    require(all(physical(a) for a in cells + rotated_cells), 'Invalid input')
    face = reconstruct(cells)
    direct = reconstruct(rotated_cells)
    transformed = rotate(face)
    require(all(physical(a) for a in [face, direct, transformed]), 'Invalid face')
    # Conservativity of the local linear reconstruction, not a PDE update proof.
    for stencil in (cells, rotated_cells):
        lo, hi = reconstruct(stencil, -1), reconstruct(stencil, 1)
        require(physical(lo), 'Invalid left face')
        require([[ (lo[i][j] + hi[i][j]) / 2 for j in range(2)]
                 for i in range(2)] == stencil[1], 'Center not preserved')
    return face, direct, transformed, subtract(direct, transformed)


def run():
    face, direct, transformed, defect = check_case(F(1))
    require(face == matrix(F(7, 20), F(7, 20)), 'Wrong original face')
    require(direct == matrix(F(13, 25), F(-3, 50)), 'Wrong rotated face')
    require(transformed == matrix(F(49, 100), F(-7, 100)), 'Wrong conjugate')
    require(defect == [[F(1, 200), F(3, 200)],
                       [F(3, 200), F(-1, 200)]], 'Wrong exact defect')
    require(frobenius2(defect) == F(1, 2000), 'Wrong defect norm')
    require(purity(direct) - purity(face) == F(29, 2000), 'Wrong purity change')
    refinements = []
    for denominator in (4, 8, 16, 32, 64, 128):
        h = F(1, denominator)
        _, _, _, difference = check_case(h)
        require(difference == [[F(3, 200) * h*h, F(3, 200) * h*h],
                               [F(3, 200) * h*h, -F(3, 200) * h*h]],
                'Smooth-grid defect formula failed')
        require(frobenius2(difference) == F(9, 10000) * h**4,
                'Defect does not have claimed quadratic norm scaling')
        refinements.append({'h': str(h), 'frobenius_squared': str(frobenius2(difference))})
    return {'passed': True, 'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'scope': 'Exact positive-real increasing line branch only; not a FLASH execution or full PDE convergence test',
            'independent_review': 'unavailable; reviewer returned a usage-limit error',
            'face_defect': [[str(v) for v in row] for row in defect],
            'frobenius_squared': '1/2000', 'purity_change': '29/2000',
            'smooth_grid_checks': refinements,
            'conclusion': 'Finite-stencil basis dependence; smooth example defect is O(h^2); no new algorithm or physical-impact claim'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    text = json.dumps(run(), indent=2) + '\n'
    if args.output:
        args.output.write_text(text)
    print(text, end='')
