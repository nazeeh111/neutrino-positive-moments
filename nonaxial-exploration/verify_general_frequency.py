"""Standalone certificate verifier: standard library, no proposer imports.

Verifies exclusions at exact rational frequency points, not measurement boxes.
Separate implementation by the same authoring agent is not independent review.
"""
from fractions import Fraction
from itertools import combinations
import hashlib
import json
import sys


def rational(value):
    if type(value) is not int and not isinstance(value, str):
        raise ValueError('Rationals must be JSON integers or strings')
    return Fraction(value)


def real_det(matrix):
    if len(matrix) == 1:
        return matrix[0][0]
    return sum((-1)**j*matrix[0][j]*real_det([row[:j]+row[j+1:] for row in matrix[1:]])
               for j in range(len(matrix)))


def parse(data):
    if not isinstance(data, dict) or set(data) != {'plus', 'minus', 'k', 'wr', 'gamma'}:
        raise ValueError('Unexpected input fields')
    out = {name: rational(data[name]) for name in ('k', 'wr', 'gamma')}
    if out['k'] < 0 or out['gamma'] <= 0:
        raise ValueError('Require k>=0 and gamma>0')
    for name in ('plus', 'minus'):
        raw = data[name]
        if not isinstance(raw, dict) or set(raw) != {'N', 'F', 'Q'}:
            raise ValueError('Invalid species fields')
        if not isinstance(raw['F'], list) or not isinstance(raw['Q'], list) or len(raw['F']) != 3 or len(raw['Q']) != 3:
            raise ValueError('Invalid species shape')
        if any(not isinstance(row, list) or len(row) != 3 for row in raw['Q']):
            raise ValueError('Invalid tensor shape')
        N = rational(raw['N'])
        f = list(map(rational, raw['F']))
        q = [list(map(rational, row)) for row in raw['Q']]
        if q != [list(row) for row in zip(*q)] or sum(q[i][i] for i in range(3)) != N:
            raise ValueError('Invalid tensor symmetry or trace')
        gram = [[N, *f]]+[[f[i], *q[i]] for i in range(3)]
        for size in range(1, 5):
            for ids in combinations(range(4), size):
                if real_det([[gram[i][j] for j in ids] for i in ids]) < 0:
                    raise ValueError('Nonpositive moment Gram matrix')
        out[name] = {'N': N, 'F': f, 'Q': q}
    return out


def canonical(value):
    if isinstance(value, Fraction):
        return str(value)
    if isinstance(value, list):
        return [canonical(v) for v in value]
    if isinstance(value, dict):
        return {k: canonical(v) for k, v in value.items()}
    return value


def digest(data):
    return hashlib.sha256(json.dumps(canonical(data), sort_keys=True, separators=(',', ':')).encode()).hexdigest()


# Exact complex pairs. These do not share arithmetic with the proposer.
def c(value):
    return value if isinstance(value, tuple) else (Fraction(value), Fraction(0))


def ca(left, right):
    a, b = c(left), c(right)
    return a[0]+b[0], a[1]+b[1]


def cn(value):
    a = c(value)
    return -a[0], -a[1]


def cm(left, right):
    a, b = c(left), c(right)
    return a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0]


def cd(left, right):
    a, b = c(left), c(right)
    norm = b[0]*b[0]+b[1]*b[1]
    if not norm:
        raise ZeroDivisionError('Complex denominator zero')
    return cm(a, (b[0]/norm, -b[1]/norm))


def complex_det(matrix):
    if len(matrix) == 1:
        return c(matrix[0][0])
    result = c(0)
    for j in range(len(matrix)):
        minor = [row[:j]+row[j+1:] for row in matrix[1:]]
        result = ca(result, cm((-1)**j, cm(matrix[0][j], complex_det(minor))))
    return result


# Rational interval pairs, with separate operations from complex pairs.
def box(value):
    return value if isinstance(value, tuple) else (Fraction(value), Fraction(value))


def ia(left, right):
    a, b = box(left), box(right)
    return a[0]+b[0], a[1]+b[1]


def inn(value):
    a = box(value)
    return -a[1], -a[0]


def im(left, right):
    a, b = box(left), box(right)
    products = [x*y for x in a for y in b]
    return min(products), max(products)


def inv(value):
    a = box(value)
    if a[0] <= 0 <= a[1]:
        raise ZeroDivisionError('Interval includes zero')
    return 1/a[1], 1/a[0]


def iq(value):
    a = box(value)
    return (Fraction(0) if a[0] <= 0 <= a[1] else min(a[0]**2, a[1]**2), max(a[0]**2, a[1]**2))


def reconstruct_lower(data, p):
    plus, minus = data['plus'], data['minus']
    k, wr, gamma = (data[name] for name in ('k', 'wr', 'gamma'))
    x, y = wr/k, gamma/k
    omega, r = (wr, gamma), (x, y)
    g0 = plus['N']-minus['N']
    g = [a-b for a, b in zip(plus['F'], minus['F'])]
    Omega = ca(omega, g0)
    Delta = ca(k*(k+g[2]), cn(cm(omega, Omega)))
    A = cd(ca(k*(k+g[2]), cm(g0, Omega)), Delta)
    q = im(-x*y, inv(p))
    h2 = ia(iq(p), x*x)
    ih2 = inv(h2)
    alpha = im(p, ia(2, inn(ih2)))
    beta = im(-y, inv(p))
    kp = im(iq(ia(p, inn(q))), ih2)
    km = im(iq(ia(p, q)), ih2)
    ep, em = ia(ia(kp, -1), ih2), ia(ia(km, -1), ih2)
    L = [[box(0) for _ in range(3)] for _ in range(3)]
    base00 = cm(ca(1, cn(cm(r, r))), A)[0]-(x*g0+g[2])/k
    L[0][0] = box(base00-(plus['N']+minus['N']-plus['Q'][2][2]-minus['Q'][2][2])/(2*gamma))
    for i in range(2):
        cross = cd(cm(Omega, g[i]), Delta)
        base0i = ia(im(q, cross[0]), im(p, cross[1]))
        weighted0i = ia(im(alpha, g[i]), im(beta, plus['Q'][2][i]+minus['Q'][2][i]))
        L[0][i+1] = L[i+1][0] = ia(base0i, im(-1/(2*gamma), weighted0i))
        for j in range(2):
            baseij = ca(1 if i == j else 0, cn(cd(g[i]*g[j], Delta)))[0]
            weightedij = ia(im(ep, plus['Q'][i][j]), im(em, minus['Q'][i][j]))
            if i == j:
                weightedij = ia(weightedij, ia(im(kp, plus['Q'][2][2]), im(km, minus['Q'][2][2])))
            L[i+1][j+1] = ia(baseij, im(-1/(2*gamma), weightedij))
    return L


def positive_pivots(L):
    # Interval LDL elimination, rather than the proposer's determinants.
    d0 = L[0][0]
    if d0[0] <= 0:
        return False
    d1 = ia(L[1][1], inn(im(iq(L[0][1]), inv(d0))))
    if d1[0] <= 0:
        return False
    off = ia(L[1][2], inn(im(im(L[0][1], L[0][2]), inv(d0))))
    d2 = ia(ia(L[2][2], inn(im(iq(L[0][2]), inv(d0)))), inn(im(iq(off), inv(d1))))
    return d2[0] > 0


def verify(raw, certificate):
    try:
        data = parse(raw)
        if certificate['input_sha256'] != digest(data):
            return False
        if data['k'] == 0:
            if set(certificate) != {'kind', 'input_sha256'} or certificate['kind'] != 'general-3d-zero-k-v1':
                return False
            plus, minus = data['plus'], data['minus']
            gp = [[plus['N'], *plus['F']]]+[[plus['F'][i], *plus['Q'][i]] for i in range(3)]
            gm = [[minus['N'], *minus['F']]]+[[minus['F'][i], *minus['Q'][i]] for i in range(3)]
            omega = data['wr'], data['gamma']
            Pi = [[ca(1 if i == j == 0 else -1 if i == j else 0, cd(gp[i][j]-gm[i][j], omega))
                   for j in range(4)] for i in range(4)]
            return complex_det(Pi) != c(0)
        if set(certificate) != {'kind', 'input_sha256', 's_bounds', 'p_bounds'}:
            return False
        if certificate['kind'] != 'general-3d-rational-point-v1':
            return False
        if not isinstance(certificate['s_bounds'], list) or not isinstance(certificate['p_bounds'], list):
            return False
        if len(certificate['s_bounds']) != 2 or len(certificate['p_bounds']) != 2:
            return False
        sl, su = map(rational, certificate['s_bounds'])
        pl, pu = map(rational, certificate['p_bounds'])
        x, y = data['wr']/data['k'], data['gamma']/data['k']
        v = 1-x*x+y*y
        discriminant = v*v+4*x*x*y*y
        if not (0 <= sl <= su and sl*sl <= discriminant <= su*su and 0 < pl <= pu):
            return False
        if pl*pl > (sl+v)/2 or pu*pu < (su+v)/2:
            return False
        return positive_pivots(reconstruct_lower(data, (pl, pu)))
    except (KeyError, TypeError, ValueError, ZeroDivisionError, IndexError):
        return False


if __name__ == '__main__':
    if len(sys.argv) != 3:
        raise SystemExit('Usage: verify_general_frequency.py INPUT.json CERTIFICATE_OR_RESULT.json')
    with open(sys.argv[1]) as f:
        data = json.load(f)
    with open(sys.argv[2]) as f:
        candidate = json.load(f)
    valid = verify(data, candidate.get('certificate', candidate))
    print(json.dumps({'valid_exclusion_certificate': valid}))
    raise SystemExit(0 if valid else 1)
