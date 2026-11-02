"""Certified sufficient exclusion at arbitrary rational frequency points.

Standard library only. Outputs excluded with a separately verified certificate,
or unresolved. Does not classify possible frequencies or measurement intervals.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from math import isqrt
from pathlib import Path
import argparse
import hashlib
import json
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from independent_verifier import C
from general_moment_certificate import validate
from verify_general_frequency import verify


@dataclass(frozen=True)
class Interval:
    lo: F
    hi: F

    def __post_init__(self):
        if self.lo > self.hi:
            raise ValueError('Reversed interval')

    def __add__(self, other):
        other = interval(other)
        return Interval(self.lo+other.lo, self.hi+other.hi)
    __radd__ = __add__

    def __neg__(self):
        return Interval(-self.hi, -self.lo)

    def __sub__(self, other):
        return self+-interval(other)

    def __rsub__(self, other):
        return interval(other)+-self

    def __mul__(self, other):
        other = interval(other)
        values = [self.lo*other.lo, self.lo*other.hi, self.hi*other.lo, self.hi*other.hi]
        return Interval(min(values), max(values))
    __rmul__ = __mul__

    def reciprocal(self):
        if self.lo <= 0 <= self.hi:
            raise ZeroDivisionError('Interval denominator includes zero')
        return Interval(1/self.hi, 1/self.lo)

    def __truediv__(self, other):
        return self*interval(other).reciprocal()

    def __rtruediv__(self, other):
        return interval(other)*self.reciprocal()

    def square(self):
        values = [self.lo*self.lo, self.hi*self.hi]
        return Interval(F(0) if self.lo <= 0 <= self.hi else min(values), max(values))


def interval(value):
    return value if isinstance(value, Interval) else Interval(F(value), F(value))


def rational(value):
    if type(value) is not int and not isinstance(value, str):
        raise ValueError('Use JSON integer or rational-string inputs, not floats')
    return F(value)


def parse(raw):
    if not isinstance(raw, dict) or set(raw) != {'plus', 'minus', 'k', 'wr', 'gamma'}:
        raise ValueError('Require plus, minus, k, wr, gamma')
    out = {key: rational(raw[key]) for key in ('k', 'wr', 'gamma')}
    if out['k'] < 0 or out['gamma'] <= 0:
        raise ValueError('Require k>=0 and gamma>0; z is the chosen wave-vector direction')
    for species in ('plus', 'minus'):
        m = raw[species]
        if not isinstance(m, dict) or set(m) != {'N', 'F', 'Q'}:
            raise ValueError('Require N, F, Q for each species')
        if not isinstance(m['F'], list) or not isinstance(m['Q'], list) or any(not isinstance(row, list) for row in m['Q']):
            raise ValueError('Moment vectors and tensor rows must be arrays')
        out[species] = {'N': rational(m['N']), 'F': list(map(rational, m['F'])),
                        'Q': [list(map(rational, row)) for row in m['Q']]}
        validate(out[species])
    return out


def canonical(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, list):
        return [canonical(v) for v in value]
    if isinstance(value, dict):
        return {k: canonical(v) for k, v in value.items()}
    return value


def digest(data):
    return hashlib.sha256(json.dumps(canonical(data), sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def sqrt_bounds(value, bits):
    if value < 0:
        raise ValueError('Negative square-root input')
    num, den = value.numerator, value.denominator
    n, d = isqrt(num), isqrt(den)
    if n*n == num and d*d == den:
        exact = F(n, d)
        return Interval(exact, exact)
    scale = 1 << bits
    lower = isqrt((num*scale*scale)//den)
    return Interval(F(lower, scale), F(lower+1, scale))


def constants(data):
    plus, minus = data['plus'], data['minus']
    k, omega = data['k'], C(data['wr'], data['gamma'])
    g0 = plus['N']-minus['N']
    g = [a-b for a, b in zip(plus['F'], minus['F'])]
    Omega = omega+g0
    Delta = k*(k+g[2])-omega*Omega
    if Delta == C():
        return None
    base = [[C() for _ in range(3)] for _ in range(3)]
    base[0][0] = (k*(k+g[2])+g0*Omega)/Delta
    for i in range(2):
        base[0][i+1] = base[i+1][0] = Omega*g[i]/Delta
        for j in range(2):
            base[i+1][j+1] = (-1 if i == j else 0)+g[i]*g[j]/Delta
    return base, g0, g


def lower_matrix(data, base, g0, g, p):
    k, gamma = data['k'], data['gamma']
    x, y = data['wr']/k, gamma/k
    q = -x*y/p
    h2 = p.square()+x*x
    alpha, beta = p*(2-1/h2), -y/p
    means = []
    for species, sign in (('plus', 1), ('minus', -1)):
        m = data[species]
        kp = (p-sign*q).square()/h2
        ep = kp-1+1/h2
        avg = [[interval(0) for _ in range(3)] for _ in range(3)]
        avg[0][0] = interval(F(3, 2)*(m['N']-m['Q'][2][2]))
        for i in range(2):
            avg[0][i+1] = avg[i+1][0] = F(3, 2)*(alpha*sign*m['F'][i]+beta*m['Q'][2][i])
            for j in range(2):
                avg[i+1][j+1] = F(3, 2)*((kp*m['Q'][2][2] if i == j else 0)+ep*m['Q'][i][j])
        means.append(avg)
    r = C(x, y)
    B = [[interval(0) for _ in range(3)] for _ in range(3)]
    B[0][0] = interval(((1-r*r)*base[0][0]).re-(x*g0+g[2])/k)
    for i in range(2):
        B[0][i+1] = B[i+1][0] = q*base[0][i+1].re+p*base[0][i+1].im
        for j in range(2):
            B[i+1][j+1] = interval(-base[i+1][j+1].re)
    return [[B[i][j]-(means[0][i][j]+means[1][i][j])/(3*gamma) for j in range(3)] for i in range(3)]


def interval_det(matrix):
    if len(matrix) == 1:
        return matrix[0][0]
    return sum(((-1)**j*matrix[0][j]*interval_det([row[:j]+row[j+1:] for row in matrix[1:]])
                for j in range(len(matrix))), interval(0))


def propose(raw, max_bits=512):
    if type(max_bits) is not int or not 1 <= max_bits <= 4096:
        raise ValueError('max_bits must be an integer in [1,4096]')
    data = parse(raw)
    binding = digest(data)
    if data['k'] == 0:
        certificate = {'kind': 'general-3d-zero-k-v1', 'input_sha256': binding}
        if verify(raw, certificate):
            return {'status': 'excluded', 'certificate': certificate}
        return {'status': 'unresolved', 'reason': 'zero-k polarization is singular'}
    known = constants(data)
    if known is None:
        return {'status': 'unresolved', 'reason': 'singular reduction Delta=0'}
    x, y = data['wr']/data['k'], data['gamma']/data['k']
    v = 1-x*x+y*y
    disc = v*v+4*x*x*y*y
    bits = min(32, max_bits)
    attempts = []
    while True:
        attempts.append(bits)
        s = sqrt_bounds(disc, bits)
        low_square, high_square = max(F(0), (s.lo+v)/2), (s.hi+v)/2
        p = Interval(sqrt_bounds(low_square, bits).lo, sqrt_bounds(high_square, bits).hi)
        if p.lo > 0:
            L = lower_matrix(data, *known, p)
            if all(interval_det([row[:j] for row in L[:j]]).lo > 0 for j in (1, 2, 3)):
                certificate = {'kind': 'general-3d-rational-point-v1', 'input_sha256': binding,
                               's_bounds': [str(s.lo), str(s.hi)], 'p_bounds': [str(p.lo), str(p.hi)]}
                if verify(raw, certificate):
                    return {'status': 'excluded', 'attempt_bits': attempts, 'certificate': certificate}
        if bits == max_bits:
            return {'status': 'unresolved', 'attempt_bits': attempts,
                    'reason': 'positive lower matrix not certified; no compatibility conclusion'}
        bits = min(bits*2, max_bits)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path)
    parser.add_argument('--max-bits', type=int, default=512)
    args = parser.parse_args()
    print(json.dumps(propose(json.loads(args.input.read_text()), args.max_bits), indent=2))
