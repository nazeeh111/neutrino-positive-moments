"""Sharp axial growth caps with independently fixed species second moments.

Only Python's standard library is needed. All returned quantities are exact
Fractions, including gamma_squared; gamma itself can be irrational. Positive
measure atoms encode number measures, not their (1-u**2)-weighted versions.
See independent-second-moment.md for the assumptions and proof.
"""

from dataclasses import asdict, dataclass
from fractions import Fraction
import argparse
import json


@dataclass(frozen=True)
class Atom:
    u: Fraction
    weight: Fraction


@dataclass(frozen=True)
class Attainer:
    k: Fraction
    omega_real: Fraction
    gamma_squared: Fraction
    plus_atoms: tuple[Atom, ...]
    minus_atoms: tuple[Atom, ...]


@dataclass(frozen=True)
class BranchBound:
    gamma_squared: Fraction
    attainer: Attainer | None


@dataclass(frozen=True)
class SecondMomentBounds:
    n: Fraction
    q_plus: Fraction
    q_minus: Fraction
    r_plus: Fraction
    r_minus: Fraction
    longitudinal: BranchBound
    transverse: BranchBound
    global_gamma_squared: Fraction
    stable_reason: str | None


def _atoms(n: Fraction, q: Fraction, sign: int) -> tuple[Atom, ...]:
    """Zero-flux positive measure with mass n and second moment n*q.

    The zero-weight endpoint is omitted when q=0. This retains a single
    atom at zero rather than an artificial endpoint or duplicate node.
    """
    nodes = (Atom(sign * q, n / (1 + q)),
             Atom(Fraction(-sign), n * q / (1 + q)))
    return tuple(atom for atom in nodes if atom.weight)


def second_moment_bounds(n, q_plus, q_minus) -> SecondMomentBounds:
    """Return sharp caps and attainers for masses n, fluxes 0, moments Q±.

    Pass int, Fraction or a rational/decimal string for exact intended inputs.
    A float is interpreted as its exact binary value, as by Fraction(float).
    Stable controls have zero caps and no UHP attainer. These are actual
    stable controls, not a claim that a root at gamma=0 is attained.
    """
    n, qp, qm = map(Fraction, (n, q_plus, q_minus))
    if n <= 0:
        raise ValueError("n must be strictly positive")
    if not (0 <= qp <= n and 0 <= qm <= n):
        raise ValueError("each second number moment must lie in [0,n]")
    rp, rm = n - qp, n - qm
    reason = None
    if qp == qm == 0:
        reason = "Both species are uniquely n*delta_0, so their signed measure vanishes."
    elif rp * rm == 0:
        reason = "At least one weighted species has zero mass; no upper-half-plane root exists."
    if reason:
        zero = BranchBound(Fraction(0), None)
        return SecondMomentBounds(n, qp, qm, rp, rm, zero, zero,
                                  Fraction(0), reason)

    ap, am = qp / n, qm / n
    s = ap + am
    k = n * (2 - s) / s
    x = n * (ap - am) / s
    gamma2 = rp * rm
    long = Attainer(k, x, gamma2, _atoms(n, ap, 1), _atoms(n, am, -1))
    trans = Attainer(k / 2, -x / 2, gamma2 / 4,
                     _atoms(n, ap, -1), _atoms(n, am, 1))
    return SecondMomentBounds(n, qp, qm, rp, rm,
                              BranchBound(gamma2, long),
                              BranchBound(gamma2 / 4, trans), gamma2, None)


def independent_r_bounds(n, r_plus, r_minus) -> SecondMomentBounds:
    """Equivalent interface with independent R± = n-Q± weighted masses."""
    n, rp, rm = map(Fraction, (n, r_plus, r_minus))
    return second_moment_bounds(n, n - rp, n - rm)


def _jsonable(value):
    if isinstance(value, Fraction):
        return str(value)
    if isinstance(value, dict):
        return {key: _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("n", help="equal species number mass")
    parser.add_argument("q_plus", help="plus species second number moment")
    parser.add_argument("q_minus", help="minus species second number moment")
    args = parser.parse_args()
    result = second_moment_bounds(args.n, args.q_plus, args.q_minus)
    print(json.dumps(_jsonable(asdict(result)), indent=2))


if __name__ == "__main__":
    main()
