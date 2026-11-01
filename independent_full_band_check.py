"""Reproducible exact geometry-to-root audit of the finite-band evaluator.

The reference branch uses containment of horizontal disk diameters, rather
than the production endpoint-inversion selector. This finite audit checks
implementation behavior; the universal argument is in the written review.
Standard library only. Writes only results/independent-full-band-validation.json.
"""
from collections import Counter
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import sys

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parent
EXPECTED_SHA256 = {
    'finite_band_envelope.py':
        'e0fd9e4dce56dbb1c115690e8d7f14f4a419e6620cac2906d61886a448a9ddb5',
    'finite-k-transition-audit.md':
        'f076d0fb272d073810fe3aaafede69d331c998abf787d6a94790047760a223e3',
}


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def geometric_width(a, y):
    """Width at N=1 from the two imaginary-axis disks, with exact ties."""
    r, t = 1-a, y*y
    radius0_sq = r*r/(4*t)
    radius1_sq = a*a*(1+t)**2/(4*t*(a+t)**2)
    distance_sq = (a*a+t)**2/(4*t*(a+t)**2)
    if radius0_sq + distance_sq <= radius1_sq:
        return 'plateau', r/y
    if radius1_sq + distance_sq < radius0_sq:
        return 'cubic', a*(1+t)/(y*(a+t))
    # In the circle-intersection case both horizontal endpoints are feasible.
    # At a Gram/intersection tie this matches the evaluator's endpoint label.
    return 'endpoint', 2*a*r/(a*a+t)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def run():
    observed = {name: sha256(ROOT/name) for name in EXPECTED_SHA256}
    report = {
        'passed': False,
        'scope': 'Finite exact geometry-to-root audit; not a universal theorem proof',
        'reference_branch': 'Squared disk-radius and center-distance containment comparisons',
        'grid': {'a': 'i/100, i=1,...,99', 'y': 'j/50, j=1,...,100',
                 'mass': '7/3', 'sectors': ['longitudinal', 'transverse']},
        'enclosure_bits': 30,
        'expected_source_sha256': EXPECTED_SHA256,
        'observed_source_sha256': observed,
        'checker_sha256': sha256(Path(__file__).resolve()),
        'python_version': sys.version.split()[0],
        'completed_cases': 0,
        'branch_counts': {},
    }
    counts = Counter()
    try:
        require(observed == EXPECTED_SHA256, 'Pinned source hash mismatch')
        from finite_band_envelope import envelope
        mass = F(7, 3)
        for ia in range(1, 100):
            a = F(ia, 100)
            for jy in range(1, 101):
                y = F(jy, 50)
                branch, width = geometric_width(a, y)
                for c, sector in ((1, 'longitudinal'), (2, 'transverse')):
                    wave_number = mass*width/c
                    exact_growth = wave_number*y
                    result = envelope(mass, mass*a, wave_number, sector, bits=30)
                    context = f'a={a}, y={y}, sector={sector}'
                    require(result['branch'] == branch,
                            f'Branch disagreement at {context}: '
                            f"expected {branch}, got {result['branch']}")
                    lower, upper = F(result['growth_lower']), F(result['growth_upper'])
                    require(lower <= exact_growth <= upper,
                            f'Exact root outside enclosure at {context}')
                    counts[branch] += 1
                    report['completed_cases'] += 1
        require(report['completed_cases'] == 19800, 'Unexpected case count')
        require(dict(counts) == {'cubic': 1086, 'endpoint': 13814, 'plateau': 4900},
                'Unexpected branch coverage')
        report['passed'] = True
    except Exception as error:
        report['error'] = f'{type(error).__name__}: {error}'
    report['branch_counts'] = dict(sorted(counts.items()))
    output = ROOT/'results'/'independent-full-band-validation.json'
    output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))
    return 0 if report['passed'] else 1


if __name__ == '__main__':
    raise SystemExit(run())
