"""Exact API, compatibility-control and adversarial certificate checks."""
from copy import deepcopy
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import random
import shutil
import subprocess
import sys
import tempfile
from general_frequency_certificate import propose, parse, constants, lower_matrix, Interval, sqrt_bounds, canonical
from general_moment_certificate import moments, evaluate
from verify_general_frequency import verify, reconstruct_lower
from check_general_moment_certificate import original_polarization
from exact_nonaxial_imaginary_bound import determinant
from independent_verifier import C


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def input_data(plus, minus, k, wr, gamma):
    return canonical({'plus': moments(plus), 'minus': moments(minus), 'k': k, 'wr': wr, 'gamma': gamma})


def run():
    root = Path(__file__).resolve().parents[1]
    fixture = json.loads((root/'examples/nonaxial-general-frequency-input.json').read_text())
    result = propose(fixture)
    require(result['status'] == 'excluded' and verify(fixture, result['certificate']), 'Irrational-root fixture failed')
    require(result['certificate']['p_bounds'][0] != result['certificate']['p_bounds'][1], 'Fixture lost irrational square root')
    require(propose(fixture, max_bits=1)['status'] == 'unresolved', 'Coarse enclosure should remain unresolved')
    rejected = 0
    for field, value in [('s_bounds', ['0', '0']), ('p_bounds', ['2', '3']),
                         ('p_bounds', ['0', '2']), ('p_bounds', ['2', '1']),
                         ('kind', 'untrusted-success'), ('input_sha256', '0'*64),
                         ('s_bounds', '01')]:
        cert = deepcopy(result['certificate'])
        cert[field] = value
        require(not verify(fixture, cert), 'Forged certificate accepted: '+field)
        rejected += 1
    changed = deepcopy(fixture)
    changed['wr'] = '1/5'
    require(not verify(changed, result['certificate']), 'Certificate reused for a different input')
    rejected += 1
    equivalent = deepcopy(fixture)
    equivalent['k'] = '198/50'
    require(verify(equivalent, result['certificate']), 'Equivalent rational encoding rejected')
    for change in ('float', 'bool', 'negative_tensor', 'wrong_trace', 'negative_gamma'):
        bad = deepcopy(fixture)
        if change == 'float': bad['wr'] = 0.1
        if change == 'bool': bad['plus']['N'] = True
        if change == 'negative_tensor': bad['plus']['Q'][0][0], bad['plus']['Q'][1][1] = '-1', '8/5'
        if change == 'wrong_trace': bad['plus']['Q'][0][0] = '1'
        if change == 'negative_gamma': bad['gamma'] = '-1'
        try:
            propose(bad)
        except ValueError:
            pass
        else:
            raise RuntimeError('Invalid input accepted: '+change)
        require(not verify(bad, result['certificate']), 'Verifier accepted invalid input')
    south, north = (F(0), F(0), F(-1)), (F(0), F(0), F(1))
    plus, minus = [(south, F(1))], [(north, F(1))]
    actual = input_data(plus, minus, F(6, 5), F(0), F(8, 5))
    require(propose(actual)['status'] == 'unresolved', 'Excluded exact growing two-beam root')
    sideways = [( (F(1), F(0), F(0)), F(1))]
    oscillating = input_data(sideways, minus, F(6, 5), F(3, 5), F(4, 5))
    require(determinant(original_polarization(sideways, minus, C(F(3, 5), F(4, 5)), F(6, 5))) == C(),
            'Nonparallel oscillating fixture is not an original root')
    require(propose(oscillating)['status'] == 'unresolved', 'Excluded nonparallel oscillating root')
    singular = input_data([(south, F(25, 24))], [(north, F(25, 24))], F(3, 4), F(0), F(1))
    require(propose(singular)['reason'] == 'singular reduction Delta=0', 'Wrong singular branch')
    zero_k_root = input_data(plus, minus, F(0), F(0), F(2))
    zero_k_nonroot = input_data(plus, minus, F(0), F(0), F(3))
    require(propose(zero_k_root)['status'] == 'unresolved', 'Excluded exact zero-k root')
    zk = propose(zero_k_nonroot)
    require(zk['status'] == 'excluded' and verify(zero_k_nonroot, zk['certificate']), 'Zero-k exact determinant certificate failed')
    empty = input_data([], [], F(1), F(10), F(1, 10**30))
    refined = propose(empty)
    require(refined['status'] == 'excluded' and max(refined['attempt_bits']) > 32, 'Cancellation refinement not exercised')
    require(verify(empty, refined['certificate']), 'Refined tiny-positive-root certificate rejected')
    for scale in (F(1, 10**30), F(1), F(10**30)):
        data = deepcopy(fixture)
        for species in ('plus', 'minus'):
            data[species]['N'] = str(F(data[species]['N'])*scale)
            data[species]['F'] = [str(F(x)*scale) for x in data[species]['F']]
            data[species]['Q'] = [[str(F(x)*scale) for x in row] for row in data[species]['Q']]
        for key in ('k', 'wr', 'gamma'):
            data[key] = str(F(data[key])*scale)
        scaled = propose(data)
        require(scaled['status'] == 'excluded' and verify(data, scaled['certificate']), 'Common scaling failed')
    sqrt_cases = 0
    for numerator in range(17):
        for denominator in range(1, 18):
            for bits in (1, 7):
                value = F(numerator, denominator)
                enclosure = sqrt_bounds(value, bits)
                require(0 <= enclosure.lo <= enclosure.hi and enclosure.lo**2 <= value <= enclosure.hi**2,
                        'Invalid exact square-root enclosure')
                sqrt_cases += 1
    grid = [(F(a, 3), F(b, 3), F(c, 3)) for a, b, c in
            [(2,2,1),(-2,2,1),(2,-2,1),(-2,-2,1),(2,1,-2),(-2,1,-2),(1,2,-2),(-1,-2,2)]]
    rng = random.Random(2026110211)
    matched = 0
    outcomes = {'excluded': 0, 'unresolved': 0}
    for index in range(18):
        pp = [(v, F(rng.randint(1, 6), 17)) for v in rng.sample(grid, 5)]
        pm = [(v, F(rng.randint(1, 6), 19)) for v in rng.sample(grid, 5)]
        circle = (F(1, 2), F(1), F(2))[index % 3]
        hyper = (F(1, 10), F(1, 3), F(1, 2))[(index//3) % 3]
        a, b = (1-circle*circle)/(1+circle*circle), 2*circle/(1+circle*circle)
        h, t = (1+hyper*hyper)/(1-hyper*hyper), 2*hyper/(1-hyper*hyper)
        k = F(rng.randint(1, 9), 4)
        old = evaluate(moments(pp), moments(pm), k, a, b, h, t)
        require(old['Delta'] != C(), 'Unexpected singular parameter fixture')
        data = input_data(pp, pm, k, k*a*h, k*b*t)
        parsed = parse(data)
        L = lower_matrix(parsed, *constants(parsed), Interval(b*h, b*h))
        independent_L = reconstruct_lower(parsed, (b*h, b*h))
        require([[entry.lo for entry in row] for row in L] == old['lower_matrix']
                and all(entry.lo == entry.hi for row in L for entry in row), 'Old exact certificate mismatch')
        require([[entry[0] for entry in row] for row in independent_L] == old['lower_matrix']
                and all(entry[0] == entry[1] for row in independent_L for entry in row), 'Independent reconstruction mismatch')
        new = propose(data)
        require(new['status'] == old['status'], 'Exact-parameter status changed')
        matched += 1
        # Fresh arbitrary rational frequencies, not tied to rational square roots.
        data['wr'], data['gamma'] = str(F(rng.randint(-9, 9), 7)), str(F(rng.randint(1, 15), 5))
        answer = propose(data)
        outcomes[answer['status']] += 1
        if answer['status'] == 'excluded':
            require(verify(data, answer['certificate']), 'Independent certificate rejection')
            Pi = original_polarization(pp, pm, C(F(data['wr']), F(data['gamma'])), k)
            require(determinant(Pi) != C(), 'Excluded an original-polarization root')
    with tempfile.TemporaryDirectory(prefix='neutrino-verifier-') as directory:
        dest = Path(directory)
        shutil.copy2(root/'nonaxial-exploration/verify_general_frequency.py', dest/'verify.py')
        (dest/'input.json').write_text(json.dumps(fixture))
        (dest/'proof.json').write_text(json.dumps(result))
        command = [sys.executable, '-B', '-O', 'verify.py', 'input.json', 'proof.json']
        standalone = subprocess.run(command, cwd=dest, env={'PYTHONPATH': '', 'PYTHONDONTWRITEBYTECODE': '1'},
                                    capture_output=True, text=True, check=True)
        require(json.loads(standalone.stdout)['valid_exclusion_certificate'], 'Standalone verifier failed')
    return {'passed': True, 'sqrt_enclosure_cases': sqrt_cases, 'exact_parameter_crosschecks': matched,
            'arbitrary_frequency_cases': 18, 'statuses_not_an_accuracy_benchmark': outcomes,
            'forged_or_wrong_input_certificates_rejected': rejected, 'invalid_input_cases': 5,
            'actual_growing_root_controls': 2, 'singular_reduction_control': True, 'zero_k_controls': True,
            'common_scale_cases': 3, 'tiny_root_refinement_bits': refined['attempt_bits'],
            'isolated_stdlib_verifier_passed': True,
            'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'proposer_sha256': hashlib.sha256((root/'nonaxial-exploration/general_frequency_certificate.py').read_bytes()).hexdigest(),
            'verifier_sha256': hashlib.sha256((root/'nonaxial-exploration/verify_general_frequency.py').read_bytes()).hexdigest(),
            'scope': 'Separate implementation by the primary agent; not independent scientific review or a physical-use pass'}


if __name__ == '__main__':
    print(json.dumps(run(), indent=2))
