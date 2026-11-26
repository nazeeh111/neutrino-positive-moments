"""Read-only recipient report for selected saved exact certificates.

The local code and manifest are trusted inputs; hashes are not authentication.
No manifest commands are executed and numerical evidence is not reproduced.
"""
import argparse
import hashlib
import importlib.util
import json
import os
import stat
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from independent_verifier import verify as verify_m1
from independent_m2_verifier import verify as verify_m2, verify_m1_possible


# Fixed checker and data paths are part of the admitted integrity coverage.
DATA_FILES = frozenset({
    'results/screen.json', 'results/stress.json', 'results/continuous-baseline.json',
    'results/m2-counterflux-screen.json', 'examples/boundary-proof.json',
    'examples/m2-counterflux-input.json', 'examples/m2-counterflux-proof.json',
    'examples/nonaxial-general-frequency-input.json', 'examples/nonaxial-general-frequency-result.json',
})
CHECKER_FILES = frozenset({'verify_recipient.py', 'independent_verifier.py',
    'independent_m2_verifier.py', 'nonaxial-exploration/verify_general_frequency.py'})
MAX_FILE = 8 * 1024 * 1024
MAX_TOTAL = 64 * 1024 * 1024
MAX_MANIFEST = 1024 * 1024
SNAPSHOT = {}


def read_regular(name, limit):
    relative = Path(name)
    if relative.is_absolute() or '..' in relative.parts or not relative.parts:
        raise ValueError('unsafe manifest path: ' + name)
    path = ROOT / relative
    # Reject symlinks in every component; do not follow the final component when
    # the host exposes O_NOFOLLOW. Metadata checks are not an adversarial sandbox.
    current = ROOT
    for part in relative.parts:
        current = current / part
        if current.is_symlink():
            raise ValueError('symlink in manifest path: ' + name)
    flags = os.O_RDONLY | getattr(os, 'O_NONBLOCK', 0) | getattr(os, 'O_NOFOLLOW', 0)
    fd = os.open(path, flags)
    try:
        info = os.fstat(fd)
        if not stat.S_ISREG(info.st_mode):
            raise ValueError('not a regular file: ' + name)
        if info.st_size > limit:
            raise ValueError('file exceeds size limit: ' + name)
        chunks, total = [], 0
        while True:
            chunk = os.read(fd, min(65536, limit + 1 - total))
            if not chunk:
                break
            total += len(chunk)
            if total > limit:
                raise ValueError('file exceeds size limit: ' + name)
            chunks.append(chunk)
        return b''.join(chunks)
    finally:
        os.close(fd)


def load(name):
    return json.loads(SNAPSHOT[name])


def verify_bytes():
    manifest = json.loads(read_regular('BUNDLE_MANIFEST.json', MAX_MANIFEST))
    if manifest['format'] != 'sha256-file-manifest-v1':
        raise ValueError('unsupported manifest format')
    files = manifest['files']
    missing = (DATA_FILES | CHECKER_FILES) - files.keys()
    if missing:
        raise ValueError('required integrity entries missing: ' + ', '.join(sorted(missing)))
    if len(files) > 10000:
        raise ValueError('too many manifest entries')
    count, total = 0, 0
    captured = {}
    for name, record in files.items():
        data = read_regular(name, min(MAX_FILE, MAX_TOTAL - total))
        total += len(data)
        if len(data) != record['bytes'] or hashlib.sha256(data).hexdigest() != record['sha256']:
            raise ValueError('manifest mismatch: ' + name)
        if name in DATA_FILES:
            captured[name] = data
        count += 1
    SNAPSHOT.clear()
    SNAPSHOT.update(captured)
    return count


def exact_report():
    counts = {'axial_screen_stress': 0, 'axial_numerical_result_certificates': 0,
              'axial_continuous_certificates': 0, 'axial_example': 0,
              'm2_counterflux': 0, 'm1_counterflux': 0, 'm2_example': 0,
              'general_nonaxial_example': 0}
    skipped = {'continuous_unresolved': 0, 'screen_unresolved': 0}
    for filename in ('results/screen.json', 'results/stress.json'):
        for row in load(filename)['cases']:
            verify_m1(row['problem'], row['result'])
            counts['axial_screen_stress'] += 1
            numeric = row.get('numeric_result')
            if numeric:
                if numeric['status'] == 'unresolved':
                    skipped['screen_unresolved'] += 1
                else:
                    verify_m1(row['problem'], numeric)
                    counts['axial_numerical_result_certificates'] += 1
    example = load('examples/boundary-proof.json')
    verify_m1(example['problem'], example['result'])
    counts['axial_example'] += 1
    for row in load('results/continuous-baseline.json')['cases']:
        for method in ('exact', 'conic'):
            result = row[method]
            if result['status'] == 'unresolved':
                skipped['continuous_unresolved'] += 1
            else:
                verify_m1(row['problem'], result)
                counts['axial_continuous_certificates'] += 1
    data = load('results/m2-counterflux-screen.json')
    entries = [row for row in data['closure_modes'] if 'result' in row]
    entries += data['continuous_designated_point_controls']
    for row in data['validated_linear_bracket_endpoint_checks']:
        entries += row['endpoint_frequency_checks']
    for entry in entries:
        verify_m2(entry['data'], entry['result'])
        counts['m2_counterflux'] += 1
        if entry.get('m1_result', {}).get('status') == 'possible':
            verify_m1_possible(entry['data'], entry['m1_result'])
            counts['m1_counterflux'] += 1
    verify_m2(load('examples/m2-counterflux-input.json'), load('examples/m2-counterflux-proof.json'))
    counts['m2_example'] += 1
    path = ROOT / 'nonaxial-exploration/verify_general_frequency.py'
    spec = importlib.util.spec_from_file_location('recipient_general_frequency', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    candidate = load('examples/nonaxial-general-frequency-result.json')
    if not module.verify(load('examples/nonaxial-general-frequency-input.json'), candidate['certificate']):
        raise ValueError('general nonaxial example certificate refused')
    counts['general_nonaxial_example'] += 1
    return counts, skipped


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, help='write report exclusively to a new file')
    args = parser.parse_args()
    report = {'format': 'recipient-verification-v1',
              'scope': 'Selected saved certificates use captured verified bytes. Local checker code and manifest are trusted; hashes do not authenticate origin or make concurrent code changes safe.',
              'byte_integrity': {'status': 'not_checked', 'checked_files': 0,
                  'required_data_files': len(DATA_FILES), 'required_checker_files': len(CHECKER_FILES),
                  'max_file_bytes': MAX_FILE, 'max_total_bytes': MAX_TOTAL},
              'exact_certificates': {'status': 'not_checked', 'checked': {}, 'skipped': {}},
              'unsupported_families': [
                  'grid witness and Farkas rows',
                  'nonaxial tangent, imaginary-bound, full-frequency and general-moment theorem artifacts',
                  'finite-band, second-moment and other theorem and attainer validation artifacts',
                  'numerical root brackets, timings, benchmarks, figures, PDF and physical claims'],
              'unsupported_family_count': 4,
              'numerical_evidence': 'unchecked; no numerical benchmarks or root finding rerun'}
    try:
        report['byte_integrity']['checked_files'] = verify_bytes()
        report['byte_integrity']['status'] = 'passed'
        counts, skipped = exact_report()
        report['exact_certificates'].update(status='passed', checked=counts, skipped=skipped)
    except Exception as exc:
        section = 'byte_integrity' if report['byte_integrity']['status'] != 'passed' else 'exact_certificates'
        report[section].update(status='failed', error=f'{type(exc).__name__}: {exc}')
    rendered = json.dumps(report, indent=2, sort_keys=True) + '\n'
    if args.output:
        with args.output.open('x', encoding='utf-8') as out:
            out.write(rendered)
    else:
        print(rendered, end='')
    return 0 if report['exact_certificates']['status'] == 'passed' else 1


if __name__ == '__main__':
    raise SystemExit(main())
