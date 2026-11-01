"""Independent floating-point check of six smooth profiles and four saved roots.

Uses SciPy's Bessel i0 and direct complex dispersion integrals, importing no
production module and no Arb callback. This is a numerical crosscheck, not a
replacement for the saved rigorous ball-arithmetic signs. Run with the local
SciPy runtime or any suitable NumPy/SciPy environment. Writes only its own JSON.
"""
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import platform
import sys

import numpy as np
import scipy
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import i0

ROOT = Path(__file__).resolve().parent


def profile(u, concentration):
    # A 3/4 uniform azimuthal ring at cos(theta)=1/3 and a 1/4 south pole,
    # convolved on the sphere with a normalized vMF kernel, then azimuth-integrated.
    ring_argument = concentration * math.sqrt(max(0.0, 8.0 * (1.0-u*u) / 9.0))
    ring = math.exp(concentration*u/3.0) * i0(ring_argument)
    pole = math.exp(-concentration*u)
    return concentration * (0.75*ring+0.25*pole) / (2.0*math.sinh(concentration))


def integrate_real(callback):
    value, error = quad(callback, -1.0, 1.0, epsabs=2e-13,
                        epsrel=2e-13, limit=300)
    return value, error


def direct_determinant(gamma, concentration):
    # N=k=1; no scalar target or reweighted moment identities are reused.
    omega = complex(0.0, gamma)
    moments = []
    errors = []
    for power in range(3):
        def integrand(u):
            signed_density = profile(u, concentration)-profile(-u, concentration)
            return u**power * signed_density / (omega-u)
        real, re_error = integrate_real(lambda u: integrand(u).real)
        imag, im_error = integrate_real(lambda u: integrand(u).imag)
        moments.append(complex(real, imag))
        errors.append([re_error, im_error])
    j0, j1, j2 = moments
    return (1.0+j0)*(-1.0+j2)-j1*j1, errors


def main():
    saved_path = ROOT/'results'/'smooth.json'
    saved = json.loads(saved_path.read_text())
    rows = []
    for case in saved['cases']:
        concentration = case['kappa']
        values, errors = zip(*(integrate_real(lambda u, j=j:
                              u**j * profile(u, concentration)) for j in range(3)))
        deviation = max(abs(value-target) for value, target in
                        zip(values, (1.0, 0.0, 1.0/3.0)))
        row = {'kappa': concentration, 'plus_moments': list(values),
               'quad_error_estimates': list(errors),
               'maximum_moment_deviation': deviation,
               'minus_moments_by_reflection': [values[0], -values[1], values[2]],
               'saved_status': case['status']}
        if deviation > 2e-11:
            raise ArithmeticError('independent smooth moments disagree')
        if case['status'] == 'certified_growing_mode':
            low, high = (float(Fraction(value)) for value in case['gamma_interval'])
            left, _ = direct_determinant(low, concentration)
            right, _ = direct_determinant(high, concentration)
            if not (left.real > 0 and right.real < 0):
                raise ArithmeticError('independent endpoint signs disagree')
            root = brentq(lambda gamma: direct_determinant(gamma, concentration)[0].real,
                          low, high, xtol=1e-14, rtol=1e-14)
            residual, root_errors = direct_determinant(root, concentration)
            if not (low < root < high and abs(residual) < 2e-11):
                raise ArithmeticError('independent root/residual check failed')
            row.update(gamma=root, saved_bracket=case['gamma_interval'],
                       endpoint_real_determinants=[left.real, right.real],
                       complex_determinant_residual=[residual.real, residual.imag],
                       absolute_residual=abs(residual),
                       root_integral_error_estimates=root_errors)
        rows.append(row)
    result = {'kind': 'independent floating-point crosscheck, not interval certification',
              'python': platform.python_version(), 'executable': sys.executable,
              'numpy': np.__version__, 'scipy': scipy.__version__,
              'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'smooth_json_sha256': hashlib.sha256(saved_path.read_bytes()).hexdigest(),
              'cases': rows, 'moment_profiles_checked': len(rows),
              'roots_checked': sum('gamma' in row for row in rows), 'passed': True}
    (ROOT/'results'/'independent_scipy_smooth.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
