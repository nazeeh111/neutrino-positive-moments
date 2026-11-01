#!/bin/sh
set -eu
cd "$(dirname "$0")"
export PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
positive_python="${POSITIVE_MOMENT_PYTHON:-python3}"
"$positive_python" -m unittest test_independent_verifier test_sharp_envelope test_second_moment_bounds -v
"$positive_python" screen.py
"$positive_python" stress.py
"$positive_python" smooth_check.py
"$positive_python" independent_scipy_smooth.py
"$positive_python" -O verify_all.py
"$positive_python" -O offaxis-independent.py
"$positive_python" -O independent_continuous_replay.py
"$positive_python" -O -m unittest -v test_independent_m2_locus test_m2_independent_integration
"$positive_python" -B -O finite-k-cap-checks.py
"$positive_python" -B -O test_finite_band_envelope.py
"$positive_python" -B -O independent_full_band_check.py
"$positive_python" counterflux_modes.py
"$positive_python" verify_counterflux_reference.py
"$positive_python" m2_counterflux_screen.py
"$positive_python" -O verify_m2_screen.py
