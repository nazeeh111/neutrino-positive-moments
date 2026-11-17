# An upper bound, with its witness

The written finite-band theorem optimizes over positive axial measures whose number moments are `(N, 0, Q)`, both species equal, and every real frequency part with `0 ≤ k ≤ K`. A single witness establishes attainability, not the universal upper bound. This case is the existing README example and endpoint construction in [the retained proof](../../finite-k-transition-audit.md#6-explicit-attainers-in-every-regime), checked by the original-factor method in [the existing tests](../../test_finite_band_envelope.py).

For `N = 1`, `Q = 1/3`, `K = 1/1000`, the evaluator selects the endpoint branch and returns exactly `γ² = 1333/3000000`. The growing frequency is `ω = i√(1333/3000000)` at `k = K`.

| Positive species | Angular cosine u | Atom weight |
| --- | --- | --- |
| μ+ | −1 | 1/4 |
| μ+ | 1/3 | 3/4 |
| μ− | −1/3 | 3/4 |
| μ− | 1 | 1/4 |

Each species has mass 1, first moment 0 and second number moment 1/3. The native exact complex arithmetic evaluates the original longitudinal dispersion factor `(1+I₀)(−1+I₂)−I₁²` to zero. A separate rational endpoint calculation reproduces the squared growth. The saved rational enclosure brackets its positive square root; it is not a statistical confidence interval.

[Exact input and four atoms](case.json) · [Unmodified evaluator output](evaluator.json) · [Actual check output](checked.json) · [Wide SVG](witness.svg) · [Phone SVG](witness-mobile.svg) · [Figure source](figure.py) · [Offline readable case](index.html)

## Replay one finite result

From the repository root, with Python 3 and no scientific packages:

```sh
python3 -I -S -B -O docs/presentation/check.py
python3 -I -S -B finite_band_envelope.py --mass 1 --second 1/3 --kmax 1/1000
```

The small wrapper hashes its three original kernel files, recomputes the exact evaluator result, checks positive support and all supplied moments, checks the original factor, and rejects altered atom weights or squared growth. Its saved source binding identifies this delivered computation; it is not an author signature or independent scientific review. The code and local hashes are trusted. No scan, optimizer, new theorem, root-finding campaign or physical measurement is involved.

[Download the self-contained example](worked-example.zip). Extract into a new directory and run `python3 -I -S -B -O check.py`. The archive includes its own short guide, all three native modules, exact input/output, figure sources and notices. Reopening or checking writes no research output. No full repository or installed package is required.

## Interpret the units and quantifiers

`N` is self-interaction-weighted number mass relative to a reference density. Frequencies and wave numbers are reduced by the reference interaction frequency; `k` includes the factor of light speed. This dimensionless example assigns no astrophysical wavelength or growth time. See [the original physical scope](../../PHYSICAL_SCOPE.md).

The classifier and finite-band bound assume positive axisymmetric angular measures and propagation along the symmetry axis, local collisionless fast evolution and explicit pointwise heavy-flavor cancellation. A possible growing completion does not predict the actual background or imply every completion is unstable. Atomic attainment does not establish smooth or occupation-constrained attainment. Different inputs on the envelope may use different completions. The general unequal-moment finite-band problem remains open.

The later 3D sufficient-exclusion interface is a different question, with separate assumptions and review status. Its `unresolved` answer is not a compatibility or stability decision. The original guide retains that method, the 5.72% strongest-cap comparison, failed counter-flux usefulness gate, unknown priority and absent demonstrated operational astrophysical benefit. The PDF still says working draft, 2 November 2026, not submitted or peer reviewed. No manuscript or historical data was changed for this example.
