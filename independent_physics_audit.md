# Independent physics and scientific-value audit

2026-11-01. Read-only audit of TASK_STATE.md, both current plans, the production geometry/moment sources, independent checker, sharp-envelope tests, baselines and canonical screen/smooth results. No 504-point rerun. Applied the research skill to primary-source comparison and direct independent numerical checks. Existing files were not edited.

## Verdict

The restricted result survives this physics attack. I found no sign, azimuthal normalization, transverse factor, smoothing, or integration-callback defect that invalidates it. Independent SciPy calculations reproduce all four smooth growing modes and the stipulated number moments.

This is a plausible narrow theoretical contribution, not a demonstrated breakthrough or transport-ready diagnostic. The generic transform-locus construction is old moment theory. The potentially substantive contribution is its specific exact constrained-frequency application, the sharp fixed-k purely imaginary envelopes, and a constructive smooth example showing substantial instability hidden even by identical isotropic second moments. Priority of those particular neutrino results remains unestablished. Failure to find a prior match is not evidence of priority.

## Physical normalization and sector algebra

Use G=mu_e-mu_bar with mu defined after azimuthal integration, so integration is with respect to du and its mass is N. For k along the axial direction, I_j=integral u^j dG/(omega-ku). Then the longitudinal determinant is

`D_L=(1+I_0)(-1+I_2)-I_1^2`,

and each of the two degenerate transverse components is

`T=-1+(I_0-I_2)/2`.

The factor 1/2 is the azimuthal average of either transverse velocity square. The recurrence I_1=(omega I_0-g_0)/k, I_2=(omega I_1-g_1)/k gives the implemented affine target in both sectors. I independently expanded these expressions. There is no arbitrary branch selection in the upper half plane. Standard flavor-wave dispersion theory already separates these axial and axial-symmetry-breaking branches; the latter are not additional unrelated particle species. See [Fiorillo, Goimil-Garcia and Raffelt, Phys. Rev. D 111, 083028 (2025), Sec. II](https://doi.org/10.1103/PhysRevD.111.083028).

The normalizations need to be stated explicitly in a manuscript. N is a number-density-derived self-interaction potential once sqrt(2) G_F and the chosen reference scale are absorbed. Bare physical number density cannot be added to a frequency. Likewise omega and k are the reduced variables after the real matter and neutrino-current shifts. A shifted zero mode is not generally a physical homogeneous mode when g_1 is nonzero. The equal-zero-flux example avoids that spatial-shift ambiguity. [The zero-mode literature distinguishes these modes](https://doi.org/10.1103/PhysRevD.101.063001).

## Sharp envelopes and their limits

For equal masses N and zero species fluxes, the one-species transform lens has a horizontal difference width

`N/[y sqrt(1+y^2)] - N/(1+y^2)`, with `y=gamma/k`.

Comparing this with the longitudinal target k/(1+y^2), or transverse target 2k/(1+y^2), recovers exactly

`gamma_L^2 <= N^2 k/(k+2N)` and `gamma_T^2 <= N^2 k/[4(k+N)]`.

This checks the necessary side without appealing to an angular grid. The reflected two-ring family supplies the sufficient side through its direct beam dispersion polynomial. The claimed theorem requires gamma>0 and k>0; gamma=0 is a limiting boundary outside the upper-half-plane root test. It bounds only roots on the imaginary omega axis. Nothing in this audit establishes an arbitrary-real-part growth maximum.

These are fixed-k improvements. The longitudinal upper envelope tends to N as k tends to infinity, so a claim of a uniformly stronger global growth bound would require separate justification. The exact M0/M1 frequency map also must not be described as an exact M0/M1/M2 frequency map: the smooth construction is one special ambiguity example with matched second moments, while the four-disk classifier constrains only mass and flux.

## Smooth profiles and independent calculation

The vMF kernel normalized on S^2 is kappa exp(kappa n dot m)/(4 pi sinh kappa). Averaging it around the ring at u=1/3 and integrating the observation azimuth yields the implemented coefficient kappa/(2 sinh kappa), with the I_0 Bessel ring factor. `0F1(1;x)=I_0(2 sqrt(x))` removes the square root without changing the real profile. The atomic species weights 3/4 at 1/3 and 1/4 at -1 have mass 1, zero dipole and second axial moment 1/3. Axial azimuthal averaging gives the full second number tensor I/3. Rotational convolution preserves zero spherical-harmonic coefficients, so the dipole and quadrupole remain exactly zero for every kappa. Exponentials and the positive kernel make the profiles strictly positive and smooth on the full sphere, including the poles.

I independently used SciPy `special.i0`, `integrate.quad` and `optimize.brentq`, without importing smooth_check.py or its density/dispersion functions. The following roots lie strictly inside the saved rational brackets:

| kappa | independent gamma | absolute direct complex determinant residual |
| --- | --- | --- |
| 8 | 0.26693923381998724 | 5.6e-17 |
| 16 | 0.41343493440451845 | 6.9e-17 |
| 32 | 0.49356167106585314 | 8.4e-17 |
| 64 | 0.5352017806759694 | 2.3e-16 |

The residuals come from separately integrating all three complex I_j over [-1,1] and evaluating D_L, not from reusing the reduced scalar equation. All six moment triples agreed with (1,0,1/3) within approximately 3e-16. The saved endpoint signs were reproduced. Kappa 2 and 4 have no saved positive-gamma bracket; neither the saved run nor this check establishes their stability.

The ignored Arb `analytic` callback flag is acceptable here: the density is entire in u, while the dispersion integrand is meromorphic, with poles at plus/minus i gamma. Installed python-flint `acb.integral.__doc__` explicitly permits ignoring that flag for meromorphic functions because pole evaluation returns nonfinite balls. It requires special handling for branch functions such as sqrt and log. The 0F1 representation avoids that defect. This audit did not independently recalculate Arb's interval enclosures, but the local implementation and independent floating-point check are consistent.

`moment_enclosure.contains(target)` alone only checks numerical consistency. Exact preservation is supported by the analytic harmonic-convolution argument above, not by interval containment alone.

## Meaning of pressure and physical admissibility

The matched tensor is exactly the second **number angular moment**. Calling it actual thermodynamic radiation pressure requires an energy assumption. With a shared angle-independent energy spectrum, or a monoenergetic realization, the corresponding energy-weighted second moment also matches. That is a valid added realization, but it is not established for arbitrary angular-energy correlations. Froustey et al. explicitly flag this units distinction for their similarly named number-pressure moment. [Their Eq. (6) and footnote](https://arxiv.org/html/2311.11968v3#S2.SS1) also show why grey energy density and energy flux are not automatically the number moments required here.

The profiles satisfy positivity but impose no fermion occupation ceiling, specified energy spectrum, collisional equilibration requirement, angular regularity prior, or transport realizability constraint beyond the stated angular moments. Finite-kappa examples can be embedded in dilute phase-space distributions with an appropriate energy spectrum, but no astrophysical dataset has been shown to realize them. Smoothness removes an atomic-measure objection; it does not prove astrophysical frequency or prevalence.

Heavy-flavor cancellation is pointwise in angle in the model. Equal heavy-flavor masses and fluxes alone would not justify cancellation of their hidden angular structure.

## Prior art and what could still matter

The exact one-species transform range on a finite interval is classical. [Di Dio, arXiv:1506.06589, Sec. 5](https://arxiv.org/html/1506.06589v1#S5) describes Weyl-circle boundaries for truncated finite-interval moment problems. [Heide et al., Complex Analysis and Operator Theory 18, 82 (2024)](https://link.springer.com/article/10.1007/s11785-024-01525-1) explicitly identify truncated Hausdorff transform sets as intersections of two matrix balls. The present scalar two-disks result is a low-order special case. Intersecting two species' translated sets is a direct mathematical construction, not by itself a new general theorem.

Closure ambiguity is also established. [Nagakura and Johns, Phys. Rev. D 103, 123025](https://doi.org/10.1103/PhysRevD.103.123025) discuss the impossibility and uncertainty of recovering full angular distributions from a few moments. [Froustey et al.](https://arxiv.org/html/2311.11968v3) already perform moment-based flavor stability analysis using a chosen closure. [Quantum closures for neutrino moment transport, Phys. Rev. D 111, 063046 (2025)](https://doi.org/10.1103/PhysRevD.111.063046) studies physically constrained quantum closures. A statement that low moments alone do not determine flavor evolution is therefore inadequate novelty.

An especially direct conceptual predecessor is [Johns and Nagakura, Fast flavor instabilities and the search for neutrino angular crossings, Secs. IV-V](https://arxiv.org/html/2104.04106v2#S4). It explicitly distinguishes closure from full angular completion and compares different completions at the same zeroth and first moments, producing crossing disagreements. Its displayed examples do not establish the present sharp imaginary-frequency envelope or the particular same-isotropic-second-moment smooth growing pair, but the information-model framing and closure-versus-completion terminology are already present there.

The distinction worth pursuing is exact frequency attainability over **all positive completions**, together with sharp elementary limits and explicit witnesses. Ordinary crossing tests give some existence information, while closure-based spectra describe one selected completion. These are different quantifiers. The smooth matched-second-moment example strengthens the concrete counterexample beyond an abstract nonuniqueness observation, but should be compared directly with prior hidden-multipole examples before claiming priority.

Real utility is not yet established. The 504 points do not prove a continuous region, a fastest mode, a useful operating boundary on actual simulation data, or a better scientific decision than existing crossing/zero-mode diagnostics. The screen has no observed grid misses and no generic conic baseline. Exact arithmetic and independent replay are useful reliability properties, but do not replace evidence that the output changes a meaningful physical inference.

## Provenance and next discriminating work

I checked that every source hash in screen.json matches the current named source and that smooth.json's source/plan hashes match. Counts are 72 possible, 432 impossible, zero unresolved. This establishes local artifact consistency, not an external preregistration timestamp. Inputs are three exact synthetic/fiducial moment cases, not a newly analyzed astrophysical dataset. The disclosed prototype failure and pre-screen proof-format amendment should remain disclosed.

The strongest next discriminator is a literature-grade comparison of the sharp constrained envelopes and matched-second-moment example, then a narrowly sourced physical case with authentic per-species number moments and correctly converted wavevector. If those outputs remain scientifically informative, this could justify a focused methods/theory paper. If the useful answer is always the old generic growth bound, known crossing ambiguity, or a conventional closure comparison, the candidate should be reduced to a clear technical note rather than presented as a major discovery.
