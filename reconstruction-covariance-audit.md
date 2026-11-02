# Reconstruction covariance: exact local diagnostic, no new algorithm

Status: exact arithmetic checked by the primary authoring agent; **no independent review completed**. The requested reviewer returned a usage-limit error. This is a scoped follow-up, outside the accepted moment-envelope theorem and unchanged paper. It is not an execution or reproduction of FLASH.

## Question and source mapping

Does reconstruction commute with a constant change of flavor basis? Write the requirement as `R(U rho_j U†) = U R(rho_j) U†`. This compares representations of the same inputs, transforming all quantities consistently. It does not assert invariance under spatially varying basis changes or under rotating densities while keeping physical couplings fixed.

[Grohs et al., arXiv:2501.07540v2, Section II.4](https://arxiv.org/html/2501.07540v2#S2.SS4) uses classical minmod for diagonal entries (Eqs. 25–26). For strictly positive, increasing, real off-diagonal entries, the phase and ordering tests select line interpolation (Eqs. 27–31), also giving ordinary minmod. Every input below satisfies those restrictions in both bases; no circular interpolation, phase wrapping or signed-modulus ambiguity is involved. This isolates a documented reconstruction subroutine, not the full transport, closure or source update.

## Exact counterexample

Let `rho(x,z) = [[1+z,x],[x,1-z]]/2`, with Bloch coordinates `(x,z)`. Choose the three cells

```
left   = (1/5, 1/10)
center = (3/10, 3/10)
right  = (1/2, 2/5).
```

All three matrices have trace one and are strictly positive. Their right reconstructed face has coordinates `(7/20,7/20)`. Set

```
U = [[2,-1],[1,2]] / sqrt(5),
(x,z) -> ((3*x+4*z)/5, (-4*x+3*z)/5).
```

Reconstructing the transformed inputs gives `(13/25,-3/50)`. Transforming the previously reconstructed face gives `(49/100,-7/100)`. Consequently

```
R(U rho_j U†) - U R(rho_j) U† = [[1,3],[3,-1]] / 200.
```

Its squared Frobenius norm, the sum of squared matrix entries, is exactly `1/2000`. The reconstructed purity `Tr(rho²)` changes by `29/2000`. Both answers remain positive, so this example demonstrates basis dependence, not negative density eigenvalues. Both left and right faces preserve their center by arithmetic averaging.

The standard-library checker evaluates matrix conjugation and entrywise reconstruction directly with rational numbers, then compares them with these independently written expected coordinates. It checks input and face positivity with exact principal minors. These are two calculation routes within one author's check, not independent scientific review.

## Smooth-grid control

Take samples of the fixed smooth profile

```
x(t) = 3/10 + 3*t/20 + t²/20,
z(t) = 3/10 + 3*t/20 - t²/20
```

at `t=-h,0,h`. For `0<h<3/7`, the face difference is

```
(3*h²/200) * [[1,1],[1,-1]],
||difference||_F² = 9*h⁴/10000.
```

The retained checker verifies six exact spacings `h=1/4,...,1/128`. The norm therefore shrinks quadratically for this sampled profile. This is not a full finite-volume convergence study, but it directly defeats an inference of nonconvergence from the finite-stencil example. No observed neutrino growth-rate discrepancy has been attributed to this mechanism.

## Prior methods and novelty decision

[Hajduk, Kuzmin and Aizinger, *Frame-invariant directional vector limiters for discontinuous Galerkin methods*, June 2018 preprint](https://www.uni-due.de/imperia/md/images/mathematik/ruhr_pad/ruhrpad-2018-15.pdf), Sections 1, 4, 5 and 7.1, already analyzes nonobjective componentwise limiting and constructs frame-dependent directions that transform consistently. It includes SVD/Gram–Schmidt approaches and an explicit rotated-data comparison. Its SVD discussion warns about degeneracy; its introduction also warns that a common correction factor can be overly diffusive. These are strong precedents for the mechanism, although spatial-frame transformations in that paper are not identical to internal flavor rotations here.

[Guermond, Popov and Tomas (2019)](https://arxiv.org/abs/1807.02563) supplies a general convex-limiting framework for hyperbolic systems with sources. [Olbrant, Hauck and Frank (2012)](https://www.sciencedirect.com/science/article/pii/S0021999112001362) develops realizability-preserving radiation moment discretization with a scaling limiter and an appropriate time-step restriction. Neither citation is claimed to implement this quantum-neutrino system.

Our inference: changing from separate component limiting to a covariant vector limiter and positive-matrix scaling does **not by itself identify a novel numerical mechanism**. Existing vector limiting and convex-domain methods are necessary serious baselines. No exact prior copy of this rational stencil was sought or is needed to reject a claim that the general defect or repair idea is new. Conversely, these precedents do not prove that all useful quantum-specific extensions are solved.

## A precise restriction on a possible repair

An elementary geometric obstruction prevents simply demanding scalar monotonicity in every basis. Suppose a fixed-trace Bloch-vector linear reconstruction uses faces `rC +/- s/2`. Require, for every linear observable `a dot r`, that its right face lies between `a dot rC` and `a dot rR`, and its left face lies between `a dot rL` and `a dot rC`.

By separation of a point from a closed line segment, the first requirement implies `rC+s/2` belongs to `[rC,rR]`; the second implies `rC-s/2` belongs to `[rL,rC]`. Thus `s` must be a nonnegative multiple of both `rR-rC` and `rC-rL`. If those increments are not collinear, `s=0`. This is a two-neighbor, every-observable monotonicity restriction, not a no-go theorem for covariance, higher-order transport or less restrictive convex-hull constraints. It motivates choosing physical invariant constraints rather than importing all scalar requirements simultaneously. No priority is claimed for this elementary argument.

## Decision and reproducibility

**Do not promote this lead to the selected moonshot algorithm.** The finite-stencil diagnostic is useful as a solver regression requirement, but it establishes neither high impact nor an algorithmic invention. A larger implementation would first need a quantum-specific mechanism beyond established vector/convex limiting, full-update conservation and positivity arguments, and a predeclared physically meaningful comparison. No installation or large scan is justified by this example alone.

Run:

```sh
python3 -B -O usefulness-assessment/reconstruction_covariance.py
```

The frozen report is `usefulness-assessment/reconstruction-covariance.json`; the script reports its own SHA256. It does not overwrite results unless explicitly passed `--output PATH`. The earlier physical-use and speed gates remain failed, and the finite-band theorem's priority and practical usefulness remain unresolved.
