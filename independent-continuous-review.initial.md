# Continuous comparator: plan recorded before benchmark execution

This comparison asks whether the specialized exact geometry provides useful capability beyond a conventional continuous conic formulation. It is not a comparison against an angular grid, and no superiority is assumed.

## Frozen experiment plan

The existing 504 screen points and 90 boundary/scaling stress points are development evidence, not held-out data. Reuse all 594 without changing them. Add 128 deterministic cases from Python Random(271828), drawing masses, flux fractions, k and complex frequency from explicit finite rational lists in baseline_continuous_exploration.py. These fresh mathematical cases are synthetic, not astrophysical data. Include both dispersion sectors. The generator is fixed before examining outputs.

Use three interleaved repetitions, with order independently shuffled by Random(1729); each method is run on every case. Python imports and independent checking are outside generation timings. Include exact problem preparation, solver construction, solve and proof extraction in each method's total. Report coverage, disagreements, unresolved results, per-case median timings and medians across cases. A floating status is not an exact certificate. Do not compare a solver's internal solve time with another method's end-to-end time.

Primary comparison: exact four-disk production classifier versus a generic Clarabel SOCP using the original degree-three Hausdorff matrices. Both return the existing independently replayable proof formats when extraction succeeds; unresolved boundary cases remain unresolved. No fallback to the exact four-disk classifier is allowed inside the conic method. A secondary ordinary floating circle-intersection method checks whether the generic Weyl geometry already gives equivalent approximate results cheaply. It must be identified as classical geometry, not a new algorithm.

Settings fixed before execution: Clarabel max_iter=200, tol_gap_abs=tol_gap_rel=tol_feas=tol_infeas_abs=tol_infeas_rel=1e-10; other settings retain defaults, max_threads=1 when exposed. A normalized conic margin greater than 1e-8 suggests possible, less than -1e-8 suggests impossible, and the intervening region is numerically unresolved. Exact certificate acceptance, rather than those suggestions or solver status, defines certified output. Floating geometry uses a normalized containment tolerance of 1e-9; no exact guarantee is attached.

No practical superiority claim is allowed unless coverage and complete generation/checking costs support it. Any finite-set speed result is limited to this small axial model, software versions and machine. Exact behavior at boundaries and rescaling is a distinct guarantee, not evidence that a modern conic solver is mathematically wrong. Development timings in the prior screen remain exploratory.

## Generic SOCP construction

For each species and each sign, the established Hausdorff condition is H_j(v)=C_j+v_0 A_j+v_1 B_j positive semidefinite. Eliminate the second species' transform using its known difference from the first, and use s_1=(N_1/Im z)(v_0+i v_1). The valid bounds are -1<=v_0<=1 and -1<=v_1<=0. Scale each complete affine matrix by the same positive exact maximum of its coefficient magnitudes and one, as in the existing numeric comparison.

Maximize delta subject to H_j(v)-delta I positive semidefinite for j=1,...,4. A real symmetric 2x2 matrix [[a,b],[b,c]] is positive semidefinite iff (a+c,a-c,2b) belongs to the three-dimensional Lorentz cone. Thus the problem has three scalar variables, four SOC(3) constraints and four nonnegative bounds. Its optimal delta has the feasibility sign; delta=0 includes tangency. This formulation uses no angular discretization or Weyl-circle reduction.

Clarabel's official Python interface requires Ax+s=b with conic slack s, sparse CSC matrices, and separate primal/dual/accuracy statuses. The implementation follows [the official getting-started documentation](https://clarabel.org/stable/python/getting_started_py/) and [solver settings](https://clarabel.org/stable/api_settings/), retrieved November 1, 2026. Positive proofs are rationalized transform points followed by exact quadrature checks. Negative proofs use rational rank-one Hausdorff dual directions and exact nonnegative combinations, as in the existing comparison backend. Every accepted proof is replayed by independent_verifier.py; sharing certificate extraction does not make solver status a proof.

## Environment inspection

Before any setup, the two existing scientific environments had SciPy but no Clarabel, CVXPY, SCS or ECOS. The bundled artifact Python had none of these numeric packages. No existing environment was modified by this inspection. The parent is preparing a separate free task-local comparator environment.
