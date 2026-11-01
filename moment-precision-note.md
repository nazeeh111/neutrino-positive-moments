# Illustrative precision of the designed moment boundaries

2026-11-01. Independent scope: an exact sensitivity calculation for the already designed axial boundary examples. These perturbations are mathematical choices. No measured error bars, astrophysical uncertainty, or physical precision requirements are supplied or inferred. No production classifier, preregistered result, paper, or comparison result is changed.

## Result

The existing exact-boundary advantage is conditional on the prescribed moment data and frequency being exact. At the two designed M1 boundaries, **every positive increase in the common species mass gives an explicitly strictly feasible rational moment certificate at the original frequency**. A floating conic solver may then certify the case if its working precision resolves the resulting margin. This calculation does not establish when it will do so, and does not undo the recorded exact-data unresolved outcomes.

Symmetric uncertainty intervals require care: the nominal compatible completion remains available. Widening an existential uncertainty class cannot turn a possible frequency into an impossible one. The sensitivity described here concerns loss of exact-boundary geometry and the addition of frequencies previously excluded by the nominal data, not removal of the nominal compatible root.

## 1. Common mass perturbation at the original fixed k

Keep individual fluxes exactly zero and require the two species to have the same mass N. The sharp M1 caps already proved over all real frequency parts are

    B_L(N,k) = N² k/(k+2N),
    B_T(N,k) = N² k/[4(k+N)],

where B is squared growth. Both increase strictly with N>0. The designed examples are

    longitudinal: N=1, k=2/3, omega=i/2;
    transverse:   N=1, k=1/3, omega=i/4.

For N=1+e, e>0, the exact squared margins above the original rates are

    B_L(1+e,2/3) - 1/4  = e(5+4e)/[4(4+3e)],
    B_T(1+e,1/3) - 1/16 = e(5+4e)/[16(4+3e)].

Thus an arbitrarily small positive mass shift moves either original frequency strictly below its M1 cap. An arbitrarily small negative fixed mass shift makes it impossible. If instead N is allowed in [1-e,1+e], the original point remains possible and admits a strict completion at the upper endpoint.

A direct rational strict certificate is available, so the conclusion does not rely on a continuity assertion or a numerical search. Write z=omega/k=i3/4, S_s=integral dmu_s/(z-u), and U_s=(z²-1)S_s-Nz. For the longitudinal example choose

    U_+ = 1/3 + i2N/3,       U_- = -1/3 + i2N/3,
    S_+ = -16/75 - i68N/75, S_- = 16/75 - i68N/75.

Then U_+-U_-=2/3=k. For the transverse example interchange species; the difference becomes -2/3=-2k. These are the original root conditions with equal masses and zero individual flux.

For the plus species, reweight dnu=dmu/|z-u|². Its first four moments are

    (t0,t1,t2,t3) = (272N/225, 16/75, 8N/25, -3/25).

For the minus species the odd moments reverse sign. The degree-three interval moment matrices are

    H_sigma = [[t0+sigma*t1, t1+sigma*t2],
               [t1+sigma*t2, t2+sigma*t3]], sigma in {-1,+1}.

Their two determinants, swapped by reflection, are exactly

    16(N-1)(4N+1)/225,
    16(N+1)(4N-1)/225.

All diagonal entries and both determinants are positive for N>1. Hence the interval moment theorem gives positive supported completions with strict positive definite matrices. At N=1 one determinant is zero. In standard convex-optimization language this is a strict feasible point, often called a Slater point. Rational N gives rational certificate entries even for extremely small e. The determinant margin vanishes linearly with e; exact existence alone does not promise robust floating-point recovery at every e.

For a point above the nominal rate, set gamma=gamma0*r, r>1. Both designed examples have the same exact minimum common mass

    N_min(r) = r[3r + sqrt(9r²+16)]/8.

The frequency is impossible below N_min, possible at N_min, and has a strict moment completion above it. Thus the precise extra mass needed is N_min(r)-1. For r=1+eta with small positive eta, this is (8/5)eta+O(eta²). For example, gamma=gamma0*(1+e/2) is already strictly attainable at N=1+e for every sufficiently small positive e, although it is impossible at N=1. This is an illustrative relative tolerance relation, not a precision requirement for any physical dataset.

## 2. Why uncertainty only in Q does not erase the M1 stress example

The designed M1 attainers have Q_+=Q_-=3/8. Intervals [3/8-delta,3/8+delta] contain these same measures for every delta>=0 after clipping to [0,1]. The sharp M1 fixed-k cap is unchanged because N and flux are unchanged. Therefore uncertainty only in those second moments cannot remove the original boundary completion or move its growth beyond the M1 cap. No full fixed-k Q-interval envelope is needed to make this statement.

Separately, the established all-positive-k M2 theorem gives a simple exact sensitivity result. Allow independent Q_+,Q_- in the same clipped interval [l,h], where l=max(0,a-delta), h=min(1,a+delta). For 0<l<1 the sharp worst-case growth maxima are

    longitudinal: 1-l;    transverse: (1-l)/2.

They are attained at Q_+=Q_-=l, at k=(1-l)/l for longitudinal and half this k for transverse. These optimizing k values move with l, so these formulas must not be presented as prescribed-k compatibility tests.

If l=0<h, the corresponding suprema are 1 and 1/2. They are not attained: approaching Q_+=Q_-=0 through positive values sends the maximizing k to infinity, while Q_+=Q_-=0 itself forces identical delta_0 species and is stable. If h=0, the uncertainty set is that stable singleton and both maxima are zero. If l=1, both species are endpoint-only and stable. These exceptions prevent a spurious claim of a discontinuity at fixed finite k. The accepted finite-K theorem instead has a sharp cube-root approach to zero at Q=0 for fixed K.

At a fixed interior a and delta<a, these all-k maxima increase exactly by delta (longitudinal) and delta/2 (transverse). A rate gamma above the nominal all-k longitudinal cap 1-a cannot be admitted by this interval unless delta>=gamma-(1-a); transverse requires delta>=2gamma-(1-a). Equality provides an attaining completion when the resulting l is positive. This concerns maximum growth somewhere in k, with the displayed attainers, rather than compatibility of an arbitrary preassigned complex frequency.

## 3. Consequence for the comparison claim

The existing conic comparison proves an exact-data certificate advantage on specified designed boundary stresses; it does not establish an advantage after measurement uncertainty. The explicit mass perturbation removes the forced singularity and supplies a concrete rational strict certificate. A proportionate follow-up could test a few declared e values with the frozen conic implementation, recording exact certificate success and residual margins. It must not call the chosen e values realistic, replace the original comparison, infer a general speed advantage, or claim exact arithmetic has lost its exact-data correctness benefit.

This note establishes mathematical sensitivity and an explicit positive margin. It does not report a new conic runtime experiment or new physical utility. The failed physical M2 utility gate remains failed.

## Checks

`python3 -B -O test_moment_precision.py` checks the displayed determinants, both original dispersion factors, exact squared margins at the four declared rational e values 1/100, 1/10000, 1/10^8, and 1/10^12, and a strictly above-nominal rate at N=1+e. For both sectors at all four e values, the unchanged independent M1 checker verifies explicit two-atom quadrature certificates at N=1+e and four-disk emptiness certificates at N=1-e. Its independently reconstructed H matrices are strictly positive definite for the positive shifts. It checks the interval endpoint/corner formulas and the already accepted all-k attainer at a single interior interval. These checks corroborate the algebra; the interval and strict-feasibility conclusions use the written theorems above.
