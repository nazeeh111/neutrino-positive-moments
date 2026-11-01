# Independent all-real-frequency audit of the complete equal-M2 envelope

2026-11-01. This is a new mathematical result/review, not a revision of the previous failed utility gate or evidence of an astrophysical cutoff. Scope: positive axial species measures on [-1,1], each with mass N>0, individual first moment zero, and the same individual second number moment Na, 0<a<1. Write r=1-a, c=1 for longitudinal and c=2 for transverse, k>0, z=omega/k=x+iy with y>0. The older `finite-k-independent.md` proved a necessary cubic cap and its sharp low-Q band; this audit resolves its remaining equal-second-moment transition, subject to the already proved two-disk moment locus and original-factor reduction. Priority and scientific impact require separate assessment.

## Verdict

The imaginary-axis piecewise envelope is also the sharp envelope over every real frequency. The missing argument is an enclosing-disk construction whose radius cannot increase off the imaginary axis. This conclusion is analytic; the bounded scan described below was only a falsifier search.

All three regimes have explicit positive atomic attainers: the original endpoint-compensated pair, the previously derived two-interior-atom pair, and a three-atom plateau pair. The maximum at fixed k is nondecreasing in k up to saturation. Therefore the sharp maximum for 0<k<=K is its value at k=K. Reflection gives the same answer for |k|<=K.

## 1. The two necessary and sufficient species disks

Let S=integral dmu(u)/(z-u) and U=(z²-1)S-Nz. The positive moment locus already established in `m2_frequency.py` is the intersection of two disks. In U coordinates their centers and radii are

    C0 = i N r/(2y),                  R0=N r/(2y),
    C1 = N[-x+i((2-a)y²+a(1-x²))/(2y)]/(a+x²+y²),
    R1² = N²a² |z²-1|²/[4y²(a+x²+y²)²].

The centers/radii agree with the parent derivation. Put L²=|C1-C0|². Equal signed zeroth and first moments reduce the ORIGINAL dispersion factors to U_+-U_-=k longitudinal and U_+-U_-=-2k transverse. Consequently a root needs two points of this same locus separated horizontally by c k. No generic Euclidean diameter is substituted for that horizontal separation.

## 2. A disk containing the intersection, for every beta

For any beta in[0,1], multiply the two disk inequalities by1-beta and beta and complete the square. Every point of their intersection lies in the disk

    C_beta=(1-beta)C0+beta C1,
    R_beta²=(1-beta)R0²+beta R1²-beta(1-beta)L².

Thus the horizontal separation of two attainable points is at most2R_beta, for every beta. The enclosing disk may not itself be an attainable locus; only containment is used. Rewriting its squared radius gives

    R_beta²=(1-beta)R0²+beta² R1²
             +beta(1-beta)(R1²-L²).

All three coefficients are nonnegative. This rewrite supplies the off-axis proof.

## 3. Both variable terms are bounded by their x=0 values

The accepted Gram-radius argument gives R1²(x,y)<=R1²(0,y), strictly if x!=0. Independently, set t=y² and w=x²+y². Directly subtracting the squared center separation from R1² gives the simpler identity

    R1²-L² = N²(1-a²)(a-w)/[4t(a+w)].

To verify it without trusting a symbolic simplifier, the numerator before cancellation is

    a²[(w-1)²+4t]-[(a²+t-x²)²+4tx²]
      = (1-a²)(a²-w²).

The denominator is4t(a+w)². The function (a-w)/(a+w) has derivative -2a/(a+w)²<0, so this term cannot increase when x² increases. Since R0 does not depend on x, it follows for EVERY fixed beta in[0,1] that

    R_beta²(x,y)<=R_beta²(0,y).

For beta>0 this is strict at x!=0. Choosing beta optimally at x=0 now bounds every real frequency, without any off-axis search or assumption that a smaller single-disk radius alone controls the lens width. This is the missing step in the earlier partial theorem.

## 4. The optimal enclosing radius at x=0 is exactly the half-width

At x=0 the two disk centers are on the imaginary axis. Write

    R0=N r/(2y),
    C1=i V1,    V1=N[a+(2-a)y²]/[2y(a+y²)],
    R1=Na(1+y²)/[2y(a+y²)],
    delta=V1-R0=N(a²+y²)/[2y(a+y²)]>0.

The unconstrained minimizing beta is

    beta0=(R0²+delta²-R1²)/(2delta²)
         =r(y²-a²)(a+y²)/(a²+y²)².

Clamp beta0 to[0,1]. The result is not merely a loose enclosing circle:

- If beta0<=0, the horizontal extrema of disk0 lie in disk1, so the half-width is R0. This is exactly y<=a.
- If beta0>=1, the horizontal extrema of disk1 lie in disk0, so the half-width is R1. This is exactly y+a/y<=r.
- If0<beta0<1, the circles intersect at two horizontally aligned points; both lie in the lens. Their height and half-separation are

      V*=N r y/(a²+y²),
      H=N a r/(a²+y²).

  The beta0 enclosing circle has center iV* and radius H, so the bound is attained.

The beta1 condition follows by multiplying positive denominators:

    beta0>=1 iff y⁴-(1-4a+a²)y²+a²<=0
                  iff (y+a/y)²<=r².

All quantities are positive, giving the stated equivalence. The two special regimes are disjoint: y<=a makes a/y>=1, so it cannot satisfy y+a/y<=r<1.

Therefore the exact maximum horizontal chord at x=0, and an upper bound at every x, is

    W(y)=N r/y,                                y<=a;
    W(y)=Na(1+y²)/[y(a+y²)],                    y+a/y<=r;
    W(y)=2Na r/(a²+y²),                        otherwise.

The formulas agree at their transition boundaries. Each active piece decreases strictly with y; the entire continuous W does also. Hence the root condition c k<=W(y) gives a unique maximal y at a given k.

## 5. Complete fixed-k and finite-K envelope

Define

    k*=N r/(c a),
    E_c(k)=sqrt(2Na r k/c-a²k²),

and let Gamma_c(k) be the unique positive root of

    Gamma³-(Na/c)Gamma²+a k² Gamma-(Na/c)k²=0.

If0<a<=3-2sqrt(2), let

    y_±=[r±sqrt(1-6a+a²)]/2,
    k_c(y)=Na(1+y²)/[c y(a+y²)],
    k_low=k_c(y_+),    k_high=k_c(y_-).

Then the sharp growth F_c(k), over all positive compatible species measures AND all Re(omega), is

    F_c(k)=N r/c,             k>=k*;
    F_c(k)=Gamma_c(k),        k_low<=k<=k_high;
    F_c(k)=E_c(k),            0<k<k* outside that cubic band.

For a>3-2sqrt(2), the cubic band is absent and E_c(k) applies throughout0<k<k*. At equality a=3-2sqrt(2), the band is a single k, where the two formulas agree. In the low-a case the ordering is0<k_low<=k_high<k*: both y endpoints exceed a and k_c decreases strictly.

This is a piecewise attained formula, not the minimum of three unrelated rate formulas. In the interior cubic band E_c(k) is only the rate of a nonoptimal endpoint pair, and is smaller than the sharp cubic rate. Taking min(E_c,Gamma_c,Nr/c) would be wrong.

E_c increases strictly until k*, the cubic increases strictly by the accepted implicit-derivative proof, and values match at the transitions. Therefore F_c is nondecreasing. The sharp maximum for0<k<=K is F_c(K). At k=0 the signed moments through second order vanish, so there is no growing root. At a=0 or a=1 positivity forces identical species measures and the envelope is zero. These degenerate cases must be handled separately from formulas dividing by a. The branch scaling identity F_2(k)=F_1(2k)/2 follows from the displayed formulas and original factors.

## 6. Explicit attainers in every regime

Endpoint regime: use

    mu+ =N/(1+a) delta_a +Na/(1+a) delta_-1,
    mu- =reflection(mu+).

Then U_+-U_-=2Na r/(a²+y²) at z=iy. Choosing y=E_1(k)/k gives the longitudinal root. Swap the species for transverse and use its c=2 formula. Mass, zero flux and individual second moment Na are immediate. This construction is valid more broadly; it is optimal only in its stated regime.

Cubic regime: retain the accepted construction with pq=a, q-p=y+a/y, q<=1, mu+=Nq/(p+q)delta_p+Np/(p+q)delta_-q and mu- its reflection. Here y=Gamma_c(k)/k. Use the original ordering for longitudinal and swap for transverse. It attains the Gram disk diameter exactly.

Plateau regime: an explicit fixed-k three-atom construction also exists, although finite-K sharpness could already use the older two-atom attainer at k*. Set y=Nr/(c k)<=a and

    w=Nr/(1-y²),
    mu+ =w delta_y
          +N(a-y)/[2(1-y)] delta_1
          +N(a+y)/[2(1+y)] delta_-1,
    mu- =reflection(mu+).

All weights are nonnegative. Directly, their mass is N, first moment zero and second moment Na. More simply, (1-u²)mu_+=Nr delta_y and the reflected weighted measure is Nr delta_-y. Therefore

    U_+=Nr(1+i)/(2y),
    U_-=Nr(-1+i)/(2y),
    U_+-U_-=Nr/y=c k.

This gives gamma=k y=Nr/c. Swap species for transverse. At y=a the positive endpoint1 weight vanishes and the construction reduces to the old two-atom pair. For y<a it is generally genuinely three-atomic. No smooth exact attainment follows.

## 7. Independent falsifier and limits

Before the enclosing-disk proof, an executable bounded scan used the existing NumPy/SciPy environment, with a in{.01,.03,.06,.1,.15,.1715,.2,.3,.5,.7,.9},40 logarithmic y points from1.001a to5, and61 x points over[0,5]. For each x it maximized the horizontal intersection-interval length over the common vertical range using scalar bounded minimization. There were no off-axis widths exceeding W(0,y) by1e-6 relative. This is only an exploratory counterexample search, not the proof. Temporary script `/private/tmp/finite_k_width_attack.py` records its implementation; no existing canonical data was changed.

The general proof is Sections2-4, and the attainers give sufficiency without assuming enclosing disks are fully attainable. The old two-disk moment locus is classical, and convex combinations of disk inequalities are elementary. Any original contribution would be the specialized complete constrained envelope and its derivation, not invention of these generic tools. This audit does not establish priority, practical decision usefulness, smooth-measure exact attainment, nonlinear conversion, a justified physical ultraviolet cutoff, actual three-dimensional transport or computational speed superiority. It does not change the failed M2 utility gate or turn the negative fusion prototype into a success. Parent executable exact checks and a further independent review are appropriate before editing the paper.
