# Independent review of the finite-wave-number result

Reviewed 2026-11-01. Scope: `finite-k-independent.md` and `finite-k-cap-checks.py`, against the original dispersion factors in `paper/draft.typ`. This review concerns positive axial measures, individual mass N, individual zero first moment, and equal supplied second moments Na, with N>0 and 0<a<1. It is a mathematical review of the stated restricted result, not a priority or physical-cutoff assessment.

## Verdict

The all-real-frequency cubic upper bound is valid. The stated two-atom measures attain it exactly in the specified low-second-moment band. Monotonicity of that upper bound then proves the sharp finite-interval maximum **within that band**, and the finite-K cube-root asymptotic is sharp for fixed K>0 as a tends to zero. I found no error invalidating these results.

This does not solve the complete finite-K problem for arbitrary a and K. Outside the attained cubic band and the separately established unrestricted-cap threshold, the note supplies upper bounds and examples rather than a sharp general formula. Its own limitations accurately say this. The two branches must also remain distinguished: Section 6 is a statement about each branch's extremum. When the longitudinal cubic is attained at K, it bounds the union of the two branches as well, because the transverse cubic root is strictly smaller there.

The interval is 0<k<=K, with K>0. At k=0 equal signed moments through order two make I0=I1=I2=0 for an upper-half-plane frequency, so neither original factor has a growing root. If negative k is included, the appropriate bounded domain is |k|<=K; reflection of both measures preserves these constraints. A literal domain k<=K with arbitrarily negative k is not a finite-wave-number cutoff.

## Independent derivation checks

### Original factors and signs

Write z=omega/k, S_s=integral dmu_s/(z-u), and U_s=(z²-1)S_s-Nz. The signed zeroth and first moments vanish, so I1=z I0 and I2=z² I0. Direct expansion of the original, unmultiplied factors gives

    D_L=-1+(z²-1)I0,
    T=-1-(z²-1)I0/2,
    I0=(S_+-S_-)/k.

Thus U_+-U_-=k for D_L=0, and U_+-U_-=-2k for T=0. This agrees with the paper's azimuthal factor of one half and the signs in the new note. No real-pole cancellation or multiplication introduces extra upper-half-plane roots.

### The single Gram disk and its limitation

In L²(mu), the functions 1 and u are orthogonal with squared norms N and Na. Projecting f(u)=1/(z-u) onto their span gives

    integral |f|² dmu >= |S|²/N + |zS-N|²/(Na).

Since integral |f|² dmu=-Im(S)/Im(z), this is the stated necessary inequality. Completing its square independently gives exactly the displayed center and radius. Affine multiplication by z²-1 gives

    r_U=Na |z²-1|/[2y(a+x²+y²)].

This uses only positivity and the three supplied moments. It is an enclosing disk, and there is no assumption that every point in it is realized by a measure supported on [-1,1]. That distinction matters outside the sharp band.

For t=y² and d=x², I independently expanded the squared-radius comparison and recovered

    d {2(a+t)[1+a+(3-a)t] + (1-a)(1+a+2t)d}.

All factors that were cleared are strictly positive. For 0<a<1 and y>0 this expression is positive whenever x is nonzero. Consequently the zero-real-frequency radius bounds every real part, and equality of the final diameter bound requires x=0. The step is an analytic all-frequency argument, not an extrapolation from rational tests.

### Cubic cap and monotonicity

The two species share the same enclosing disk, so their difference has modulus at most twice its radius. With c=1 for longitudinal and c=2 for transverse, this gives

    c k <= Na(1+y²)/[y(a+y²)],  y=gamma/k.

Equivalently f(gamma)<=0, where

    f(gamma)=gamma-(Na/c)(k²+gamma²)/(a k²+gamma²).

For positive k and 0<a<1, its gamma derivative is strictly positive, f(0)=-N/c, and f tends to positive infinity. Hence it has exactly one positive zero, the stated cubic root Gamma_c(k), and every growing root satisfies gamma<=Gamma_c(k). The derivative in k stated in the note is also correct and strictly negative at fixed positive gamma. Implicit differentiation therefore proves dGamma_c/dk>0.

The root lies between Na/c and N/c for k>0. Its positive limit Na/c as k tends to zero is only the limit of a loose upper bound; it does not assert instability at k=0.

## Attainment and exact sharp domain

Let d=y+a/y, pq=a, and q-p=d. The quadratic definitions imply p>0, q>0, and q>=p. The condition q<=1 is equivalent to d<=1-a: q solves q²-dq-a=0, and evaluating this increasing positive-root threshold at q=1 gives the equivalence. It also implies p<=1.

For the displayed weights, each species has mass N, first moment zero, and second moment Npq=Na. All weights are positive. Direct summation, without a disk sufficiency claim, gives

    I0=-2N k³ a d/[(omega²-k²p²)(omega²-k²q²)].

At omega=iky,

    (y²+p²)(y²+q²)=(y²+a)²+y²d²=2(y²+a)².

Thus I0=-Na/[k y(a+y²)], and the original longitudinal factor vanishes precisely at the stated k_L. Swapping species and halving k and omega reverses and doubles I0, so the original transverse factor vanishes as claimed. This independently verifies attainment with support on [-1,1].

The support inequality is y²-(1-a)y+a<=0. It has positive solutions exactly when a<=3-2sqrt(2), and the roots y_- and y_+ in the note are correct. At equality in this restriction the two roots coincide and the sharp band is a single wave number, not an open interval.

For explicit translation into supplied K, define

    k_c(y)=(Na/c)(1+y²)/[y(a+y²)].

Its derivative has the sign of

    -[y⁴+(3-a)y²+a],

so it decreases strictly. The cubic-attainment band is therefore exactly

    k_c(y_+) <= K <= k_c(y_-),

provided 0<a<=3-2sqrt(2). Equivalently, substitute y=Gamma_c(K)/K into y+a/y<=1-a. This makes the note's implicit condition checkable from the inputs and prevents treating “low second moment” alone as sufficient for arbitrary K.

Within this band, every smaller positive k obeys gamma<=Gamma_c(k)<=Gamma_c(K), and the explicit measure at K realizes equality. Sharpness at all smaller k is unnecessary. Outside this band the proof still gives the cubic upper bound, but it gives no attaining pair for that bound.

For equal second moments the separately proved unrestricted cap N(1-a)/c is attained at k=N(1-a)/(c a). Consequently the exact finite-interval maximum equals that cap whenever K is at least this threshold. The cubic and unrestricted caps can be combined as upper bounds elsewhere; the minimum of two bounds is not automatically an attained envelope.

## Degenerate finite-K limit

For fixed N>0 and K>0, put A=Na/c and Gamma=(A K²)^(1/3)h. The cubic becomes

    h³-(A/K)^(2/3)h²+a(K/A)^(2/3)h-1=0.

Both middle coefficients tend to zero as a tends to zero. Evaluating this expression at h=1-epsilon and h=1+epsilon and using uniqueness of the positive root proves h tends to one. This supplies the needed root-control argument behind the asymptotic, rather than relying only on formal coefficient limits.

It follows that y=Gamma/K is of order a^(1/3), while a/y is of order a^(2/3). Their sum tends to zero and is eventually below 1-a. Thus K is eventually inside the sharp cubic-attainment band for each branch. The sharp finite-K maximum really is asymptotic to (Q K²/c)^(1/3), and the transverse-to-longitudinal ratio tends to 2^(-1/3). At a=0 positivity forces both measures to be N delta_0, so there is no growing root. Defining the worst-case growth as zero in that stable case gives the stated continuity. The incompatible unrestricted-k limit is explained by its attaining wave number diverging like 1/a.

## Executed checks and evidence limits

I ran all three test methods under `python3 -B -O`, importing the file through `runpy` with a non-main name so the test runner did not rewrite its saved results JSON. All three passed. The checks cover 1,000 rational radius identities, 400 Gram-disk square completions, and 72 direct original-factor attainers, including moment and support checks. The construction count and 60 strict endpoint-family improvements match the source claims.

These finite computations corroborate algebra and examples. They do not establish the general inequalities, arbitrary-real-frequency optimization, support equivalence, monotonicity, or limit theorem; those were checked above from the written formulas. I did not run a continuous optimization, establish a complete remaining transition envelope, or assess novelty or smooth exact attainment.

Reviewed source SHA-256 values:

    finite-k-independent.md  5477181a3b56de3f465dddc36967e2206957d9016c4722e6847c13d08935fc4d
    finite-k-cap-checks.py   d030422a6beb95303f7c01f9769fcccb4d9cbb384ae9f2bff3c4b1e2e9b7738b
    independent_verifier.py 1f79a4dfbcfe9ccc180ff1bdd01579c25b0902fa6960ae19b7448ec10a347bce
    paper/draft.typ         788fab4f4f64ebc4f6dca523be64aa9f5b2a69e1d35f6f0a7b7acaa264810c8c

## Paper integration review

I subsequently checked the integrated section “Bounded wave numbers: a sharp subregime and a failed extrapolation” in `paper/draft.typ`, SHA-256 `562a4e79a60960def2f4cfb3f933dba1eba1dccba4915eac4e67d67f07be9b5a`. The cubic coefficients, longitudinal/transverse factor, counterexample values, support inequality, y endpoints, branch-specific K band, unrestricted saturation threshold, and ratio-to-one limit notation agree with the accepted proof. Earlier definitions correctly identify z=omega/k=x+iy. The section retains the finite mathematical interval, restricted sharpness, unresolved general transition, and physical/smooth-attainment limitations. I found no mathematical integration gap. No tests were repeated and no paper source was changed for this follow-up review.
