# Fixed wave number with supplied equal second moments

This note closes part of the finite-k question and rejects a tempting wrong generalization. It does not assert that a physical ultraviolet cutoff is known, and it does not give the full fixed-k envelope for every second moment. No original screen, preregistered gate, core implementation or paper is changed here. The original unrestricted-k theorem remains valid.

## 1. The endpoint-compensated family is not generally sharp at fixed k

For equal mass N and zero flux, with Q+=Q-=Na, the previously used endpoint-compensated pair has longitudinal growing rate

    gamma_family² = 2N(1-a)ka - k²a².

This pair attains the unrestricted-k maximum at its optimizing k. Its formula must not be called a sharp envelope at other k.

A compact exact counterexample has N=1, a=3/50 and k=109/75. Set

    mu+ = (6/7) delta_(1/10) + (1/7) delta_(-3/5),
    mu- = reflection of mu+.

Both species have mass one, flux zero and second number moment 3/50. Direct substitution into the ORIGINAL longitudinal dispersion factor gives a root

    omega = i*109/250.

Its squared growth is 11881/62500 = 297025/1562500, strictly above the endpoint-family value 244269/1562500 at exactly the same k and supplied moments. The new root is not an approximate grid artifact. It also saturates the global fixed-k upper bound below, so its maximum over all real frequencies is proved.

The exploratory file finite-k-exploration.py first tested 20 explicitly listed LP cases. Three finite-grid feasible outputs exceeded the endpoint-family rate. Grid infeasibility was not treated as continuous infeasibility. A separate rational bracket counterexample at a=1/10,k=1/2 is retained in finite-k-counterexample.py/results. The compact rational root above is stronger and simpler.

## 2. A necessary Gram disk with the second moment supplied

Let 0<a<1, N>0 and mu be any positive measure on [-1,1] with moments N,0,Na. For z=x+i y, y>0, write

    S = integral dmu(u)/(z-u),
    U = (z²-1)S - Nz.

The Gram matrix of the functions 1,u,1/(z-u) is positive semidefinite. Its Schur complement gives the elementary necessary inequality

    -Im S/y >= |S|²/N + |zS-N|²/(Na).

Complete the square. S lies in a disk with center

    C_S = N[x-i(y+a/(2y))]/(a+x²+y²)

and squared radius

    r_S² = N²a²/[4y²(a+x²+y²)²].

The affine U map therefore puts every attainable U in a disk with radius

    r_U(x,y) = Na |z²-1| / [2y(a+x²+y²)].

Only necessity is needed; no sufficiency theorem for a two-disk locus is assumed in this argument.

## 3. The radius is largest at zero real frequency

For every x real and y>0,

    r_U(x,y) <= Na(1+y²)/[2y(a+y²)] = r_U(0,y).

Here is an explicit algebraic proof. Put t=y² and d=x². After multiplying by positive squared denominators, the difference between the squared right and left sides is proportional to

    (1+t)²(a+t+d)² - (a+t)²[(t+d-1)²+4t]
      = d { 2(a+t)[1+a+(3-a)t] + (1-a)(1+a+2t)d } >= 0.

For 0<a<1 it is strictly positive when x is nonzero. This argument bounds every real frequency, without a frequency scan.

## 4. A global fixed-k cubic bound

The signed zeroth and first moments vanish. Consequently the longitudinal root condition requires U+ - U- = k and the transverse condition requires -2k. Both species have the same necessary Gram disk. A difference of two disk points has magnitude at most its diameter. Writing c=1 for longitudinal and c=2 for transverse gives

    c*k <= Na(1+y²)/[y(a+y²)],       y=gamma/k,

or

    gamma <= (Na/c)(k²+gamma²)/(a k²+gamma²).

Let Gamma_c(k) be the unique positive solution of

    Gamma³ - (Na/c) Gamma² + a k² Gamma - (Na/c) k² = 0.

Then EVERY growing root, at any real frequency, obeys gamma<=Gamma_c(k). Uniqueness is immediate from the strictly increasing function

    f(gamma)=gamma-(Na/c)(k²+gamma²)/(a k²+gamma²),
    f'(gamma)=1+2(Na/c)gamma k²(1-a)/(a k²+gamma²)² > 0.

Its value at zero is negative and it tends to infinity. The bound does not rely on trusting a cubic root finder. Other independently proved bounds, such as gamma_L<=N(1-a), can still be tighter outside the sharp regime below.

## 5. Exact attainment in a low-second-moment regime

Take a y>0 satisfying

    y + a/y <= 1-a.

Put d=y+a/y. Let p,q be the positive real numbers with pq=a and q-p=d:

    p = [sqrt(d²+4a)-d]/2,
    q = [sqrt(d²+4a)+d]/2.

The condition is exactly q<=1; hence both p and q are valid polar cosines. The positive measures

    mu+ = N*q/(p+q) delta_p + N*p/(p+q) delta_(-q),
    mu- = reflection of mu+

have mass N, zero flux and second moment Na for each species. Set

    k_L = Na(1+y²)/[y(a+y²)],
    gamma_L = k_L*y.

Direct original-dispersion substitution gives D_L(i gamma_L,k_L)=0. For example, the exact integral is

    I0 = -2N k³ a d / [(omega²-k²p²)(omega²-k²q²)].

Substituting pq=a, q-p=y+a/y and omega=i*k*y verifies the root. This growth saturates the global disk-diameter bound, so it is the maximum over ALL real frequencies and compatible positive measures at that fixed k.

For the transverse branch interchange the species and take k_T=k_L/2 and gamma_T=gamma_L/2. It saturates the c=2 bound and directly gives T=0. The species remain in the same equal-Q moment class.

The support condition has solutions only when a<=3-2sqrt(2). More explicitly, y must lie between the roots

    y± = [(1-a) ± sqrt(1-6a+a²)]/2.

Thus this is a sharp subregime, not a solution for every a and k. The 3/50 rational example in Section 1 has p=1/10,q=3/5,y=3/10. The same atoms also give another exact sharp example at y=1/5,k=78/25,gamma=78/125.

## 6. A genuinely sharp finite-K statement where attainment holds

The cubic upper cap Gamma_c(k) increases strictly with k. Indeed, at fixed gamma>0,

    partial f/partial k = -2(Na/c)k gamma²(1-a)/(a k²+gamma²)² < 0,

while partial f/partial gamma>0. Thus dGamma_c/dk>0.

If K is a wave number in the attained subregime from Section 5, then every root with 0<k<=K has gamma<=Gamma_c(k)<=Gamma_c(K), and the constructed distribution at k=K attains equality. Therefore Gamma_c(K) is the sharp maximum over that finite interval. This conclusion does not require the bound to be sharp at every smaller k.

For general K outside this subregime, the complete finite-K answer remains open here. When K is at least the previously proved unrestricted-k attaining threshold, the unrestricted-k cap is of course attainable within the interval. No full transition formula or universal physical cutoff is claimed.

## 7. Exact independent checks and claim limits

Run `python3 -B -O finite-k-cap-checks.py`. Three test methods pass:

- 1,000 rational checks of the displayed factored all-real-frequency radius inequality;
- 400 exact Gram-Schur/disk identities;
- 72 exact sharp attainers, evaluating the ORIGINAL I0,I1,I2 and D_L/T, with all positive weights, support, mass, flux and second moments checked. Sixty strictly exceed the endpoint family at the same k.

These use only standard-library Fraction arithmetic plus the existing independent complex arithmetic helper, with no production or numerical imports. The rational examples use p=1/m,q=6/m,y=2/m or3/m, m=6,7,8,10,20,100, at three masses and both sectors. `results/finite-k-cap-checks.json` records their scope and source hash. Finite checks corroborate the written proof; they do not replace it.

The generic Gram disk, its diameter inequality and cubic root algebra are established mathematical tools. Priority or high-impact novelty of this restricted sharp extremum is not established by this note. Positive atomic attainment is in the measure class. No smooth-density exact attainment, realistic transport snapshot, nonlinear conversion, arbitrary three-dimensional propagation or astrophysical wave-number cutoff is inferred.

## 8. The finite-K degenerate limit is continuous, with a cube-root scale

Keep finite K,N fixed and let a=Q/N tend to zero through positive values. The unique cubic root satisfies

    Gamma_c(K) ~ (Na K²/c)^(1/3) = (Q K²/c)^(1/3).

For a direct derivation, set A=Na/c and divide the cubic by A K², writing Gamma=(A K²)^(1/3) h. Its equation becomes h³ plus terms whose coefficients tend to zero, minus one, equal to zero. The positive root therefore has h tending to one. Equivalently one can sandwich the monotone function f at (1±epsilon)(A K²)^(1/3) for any fixed epsilon>0.

Then y=Gamma/K tends to zero as a^(1/3), and a/y tends to zero as a^(2/3). Consequently the support condition y+a/y<=1-a eventually holds. Thus the asymptotic is the sharp finite-K maximum, not merely an upper-bound scaling.

This maximum tends to zero continuously at the forced-stable a=0 corner. The previous unrestricted-k optimization can instead approach N because its maximizing wave number diverges. The endpoint-compensated family at fixed K has only square-root-in-Q growth, so it misses this larger cube-root regime. At the same finite K the transverse/longitudinal maximum ratio tends to 2^(-1/3), rather than the one-half ratio of the unrestricted-k maxima. These statements concern the mathematical uncertainty class with an imposed finite interval; they do not assert an astrophysical cutoff or smooth-profile realization.
