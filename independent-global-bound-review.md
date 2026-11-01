# Independent mathematical review of the global growth cap

Reviewed 1 November 2026. Scope: correctness of `offaxis-independent.md`, its checking script and saved report, against the tensor definitions in `paper/draft.typ`. This review does not assess novelty, physical importance or publication readiness. No reviewed source file was changed.

## Verdict

**Accept the global sharp-bound theorem in the stated measure class.** For each fixed real `k>0` and `N>0`, consider all pairs of finite positive measures on `[-1,1]` with individual mass `N` and individual first moment zero. In the specified axial, collisionless, fast model, every upper-half-plane longitudinal root satisfies

    (Im omega)^2 <= N^2 k/(k+2N),

and every transverse root satisfies

    (Im omega)^2 <= N^2 k/[4(k+N)].

These are maxima over both the compatible pair of measures and arbitrary `Re omega`. They are attained at `Re omega=0` by the displayed atomic pairs. Every positive imaginary rate up to either maximum is attained by some compatible pair at zero real frequency. Conversely, attaining either maximum forces zero real frequency. The argument is a proof, not an extrapolation from the finite tests.

The paper source currently states only the previous imaginary-axis result and explicitly disclaims a global cap. That older conservative scope is consistent with its earlier proof, but it does not yet incorporate this extension.

## Algebra checked independently

For `z=x+i y`, `y>0`, the proposed map

    U=(z^2-1)S-Nz

is invertible: `z` cannot equal either real number `+1` or `-1`. Zero first moment also gives the useful integral identity

    U = - integral (1-u^2)/(z-u) dmu(u).

This identity fixes both the normalization and transform sign. It follows by writing

    (z^2-1)/(z-u)-z = u + (u^2-1)/(z-u)

and integrating the `u` term to zero.

The four proposed moments are exactly the moments of `dnu=dmu/|z-u|^2`. Each `H_sigma` is therefore the Gram matrix of `(1,u)` for the positive measure `(1+sigma u)dnu`; in particular its determinant is nonnegative. The circle identity is correct:

    det H_sigma = (R_sigma^2-|U-C_sigma|^2)/|z-sigma|^2,
    C_sigma = -sigma N/2 + i N(1-sigma x)/(2y),
    R_sigma^2 = N^2[(x-sigma)^2+y^2]/(4y^2).

I checked the identity without importing the supplied verifier or geometry. Let `A=x^2-y^2-1`, `B=2xy`, `d=A^2+B^2`. Substituting

    Re S = [(p+Nx)A+(v+Ny)B]/d,
    Im S = [(v+Ny)A-(p+Nx)B]/d

into the four moment definitions and clearing denominators produces an identically zero polynomial for both signs. An independent standard-library multivariate polynomial expansion checked every coefficient, with zero residual coefficients in both cases. This is an algebraic identity check over formal variables, not a random parameter sample. SymPy was unavailable; it was not installed or needed.

The tensor reduction is also correct. Since signed zeroth and first moments vanish, `I1=z I0` and `I2=z^2 I0`. Thus

    D_L = -1+(z^2-1)I0,
    T   = -1-(z^2-1)I0/2,
    I0=(S_plus-S_minus)/k.

Consequently the root conditions are respectively

    U_plus-U_minus=k,       U_plus-U_minus=-2k.

No complex denominator has been canceled at a possible upper-half-plane zero. Because the individual moment constraints agree, the two species have the same necessary lens. Their root values must be the endpoints of a horizontal segment of length `t=k` or `t=2k`.

## The global geometric estimate is valid

Let the segment endpoints be `m-h+i v` and `m+h+i v`, where `h=t/2`, and put `q=v-N/(2y)`. In a disk centered at real coordinate `-sigma N/2`, the farther endpoint has squared horizontal displacement `(|m+sigma N/2|+h)^2`. This handles both endpoints exactly, even if the midpoint is outside the interval between the centers.

Expanding the two necessary disk inequalities and averaging cancels both terms involving `x`. The resulting necessary inequality is

    [ (|m+N/2|+h)^2 + (|m-N/2|+h)^2 ]/2 + q^2
        <= N^2(1+y^2)/(4y^2).

The bracketed average is at least `(N/2+h)^2`: the sum of the two absolute values is at least `N`, and the mean of two squares is at least their squared mean. All subsequent squaring and division use positive quantities. Therefore

    y^2 <= N^2/[t(t+2N)].

With `y=Im(omega)/k`, this gives exactly the two claimed growth caps, for every real `x`.

Only necessity of the two disks is used. A possible gap between determinant disks and positive-semidefinite Gram matrices would not invalidate this upper bound, since every actual positive measure is already known to satisfy the determinant inequalities. Likewise, sufficiency of a general horizontal chord is not assumed. Sharpness and all-sub-cap attainability are provided separately by actual measures.

## Equality and attainability

Equality in the cap forces `q=0`, equality of the two squared terms, and equality in their sum bound. Since `N>0`, equality of `|m+N/2|` and `|m-N/2|` forces `m=0`. At this midpoint and height, the two least endpoint slacks are `-N^2 x/(2y^2)` and `+N^2 x/(2y^2)`. Their simultaneous nonnegativity forces `x=0`. This is a necessity for every maximizing root, not merely a property of the constructed family.

For the displayed family, the positive weights sum to `N` and the first moment vanishes separately in each species. Direct beam summation gives the stated longitudinal determinant, and after swapping the species the stated transverse factor. Their growing roots obey

    gamma_L^2=2Nka-k(k+2N)a^2,
    gamma_T^2=Nka-k(k+N)a^2.

The maximizing values `a=N/(k+2N)` and `a=N/[2(k+N)]` lie strictly inside `(0,1)`. Each parabola increases continuously from zero on the interval from `a=0` to its optimizer. Hence every strictly positive rate at or below the cap has a legitimate attaining pair. The zero-rate endpoint is not being counted as an unstable upper-half-plane root.

There is a minor wording problem in the source: “The measures attaining it need not be unique” could be read as asserting that nonunique maximizing pairs exist. The proof there does not establish that. Replace it with “No uniqueness statement is made.” This does not affect the accepted theorem. In fact a further argument is available: each maximizing transform saturates the appropriate one of the two Gram determinants; its null linear polynomial confines the measure to one endpoint and one interior atom. Together with mass and first moment, this appears to pin down precisely the displayed pair in this fixed model. Such a uniqueness result should be stated and proved separately if desired.

Exact attainment is established for positive measures including atoms. Smooth approximation can preserve mass and zero flux by rotationally symmetric smoothing on the sphere; the atomic growing roots are simple, so sufficiently small smoothing preserves nearby upper-half-plane roots. Thus the same values are suprema over that smooth class. The current argument does not establish smooth exact attainment at the cap. It also does not give the full attainable two-dimensional frequency region.

## Verification and limits

I independently loaded the checking module through `importlib` and ran its four tests with `python3 -B -O`, bypassing its report-writing entry point. All four tests passed. The saved JSON accurately describes that suite. Those 1,200 determinant evaluations, 1,000 chord checks and 192 off-axis rejections are useful corroboration; none is the reason the universal theorem is accepted. Acceptance rests on the Gram necessity, exact identities, global inequality, equality analysis and explicit attaining measures above.

The result remains restricted to the prescribed axial model, equal individual masses, zero individual first moments, fixed positive axial wave number and positive measure uncertainty class. It does not predict the distribution selected by a physical system, guarantee instability for every completion, cover additional independently fixed moments, or establish a theorem for arbitrary three-dimensional angular distributions or nonzero common flux.

## Addendum: independently fixed second moments

After the initial review, I also reviewed `independent-second-moment.md`, `second_moment_bounds.py`, its exact test file, and the newly integrated mathematical sections of `paper/draft.typ` (the two sharp-bound sections). **Accept this second theorem with its explicit optimization over all positive wave numbers.** No material mathematical error was found in those paper sections. They now incorporate the first theorem and remove the unsupported nonuniqueness wording discussed above.

With `Q_s=integral u^2 dmu_s` and `R_s=N-Q_s`, the positive weighted measures `lambda_s=(1-u^2)mu_s` have masses exactly `R_s`. Their transform difference `F=integral d(lambda_plus-lambda_minus)/(omega-ku)` satisfies the original, unmultiplied identities

    D_L=-1-F,       T=-1+F/2.

For physical frequency `omega=Omega+i gamma`, the elementary circle identity is

    |1/(Omega+i gamma-ku)+i/(2gamma)| = 1/(2gamma).

Thus the individual transform disks use **gamma**, not the dimensionless `y=gamma/k` from the earlier proof. Averaging over a positive measure gives center `-iR_s/(2gamma)` and radius `R_s/(2gamma)`. The difference disk contains a real value `c` only if

    c^2 + (R_plus-R_minus)^2/(4gamma^2)
        <= (R_plus+R_minus)^2/(4gamma^2),

which is precisely `gamma^2 <= R_plus R_minus/c^2`. The targets `c=-1` and `c=2` give the stated longitudinal and transverse bounds at every wave number and every real frequency.

For the proposed attainer let `a=Q_plus/N`, `b=Q_minus/N`, `s=a+b>0`, with `a,b<1`. The actual number measures have the correct mass, first moment zero, and second moments `Na,Nb`. Multiplying by `1-u^2` leaves only the central atoms, of masses `R_plus,R_minus`, at `a,-b`. At

    k_L=N(2-s)/s,
    Omega_L=N(a-b)/s,
    gamma_L^2=R_plus R_minus,

the real denominator parts are `Omega_L-k_L a=-R_plus` and `Omega_L+k_L b=R_minus`. Therefore direct substitution gives `F=-1`. Reflecting the two measures and using `(k_T,Omega_T,gamma_T)=(k_L/2,-Omega_L/2,gamma_L/2)` yields denominator real parts `R_plus/2,-R_minus/2` and gives `F=2`. These checks prove actual attainment without relying on sufficiency of the relaxed transform disks.

The edge quantifiers are essential and handled correctly. If either `R_s=0`, the disk inequality forbids a nonzero real target at positive gamma. If both `Q_s=0`, positivity forces both number measures to be `N delta_0`, so their signed measure vanishes and there is no instability. In every other admissible case, the displayed finite positive wave number and positive rate exist. A zero `Q` for only one species causes no problem: its zero-weight endpoint is omitted. The discontinuity at the both-zero corner is legitimate because optimizing wave numbers diverge as `a+b` approaches zero.

The theorem gives maxima over all compatible measures, all real frequencies and all `k>0`. It does not give a generally sharp formula for each prescribed wave number. To avoid suggesting the product cap is never attained at fixed k, the paper's sentence “It is not a sharp fixed-k formula” can be written more precisely as “It does not supply the sharp bound for each prescribed k.” This is a wording refinement, not an objection to the theorem.

I reran `python3 -B -O -m unittest -v test_second_moment_bounds.py`: all five tests passed, including 1,170 exact substitutions in the original dispersion factors and independent mass/flux/second-moment checks. As before, the written proof establishes the universal claim; the finite suite checks the implementation and representative attainers.

## Integration record from parent

The earlier statement about an imaginary-axis-only paper describes the source at the first review. The current seven-page draft incorporates both accepted extensions, applies the precise prescribed-k wording above, and makes no distribution-uniqueness assertion. Its current source/PDF hashes are in results/paper-build.json. This integration note is a parent record, not a further independent review.
