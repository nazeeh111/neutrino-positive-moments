# Independent extension: sharp growth cap over every real frequency

The earlier purely imaginary attainable envelopes extend to **all real parts of frequency**, under exactly the same equal-mass, zero-flux, axial model. This is a written proof, corroborated by `offaxis-independent.py`; it is not inferred from a frequency scan. Priority and physical utility remain separate unresolved questions.

## Statement and scope

Let each positive angular measure have mass N>0 and first moment zero, with support u in [-1,1]. Assume the already specified local collisionless axial fast-flavor model and k>0. A longitudinal root omega=Omega+i gamma, gamma>0, for some compatible pair obeys

    gamma^2 <= N^2 k/(k+2N).

A transverse-sector root obeys

    gamma^2 <= N^2 k/[4(k+N)].

Both are sharp over **all** compatible measures and **all** real Omega. Equality requires Omega=0, and the explicit positive two-atom measures below attain equality. Every positive rate below each cap is attainable with Omega=0. The result is existential over completions. It neither predicts the actual angular distribution nor says that every completion is unstable. It does not supply the full two-dimensional frequency region, arbitrary three-dimensional distributions or nonzero common flux.

## 1. A transformed species lens from positive-moment determinants

Write z=omega/k=x+i y, y>0, and let

    S = integral dmu(u)/(z-u),       U=(z^2-1) S - N z.

The complex affine map is invertible because z is in the upper half-plane. Set U=p+i v. Reweight dnu=dmu/|z-u|^2. Its four real moments are

    t0=-Im(S)/y,
    t1=x t0-Re(S),
    t2=N+2x t1-(x^2+y^2)t0,
    t3=2x t2-(x^2+y^2)t1.

For sigma in {+1,-1}, positivity of dmu on [-1,1] makes the Gram matrix

    H_sigma = [[t0+sigma t1, t1+sigma t2],
               [t1+sigma t2, t2+sigma t3]]

positive semidefinite. Direct substitution, without the production disk formulas, gives the identity

    det H_sigma = [R_sigma^2-|U-C_sigma|^2] / |z-sigma|^2,

where

    C_sigma = -sigma N/2 + i N(1-sigma x)/(2y),
    R_sigma^2 = N^2[(x-sigma)^2+y^2]/(4y^2).

The denominator is strictly positive. Consequently each attainable U lies in both closed disks. Sufficiency of these disks is not needed for the upper bound; only the elementary Gram-matrix necessity is used. The tests evaluate the original moment Gram determinants and this formula at 1,200 independently generated rational choices, checking exact equality with standard-library Fraction arithmetic.

## 2. The root condition is a horizontal chord

With equal species masses and zero first moments, g0=g1=0. Direct recurrence reduction of the longitudinal determinant gives

    (z^2-1)(S_plus-S_minus) = k.

The transverse condition gives the same left-hand side equal to -2k. The subtraction of Nz in U cancels between species. Thus two points of the **same** species lens must differ horizontally by t=k for longitudinal roots or t=2k for transverse roots. The two points can be ordered either way, so only t>0 matters.

## 3. Exact bound for every horizontal chord, including unequal disk radii

Suppose the two chord endpoints are

    U_1=m-h+i v,       U_2=m+h+i v,       h=t/2.

Write q=v-N/(2y). For a disk centered at C_sigma, the larger of the two squared horizontal offsets is (|m+sigma N/2|+h)^2. Requiring **both** chord endpoints to lie inside **each** disk therefore gives

    (|m+sigma N/2|+h)^2 + (q+sigma N x/(2y))^2
       <= N^2[(x-sigma)^2+y^2]/(4y^2).

After expansion, this is

    (|m+sigma N/2|+h)^2 + q^2
       + sigma N x q/y + sigma N^2 x/(2y^2)
       <= N^2(1+y^2)/(4y^2).

Average the sigma=+1 and sigma=-1 inequalities. Both terms involving x cancel exactly, even though the original radii are unequal:

    {(|m+N/2|+h)^2+(|m-N/2|+h)^2}/2 + q^2
       <= N^2(1+y^2)/(4y^2).

Let a=|m+N/2|+h and b=|m-N/2|+h. The triangle inequality gives a+b>=N+2h, and (a^2+b^2)/2 >= ((a+b)/2)^2. Therefore

    (N/2+h)^2 <= N^2(1+y^2)/(4y^2),
    t <= N[sqrt(1+y^2)/y-1],
    y^2 <= N^2/[t(t+2N)].

Substitute y=gamma/k and t=k or 2k to obtain the stated caps. No optimization, angular-grid approximation, frequency sampling or assumption Omega=0 occurs in this proof.

## 4. Equality forces zero real frequency

At equality, q=0, a=b and a+b=N+2h. Since N>0, a=b implies |m+N/2|=|m-N/2| and hence m=0. At m=q=0, equality in the averaged bound makes the two individual disk slacks equal to

    -N^2 x/(2y^2),       +N^2 x/(2y^2).

Both can be nonnegative only when x=0. Therefore every root attaining the global cap has Re(omega)=kx=0. No uniqueness statement about the distributions is made.

## 5. Explicit attainers and all sub-cap rates

For a in (0,1), set

    mu_plus = N/(1+a) [delta_a + a delta_{-1}],
    mu_minus = N/(1+a) [delta_{-a} + a delta_1].

Each measure has mass N and zero first moment. Direct substitution into the dispersion tensor gives

    D_L = -1 - 2Nka(1-a)/(omega^2-k^2 a^2),
    gamma_L^2 = 2Nka-k(k+2N)a^2.

The parabola has maximum N^2 k/(k+2N) at a=N/(k+2N). For transverse roots swap the two species. Direct substitution then gives

    T = -1 - Nka(1-a)/(omega^2-k^2 a^2),
    gamma_T^2 = Nka-k(k+N)a^2.

The maximum is N^2 k/[4(k+N)] at a=N/[2(k+N)]. Varying a continuously from zero to the maximizing value realizes every positive rate up to the maximum. These positive atomic distributions belong to the stated uncertainty class. Smooth distributions approach the bound by convolution and root continuity; the present theorem claims exact attainment only in the measure class. Existing separately validated smooth examples remain valid evidence of regularity and low-moment ambiguity, not new global-cap proofs.

## Exact checks and what they establish

Run `python3 -O offaxis-independent.py` from this task directory. Four tests pass:

- 1,200 exact determinant/circle identities from the independent moment Gram matrices;
- 1,000 exact horizontal-chord cancellation and elementary lower-bound checks;
- sharp-bound equality residuals at three rational configurations;
- 192 exact off-axis impossibility certificates checked by `independent_verifier.py`, with real parts as small as 1e-5, growth at the bound or above it, and common scales 1e-6,1,2,1e6.

The checker uses no production geometry, optimizer, NumPy or SciPy. The 192 examples corroborate the general proof; they do not replace it. Existing independent attainer tests check the direct beam dispersion and all-subthreshold construction. This extension does not alter any earlier predeclared screen, gate or failure. No paper or production core was modified for this derivation.
