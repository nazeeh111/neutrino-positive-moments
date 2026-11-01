# Exact axial envelope with separately fixed second number moments

This independent derivation uses no proposed numerical optimum. Its scope is
the two stated axial dispersion factors, positive species measures on [-1,1],
equal number mass N>0, zero flux for each species, and separate second number
moments Q+ and Q- in [0,N]. It does not optimize a nonaxial wave vector.

## Result

Write R±=N-Q±. If Q++Q->0 and R+R->0, the exact attained maximal growth is

\[
\gamma_{L,\max}=\sqrt{R_+R_-},\qquad
\gamma_{T,\max}=\frac12\sqrt{R_+R_-}.
\]

Thus the maximum over the union of the two branches is the longitudinal
value. If both Q values vanish, or if either R value vanishes, there are no
upper-half-plane roots in either factor. For those stable controls the
reported growth cap is zero; a maximum over a nonempty unstable-root set
is not claimed. In particular, the first formula cannot be extended to the
single corner Q+=Q-=0, where it would incorrectly give N.

## Reduction from the original factors

Let G=mu+ - mu-, w=x+iy with y>0, k>0, and

\[
I_j=\int\frac{u^j\,dG(u)}{w-ku},\qquad
H_\pm=(1-u^2)\mu_\pm,\qquad
J=\int\frac{d(H_+-H_-)(u)}{w-ku}=I_0-I_2.
\]

Signed mass and signed flux vanish, so polynomial division gives
I1=(w/k)I0 and I2=(w/k)I1. Consequently the original, unmultiplied factors obey

\[
D=(1+I_0)(-1+I_2)-I_1^2=-1-J,
\qquad T=-1+J/2.
\]

For y>0 none of these integrals has a pole on the integration support.
No denominator clearing or extraneous-root argument is needed. Longitudinal
roots solve J=-1 and transverse roots solve J=2. The positive measures H±
have masses exactly R±, regardless of any unconstrained higher moments.

## Global upper bound at arbitrary real frequency and wave number

The following scalar bound uses positive measures only. If positive H± have
masses R± and

\[
F_+-F_-=c\ne0,\qquad
F_\pm=\int\frac{dH_\pm(u)}{x+iy-ku},
\]

where c is real, then y <= sqrt(R+R-)/|c|. Indeed, for every real u,

\[
\left|\frac1{x+iy-ku}+\frac{i}{2y}\right|=\frac1{2y}.
\]

Integrating and using the triangle inequality gives

\[
\left|F_\pm+\frac{iR_\pm}{2y}\right|\le\frac{R_\pm}{2y},
\quad
\left|c+\frac{i(R_+-R_-)}{2y}\right|
\le\frac{R_++R_-}{2y}.
\]

Squaring yields c² <= R+R-/y². This proves the two bounds for all admissible
measures, all real x, and all k>0, including measures with continuous or
singular parts. It also proves absence of unstable roots if R+R-=0.
When Q+=Q-=0, positivity forces mu+=mu-=N delta_0, so G=0 and D=T=-1.
These arguments exhaust the stable controls listed above.

## Explicit feasible longitudinal attainer

Set a=Q+/N, b=Q-/N and s=a+b. Suppose s>0 and a,b<1. Choose

\[
\mu_+=\frac{N}{1+a}\delta_a+\frac{Na}{1+a}\delta_{-1},
\qquad
\mu_-=\frac{N}{1+b}\delta_{-b}+\frac{Nb}{1+b}\delta_{+1}.
\]

Both measures have positive mass N, zero flux, and their required separate
second moments Na and Nb. Their weighted versions are simply
H+=R+ delta_a and H-=R- delta_-b. In particular, the endpoint masses are
number masses and do not survive multiplication by 1-u².

The following finite k and complex frequency attain the upper bound:

\[
k_L=\frac{N(2-s)}s>0,\qquad
x_L=\frac{N(a-b)}s,\qquad
w_L=x_L+i\sqrt{R_+R_-}.
\]

One can verify the scalar equation directly, since
xL-kL*a=-R+ and xL+kL*b=R-. With y²=R+R-,

\[
J=\frac{R_+}{-R_++iy}-\frac{R_-}{R_-+iy}=-1.
\]

This is an original longitudinal root because of the exact identity above.
For a=0 or b=0 the corresponding zero-weight endpoint is omitted and that
species is precisely N delta_0; the construction remains valid. Values
a=1 or b=1 belong to the stable controls, so the code does not output
spurious real-pole 'attainers' at those boundaries.

## Explicit transverse attainer

Reflect both longitudinal measures, replacing each node u by -u. The number
moments and flux constraints are preserved. Set

\[
k_T=k_L/2,\qquad x_T=-x_L/2,\qquad
w_T=x_T+\frac{i}{2}\sqrt{R_+R_-}.
\]

The interior denominators are now (R++i sqrt(R+R-))/2 and
(-R-+i sqrt(R+R-))/2, respectively. Direct substitution gives J=2 and T=0.
The measures attaining the two branch maxima need not be the same.

## Interpretation of the corner

For arbitrarily small positive s the optimum tends to N while kL diverges
like 2N/s. At s=0 the signed measure vanishes exactly and no instability
exists. This discontinuity is consistent with optimization over unbounded
k. A bounded-k problem is different and is not covered by the stated cap.

## Reusable exact implementation and checks

`second_moment_bounds.py` provides `second_moment_bounds(N,Qplus,Qminus)`
and the equivalent `independent_r_bounds(N,Rplus,Rminus)`. Inputs may be
rational strings or `fractions.Fraction`. All outputs are exact rational
quantities: squared growth caps, k, real frequency, and positive atom
weights/nodes. Irrational growth is represented by its rational square.
Stable controls return a reason and no UHP attainer.

The stdlib-only test file `test_second_moment_bounds.py` checks 585 unstable
parameter triples at three mass scales and 15 rational normalized moment
nodes. Both branch attainers are checked, giving 1,170 original-factor root
substitutions. Each substitution evaluates I0,I1,I2 from the returned NUMBER
atoms and computes `(1+I0)*(-1+I2)-I1*I1` or `-1+(I0-I2)/2` directly.
Arithmetic is exact in Q[t]/(t²+gamma²), so every residual is exactly zero,
even when gamma is irrational. The same cases check positive weights,
support, number mass, zero flux, separate second moments and positive k.
Additional tests cover stable controls, the zero-central-node case,
equivalence of the Q and R interfaces, rational strings and invalid inputs.
The tests verify the implementation and finite examples; the global
optimality and stable-case statements follow from the proof above.

Run with `python3 -m unittest -v test_second_moment_bounds.py`, or invoke
`python3 second_moment_bounds.py 1 1/3 2/3` for a JSON attainer. Source hashes
and executed test counts are recorded in
`results/second_moment_bounds_validation.json`.
