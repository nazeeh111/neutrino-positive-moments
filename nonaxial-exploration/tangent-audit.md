# First variation at the isotropic axial attainer

Primary-derived, exact-arithmetic checked, not independently reviewed. This is a local result about a specified root branch. It does not prove the global nonaxial envelope or establish scientific priority.

## Target and exploratory result

TANGENT_PLAN.md fixes the 30 rational sphere nodes, k=2, unit masses and full isotropic moments before the tangent calculation. A numerical eigenvector/linear-program calculation found maximum first-order gain about 4.69e-17, with moment residual 2.78e-17. Its optimizer returned the starting axial weights. Those rounding-scale numbers are not a proof. The exact calculation below supersedes the finite-node restriction and does not depend on the optimizer, its tolerance or a grid.

Let mu+ give weight 3/16 to each `(±2/3,±2/3,1/3)` and weight 1/4 to `(0,0,-1)`. Let mu- be its inversion. Both have mass one, zero flux vector and second tensor I/3. Their conditional first and second transverse moments on the latitude ring equal an azimuthally uniform ring, so their **full 4x4 polarization** is the same as that axial reference at every frequency. This statement does not azimuthally average arbitrary inputs.

For any two admissible positive measures nu+,nu- with those same moments, take the path mu_s(t)=(1-t)mu_s+t nu_s. The root near omega0=2i/3 is simple and hence has an analytic continuation for sufficiently small t. The measure dependence of the polarization entries is affine, and the upper-half-plane denominators are bounded away from zero near omega0.

## Exact derivative

With z=v_z and k=2, put I_j=integral z^j d(mu+-mu-)/(omega-2z). At omega0,

```
(I0,I1,I2)=(-9/10,-3i/10,1/10),
D=(1+I0)(-1+I2)-I1²=0,
partial_omega D=3i/2.
```

Both transverse polarization factors equal -3/2. Thus the derivative of the full determinant is 27i/8 and is nonzero. Its first variation contains only the longitudinal 2x2 block: off-block couplings cannot contribute to first order at this block-diagonal matrix. This is the step allowing arbitrary nonaxial perturbations rather than silently assuming they stay axial.

For a signed perturbation deltaG=delta mu+ - delta mu-, the eigenfrequency derivative is the integral of

```
d(z) = (z²-9+6i*z)/(10*(1+3i*z)),
h(z) = Im d(z) = 3*z*(11-z²)/(10*(1+9*z²)).
```

Therefore gamma'(0)=integral h d(nu+-nu-) - integral h d(mu+-mu-). The latter integral is 2/3.

## Continuous dual certificate

The following exact identity holds on the whole interval:

```
q(z) = (1+3*z/5-z²)/2,
q(z)-h(z) = (1-3*z)²*(1-z²)/(2*(1+9*z²)) >= 0.
```

Every positive unit measure with E[z]=0 and E[z²]=1/3 has E[q]=1/3. Hence E_nu+[h] <= 1/3. Since h is odd, applying the reflected inequality to nu- gives E_nu-[-h] <= 1/3. Adding proves **gamma'(0)<=0**. The two species may vary independently. Full tensor constraints are stronger than the projected constraints used by this certificate.

This controls the continuation of the specified simple mode at k=2; it does not bound another mode or finite changes in the distribution. It is consistent with the accepted axial envelope and should not be counted as a separate novel high-impact algorithm.

## Equality directions and what is still missing

Equality for nu+ requires its z support to lie in `{-1,1/3,1}`. The three moment equations uniquely set their masses to `(1/4,3/4,0)`. Equality for nu- is the reflected statement. For either species, the full zero-flux and isotropic-tensor requirements then fix all first and second transverse moments on the only nonpolar ring. The numerator of each polarization entry has degree at most two in the velocity components, while its denominator depends only on z. Consequently every such equality case has exactly the same full polarization matrix as the base, despite possible higher azimuthal structure. Its dispersion relation is unchanged along the mixture path.

Every fixed target with a different z marginal therefore initially decreases this branch's growth; equality targets leave the full polarization unchanged. This still does not supply a uniform neighborhood of global optimality over all positive measures: dual slack approaches zero near the contact latitudes, and higher-order effects would need uniform control. It also says nothing about a distant completion with growth above 2/3. Do not replace those missing arguments with a larger random scan.

The remaining universal question is still whether gamma<=2/3 for all nonaxial isotropic completions. The accepted safe bound remains gamma<=1. A useful next proof must control mixed longitudinal/transverse polarization at finite amplitude, or a counterexample must depart from this first-order route. Independent review is unavailable following the recorded reviewer usage-limit error.

Run `python3 -B -O nonaxial-exploration/exact_tangent_certificate.py`. It checks the root and full-determinant derivative, the derivative-polynomial coefficients, the global slack factorization and its moment contact exactly. Finite rational evaluations are implementation cross-checks, not the basis of the universal inequality.
