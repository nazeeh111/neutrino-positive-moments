# Sharp nonaxial cap for purely growing modes

Primary derivation and exact algebraic checks. Independent review has not been completed because the reviewer returned a usage-limit error. The later full-frequency extension is now included in Section7 of the research draft; this earlier imaginary-axis derivation is retained as a historical step. It establishes a restricted mathematical extension, not scientific priority or practical astrophysical value.

## Statement and limits

Let each positive species measure on the unit sphere have mass N>0, zero vector flux and second tensor `(N/3) I`. In the collisionless fast model and the shifted frequency convention of this repository, any **purely imaginary** growing root `omega=i gamma`, gamma>0, at any real wave vector obeys

```
gamma <= 2N/3.
```

No axial symmetry or inversion relation between species is required. The cap is sharp: the existing axial isotropic attainers are included in this larger class. This is an all-wave-number cap, not the complete fixed-k envelope. It does not cover roots with nonzero real frequency; the previously proved safe cap N remains the accepted general all-real-frequency bound.

## Reduce the full polarization without discarding angular structure

Choose z along the wave vector, k>0; reflect the axis for k<0. Write `u=(v_x,v_y)`, `r=omega/k`, and `G=mu+ - mu-`. Equal species moments give zero signed moments through degree two. Define

```
S = integral dG/(omega-kz),
T = integral u dG/(omega-kz),
R = integral u u^T dG/(omega-kz).
```

In coordinates `(0,z,x,y)`, the original polarization matrix is

```
Pi = [[1+S, r*S, T^T],
      [r*S, -1+r²*S, r*T^T],
      [T, r*T, -I+R]].
```

The identities in the z row use the vanishing signed mass and vector flux, not axial symmetry. If Pi has a null vector `(a0,az,b)`, subtracting r times its first row from its second gives `az=-r*a0`. Setting `q=(1-r²)*a0` yields the equivalent reduced system

```
M [q,b]^T = 0,
M = [[1/(1-r²)+S, T^T], [T, -I+R]].
```

For r=i y, y=gamma/k>0, this change is invertible. Put `J=diag(i,1,1)` and `H=-J M J`. Congruence preserves singularity. Because M is complex symmetric, the Hermitian part of H is its entrywise real part. Define the positive unit measure

```
nu = (mu+ + reflected(mu-))/(2N),
reflected(mu-)(A) = mu-(-A).
```

It has zero mean and second tensor I/3 even when the species are independent. Each kernel contributing to the Hermitian part is odd under full inversion. With `d=z²+y²`, `a=1/(1+y²)` and `t=2N/(3 gamma)`, direct calculation gives

```
Herm(H) = (1-t)*diag(a,I) + t*E_nu[P(v)],
```

where the pointwise positive matrix P is given next. This averaging identity concerns the Hermitian part used for exclusion; it does not identify the original dispersion relation with that of an inversion-paired background.

## Pointwise positive matrix certificate

Define

```
c = 3/(2*(1+y²)),
alpha = c*(1+2*y²), beta = -c*y,
q(z) = c*(1-z²+2*y*z),

P00 = q(z) - 3*y*z/d,
P0u = (-3*y²/d + alpha + beta*z)*u^T,
Puu = (3/2)*z²*I + (3/2 + 3*y*z/d)*u*u^T.
```

The supplied moments imply `E[q]=a`, `E[u]=E[z*u]=0`, `E[z²]=1/3` and `E[u*u^T]=I/3`. These establish the averaging identity above with every mixed entry retained.

To prove P is positive semidefinite, let `s=1-z²=|u|²` and `B=y+(1+y²)*z-y*z²`. Algebra gives

```
P00 = c*(z-y)²*s/d,
P0u = c*(z-y)*B*u^T/d.
```

In the transverse direction perpendicular to u, the eigenvalue is `(3/2)z² >= 0`. In the direction parallel to u, the lower diagonal entry is

```
L = 3/2 + 3*y*z*s/d.
```

It is nonnegative because

```
(2/3)*d*L = d+2*y*z*s
             = (z+y*s)² + y²*z²*(2-z²) >= 0.
```

The determinant of the remaining 2x2 block is

```
P00*L - c²*(z-y)²*B²*s/d²
    = c²*y²*(z-y)²*s²/d >= 0.
```

Its key polynomial identity is

```
(1+y²)*(d+2*y*z*s) - B² = y²*s*d.
```

Thus P is positive semidefinite everywhere on the sphere, including the poles by direct substitution. No angular grid, optimization tolerance or first-order approximation enters this proof.

For gamma>2N/3, `0<t<1`. The positive definite term `(1-t)*diag(a,I)` then makes Herm(H) positive definite. A matrix with positive definite Hermitian part cannot have a null vector, since `Re(x† H x)>0` for nonzero x. Hence Pi cannot be singular above the claimed cap. At k=0, the zero signed moments give Pi=eta directly, so there is no growing root there either.

## Sharpness and verification boundary

For N=1, the five-atom plus distribution and its inverted minus distribution in tangent-audit.md have k=2 and omega=2i/3. Their exact full polarization and null determinant are already checked in exact_tangent_certificate.py. Scaling both masses, k and omega by N proves sharpness for every N>0. The broader axial plateau family also attains this cap at k>=2N, but is not needed for sharpness here.

`exact_nonaxial_imaginary_bound.py` checks the polynomial identities as bivariate coefficient equalities, constructs full polarization matrices from rational isotropic atomic measures with independently chosen species, and verifies the reduction and Hermitian-part identity using exact complex arithmetic. These finite construction checks support implementation correctness; the pointwise proof supplies the universal quantifier.

The earlier local tangent certificate remains correct but is no longer the strongest evidence on the imaginary axis. The remaining global question is whether mixed modes with **nonzero real frequency** can exceed 2N/3. The present proof uses r=i y in its positive-definite reduction and cannot be extended to those modes just by substituting a complex y. No numerical absence claim, novel algorithm, independent review or operational physics gain is asserted.
