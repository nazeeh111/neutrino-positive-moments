# Sharp full-frequency nonaxial isotropic growth cap

**Review status:** primary derivation and exact algebraic checks only. Independent review is unavailable following the recorded reviewer usage-limit failure. This extension is included in Section7 of the November2 research draft; the earlier independent reviews do not cover it. No priority, publication-readiness or practical-impact claim is made.

## Statement

In the repository's local collisionless fast-flavor model, let each positive species measure on the unit sphere have mass N>0, zero vector flux and second tensor `(N/3) I`. For every real wave vector and every growing frequency

```
omega = omega_R + i gamma,  gamma>0,
```

the full, nonaxial polarization determinant can vanish only if

```
gamma <= 2N/3.
```

The cap is sharp. It permits arbitrary angular measures and independently chosen species. It does not require axial symmetry, inversion pairing or zero real frequency. It is an all-wave-number growth cap; the exact fixed-k envelope for this full nonaxial class is still unresolved. This is a mathematical information bound, not a prediction that the actual background attains it.

## Full polarization reduction

For k>0 choose the z axis along the wave vector, write `u=(v_x,v_y)`, `r=omega/k=x+i y`, and use the reduced matrix from imaginary-axis-bound.md:

```
M = [[1/(1-r²)+S, T^T], [T, -I+R]],
S = integral dG/(omega-kz),
T = integral u dG/(omega-kz),
R = integral u u^T dG/(omega-kz),  G=mu+ - mu-.
```

Zero signed mass and flux make singularity of this matrix equivalent to singularity of the original 4x4 polarization. For Im(r)>0, `1-r²` is nonzero. All mixed entries are retained. At k=0, zero signed moments through degree two give the original polarization equal to its nonsingular Minkowski metric, so no growing root occurs.

Choose `p+i q=sqrt(1-r²)` with p>0, and define real parameters

```
h = sqrt(p²+x²),  t = sqrt(h²-1),
a = x/h,         b = p/h.
```

Then `a²+b²=1`, `h²-t²=1`, `b>0`, `t>0`, and

```
x=a*h, y=b*t, p=b*h, q=-a*t.
```

These identities follow from `p²-q²=1-x²+y²` and `p*q=-x*y`. Thus this parameterization covers every upper-half-plane r, not just a discrete set or an imaginary-axis subset.

Put `J=diag(i*(p+i q),1,1)` and `H=-J M J`. Complex congruence preserves singularity. The matrix is complex symmetric, so its Hermitian part is its entrywise real part. With

```
s=1-z²=|u|²,
D=(z-a*h)²+b²*t²,
F_a(v) = 1/D * [[-(z-a*h)*s, -t*(h-a*z)*u^T],
                [-t*(h-a*z)*u, (z-a*h)*u*u^T]],
```

direct expansion gives

```
Herm(H) = I + (1/k) integral F_a(v) dG(v).
```

For the upper-left entry, polynomial division uses
`(1-r²)/(r-z)=(1-z²)/(r-z)-(r+z)`; the last term integrates to zero by the signed moment constraints. Also `F_a(-v)=-F_-a(v)`. Hence, if `nu+=mu+/N` and `nu-=reflected(mu-)/N`, both positive unit measures with isotropic moments,

```
Herm(H) = I + (N/k)*(E_nu+[F_a] + E_nu-[F_-a]).
```

This averaging is only an exclusion argument for the Hermitian part; the original spectra are not identified with those of an inversion-paired background.

## A positive matrix for each species

For each sign of a use its corresponding F and set

```
K = b*h+a*t,
e = (K²-t²)/h²,
alpha = b*(1+2*t²)/h,  beta=-t/h,

Q_a(v) = (3/2) * [[s, (alpha+beta*z)*u^T],
                  [(alpha+beta*z)*u, (K²/h²)*z²*I + e*u*u^T]],
P_a(v) = Q_a(v) + 3*y*F_a(v),  y=b*t.
```

Although e may be negative, P is positive semidefinite. The proof concerns the complete matrix, not the signs of its individual coefficients. Define

```
ell = -t*D + h*((b*h-a*t)*(z-a*h) + (b*h+a*t)*b*t),
A = z-a*h-b*t.
```

Exact polynomial identities, using only `a²+b²=1` and `h²-t²=1`, give

```
P00 = (3/2)*s*A²/D,
P0u = (3/2)*A*ell*u^T/(h*D).
```

In the transverse direction perpendicular to u, the eigenvalue is
`(3/2)*(K²/h²)*z² >= 0`. In the parallel direction, the remaining 2x2 block is

```
(3/(2*D)) * [[s*A², sqrt(s)*A*ell/h],
             [sqrt(s)*A*ell/h, ell²/h²]].
```

It is an outer product with itself, hence positive semidefinite. The key identity for its lower entry is

```
ell² = D*(K²*z² + (K²-t²)*s)
       + 2*h²*b*t*(z-a*h)*s.
```

The pole cases `s=0` follow directly from Q and F, with nonnegative transverse diagonal entries. There are no angular discretization or numerical-tolerance assumptions in this proof.

## Average the certificate and exclude larger growth

For either isotropic unit measure, its moment constraints give

```
E[Q_a] = diag(1, lambda_a, lambda_a),
lambda_a = (K²-t²/2)/h².
```

The two transverse constants need not be equal or separately positive. What matters is their average:

```
lambda_bar = (lambda_a+lambda_-a)/2
           = (b²+t²/2)/h²,
1-lambda_bar = (a²+t²/2)/h² >= 0.
```

With `tau=2N/(3 gamma)` and `gamma=k*y`, the exact decomposition is

```
Herm(H) = (1-tau)*I
          + (tau/2)*(E_nu+[P_a] + E_nu-[P_-a])
          + tau*diag(0,1-lambda_bar,1-lambda_bar).
```

If gamma>2N/3, then `0<tau<1`. The first term is positive definite and the other terms are positive semidefinite. H therefore has positive definite Hermitian part and cannot be singular: a nonzero null vector would contradict `Re(w† H w)>0`. This excludes the original polarization root and proves the cap for every real frequency part.

The imaginary-axis proof is a special case, a=0. The new argument needs separate species certificates and their averaged transverse constants; merely substituting a complex y in the earlier certificate is invalid.

## Exact oscillating attainers

Let gamma=2N/3 and choose any real omega_R and positive k satisfying

```
k >= 3*abs(omega_R)+2N.
```

Set `u+=(omega_R+gamma)/k`, `u-=(omega_R-gamma)/k`; both lie in `[-1/3,1/3]`. For species s, place mass

```
m_s = (2N/3)/(1-u_s²)
```

uniformly on the latitude z=u_s, or on four azimuthal points with the same conditional first and second transverse moments. Put the remaining endpoint masses

```
N*(1/3-u_s)/(2*(1-u_s)) at z=+1,
N*(1/3+u_s)/(2*(1+u_s)) at z=-1.
```

They are nonnegative and give exactly the stated full isotropic moments. For the axial longitudinal transform, only the interior atom contributes to the `(z²-1)`-weighted difference. At `omega=omega_R+i gamma`, the two interior denominators have offsets `-gamma+i gamma` and `+gamma+i gamma`. Their difference gives the exact longitudinal root condition. The full transverse factors remain nonsingular for the explicit checked example below.

For N=1, k=4 and omega_R=2/3, one species uses latitude1/3 with mass3/4 and south-pole mass1/4; the other is the six-axis octahedral measure with weights1/6. The original full polarization has the exact root `omega=2/3+2i/3`. Consequently a proposed universal cap strictly below2N/3 at every nonzero real frequency would be false. This attainer supplements, rather than replaces, the earlier omega_R=0 witness.

## Verification and remaining requirements

`exact_nonaxial_full_bound.py` checks the polynomial identities after exact reduction by the two parameter constraints, directly assembles the original polarization for independent rationally rotated species, verifies the full Hermitian decomposition and strict positivity above the cap, and replays an exact oscillating attainer. Finite cases verify implementation details; the outer-product and averaging argument carries the universal quantifier.

Independent review, a refreshed priority search for this full nonaxial theorem, incorporation into a reviewed paper, and a useful physical decision remain outstanding. The result does not establish the complete finite-k nonaxial envelope, general unequal/nonzero-flux moment bounds, universal instability of a background, or exceptional practical advantage. The original high-impact research goal is not complete.
