# Sufficient nonaxial frequency exclusion with general species moments

Primary derivation and feasibility checks, not independent review. This extends the pointwise matrix certificate's input scope; it does not establish priority, practical astrophysical utility or a general exact compatibility decision.

## Inputs and conclusion

Each positive measure on the unit sphere supplies its mass N_s, vector flux f_s and full second tensor Q_s. Densities and fluxes may differ and need not be collinear. Work in the local collisionless fast model, with number-weighted moments and the repository's explicit heavy-flavor cancellation assumption. Frequencies and k are reduced quantities; known real matter/current shifts must be treated consistently.

For k>0 along z, define signed moments `g0=N+-N-`, `g1=f+-f-`, `Omega=omega+g0`, `Kz=k+g1z`, and

```
Delta = k*Kz - omega*Omega.
```

If Delta is nonzero, a computable real symmetric3x3 matrix L below is a rigorous sufficient exclusion certificate: **L positive definite implies no positive angular completion of those moments has that frequency.** Failure of this condition means unresolved, not possible or stable. Delta=0 is also unresolved by this reduction. The initial executable uses exact rational parameters for the frequency square root; arbitrary floating frequency inputs and interval-valued measurement data are not implemented.

## Correct reduction for nonzero signed density/current

Let `r=omega/k`, `u=(v_x,v_y)`, and define S,T,R as in the full isotropic proof. The original polarization, in coordinates `(0,z,x,y)`, is

```
Pi = [[1+S, r*S-g0/k, T^T],
      [r*S-g0/k, -1+r²*S-r*g0/k-g1z/k, r*T^T-g1perp^T/k],
      [T, r*T-g1perp/k, -I+R]].
```

Subtracting r times its first row equation from its z row gives the exact Ward relation

```
Omega*a0 + Kz*az + g1perp dot b = 0.
```

Set `q0=a0+r*az`. If Delta is nonzero,

```
a0 = (k*Kz*q0 + omega*(g1perp dot b))/Delta,
az = -k*(Omega*q0 + g1perp dot b)/Delta.
```

Substitution into the remaining equations gives `M(q0,b)^T=0`, where

```
M = M0 + [[S,T^T],[T,R]],
M0 = [[(k*Kz+g0*Omega)/Delta, Omega*g1perp^T/Delta],
      [Omega*g1perp/Delta, -I+g1perp*g1perp^T/Delta]].
```

This symmetric reduction is valid even if Kz=0; only Delta is divided by. The exact determinant relation is `det(Pi)=-Delta*det(M)/k²`. It recovers the earlier isotropic reduction when g0 and g1 vanish. Omitting these correction terms outside that case is invalid.

## Retain the universal pointwise positive matrix

For r=x+iy, use the real parameters a,b,h,t and the pointwise matrices F_a,Q_a,P_a from full-frequency-bound.md. The proof that `P_a=Q_a+3y F_a` is positive semidefinite holds at every sphere point, independent of any moment assumptions. Let `J=diag(i*sqrt(1-r²),1,1)`, `H=-J M J`, and define the known real symmetric matrix

```
B0 = Herm(-J*M0*J) - diag((x*g0+g1z)/k,0,0).
```

The subtraction is essential: polynomial division of `(1-r²)/(r-z)` now leaves a nonzero integral of `-(r+z)`. Direct expansion gives

```
Herm(H) = B0 + (1/k)*integral F_a d(mu+-mu-).
```

Let the reflected minus measure be mu- evaluated at -v. It stays positive, its flux changes sign and its second tensor is unchanged. Since `F_a(-v)=-F_-a(v)`, define

```
L = B0 - [integral Q_a dmu+ + integral Q_-a d(reflected mu-)]/(3*gamma).
```

Then

```
Herm(H) = L + [integral P_a dmu+ + integral P_-a d(reflected mu-)]/(3*gamma).
```

The second term is positive semidefinite. Positive definite L therefore excludes a null vector of H, M and the original polarization. This argument uses the original equations, not an assumed angular closure.

Every entry of the integrals of Q is available from N,f,Q. For `K_a=bh+at`, `e_a=(K_a²-t²)/h²`, `alpha=b(1+2t²)/h`, and `beta=-t/h`,

```
(integral Q_a)00 = (3/2)*(N-Qzz),
(integral Q_a)0i = (3/2)*(alpha*f_i + beta*Qzi),
(integral Q_a)ij = (3/2)*[(K_a²/h²)*Qzz*delta_ij + e_a*Qij],
```

where i,j run over x,y. Reflection flips f but not Q, and the sign of a is changed for the minus term. This yields a fixed-size calculation with no angular grid or numerical optimizer. Its bit cost still depends on rational input sizes; no timing superiority is claimed.

## Fixed feasibility comparison and adversarial controls

GENERAL_CERTIFICATE_PLAN.md fixed the coordinate-axis measures, k=99/25 and omega=4i/5 before evaluation. Their masses are1 and4/5, fluxes `(1/10,1/10,0)` and `(-1/10,0,1/10)`, and second tensors `diag(3/10,3/10,2/5)` and `diag(3/10,1/5,3/10)`. The certificate's L is exactly positive definite, so the specified frequency is excluded for **every** completion with those moments.

The serious norm baseline bounds growth by the norm of the cross-species integral operator T with kernel `1-v dot v'`. Its norm depends only on the two full moment Gram matrices. Applying T to the constant function gives the exact Rayleigh lower bound

```
||T||² >= N+*N- - 2*f+ dot f- + (f-^T Q+ f-)/N- = 663/800.
```

Here gamma²=16/25 is smaller, so that norm bound cannot exclude the frequency. This is a strict synthetic exclusion improvement at the fixed fixture, not a uniform dominance or physical-use result.

The two-beam control has unit species at opposite z poles, k=6/5 and omega=8i/5. It is an exact growing mode. The correct general certificate leaves it unresolved, whereas incorrectly dropping the signed-current terms produces a positive-definite lower matrix and would falsely exclude it. A separate Delta=0 case must return unresolved without division. Random rational positive atomic backgrounds test the exact determinant relation and full Hermitian decomposition with unequal densities and nonparallel fluxes.

## Status and boundary

`general_moment_certificate.py` provides the exact research helper, and `check_general_moment_certificate.py` retains the algebra and fixed-fixture checks. Its rational frequency parameters satisfy a²+b²=1,h²-t²=1,b,t>0; it does not yet accept every rational complex frequency directly. It validates symmetric moment tensors, trace Q=N and positive moment Gram matrices, and the derivation is conditional on positive measures. It is not a complete realizability solver or a complete frequency classifier.

This passes only the declared mathematical feasibility gate. Independent scrutiny, strongest applicable literature comparisons, robust general input support and a provenance-supported physical decision remain required before treating it as an operational algorithm. The prior failed counterflux experiment and conic timing result are unchanged. No new astrophysical-use pass, practical speedup, priority finding or high-impact breakthrough is claimed.
