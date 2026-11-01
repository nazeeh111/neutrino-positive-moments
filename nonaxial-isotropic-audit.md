# Nonaxial isotropic-moment audit

2026-11-01. Independent specialist audit, separate from the parent's rotated-tetrahedron search. Scope: positive species measures on the unit sphere, each mass1, mean velocity zero and second tensorI/3; wave vector k alongz. This removes axial symmetry while retaining the full low-moment tensor. No installed environment, published source, paper or existing result was changed by this audit.

## Result

The proposed discrete matrix has the correct sign and reproduces the full polarization determinant for upper-half-plane frequencies. A concrete eight-beam construction has an exactly verified growing mode omega=i/3 at k=1/sqrt3, while azimuthally averaging either species makes their signed angular distribution identically zero. This proves that averaging away angular structure can erase an instability in this information class. It does **not** violate the published axial worst-case growth envelope, establish a full nonaxial extension, or establish novelty of azimuthal instabilities.

A safe general all-wave-number bound is gamma<=1 (or N for common species massN), obtained directly from the species block operator. It is weaker than the axial unrestricted cap2/3. Thus neither that general bound nor a failed counterexample search justifies carrying the axial sharp envelope into the nonaxial setting.

## Original matrix and polarization relation

For directions v_i and signed masses g_i, define K_ij=1−v_i dot v_j and

    A_ij=k*v_i,z*delta_ij−K_ij*g_j.

Put V_i=(1,v_i), eta=diag(1,−1,−1,−1), and D=diag(omega−k*v_i,z). Then K=V*eta*V^T and

    det(omega I−A)=det(D)*det(I+eta*V^T*diag(g)*D^(-1)*V).

Multiplication by eta changes only a nonzero determinant factor. Consequently, for Im omega>0, the matrix roots are exactly roots of

    Pi(omega,k)=eta+sum_i g_i V_i V_i^T/(omega−k*v_i,z).

All entries of D are nonzero in that domain. Real streaming-pole eigenvalues require separate interpretation; the identity is not used there. Equal per-species mass and first moments give zero signed density/current shifts. This verifies the matrix sign against the original dispersion convention, rather than only comparing two numerical implementations of the same matrix.

## Explicit nonaxial example with identical azimuthal averages

Let a=1/sqrt3 and b=sqrt(2/3). Give each of the following four directions mass1/4 in mu_plus:

    (+b,0,+a), (−b,0,+a), (0,+b,−a), (0,−b,−a).

Let mu_minus be its reflection z -> −z. Each direction has unit length. The opposite x/y pairs cancel those flux components; the two z levels cancel z flux. Each diagonal second moment is1/3, and each off-diagonal second moment vanishes. Thus both species have exactly the stated mass, zero vector flux and full isotropic second tensor.

Azimuthal averaging puts mass1/2 on each of the rings z=+a andz=−a for BOTH species. Their averaged difference is therefore identically zero. The full distribution is nonaxial: different azimuthal pairs occupy the two z levels.

In the unaveraged polarization tensor all signed entries except xx and yy vanish. Define

    h=2ka/[3(omega²−k²a²)].

Then Pi=diag(1,−1+h,−1−h,−1). The yy factor has an upper-half-plane root whenever

    gamma²=2ka/3−k²a²>0.

Its maximum is gamma=1/3 at k=a. The parent-provided axial total envelope at that wave number is sqrt(4a/9−a²/9), about0.468, so this example does not disprove that upper bound. It does disprove any assumption that azimuthal projection preserves the actual instability spectrum.

At k=a, every entry of the original8x8 matrix is rational: k*v_z=±1/3 and every dot product is rational. Hence an exact rational-complex eigenvector verifies omega=i/3 without rounding the irrational direction coordinates. The following standard-library reproduction was executed successfully. It imports only the existing independent Fraction-complex helper, not a production stability implementation.

```python
from fractions import Fraction as F
from itertools import product
from independent_verifier import C

labels = list(product(('x','y'), (-1,1), (-1,1)))
g = [F(z if axis == 'x' else -z,4) for axis,s,z in labels]
A = []
for i,(axis,s,z) in enumerate(labels):
    row = []
    for j,(axis2,s2,z2) in enumerate(labels):
        dot = F(2,3)*s*s2*(axis == axis2)+F(1,3)*z*z2
        row.append((F(z,3) if i == j else 0)-(1-dot)*g[j])
    A.append(row)
omega = C(0,F(1,3))
Q = [C(s if axis == 'y' else 0)/(omega-F(z,3))
     for axis,s,z in labels]
for i in range(8):
    residual = sum((A[i][j]*Q[j] for j in range(8)),C())-omega*Q[i]
    if residual != C():
        raise RuntimeError('Nonzero original-matrix eigenvector residual')
```

Q is nonzero. Its spatial y factor b has been divided out as a common scalar. The analytic diagonal Pi calculation independently explains the root and its finite-k interval. The original matrix has characteristic polynomial

    (omega²−k²a²)² * [(omega²−k²a²)²−(2ka/3)²].

Only the upper-half-plane root is needed for the audit; no claim about the physical relevance of the real streaming-pole factors is made.

## Parent transverse-third-moment certificate checked

The parent's separate sharp third-moment claim is correct. For an isotropic zero-mean probability measure, consider f(v)=v_x²*v_z. With a=1/sqrt3, the pointwise dual slack is

    P=a/12+v_z/6+(a/2)*v_x²+(a/4)*v_z²−v_x²*v_z.

For fixed z, it is affine in x² on[0,1−z²]. At its two endpoints,

    P|x²=0=(a/4)(z+a)²,
    P|x²=1−z²=(z−a)²(z+7a/4).

Both are nonnegative for −1<=z<=1 because7a/4>1. Taking expectations uses only the supplied moments and gives E[f]<=a/3=1/(3sqrt3). The four-atom mu_plus above attains equality. This is a checked moment extremum; its coefficient is below the axial longitudinal coefficient and does not settle the complete nonaxial finite-k problem.

## A safe general bound

Represent the linearized problem on L²(mu_plus) direct-sum L²(mu_minus), with ordinary positive inner products. Its diagonal blocks are self-adjoint. The off-diagonal blocks are T and−T*, where

    (T f)(v)=integral (1−v dot v') f(v') dmu_minus(v').

Under the common isotropic moments, the feature maps with columns1,sqrt3*v_x,sqrt3*v_y,sqrt3*v_z are isometries fromR4 into each species space. Therefore T has nonzero singular values1,1/3,1/3,1/3. The imaginary Hermitian part of the full operator has norm||T||=1, giving |Im omega|<=1 for every eigenmode, independent of k. Common massN rescales this bound toN. This standard operator-norm argument supplies only a safe envelope; sharpness and novelty are not asserted.

## Limited numerical exploration

An additional independent finite-grid moment LP used8,006 directions (4,000 seeded random sphere points, their antipodes and six coordinate axes), exact intended constraints in floating arithmetic, and objective E[(e dot v)²*v_z] for seven prescribed polarization angles between e andz. Its maximum ranged from approximately.222222 for parallel polarization to.192337 for perpendicular polarization, with reported moment residuals below2e−13. This was an exploratory lower approximation to a moment extremum, not an upper bound, and supplied no violation of the axial small-k coefficient. No canonical scientific claim relies on this LP and no large campaign is recommended from it.

## Prior art and next-action boundary

Azimuthal structure and symmetry-breaking flavor instabilities are established topics. A primary antecedent is Raffelt, Sarikas and deSousaSeixas, *Axial symmetry breaking in self-induced flavor conversion of supernova neutrino fluxes*, PRL111,091101(2013), https://arxiv.org/abs/1305.7140 . That paper's regime and instability mechanism are not identified with this exact moment-constrained fast-limit example; it establishes why the broad claim “azimuthal structure can matter” is not new.

The useful evidence here is an exact scope counterexample to replacing a nonaxial distribution by its azimuthal average, plus a verified discrete-matrix route for further falsifiers. The published axial theorem remains valid in its stated class. Its extension to all nonaxial isotropic completions remains unproved and has not been disproved by this construction. Any new global theorem must control the full matrix-valued transform, not silently use the scalar axial factors. Do not promote the bounded scan, this projection counterexample or the safegamma<=N estimate into the original high-impact goal's completion.


## Retained parent implementation review

Read nonaxial-exploration/exact_projection_counterexample.py, SHA256 e4867630a02b135913b35735ce13a26c8addf382e3e8a5eeffd1abfaa484328c. Its encoded coordinate scales, unit-sphere identity, per-species means/tensors, equal azimuthal projection test, original8x8 eigen-equations and full polarization diagonal all match the independent derivation above. The off-diagonal checks correctly omit nonzero common coordinate scale factors because they test exact zero before scaling. No consequential implementation error was found. The independently executed Fraction-complex witness above supplies a separate construction/order of the same exact eigen-equation. The retained parent script is the preferred executable reproduction artifact. No larger search or paper change is justified merely to repeat this validated scope example.
