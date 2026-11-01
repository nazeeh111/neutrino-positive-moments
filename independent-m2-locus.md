# Independent audit: exact Stieltjes locus with a supplied second moment

This derivation is independent of the production implementation. It concerns
one positive finite measure, not a neutrino dispersion relation. No numerical
optimizer, angular discretization, or unproved moment-completion assertion is
used. The associated independent test is `test_independent_m2_locus.py`.

## Statement and admissible input

Let N>0, M,Q be real, Delta=NQ-M²>=0, Q<=N, z=x+iy with y>0. These conditions
already imply 0<=Q<=N and |M|<=N. A complex number S is the transform

    S = integral dmu(u)/(z-u)

of a positive measure supported on [-1,1] with moments (N,M,Q) if and only if
it belongs to BOTH of the following closed disks. Write D=Q-2Mx+N|z|²>0.

    Disk A:
    |S - ( N(Nx-M)/D - i(2N²y+Delta/y)/(2D) )|²
        <= Delta²/(4y²D²).

    Disk B:
    U=(z²-1)S-Nz-M,
    |U-i(N-Q)/(2y)|² <= (N-Q)²/(4y²).

The map from S to U is invertible: z²-1 is nonzero in the upper half-plane.
The disks may collapse to points. They remain exact in those cases.

For N=0, positivity instead requires M=Q=0 and S=0. Disk A is undefined
and must not be evaluated. Merely keeping Delta>=0 and Q<=N at N=0 is
insufficient: it would admit M=0,Q<0, which is not a moment sequence.
Negative N is infeasible.

## Reweighting is a bijection

Put d rho=dmu/|z-u|². Its degree-four moments t_j satisfy

    t0 = -Im(S)/y,
    t1 = x*t0-Re(S),
    t2 = N+2*x*t1-|z|²*t0,
    t3 = M+2*x*t2-|z|²*t1,
    t4 = Q+2*x*t3-|z|²*t2.

Conversely any positive rho on [-1,1] with these moments gives the desired
mu by dmu=|z-u|² d rho. Its first three moments are (N,M,Q), and its transform
is (bar z)*t0-t1=S. Thus it is enough to characterize this degree-four
truncated moment sequence exactly.

## Necessity and the precise congruences

Define the real symmetric matrices

    H = [[t0,t1,t2], [t1,t2,t3], [t2,t3,t4]],
    L = [[t0-t2,t1-t3], [t1-t3,t2-t4]].

They are positive semidefinite because they integrate |p(u)|² and
(1-u²)|q(u)|² for polynomials p of degree<=2 and q of degree<=1.

Use the convention G_ij=integral p_i*conjugate(p_j) d rho. Changing the
polynomial basis in H to (1,z-u,u(z-u)) is invertible and gives

    G = [[t0, S, z*S-N],
         [bar S, N, M],
         [bar(z*S-N), M, Q]].

The basis-change determinant is 1. In particular det G=det H, exactly.
For Delta>0 the lower 2x2 block is positive definite. Its Schur complement
is nonnegative if and only if

    F_A = D*|S|²-2N(Nx-M)*Re(S)
          +(2N²y+Delta/y)*Im(S)+N³ <= 0.

Completing the square gives precisely disk A, since

    F_A = D*( |S-center_A|²-radius_A² ),
    det H = -F_A.

The Schur-complement argument proves the ENTIRE positive-semidefinite
condition, not merely the determinant inequality. This is why no diagonal
conditions are missing when Delta>0.

For L use the basis (1,z-u), with the positive weight (1-u²)d rho. The
invertible change gives

    K = [[Im(U)/y, -U], [-bar U, R]],  R=N-Q>=0.

Hence det L=R*Im(U)/y-|U|². When R>0 its nonnegative determinant implies
the other diagonal is nonnegative, so this condition is equivalent to
K>=0 and exactly to disk B. When R=0, disk B forces U=0, which also forces
the other diagonal to vanish. No extra diagonal condition is needed there
either.

## Collapsed first disk: Delta=0

Set m=M/N. Then Q=Nm² and |m|<=1. Positivity and the supplied moments force

    integral (u-m)² dmu=Q-M²/N=0,
    mu=N*delta_m,  S=N/(z-m).

Disk A has radius zero and its center is exactly this S. Its H is positive
semidefinite, and disk B holds because this actual positive measure exists.
This separately proves both directions in the singular lower-block case;
one must not use the inverse of that block when Delta=0.

At Q=N, positivity forces support {-1,+1}. The unique measure has weights
(N-M)/2 at -1 and (N+M)/2 at +1. Disk B becomes U=0 and therefore gives
exactly S=(Nz+M)/(z²-1). Disk A then holds automatically. The overlap with
Delta=0 is an endpoint atom and is consistent with both formulas.

## Constructive sufficiency of H>=0 and L>=0

Here is a direct proof of the needed degree-four moment theorem, including
singular cases. Let ell(u^j)=t_j. These conditions say ell(p²)>=0 for
degree(p)<=2 and ell((1-u²)q²)>=0 for degree(q)<=1. The construction below
uses at most three nodes. It supplies an exact witness, rather than a
numerical positive-semidefinite status.

If t0=0, H>=0 and L>=0 force all t_j=0: the zero diagonal in H first kills
t1,t2, then its other zero diagonal kills t3, and t4 is both >=0 and <=0.
The zero measure is the witness. For the present N>0 problem this case
cannot occur, since then N=t2-2xt1+|z|²t0 would vanish.

### Positive definite H: three-node Gauss-Radau construction

Define v_j=t_j-t_(j+1), j=0,...,3, and

    A = [[v0,v1],[v1,v2]],   B = [[v1,v2],[v2,v3]].

These are the degree-three moment and multiplication matrices for the
formal functional ell((1-u)*f). A is positive definite, because for any
nonzero linear p,

    2*ell((1-u)p²)=ell(((1-u)p)²)+ell((1-u²)p²)>0.

Moreover A-B is positive definite by H>0, and A+B=L>=0. Consequently the
two distinct generalized eigenvalues lambda_1<lambda_2 of B w=lambda A w
lie in [-1,1). Distinctness follows because the associated multiplication
operator maps 1 to u and A is positive definite.

Equivalently they are the roots of q(u)=u²+a*u+b, where

    A*[b,a]^T = -[v2,v3]^T.

Let

    alpha_1=(lambda_2*v0-v1)/(lambda_2-lambda_1),
    alpha_2=(v1-lambda_1*v0)/(lambda_2-lambda_1).

Both are strictly positive: they are the spectral weights of the vector 1
for that self-adjoint multiplication operator. They integrate v0,...,v3
exactly. One can also verify the last two moments directly from q's two
orthogonality equations.

Set beta_i=alpha_i/(1-lambda_i) and w=t0-beta_1-beta_2. The signed measure

    rho_* = beta_1*delta_(lambda_1)+beta_2*delta_(lambda_2)+w*delta_1

has moments t0,...,t4: its mass is t0 and the differences of successive
moments are v0,...,v3. Since q(lambda_i)=0, exact moment matching gives

    w*q(1)² = ell(q²)>0.

Thus w>0 and rho_* is positive. Its three nodes are in [-1,1]. This proves
sufficiency constructively, even when L is singular. The strict upper
endpoint exclusion for lambda_i ensures all divisions are legitimate.

### Singular H with positive variance: two-node construction

Suppose d=t0*t2-t1²>0 and H is singular. Then H has rank two. Define

    A0=[[t0,t1],[t1,t2]], B0=[[t1,t2],[t2,t3]], T=A0^(-1)*B0.

The matrix T represents multiplication on the basis (1,u), is self-adjoint
for A0, and maps 1 to u. The rank-two condition supplies the quadratic
recurrence and makes its square reproduce t4 exactly. More explicitly,
the Schur complement gives

    t4=[t2,t3]*A0^(-1)*[t2,t3]^T,
    A0-B0*A0^(-1)*B0=L>=0.

It follows that T is a contraction, so its distinct real eigenvalues lie
in [-1,1]. Their positive spectral weights for 1 give a two-node rho with
moments through degree four. Nodes and weights can be obtained with the
same quadratic formulas as above, replacing v by t. This includes endpoint
nodes. No inversion of the singular 3x3 matrix is required.

### Zero variance: one node

If t0>0 and d=0, write m=t1/t0. Positive semidefiniteness forces
t2=m²t0 and t3=m³t0. L's top diagonal gives |m|<=1. H gives t4>=m⁴t0,
while L's determinant gives the reverse inequality for |m|<1; when
|m|=1, L's bottom diagonal gives that reverse inequality directly.
Thus t4=m⁴t0 and rho=t0*delta_m is the witness.

These cases exhaust H>=0. Finally multiply each rho weight by |z-u|² to
obtain mu and verify its mass, flux, second moment and transform. This
completes the necessity AND sufficiency proof for the two disks.

## Exact witness arithmetic and a necessary three-atom example

For rational x,y,N,M,Q,Re(S),Im(S), all t_j and quadratic coefficients are
rational. The construction therefore requires only rational arithmetic
and one square root. Its signs can be checked exactly by comparing rational
squares after accounting for operand signs. A verifier can check nonnegative
weights, support, three supplied moments, and the real/imaginary parts of
the transform directly; it need not trust an optimizer or this construction.

A universal two-atom witness is false once Q is supplied. For example

    N=1, M=0, Q=2/3, z=i, S=-2i/3,
    mu=(delta_-1+delta_0+delta_1)/3.

The uniquely implied reweighted moments are (2/3,0,1/3,0,1/3), and det H=1/27.
Any measure on at most two atoms would have rank(H)<=2, a contradiction.
Three atoms are both sufficient in general and sometimes necessary.

## Check scope

The independent standard-library test checks exact rational congruence
determinants and full principal-minor positive semidefiniteness against the
two disk inequalities for arbitrary proposed S, including infeasible points,
Delta=0, Q=N and endpoint intersections. It constructs exact three-node
Radau witnesses in quadratic fields from independent rational atomic inputs
and checks every recovered moment and transform. Tests establish the stated
finite checks; the argument above establishes all real input cases.
