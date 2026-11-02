# Quantum maximum-entropy closure: a conditional exclusion obstruction

2026-11-02. Primary derivation and exact algebra checks, not independently reviewed. This is a decision about a proposed benchmark, not a new successful utility experiment. The previous counter-flux gate remains failed.

## Decision and source boundary

A broad scan of the quantum maximum-entropy closure has not earned implementation. In the ordered-flavor class specified below, every longitudinal growing closure mode already has a positive angular completion matching the electron species' first three moments. A sound existential moment diagnostic therefore cannot exclude it. Cases outside that class remain open; this note does not classify the published merger snapshot.

[Froustey, Kneller and McLaughlin, arXiv2409.05807v2](https://arxiv.org/html/2409.05807v2), Eqs.13–17, define a positive logarithmic-mean weight h and an off-diagonal closure P=kappa*N+eta*F from its first four moments. Eq.24 gives the two characteristic speeds; SectionIII checks their reality and causality numerically. SectionIV uses axisymmetric longitudinal perturbations and removes transverse background fluxes after rotation. That projected calculation is not full three-dimensional snapshot data. The following quadrature argument and conditional completion are our deductions from these equations, not claims made by the source.

## Two-node representation

Absorb the angular normalization into a positive measure dh on [-1,1]. Write I_j=integral u^j dh and assume D=I_0*I_2-I_1^2>0. Define

    kappa=(I_2^2-I_1*I_3)/D,
    eta=(I_0*I_3-I_1*I_2)/D,
    p_2(u)=u^2-eta*u-kappa.

Direct substitution gives integral p_2 dh=integral u*p_2 dh=0. Thus p_2 is the degree-two orthogonal polynomial of dh. Its two distinct real roots lie in [-1,1], and its Gaussian quadrature weights w_1,w_2 are positive. Quadrature is exact through degree three. These are classical Gaussian quadrature facts, not an invented algorithm. For a strictly positive density throughout the interval, both roots are interior. Rank-one measures D=0 are outside this proof and require their limiting or one-node treatment.

A separate proof of real, causal speeds avoids locating roots. Let H=[[I_0,I_1],[I_1,I_2]] and K=[[I_1,I_2],[I_2,I_3]]. H is positive definite, K is symmetric, and -H<=K<=H because integral (1+/-u)*(c_0+c_1*u)^2 dh>=0. The generalized eigenvalues of K v=lambda H v are real and belong to [-1,1]. They are exactly the roots of p_2. This is a standard positive-weight moment argument, not a claim of new general hyperbolicity theory.

## Additional condition needed for positive source weights

For each neutrino or antineutrino sector s, let the two diagonal Maxwell–Boltzmann distributions be strictly positive exponentials. Their logarithmic ratio is affine:

    log(rho_e,s/rho_x,s)=A_s+B_s*u,
    g_s=rho_e,s-rho_x,s=(A_s+B_s*u)*h_s.

The logarithmic mean has its continuous value at equal diagonal densities. Gaussian quadrature for h_s defines source weights

    d_s,j=w_s,j*(A_s+B_s*u_s,j).

For n=0,1,2, sum_j d_s,j*u_s,j^n=integral u^n*g_s du, since the integrand against h_s has degree at most three. Consequently Q_g,s=kappa_s*N_g,s+eta_s*F_g,s. This background identity concerns the flavor difference g_s, not either electron distribution separately.

If A_s>=abs(B_s), then g_s is nonnegative throughout [-1,1], so d_s,j>=0. This simple sufficient condition is the ordered-flavor class used here. Nonnegative values at the two nodes are enough for the quadrature construction, but no claim about a whole snapshot follows without checking them. If the weights have either sign, the positive-completion conclusion below is unavailable. Positivity of h alone does not imply positivity of g or of d.

## Exact equivalence of the linear operators

Use the reduced real wave number k and reduced frequency omega, with the same real density/current shifts as PHYSICAL_SCOPE.md. Let sigma_plus=1,sigma_minus=-1. The two-node angular density-amplitude operator is

    L_(s,j),(t,l)=k*u_s,j*delta_(s,j),(t,l)
                   -d_s,j*(1-u_s,j*u_t,l)*sigma_t.

Define T_s=[[1,1],[u_s,1,u_s,2]] and T=diag(T_plus,T_minus). Distinct nodes make T invertible. The moment operator M has streaming blocks k*[[0,1],[kappa_s,eta_s]] and interaction blocks

    M_source[s,t]=sigma_t*[[-N_g,s,F_g,s],[-F_g,s,Q_g,s]].

Then M*T=T*L. Streaming follows from p_2(u_s,j)=0; the interaction follows from exact matching of the first three g moments. This proves equality of the entire four-variable longitudinal spectra. For a zero source weight, its angular row is free streaming and contributes no growing eigenmode. No finite numerical scan is needed for the general identity.

Now assume the heavy-flavor angular backgrounds are the SAME nonnegative measure xi in both sectors, as required by the pointwise heavy-flavor cancellation model. The positive electron completions

    mu_s=xi+sum_j d_s,j*delta_(u_s,j)

have the supplied electron N,F,Q, including the actual diagonal exponential second moments. Their signed difference is exactly the quadrature source because xi cancels. Hence each growing closure root is a root of an admissible completion's longitudinal dispersion relation. For the three-dimensional positive-measure model, lift each node to a uniform azimuthal ring; the corresponding axisymmetric full second tensor also matches. This only represents longitudinal perturbations.

Therefore a sound M1/M2 existential compatibility test cannot reject these roots. They can still differ from the growth of the particular continuous exponential completion. Compatibility with some completion is not accuracy for that chosen completion.

## Checks and scientific implications

`qme_quadrature_gate.py` checks the operator identity using exact fractions for 96 two-sector rational moment fixtures, including nonzero eta and unequal sector parameters. The positive h fixtures have three support points and known two-node Gaussian rules. They exercise the general algebra; they are not physical merger data or numerical reconstructions of exponential h. Each checks four h moments, three g moments, positive source weights and M*T=T*L. An exact growing root at reduced k=16/25, omega=102i/125 is checked directly in the original longitudinal polarization and accepted by the existing M2 classifier plus its separate verifier. A sign-changing g control demonstrates why h positivity alone is insufficient. A common heavy component checks electron moment preservation and exact cancellation.

This closes one tempting but invalid route to claiming utility. It does not establish priority of our sharp cap, improve the failed physical gate, or prove that every quantum closure is compatible. Unequal heavy backgrounds, sign-changing flavor differences, transverse modes, full tensor closures, vacuum terms and collisions remain outside the conclusion. The full nonaxial certificate is sufficient rather than complete and cannot overcome this information-theoretic obstruction.

The next physical test must identify a mode and independently justified input constraints that are not already guaranteed a positive completion by the closure construction. Printed growth summaries and closure-derived pressures alone do not establish that. Do not build an uncertainty-box API solely to run this covered class, and do not relabel this classical quadrature deduction as the high-impact algorithm sought by the goal.

Reproduce from the repository root:

    python3 -B -O usefulness-assessment/qme_quadrature_gate.py
