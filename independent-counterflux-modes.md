# Independent counter-flux mode reference

This work reproduces a published controlled input and its axial linear modes. It is not a neutron-star-merger snapshot replay, nonlinear simulation or all-root census. No production moment classifier, Weyl disks or exact feasibility backend is imported by counterflux_modes.py.

## Source anchors and normalization

The prescribed continuum input is the azimuth-integrated version of Richers et al. (2021), Sec. III A Eq. (8): positive number measures dmu±=(1±u)du/2, with masses one and fluxes ±1/3. Their second number moments are exactly 1/3 by direct polynomial integration. The dimensional benchmark uses N=4.89e32 cm^-3; only its dimensionless form is used here. [Primary input](https://arxiv.org/html/2101.02745#S3.SS1).

Froustey et al. (2024, corrected arXiv v3) Eq. (23) fixes the sign of the frame shift: reduced k = physical K - 2/3, and omega is unchanged because signed mass vanishes. Equations (12), (19), (24)-(27) specify the reproduced polynomial Minerbo/FT quantum closure. The corrected matrix convention is checked through those component equations rather than relying on an older Eq. (34). [Primary source](https://arxiv.org/html/2311.11968v3#S3.SS1).

Units set mu0=sqrt(2) G_F N=1 and c=hbar=1. The fixed physical analysis grid is K=j/8 for j=-8,...,8. This is a chosen analysis band, not a physical ultraviolet cutoff. In particular, replacing it with reduced |k|<=1 would be a different experiment.

## Three distinct descriptions

1. The full-angular linear input has signed distribution G(u)=u and second moments Q±=1/3.
2. The exact maximum-entropy angular completion has dmu±=beta exp(±beta u)du/[2 sinh(beta)], where coth(beta)-1/beta=1/3. Numerical beta is 1.0745628999535313, and Q±=0.37959270071998935. Its full-angular dispersion is evaluated independently.
3. The published quantum-moment closure uses the approximate polynomial Minerbo value chi=17/45, not that exact entropy second moment. Its evolving off-diagonal closure is also not the full-angular dispersion of the reconstructed background. The source itself identifies the polynomial as an approximation.

These descriptions share M0/M1, but they must not be conflated. Imposing supplied Q=1/3 on a proposed closure mode is an intentional consistency audit against the linear-input information class, not proof that the published method used that Q.

## Reproduced closure formulas

From the component equations, set f=1/3 and chi=17/45. For longitudinal amplitudes A,B,Abar,Bbar, the omega matrix is

    [-1,       k+f,        1,       -f]
    [chi*k-f,  chi,        f,       -chi]
    [-1,       -f,         1,       k+f]
    [f,        chi,        chi*k-f, -chi].

Subtracting the species equations closes the two variables A-Abar and B-Bbar. It yields

    omega² = (k+2f)(chi*k-2f),
    gamma² = (k+2f)(2f-chi*k)
           = K[2f(1+chi)-chi*K].

The independent four-by-four eigenvalue replay agrees on every growing grid point to 6.11e-16. Its characteristic-polynomial coefficients match the factored analytic polynomial to 2.78e-15. The formula reproduces the source's reported global peak: reduced k=28/51 and gamma=62/(9 sqrt(85)). Its physical K=62/51 lies outside this chosen band and is not silently added to the grid.

For either axial transverse direction, Pi_ex[x,z]=DeltaF_x=0. Equations (26)-(27) leave the exact block qT[[1,-1],[1,-1]], qT=(1-chi)/2=14/45. Its square is zero and both eigenvalues are zero for every K. This statement is specific to the reproduced quantum closure; the full-angular reference does have transverse growing modes.

## Continuous angular reference and certified controls

The fixed imaginary-frequency search uses 160 logarithmic gamma points from 1e-7 to 3. At each K, each sign-changing bracket is refined with SciPy brentq. Pairing u and -u analytically leaves smooth real integrals. The value of G(u)/u at zero is subtracted before adaptive integration, avoiding a narrow low-gamma numerical layer. Every retained numerical root is independently checked using direct 256- and 512-point Gauss-Legendre evaluations of the ORIGINAL I0,I1,I2 factors. No recurrence-derived determinant replaces that independent check. The maximum 512-point factor residual over the retained roots is 2.67e-13.

For the linear input, the same integrals have elementary formulas. With A=integral_0^1 u²/(gamma²+k²u²)du and B=integral_0^1 u⁴/(gamma²+k²u²)du,

    A=[1-gamma/|k| atan(|k|/gamma)]/k²,
    B=[1/3-gamma² A]/k²,
    I0=-2k A, I1=-2i gamma A, I2=-2k B.

At k=0 the continuous limits A=1/(3 gamma²), B=1/(5 gamma²) are used. Arb at 160 bits evaluates these ORIGINAL dispersion factors at rational endpoints, with strict opposite signs and interval width 1e-8. Continuity then certifies existence of all 13 retained linear-input growing roots: eight longitudinal and five transverse. This is an IVT existence proof, not uniqueness, a count of every complex root, or a stability proof when the search finds nothing. The exact-entropy roots remain numerical references, not interval-certified roots.

At K=1/8,1/2,7/8,1 the linear longitudinal rates are approximately 0.3441151,0.6188162,0.6892387,0.6866497, while the quantum-closure rates are 0.3300182,0.6039990,0.7172640,0.7353508. These differences alone do not prove incompatibility: another positive completion with the same supplied moments could realize the proposed frequency. The parent moment audit must decide that question independently.

## Development failures and verification limits

The first matrix replay stopped before producing a result because my initial assembly used incorrect signs for the antineutrino flux entries. The original source and failed logs are preserved. Correcting those entries from Eqs. (25),(27) fixed the growing-branch disagreement. At physical K=0, floating eigvals split a defective exact-zero eigenvalue into a tiny apparent complex pair, about 2.7e-8. That numerical splitting is retained as a diagnostic, never counted as a physical growing mode; the exact analytic formula and independent characteristic-polynomial check establish the zero value. No input grid, search range or scientific success gate was changed to hide a mode discrepancy.

Current outputs are results/counterflux-modes.json and its log. The JSON includes exact rational closure gamma² values, numerical root residuals, 13 validated linear brackets, the exact transverse block, source URLs, runtime information and the current source SHA256. The earlier completed numerical-only output is also retained. The initial source is historical failed development code, not a production alternative.

Reproduce with the existing free scientific dependencies (NumPy, SciPy, python-flint):

    PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python counterflux_modes.py

No formal claim about the parent M2 feasibility method, priority, astrophysical utility or publication readiness follows from these reference calculations alone.
