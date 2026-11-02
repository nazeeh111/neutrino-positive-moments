# General-moment certificate: bounded feasibility gate

The current full-frequency isotropic theorem is proved and primary checked. This new lead asks whether its pointwise positive matrix can give a sufficient frequency-exclusion certificate from arbitrary three-dimensional number moments `(N,F,Q)`, including unequal species masses and vector fluxes. It must not reuse the zero-signed-current reduction outside its assumptions.

Before a larger implementation, require:

1. A written reduction from the original4x4 polarization, with all known signed density/current terms retained and a correct exceptional-denominator branch.
2. Exact original-versus-reduced determinant and Hermitian-decomposition checks on positive rational atomic backgrounds with unequal masses and nonparallel fluxes. A singular reduction must return unresolved, never excluded.
3. No exclusion of the exact two-beam root with plus at-z, minus at+z, unit masses, k=6/5 and omega=8i/5.
4. At least one certified exclusion not already supplied by the strongest species-block operator-norm bound. Freeze this first fixture: plus coordinate-axis weights `(1/5,1/10,1/5,1/10,1/5,1/5)` and minus weights `(1/10,1/5,1/10,1/10,1/5,1/10)`, in order `(+x,-x,+y,-y,+z,-z)`. Test reduced k=99/25 and omega=4i/5. The full moments are computed from these specified positive measures. No alternate fixture will replace a failure of this particular test.

For the norm baseline, the cross-species operator with kernel `1-v dot v'` has a norm determined by the two full moment Gram matrices. An exact Rayleigh lower bound above the tested gamma is enough to prove that this baseline cannot exclude the frequency; a looser mass-only bound is not an acceptable substitute.

The first implementation may use rational parameters for the complex square root in the positive certificate. It must state this interface restriction. Passing this mathematical feasibility gate is not an astrophysical-use pass, a speed result, priority evidence or independent validation. General frequency input, error enclosures, mode-resolved physical data and an actual decision test remain separate requirements. The earlier failed counterflux gate is not reopened. No installation, paid API, external communication or usage-limit bypass is authorized by this plan.
