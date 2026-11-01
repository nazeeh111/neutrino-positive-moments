# Gate review: published advection discrepancies as a diagnostic use case

2026-11-01. Independent review of one lead: Grohs et al., *Advection Algorithms for Quantum Neutrino Moment Transport*, arXiv2501.07540v2. This does not rerun or change the earlier fixed counter-flux utility experiment. No simulation, package installation, classifier scan, external message or publication was performed.

## Decision

**The published summary tables do not yet define a valid exact exclusion experiment for the current axial diagnostic.** They identify a real numerical discrepancy, but extracting a classifier input requires information the tables do not supply. Do not turn a rounded growth/peak pair into an exact complex mode, assume a missing phase frequency or wave-vector sign, or average a nonaxial background into the axial model. This is a bounded input-adequacy decision, not proof that the underlying numerical error could never be excluded with more data.

## What the primary source actually supplies

[The paper](https://arxiv.org/html/2501.07540v2), Tables1/3/6/7, supplies number densities, vector fluxes, fitted growth and Fourier peaks. Axial Fiducial has equal densities4.89e32 cm^-3 and flux factors±1/3; old FLASH(ri) gives7.1e10 s^-1 and3.9(4) cm^-1. AppendixA describes peak-width uncertainty. SectionIV.C reports multiple phases and harmonics for axial TwoThirds. NSM vectors are generally nonparallel;90Degree is explicitly nonaxial. AppendixB includes vacuum mixing in FLASH but not Emu. SectionIV.A.3 supplies modeled closure pressures, not independent measured Q. These are source reports, not local replications.

## Why table-center classification would not prove usefulness

The frequency diagnostic takes a specified complex frequency and a signed reduced wave number. A fit to a domain-averaged coherence magnitude supplies an estimated growth rate, not its real oscillation frequency. A peak in a three-dimensional spatial power spectrum does not by itself supply the signed wave vector of that same temporal mode. Superposed modes, time-window effects and harmonics can matter, as the source itself demonstrates. No printed uncertainty for the fitted growth rate establishes an exact interval suitable for a proof.

The known real density/current shifts must be applied consistently. Neither taking the displayed positive peak as signed reduced k nor setting Re(omega)=0 by convenience is justified by the table alone. A symmetry argument for a particular continuum operator is not automatically a measurement of the phase of a discrete nonlinear-limiter calculation.

The broadest zero-flux sharp envelope is also unavailable for Fiducial: its individual fluxes are nonzero. The general pointwise M1/M2 classifier does support those fluxes, but still needs the frequency/shift inputs. A valid all-real-frequency exclusion could in principle avoid the missing phase frequency; that requires an applicable all-real bound for the actual moment class. The presently published equal-zero-flux theorem is not such a bound for these cases.

The current M2 question additionally needs independently justified second number moments. Reusing the closure's own pressure as the independent truth tests a different information class. A prescribed continuum completion can provide a theoretical Q, but it must be derived from that actual completion; the polynomial Minerbo coefficient and exact maximum-entropy second moment must not be silently identified.

Finally, the nonzero vacuum term means a literal exact fast-limit certificate is not automatically a certificate for the reported full simulation. A controlled comparison may be possible because the scales are very different, but a margin/error argument must be supplied rather than assumed. The tabulated discrepancy does not remove these quantifier problems.

## Structural obstruction for fixed scalar longitudinal closures

A useful no-go statement can be proved without fitting any table point. Consider the one-dimensional longitudinal linearized fast problem, with each species s having N_s>0, first moment F_s, and a constant scalar closure

    Q_s=chi_s*N_s,
    deltaP_s=chi_s*deltaN_s,
    0<chi_s<=1,  F_s²<=chi_s*N_s².

The coefficients may differ between species. Put v_s=sqrt(chi_s) and define positive two-node background weights

    w_s,+=(N_s+F_s/v_s)/2,
    w_s,-=(N_s-F_s/v_s)/2,
    mu_s=w_s,+ delta_(v_s)+w_s,- delta_(-v_s).

They have exactly the supplied N_s,F_s,Q_s and support[-1,1]. Let p_s,± be their angular perturbation amplitudes. The invertible map

    (deltaN_s,deltaF_s)^T = T_s (p_s,+,p_s,-)^T,
    T_s=[[1,1],[v_s,-v_s]]

makes deltaP_s=chi_s*deltaN_s an identity. Streaming transforms from diag(k*v_s,-k*v_s) to [[0,k],[k*chi_s,0]]. The axial interaction kernel1-u*u' is affine in either angle. Thus its zeroth and first moments involve only (N_s,F_s,Q_s), and taking moments on these two nodes closes exactly.

More explicitly, use the reduced operator and species signs sigma_s=+1,-1. The moment-space source block from species t into species s is

    sigma_t * [[-N_s,F_s],[-F_s,Q_s]].

The angular density-amplitude matrix has entries

    A_ij = k*u_i*delta_ij - w_i*(1-u_i*u_j)*sigma_j.

With T=diag(T_plus,T_minus), the closed moment matrix M obeys

    M*T = T*A.

Hence every growing eigenvalue of this constant-scalar longitudinal closure is an exact growing mode of an admissible positive atomic completion. The existential M2 diagnostic must accept it. This is the previously established admissible-quadrature obstruction applied at the operator level. It also shows why merely comparing a closure to its own supplied Q cannot automatically create an exclusion test.

I independently checked240 exact rational MT=TA identities with the standard library, seed250107540. All background weights were nonnegative and all matrix identities passed exactly. This inline exploration is supplementary; the retained reproduction artifact is the parent's `usefulness-assessment/closure_equivalence.py`, reviewed below. Finite checks corroborate the displayed general similarity proof; they are not a new scientific utility benchmark.

**Boundary of this statement:** it does not assert that the full FLASH tensor closure or nonlinear modulus/phase advection scheme is this constant-scalar operator. Different closures for background and perturbations, varying orientation tensors, additional transverse modes, nonlinear reconstruction, collisions or vacuum mixing require their own derivation. In particular, it does not imply that every FLASH(ri) fitted rate is compatible. It identifies a class for which exclusion is structurally impossible.

## What would make this lead executable

Obtain a mode-resolved time series or a reproducible linearized discrete operator for one axial case, with the signed wave vector, complex frequency, actual background number moments and a defensible error interval. Keep the independent supplied-Q class distinct from the closure's model. If using the full simulation, justify the fast-limit comparison or extend the operator. First compare with direct original-dispersion residuals, the paper's own analytic/LSA reference, and applicable generic bounds. A proof of exclusion with margin must add a decision those alternatives do not already provide.

The beam tests in AppendixA already have analytic reference solutions and show agreement in their reported growth rates. Repackaging those tests through the moment diagnostic would not establish a new capability. The published advection correction remains real and important, but the summary tables alone do not give this candidate a successful practical-use gate.

This review supplies a narrower next requirement, not permission for a large FLASH installation or an expanded parameter scan. The earlier failed counter-flux gate remains failed. No artificial mode or invented measurement precision was introduced.


## Review of the retained executable and degenerate weights

I read `usefulness-assessment/closure_equivalence.py` against the independent derivation above. Its density-amplitude angular operator, moment source blocks and block-diagonal T are correct. Its domain has strictly positive v_s, hence invertible T_s. The theorem and executable do not claim that this similarity works at chi_s=0.

Zero background weights are allowed. If w_i=0, the corresponding angular row is free streaming: (omega-k*u_i)p_i=0. For Im(omega)>0, this forces p_i=0. Thus any extra zero-weight coordinate contributes only a real streaming eigenvalue and cannot create a spurious growing completion. This addresses the positivity boundary without silently discarding it.

The retained result reports200 exact identities, including61 cases with a zero weight; I checked that its saved source hash matches the reviewed executable. The parent also supplies a concrete independently checked M2 root: N_plus=N_minus=1, F_plus=1/3, F_minus=-1/3, Q_plus=Q_minus=1/4, reduced k=1 and omega=i*5/6. The positive quadrature weights are5/6 and1/6 on u=±1/2, with reflection for the other species. Direct summation gives I0=-12/17, I1=-10i/17 and I2=-3/17, so the original longitudinal factor is exactly zero. The saved certificate is accepted by the unchanged production M2 classifier and separate Hausdorff verifier. This is a consistency witness, not a naturally observed erroneous simulation mode.

Reproduce from the research root:

    python3 -B -O usefulness-assessment/closure_equivalence.py

Reviewed executable SHA256: fd49749cc1cf499b6fe07a4caf3674f0275a2b4dbeac39fd4660273f51697e9c.

Do not generalize this result to every quantum maximum-entropy closure: a relation deltaP=kappa*deltaN+eta*deltaF need not obey the corresponding background identity Q=kappa*N+eta*F. Likewise, this scalar proof is not a derivation of the tensor closure or nonlinear modulus/phase reconstruction used by FLASH. Those distinctions preserve the possibility of a different admissible test while explaining why the present source-summary lead has not passed its input gate.
