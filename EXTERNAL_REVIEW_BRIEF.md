# Scoped scientific review request, prepared locally only

The research repository is public. No expert message or paper submission has been sent. This is a review brief for the verified research package. The full moonshot objective is incomplete; high impact, practical superiority and priority are not established.

## Current priority: the full nonaxial theorem

The November2 extension is [the sharp full-frequency nonaxial isotropic cap](nonaxial-exploration/full-frequency-bound.md), now Section7 of the13-page research draft. Independent positive sphere measures with common massN, zero vector flux and full second tensorNI/3 have maximum growth at most2N/3 over every real wave vector and every real frequency part. An exact oscillating mode attains the cap. The proof uses a full polarization reduction, separate positive matrix certificates for the species, and an outer-product identity. The exact checker is `nonaxial-exploration/exact_nonaxial_full_bound.py`.

This extension has primary checks only. The earlier independent model-agent reviews do not cover it, and the requested new reviewer returned a usage-limit error. Please distinguish mathematical correctness, priority and usefulness in any assessment.

The most discriminating checks are whether the upper-half-plane parameterization is complete; whether the original polarization and its reduced matrix have equivalent nullspaces; whether the reflection identity and the two species' different transverse constants are handled correctly; and whether the averaged positive-definite term excludes every claimed frequency. A simpler known operator inequality giving the same sharp2N/3 cap would materially change the novelty assessment.

The [information-limit corollary](usefulness-assessment/nonaxial-information-limit.md) shows that identical supplied moments admit actual maximum growth spanning all of `[0,2N/3]`. Thus the cap does not predict the actual rate, and the optimal worst-case absolute error of a deterministic point estimate from these moments alone isN/3. Is there a concrete scientific workflow that benefits from this exact exclusion guarantee despite that limitation? No operational-use pass is claimed.

A later [general-moment feasibility extension](nonaxial-exploration/general-moment-certificate.md) derives the nonzero signed-density/current reduction and a sufficient exclusion matrix for general species moments. Please check its Ward relation, determinant factor and Delta=0 branch before treating it as an algorithm. Its fixed synthetic improvement over the species-block norm bound is not a physical benchmark; general-frequency input and independent review remain incomplete.

The earlier axial claims and their separate reviews follow below. This updated brief is local/public documentation, not authorization to send a message. The existing exact email drafts remain unsent and unchanged.

## Candidate claim, with direct predecessors acknowledged

For a local collisionless axial fast-flavor model, positive per-species angular moments define an exact frequency-compatibility question: SOME completion has the specified unstable temporal frequency, or NONE does. Classical Weyl/Hausdorff theory gives a constant-size disk diagnostic. That generic theory, the angular weight (1-u²), beam growth formulas and geometric-mean estimates are not claimed new.

The possible new theoretical content is the exact optimization over all compatible positive species measures, including sharp all-real-frequency fixed-k caps, supplied-second-moment all-k caps with essential stable corners, and the complete equal-second-moment finite-band envelope over all real frequencies, with endpoint/cubic/plateau transitions and explicit attainers. The paper supplies proofs, bounds their scope and preserves an exact counterexample to an initially tempting fixed-family optimum.

## Most important falsifying questions

1. Are these exact moment-constrained sharp extrema or equivalent statements already known? Inspect the quantifier over ALL positive species completions and the moment-preserving attainers, rather than the generic disk or beam algebra.
2. Is this uncertainty class useful for a real transport-validation decision? The diagnostic can only certify information compatibility, not actual growth, nonlinear conversion or astrophysical stability. Which independently supplied number moments, angular regularity or occupation constraints are needed before it changes a scientific decision?
3. Does any written argument miss a sign, azimuthal factor, current shift, positivity condition, singular moment case or distinction between number second moments and thermodynamic pressure?
4. Is an exact boundary-resolution guarantee useful enough to justify the implementation? The stronger generic conic comparison is faster end to end on typical jointly certified cases, so there is no general runtime superiority claim.

## Evidence that cuts against the candidate

The fixed counter-flux test did not meet its declared extra-information gate. All eight growing reproduced closure modes remain M1/M2 compatible; no different band or closure was substituted afterward. The idealized axial model omits collisions, arbitrary directions, nonlinear evolution and actual merger-snapshot validation. Atomic sharpness does not prove smooth exact attainment. No physical upper-k cutoff was verified. Priority searches are bounded negative evidence and independent reviewers so far were model agents, not external human referees.

## Review entry and minimal reproduction

Start with README.md, paper/research-draft.pdf, NOVELTY.md and CLAIMS.json. Then inspect offaxis-independent.md, independent-second-moment.md, finite-k-independent.md, independent-finite-k-review.md and independent-m2-locus.md. The fixed utility test is specified in M2_COMPATIBILITY_PLAN.md.

Standard-library-only checks:

    python3 -O verify_all.py
    python3 -O independent_continuous_replay.py
    python3 -O offaxis-independent.py
    python3 -O -m unittest test_second_moment_bounds test_independent_m2_locus
    python3 -B -O finite-k-cap-checks.py
    python3 -O verify_m2_screen.py

The separate old-M1 integration comparison needs the pinned NumPy/SciPy environment:

    .venv/bin/python -O -m unittest test_m2_independent_integration

Arb reproduction requires the pinned free python-flint dependency:

    python3 verify_counterflux_reference.py

The package manifest pins source and data. Its external verification record identifies the checked ZIP and proof replays from a fresh extraction. An exact certificate checks the mathematical model; it is not evidence of novelty, usefulness or publication readiness. A specialist in neutrino linear stability and a separate truncated-moment expert are appropriate reviewers. Any outreach, submission or public publication still requires the corresponding authorization.


The complete finite-band extension is documented in finite-k-transition-audit.md and independently reviewed in independent-full-band-review.md. Please assess priority of this exact constrained optimization, including its all-real enclosing-disk argument and fixed-k three-atom plateau attainers. The reduced-zero-mode moment diagnostic of Dasgupta/Mirizzi/Sen2018 is an explicit predecessor; its known moment-only growth calculation is not claimed here. General unequal supplied second moments at fixedK remain unsolved.
