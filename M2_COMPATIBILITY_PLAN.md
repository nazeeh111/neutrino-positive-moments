# Fixed-second-moment compatibility: experiment fixed before classifier results

Date: 2026-11-01. This extends the information supplied to the same axial fast-flavor dispersion problem. The original M1 screen, stress gates, proofs and failures remain unchanged. No speedup or astrophysical stability claim is proposed.

## Mathematical question

For positive species measures on [-1,1], supply each number triple (N,M,Q), with N >= 0, M² <= NQ and 0 <= Q <= N. Decide whether SOME completion has a specified unstable frequency and real nonzero reduced wave number. Extend the exact frequency representation to omega = wr + i sqrt(gamma2), with rational wr,gamma2>0 and rational k>0. This includes exact rational squared growth rates from a real closure matrix. Represent S=p+i sqrt(gamma2)q. All disk centers in (p,q) and squared radii are rational, with metric diag(1,gamma2). This is an algebraic representation, not rationalization of a floating root.

Per-species proposed locus: one disk from the Gram matrix of (1,u,1/(z-u)), and one from the positive weighted measure (1-u²)mu. Prove necessity AND sufficiency through the degree-four Hausdorff moment conditions, including zero mass, zero variance and endpoint-supported cases. These are direct classical moment-theory specializations, not new generic mathematics. General fixed-second-moment values can require three atoms; never retain the M1 two-atom guarantee here.

Every possible decision must pass a separate standard-library checker using the original degree-four moment matrices and original sector equation. Every impossible decision must pass independent disk construction and intersection enumeration. The checker may reuse its old independent exact real-field arithmetic but must not import production M2 formulas or a numerical optimizer. Include malformed, false-status and moment-invalid controls, including optimized Python. Compare with known positive atomic measures and the unchanged M1 information class. Necessary invariants: M2-possible implies M1-possible at identical inputs; forced single species has no growth.

## Physics validation target and fixed points

Use the controlled continuum input in Richers et al. 2021, not a claimed measured merger snapshot: dmu±=(1±u)du/2 after density normalization. Its exact number triples are (1,+1/3,1/3) and (1,-1/3,1/3). This is distinct from a maximum-entropy reconstruction with the same M0/M1 but a different M2. Never call that different completion itself a physical error.

Evaluate physical axial K=j/8 for j=-8,...,8, both longitudinal and transverse sectors. This finite analysis band is chosen for reproducible code validation, not a physical ultraviolet cutoff. Use the verified real current shift reduced k=K-2/3 and omega=physicalOmega. Negative reduced k is handled by reflecting the axis and each species flux; include the special k=0 limit only if it occurs, without dividing by zero. Do not silently treat physical and reduced wave numbers as equal.

A separate implementation reproduces the published Minerbo quantum-moment closure with f=1/3 and chi=17/45, plus the full angular linear input and a maximum-entropy completion as continuous numerical references. It checks the analytic closure frequency against its independently assembled source matrix. Preserve failed reference checks before correcting implementation mistakes. Store exact rational closure gamma2 where available; approximate full-angular roots remain numerical estimates unless separately validated.

## Decision criteria

A useful positive outcome requires at least one faithfully reproduced growing closure mode on the FIXED physical-K grid that is exactly impossible for the independently supplied (N,M,Q) class, while its frequency is possible in the M0/M1-only class. That would show a concrete extra decision enabled by supplied second moments. It does not condemn a closure for a different information class or establish nonlinear astrophysical consequences.

If no such exclusion occurs, report the failed utility gate. Do not widen the K band, change moments, switch closure after inspecting results, or count artificially overlarge-growth controls as utility success. Full-angular and maximum-entropy modes provide distinct consistency controls, not exact truth inferred from an optimizer status. No runtime comparison is needed for this information-capability question; the prior conic timing loss remains a loss.

## Parallel theoretical check

A finite-k zero-flux extension is exploratory. A simple two-endpoint/interior family was conjectured to be sharp but a preliminary continuous-grid check challenges it. Preserve that counterevidence. An imaginary-axis or finite-grid observation cannot establish the maximum over all real frequencies. This question is secondary to the counterflux validation target and supplies no justified physical cutoff.

Sources: https://arxiv.org/html/2101.02745 and https://arxiv.org/html/2311.11968v3 . The latter contains the 2025 corrected closure equations; generated HTML header dates do not change publication/version identity.
