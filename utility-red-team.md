# Scientific usefulness red team

2026-11-01. Independent review of README, PHYSICAL_SCOPE, M2_COMPATIBILITY_PLAN, physical-usecase-review, uncertainty-reality-review, the declared failed experiment and the current complete equal-second-moment claim. Primary literature and public dataset metadata were refreshed. No new benchmark, simulation, archive download, installation, expert contact or public write was performed. This note assesses utility; it does not determine priority or re-review the universal proof.

## Verdict

The result has demonstrated **mathematical information-limit value** in its precisely stated axial class. It supplies attained worst-case rates and a reproducible evaluator for a completion problem that the existing moment data do not resolve uniquely. That is a defensible theoretical result to assess for publication if priority survives.

It has **not demonstrated operational astrophysical usefulness or a high-impact numerical advantage**. The actual fixed closure-exclusion test failed. The generic continuous comparator is typically faster on jointly certified frequency points. The full equal-M2 finite-band theorem is a different task from that pointwise comparator, so the timing loss must not be misapplied as a direct runtime comparison for the new closed-form bound. Conversely, the new formula cannot erase the old loss or count as a demonstrated transport improvement.

The strongest presently defensible practical use is a supplementary mathematical regression oracle. Even that needs honest competitors: checking the original dispersion residual when a full angular background is known, standard analytic beam tests, realizability checks, and the combination of already available growth caps. It is not yet a justified reason to build a larger simulation or uncertainty pipeline.

## 1. The claimed isotropic improvement needs the strongest existing cap

The parent identified an omitted internal baseline: the finite-band figure's older upper bound combines the Gram cubic and the all-k second-moment cap but omits the previously proved M1 fixed-k cap. I independently checked the resulting comparison.

At N=1,Q=1/3, longitudinal growth for0<k<2 has sharp square

    F(k)^2=k(4-k)/9,

and is2/3 for k>=2. The M1 cap is sqrt(k/(k+2)); the all-k Q cap is2/3. Their minimum changes at k=8/5. For0<k<=8/5,

    min(M1,2/3)^2/F(k)^2 = 9/[(k+2)(4-k)] <= 9/8.

For8/5<=k<=2,

    min(M1,2/3)^2/F(k)^2 = 4/[k(4-k)] <= 25/24 < 9/8.

Above2 both plateau. Adding the Gram cubic only tightens the older comparator. Thus the relative improvement of the new cap over the strongest combination is at most

    1-sqrt(8/9) = 5.7191... percent.

The limiting equality occurs as k approaches zero. A comparison suggesting an order-of-magnitude improvement for isotropic second moments at small k is therefore not supported. The complete theorem still adds sharp attainment and resolves the transitions, but its practical numerical increment in this particular important idealized class is modest. Extreme second moments can behave differently; choosing them after this finding would not establish a physical use case without provenance.

## 2. An impossibility result for a tempting spurious-mode application

Suppose a code replaces a continuous positive species background by positive atomic quadratures nu_plus,nu_minus that preserve every moment supplied to the diagnostic. Let omega be an exact growing root of the quadrature's dispersion equation, under the same axial model and real current shift.

Those two quadratures are themselves admissible members of the diagnostic's positive-measure completion class. They therefore witness that omega is possible. An exact existential moment diagnostic MUST accept this mode, even if it is absent from the intended continuous background.

This is a set-membership proof, with no numerical experiment required. It rules out using the current diagnostic as a general filter for discretization-induced spurious modes when the quadrature preserves the supplied moments. If quadrature moments differ, rejecting its mode may diagnose that mismatch; it does not by itself identify the continuum error. Adding moment uncertainty containing the quadrature makes this obstruction persist.

The practical problem is real: Morinaga and Yamada, PRD97,023024(2018), SectionsIII-IV andV, explain how angular discretization replaces branch-cut structure by artificial poles and demonstrate analytic angular integration as a remedy. Their method uses more distributional information. Primary source: https://arxiv.org/html/1803.05913 ; publisher identity/date: https://doi.org/10.1103/PhysRevD.97.023024 . A generated HTML header date must not replace the2018 publication identity.

This argument also applies to any closure whose mode can be shown to arise from positive quadratures preserving those inputs. Whether a specific quantum closure has that representation requires a separate proof; it is not assumed here.

## 3. What the theorem can and cannot decide

Within the equal per-species(N,0,Q) class, choosing identical species gives a stable completion. For every0<Q<N and positive finiteK, the accepted theorem also gives an unstable completion with positive growth. Thus these supplied moments alone cannot decide whether the actual member is stable. The exact maximum quantifies that ambiguity. Its positive value is not a prediction that flavor conversion occurs.

A sharp upper cap can rule out a proposed rate above it in the stated model. No current retained experiment shows that this rejection changes a real transport decision. The fixed naturally sourced counter-flux closure test gave zero added exclusions. An artificially inflated mode, an injected sign error or an invented threshold is a sanity test, not a demonstration of practical scientific utility.

FiniteK is a legitimate analysis band, but no physical ultraviolet cutoff has been justified. A numerical resolution limit does not prove stability of unresolved modes. Nor does F(K) being close to the unrestricted maximum prove that a numerical band captures a similar fraction of the fastest growth for each background: the maxima over background and wave number may be attained by different backgrounds. That change of quantifiers would need a new theorem.

The instability-growth cap also does not by itself set a stable numerical time step. Real streaming frequencies, the full discretized operator and possible nonnormal transient behavior remain relevant. No CFL or nonlinear-conversion claim follows merely from an upper bound on imaginary eigenfrequency.

## 4. Physics and data requirements remain binding

The current bound assumes equal species density, individual zero flux, equal supplied number second moments and axial propagation. Published transport workflows generally evolve density and flux and close the next moment. Froustey, Kneller and McLaughlin's quantum maximum-entropy closure explicitly evolves number moments and derives a closure in the small-coherence regime; SectionI states the axial presentation, and SectionIV compares to multiangle stability calculations. This is a serious conventional competitor and evidence of the actual modeling task, not evidence that its closure-estimated Q is independent truth. Primary source: https://arxiv.org/html/2409.05807v2 . The version is dated10March2025.

Froustey etal's earlier moment linear-stability framework evaluates unstable rates and wavelengths in a neutron-star-merger model and discusses successful and unsuccessful predictions. The candidate must improve a specified decision in that setting, rather than only reproduce its general motivation. Primary source: https://arxiv.org/html/2311.11968v3 , SectionsIII-V and AppendixB; v3 dated14January2025 includes the corrected equations.

The refreshed Zenodo13679313 metadata lists initial number four-fluxes and measured growth rates, but not independently supplied initial second number tensors in the exposed schema. The final-state standard deviation is a variation of the time series being averaged, not initial-moment uncertainty. Its simulation fields include general three-dimensional fluxes. These facts do not justify applying the equal-zero-flux theorem to a snapshot, treating a closure-derived second moment as measured, or projecting away nonaxial structure. Source: https://zenodo.org/records/13679313 . No file bytes were downloaded, so this is a schema assessment, not an exhaustive assertion about every archived field.

The new exact nonaxial projection-loss example gives an especially concrete warning: equal isotropic second tensors do not license azimuthal averaging of an arbitrary background. The current theorem remains valid in its own class; the example does not show that its numerical cap is false in every broader class.

Atomic measures and smooth approximations are useful mathematical completions. They are not automatically consistent with energy-resolved fermionic occupancy, transport smoothness, collisions or a calibrated uncertainty model. None of those constraints should be invented solely to create a successful benchmark.

## 5. What would earn another executable utility experiment

A meaningful operational test needs an exogenous candidate mode and a decision defined before running this classifier. The strongest remaining narrow route is a closure audit with independently supplied number moments and a justified axial setup, using a published discrepancy not selected after scanning for exclusions. It must satisfy all of the following:

1. The background moments and any uncertainty are independently known from the prescribed physical setup or data, rather than borrowed from the closure being tested. The relevant symmetry, species/current shifts and number weighting are verified.
2. The decision is explicit: reject a claimed linear mode, require a more detailed angular calculation, or prove a specified growth threshold impossible in the valid band. Threshold and band come from the problem; they are not tuned after results.
3. Compare against the minimum of every applicable existing cap, the zero-mode and realizability diagnostics, the strongest continuous moment solver, and a direct original-dispersion check when the full angular distribution is known.
4. Require a naturally occurring discrepancy that the new information excludes with margin larger than justified input uncertainty, while the existing moment-only alternatives do not provide the same decision at comparable cost. If a generic conic certificate provides the same exclusion, describe the contribution as proof representation or exactness, not a unique physical capability.
5. Preserve all planned cases and a negative outcome. Do not reuse the failed counter-flux test under a wider grid or changed closure and present it as the same gate passing.

I did not locate such a fully specified new case in the bounded primary-source review. Therefore there is **no currently earned large executable campaign**. The existence of current closure-development work establishes demand for better models, not demand for this particular information bound.

A smaller code-verification study can be honest, but should require a documented defect or an independently maintained implementation and should compare with ordinary residual checks and analytic tests. Its result would be software-validation usefulness, not astrophysical inference or high impact. A collection of self-injected errors alone does not satisfy the requested scientific usefulness claim.

## 6. Recommended decision

Preserve the complete theorem as a candidate theoretical information-limit contribution. Correct the omitted-cap comparison before discussing its improvement. Do not market exact boundary stress cases, projection-loss examples or the tested closure inclusions as demonstrated astrophysical benefit.

The next material evidence is focused expert assessment of whether the exact constrained extremum is new and whether its uncertainty class answers a worthwhile scientific question, or a genuinely new provenance-supported physical constraint/decision. No outreach was authorized or sent by this reviewer. Absence of external feedback does not prove uselessness, but it also does not turn a narrow mathematical success into the requested high-impact outcome.
