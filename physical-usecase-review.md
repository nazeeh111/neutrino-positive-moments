# Physical use-case review: finite axial wave-number bands

Reviewed 2026-11-01. Scope: targeted primary-source discovery plus inspection of the current README, PHYSICAL_SCOPE, NEXT_RESEARCH, moments.py and second_moment_bounds.py. No simulation or archive download was performed. This review owns only this file.

## Decision

A legitimate next experiment exists: use exact number moments to audit whether a proposed *linear axial mode* from a moment closure is compatible with any positive angular background in a specified wavelength band. The strongest immediately reproducible source is a published counterflux code benchmark, not an astrophysical snapshot. It requires extending the second-moment treatment to nonzero species fluxes. Finite k alone in the current equal-mass/zero-flux class does not resolve this utility bottleneck.

I found no source establishing a physical upper wave-number cutoff for these collisionless fast-flavor examples. A simulation's maximum represented wave number is a numerical limit. An imposed wavelength band is a legitimate scope restriction for a validation question, but excluding growth there does not establish astrophysical stability. The bounded search does not establish that a physical cutoff cannot be justified in another model.

## Verified source facts

**Prescribed axial counterflux input.** Richers et al., *Particle-in-cell simulation of the neutrino fast flavor instability*, PRD 103, 083013 (2021), Sec. III A, Eq. (8), prescribes equal electron-neutrino and antineutrino number densities N=4.89e32 cm^-3, opposite number flux factors ±1/3, and no initial heavy flavors. Its continuum angular prescription is linear in direction cosine. This is controlled input, not a fitted measurement. The work studies axial symmetry preserving and breaking modes and provides a resolution study. It identifies open code and a companion data DOI. [Primary paper](https://arxiv.org/pdf/2101.02745), [published article](https://doi.org/10.1103/PhysRevD.103.083013), [companion archive DOI](https://doi.org/10.5281/zenodo.4423095).

**Same low moments, different completion.** Froustey et al., PRD 109, 043046 (2024), Sec. IV/Table 1, uses those densities and fluxes for its Fiducial counterbeam test. Its number-moment definitions are Eq. (6); the authors explicitly distinguish their second moment from actual pressure. Its maximum-entropy completion differs from the earlier linear distribution. The paper explicitly attributes the difference in fastest-mode predictions to that difference. Its closure analysis also discusses spurious modes and potentially missed instabilities. These observations support testing compatibility, not an assertion that a reported closure mode is already proved impossible. [Primary full text](https://arxiv.org/html/2311.11968), especially Secs. II.1, III.2 and IV. The publisher records a 2025 correction to Eq. (34): [version/status](https://doi.org/10.1103/PhysRevD.109.043046).

**Available astrophysical data do not yet fill the M2 requirement.** The public 2024 Richers et al. archive documents initial/final *number* four-fluxes F4, measured linear growth rates, and merger-snapshot provenance. Its listed data schema does not supply independent species second angular number tensors. That is a real M0/M1 resource, but not verified independent M0/M1/M2 truth; it also contains general three-dimensional fluxes. Metadata availability is verified, file bytes and a usable axial subset are not. [Dataset description](https://zenodo.org/records/13679313), Emu_data/F4_initial and SpEC_data sections.

The 2021 companion archive landing page returned an internal accessibility error through the browsing tool. No alternate download path was attempted. Its public availability is a paper report, not a verified local data import. The continuum benchmark input is fully specified by the accessible paper and does not require that archive.

## Exact benchmark moments, derived here

After integrating azimuth, the unperturbed continuum input is

    dmu_plus(u)  = N(1+u) du/2,
    dmu_minus(u) = N(1-u) du/2,       -1 <= u <= 1.

Direct integration gives

    (M0, M1, M2)_plus  = (N, +N/3, N/3),
    (M0, M1, M2)_minus = (N, -N/3, N/3).

Both measures are positive. Their second moments are known from the prescribed angular input; no closure-derived pressure is promoted into supplied truth. These are continuum benchmark values. A finite particle quadrature and the small flavor seed can change realized moments and must be measured separately in any code replay.

The existing sharp M2 cap assumes **each species flux is zero**. The two nonzero fluxes here do not cancel that requirement. Equal total species masses and zero *combined* flux are insufficient. Nor can a real phase/current shift erase the species moment constraints. Thus feeding Q_plus=Q_minus=N/3 into second_moment_bounds.py would not produce a valid cap for this benchmark.

## One reproducible, decision-relevant experiment

1. Normalize to N=1 and use the two rational moment triples above. Extend the exact positive-measure frequency-feasibility diagnostic to supplied M2 with nonzero M1. Keep the original axial longitudinal and transverse factors and all real frequency parts. Track the real current shift between physical and reduced k explicitly. This is a substantive requirement, not an inference from the current zero-flux proof.
2. Fix a wavelength *analysis band* before searching. A simple reproducible choice is physical |K_z| <= mu0/c, with mu0=sqrt(2) G_F N/hbar in consistent units; this is dimensionless |c K_z/mu0| <= 1. It is a chosen band, not a physical ultraviolet cutoff. Also test the full-k cap where provable. Do not silently impose |reduced k| <= 1 because the signed current is nonzero here.
3. Independently solve the full angular linear problem for the specified linear distributions, then compare with the published maximum-entropy reconstruction and the quantum-moment closure. For the latter two, distinguish M0/M1-only compatibility from compatibility with the independently prescribed M2=N/3. If an alternative reconstruction has a different M2, that is an intentionally different information class, not an error.
4. Generate exact or validated inclusion certificates for full-angular control modes. Ask whether the closure's reported modes lie outside the possible region for the supplied moments in the analysis band. An exclusion is a concrete result: that mode cannot represent any positive angular completion of those inputs under the local model. An inclusion only establishes compatibility, and cannot validate the closure's actual growth prediction.
5. Include intentionally overlarge-growth controls, but count those as checker sanity tests. The utility gate is a naturally occurring closure discrepancy excluded by the independent moment constraints, or a quantitatively useful strict cap unavailable from the existing all-k M2 theorem. Report the negative outcome if no such case occurs. A smaller bound merely because the band was made tiny is not a successful scientific demonstration.

This experiment can change a code-validation decision: reject an incompatible closure-generated linear mode or identify where a more detailed angular calculation is necessary. It cannot decide nonlinear flavor equilibration, merger nucleosynthesis or macroscopic evolution. A small temporal periodic box and its late-time averages do not provide those conclusions.

## Why finite k is currently a secondary extension

Finite k regularizes the current theorem's near-Q=0 singular limit, where its optimizing k diverges. That is mathematically informative. In the physically motivated public counterflux input, however, the unsatisfied condition is nonzero species flux, not a missing finite-k theorem in the zero-flux class. A transport validation use case also needs uncertainty in supplied moments and eventually nonaxial geometry. Those requirements should not be hidden by calling the resolved simulation band physical.

Recommendation: retain the bounded-k theorem as a useful analytical extension if inexpensive, but do not sell it as the missing astrophysical use case. Prioritize the rational counterflux M0/M1/M2 compatibility experiment above. If that experiment produces only inclusions and no useful new exclusion, record that the physical-impact case remains unestablished rather than expanding a catalog of examples.
