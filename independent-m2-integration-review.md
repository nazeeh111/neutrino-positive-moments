# Independent integrated-M2 review addendum

Review date: 2026-11-01. Scope: only the integrated general-M2 locus section
and counterflux utility paragraphs, not the finite-band theorem, external
priority, closure derivation, PDF layout, or an all-root census.

## Versions actually inspected

The requested `paper/stagedraft.typ` does not exist. The integrated source is
`paper/draft.typ`, reviewed at lines 139-152 and 190-198, with the earlier
transform definitions and dispersion target checked for notation/signs.

SHA256 snapshots:

- `paper/draft.typ`: `562a4e79a60960def2f4cfb3f933dba1eba1dccba4915eac4e67d67f07be9b5a`
- `independent-m2-locus.md`: `257ffd32d93462ac8ce7238b7ddf994a8a4090911498224ca40bf3734bea584c`
- `m2_frequency.py`: `cdf04746136abad80ed83eaac13f15dcd06ed8a83b1ef3b301ee53cf51cd7da7`
- `m2_counterflux_screen.py`: `a69f405bce207ad8ba0eff89eaa11c2e6e20a3351d4fb8e43fe28542f24ffce1`
- `results/m2-counterflux-screen.json`: `2f8435526cdee3c6a751ff3a71ed8b31097762fa5dc6e3da95b800d5fda2b81c`
- `results/counterflux-modes.json`: `8e7785bfc7e580b158815710f9194c20f60dda2d7357f64c29380af120d7ae8d`

All seven source hashes saved inside the M2 screen and its reference-result
hash match the corresponding current files. This audit read those artifacts;
it did not rerun the already-passing nine tests or numerical reference search.

## One mathematical wording correction required

At reviewed line 150, the sentence "singular cases reduce to one or two"
must specify **singular H**. Singularity of L does not permit that claim.
For example N=1,M=0,Q=2/3,z=i,S=-2i/3 has implied moments
(2/3,0,1/3,0,1/3), det(H)=1/27, and L=diag(1/3,0). Its L is singular,
yet any matching measure requires at least three atoms. The positive
three-atom measure is (delta_-1+delta_0+delta_1)/3.

Suggested replacement: "A constructive Gauss-Radau argument supplies at
most three atoms from H,L; singular H cases reduce to one or two."

Also clarify that the k in the displayed closure growth formula is the
**signed reduced wave number before reflection**. The formula is correct
as written with k=K-2/3. The positive-k classifier instead uses |k| and
reflects both species fluxes when k<0; substituting |k| into the displayed
unreflected closure formula would give the wrong rate. The present preceding
definition supplies the right meaning, but an explicit label removes this
sign ambiguity.

## Findings that pass

The two disk formulas have the same centers, radii, transform convention,
and signs as the independent derivation. For N>0, Delta>=0 and Q<=N imply
the other realizability inequalities. The strictly positive d is justified
by y>0. The Delta=0, Q=N and N=0 cases are explicitly separated. The text
limits the nonsingular Schur-complement equivalence to Delta>0, and does
not silently replace positive semidefiniteness by a determinant test.

The reweighted degree-four H and L matrices and t4 recurrence are correct.
At most three atoms is the correct general witness bound; three can be
necessary. The proof objects may establish existence through those exact
moment matrices without listing reconstructed atoms. This does not claim
that the current interface emits a quadrature witness.

For rational supplied moments and rational wr,g,k, the coordinates
S=p+i*sqrt(g)*q have rational centers and metric diag(1,g). Pair candidates
need one real quadratic field for their coefficients p,q. This is **not**
a claim that the entire complex transform, its real and imaginary parts
in the ordinary basis, or recovered atom nodes have degree at most two
over the rationals. The draft's explicit phrase "coefficient field" is
accurate, and the earlier atom-recovery section allows another square root.
The negative checker independently obtains its necessary ellipses from
determinants; positive checking uses every principal minor and the original
dispersion factor. No production geometry is imported by that checker.

The fixed screen contains eight growing longitudinal closure modes, all
eight M2-possible and M0/M1-possible, zero M2-impossible, zero added-information
exclusions, and utility_gate_passed=false. The draft reports **failure** of
the declared utility criterion without changing the fixed physical-K band.
Transverse closure nilpotence is not promoted to absence of angular
transverse instability.

All 26 designated rational points near linear/entropy numerical references
are stored as distinct exact point checks, and all are M2-possible. The
13 independent linear-input brackets contribute 26 separately checked
endpoints. The draft correctly says that these point certificates neither
identify a rounded reference as an exact root nor certify a full interval.
Existence of actual linear-input roots is attributed to the separate Arb
sign/continuity argument. Entropy references remain numerical and have a
different second moment. Compatibility is not promoted to correctness of
the closure's physical growth prediction or practical astrophysical utility.

Verdict for the inspected version: the mathematical content and negative
utility interpretation pass, subject to the singular-H wording correction.
The signed-k clarification is recommended. No core, proof, test or paper
source was edited in this review.

## Resolution checked against current source

The current `paper/draft.typ` was read again after integration. Line 150 now
says "singular H cases reduce to one or two," and line 192 explicitly labels
the closure formula's k as the signed reduced wave number before axis
reflection. Both requested clarifications are resolved. Current source
SHA256: `17f8c55eee161a4806ca95ab75b2a9ba33fa7af9212c617a62203002b03a2cec`.
The scoped integrated-M2 review therefore passes without an outstanding
mathematical or interpretation issue. No tests were repeated and only this
review addendum was edited during the resolution check.
