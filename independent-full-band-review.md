# Independent review of the complete equal-second-moment finite-band envelope

Date: 2026-11-01. This review independently checked the supplied derivation and executable. The reviewer wrote this report, `independent_full_band_check.py` and its new result artifact. The proof, algorithm, supplied tests, paper and preexisting result artifacts were not modified. No numerical optimizer, external research, installation or credential access was used.

## Verdict

**No blocking mathematical or executable defect found.** The all-real-frequency upper bound, three attained regimes, finite-band maximization and exact rational branch selector are consistent. The proof supplies a universal argument; the successful finite tests below are corroboration and were not used as a replacement for it.

This verdict is for the exact axial, positive-measure uncertainty class in the supplied original dispersion model, with separate species moments `(N,0,Na)`, `N>0`, `0<=a<=1`. It does not endorse novelty, practical physics, computational superiority or publication readiness. It does not establish exact attainment in a smooth-density class.

## 1. Original-factor reduction and quantifiers

With signed moments `g0=g1=0`, the original identities give `I1=z I0`, `I2=z² I0`. Consequently

```
D_L = -1 + (z²-1) I0 = -1 + (U_plus-U_minus)/k,
T   = -1 - (U_plus-U_minus)/(2k).
```

Here `I0=(S_plus-S_minus)/k`, and the equal individual masses cancel the `-Nz` term in the difference of U values. Thus the targets are exactly `+k` longitudinal and `-2k` transverse. Since these targets are real, the needed separation is horizontal. The disk-containment argument correctly bounds that separation without replacing the lens by an asserted attainable enclosing disk.

For each fixed parameter choice and sector, the maximizing pair can differ. The theorem is an extremum over all compatible pairs, not a lower bound on any fixed background's growth. The two sectors are treated separately with `c=1,2`; species exchange produces the transverse sign.

## 2. Universal off-axis bound

I checked the transformed species disks against the supplied exact locus and rederived their squared center separation:

```
L²/N² = [(a²+y²-x²)² + 4x²y²] /
         [4y²(a+x²+y²)²].
```

Subtracting this expression from `R1²` gives exactly the identity in the proposed proof. For fixed `y>0`, `(a-w)/(a+w)` strictly decreases as `w=x²+y²` increases. The other needed inequality is also valid, although `R1²` need not itself be monotone in `x²` all the way to infinity. The correct pointwise comparison with its x=0 value follows from

```
(1+t)²(a+t+d)² - (a+t)²[(t+d-1)²+4t]
= d[2(a+t)(1+a+(3-a)t)+(1-a)(1+a+2t)d] > 0
```

for `d=x²>0`, `t=y²>0`, `0<a<1`. This distinction causes no error in the proof, which claims the pointwise comparison, not global monotonicity of R1.

The rewritten mixed radius has coefficients `1-beta`, `beta²`, and `beta(1-beta)`, all nonnegative. Hence every one of these enclosing radii at arbitrary real x is bounded by the corresponding radius at x=0. Choosing the x=0 minimizer proves the advertised upper bound for every real frequency part. The mixture radius is well defined: the true species locus is nonempty for these realizable moments, so any containing disk obtained by the completed-square inequality has nonnegative squared radius.

## 3. Sharp lens width and transitions

For x=0 the centers lie on the imaginary axis, with strictly positive separation delta. Minimizing the quadratic mixed radius gives the stated

```
beta0 = (1-a)(y²-a²)(a+y²)/(a²+y²)².
```

The three geometric cases are exhaustive:

* `beta0<=0` means `R0²+delta²<=R1²`, so the horizontal diameter of disk 0 is contained in disk 1.
* `beta0>=1` means `R1²+delta²<=R0²`, so the horizontal diameter of disk 1 is contained in disk 0.
* Otherwise the two circle-intersection points have height `Nr y/(a²+y²)` and half-separation `Na r/(a²+y²)` and attain the mixed-disk diameter.

These establish equality in the width bound, including both ties. The second case is equivalent to `y+a/y<=r`; the first is `y<=a`; they cannot overlap for `0<a<1`. The intermediate-circle formulas and the stated three width expressions agree at the boundaries.

Each active width decreases strictly: for the Gram-disk expression the derivative's sign is the sign of `-a+(a-3)y²-y⁴`, which is negative. The other two derivatives are immediately negative. Thus each prescribed positive k has a unique maximal y, and inversion gives the displayed fixed-k envelope. The cubic window exists exactly when `a<=3-2sqrt(2)`. Its endpoints exceed a, and the inverse width is strictly decreasing in y, giving the proposed k ordering and degenerate single-point window.

The endpoint growth increases until `k*=Nr/(ca)`. Along the cubic branch, `k(y)` decreases while `gamma(y)=(Na/c)(1+y²)/(a+y²)` decreases, so growth increases with k. Continuity at transitions therefore proves that the finite-band maximum is the fixed-k value at K.

## 4. Positive attaining measures

All three constructions satisfy support, positivity and the three separate moments.

* **Endpoint:** the two weights give mass N, first moment zero and second moment Na. The weighted measure `(1-u²)mu_plus` is `Nr delta_a`, so the transform difference is exactly `2Na r/(a²+y²)`.
* **Cubic:** `p,q>0`, `pq=a`, `q-p=y+a/y<=1-a` imply `q<=1` and `p<=1`. The displayed weights are positive and give moments `(N,0,Na)`. Direct evaluation at `z=iy` gives the Gram-disk rightmost value `R1+iV1`; reflection gives its leftmost counterpart. Thus their difference is exactly the claimed width. No assumption that arbitrary disk points have two-node witnesses is needed.
* **Plateau:** `0<y<=a<1` makes every denominator positive and each weight nonnegative. Their mass is N, first moment zero and second moment Na; `(1-u²)mu_plus=Nr delta_y` verifies the stated U values directly. At `y=a` the positive endpoint weight vanishes, agreeing with the endpoint pair. Three nodes are legitimate in the allowed measure class.

Exchanging the species changes the sign of U separation and verifies the transverse root in the original factors.

## 5. Exact selector and root enclosure

The executable normalizes `kappa=ck/N` and growth `g=c gamma/N`. Below the plateau threshold it defines `t=2ar/kappa-a²`, precisely the squared y obtained by inverting the endpoint-circle formula. The strict selector

```
(t+a)² < r²t
```

is the squared form of `sqrt(t)+a/sqrt(t)<r`, since t is positive. This is equivalent to t lying strictly between the cubic-window y² endpoints. At either boundary the endpoint and Gram inverse widths agree, and both inverse widths decrease strictly. Therefore this test selects exactly the cubic k interval. It requires no irrational threshold or transition calculation. Ties correctly use the endpoint branch and return the same exact value.

For the normalized cubic

```
p(g)=g³-a g²+a kappa² g-a kappa²,
```

`p(g)<0` for `0<g<=a`, because it equals `g²(g-a)+a kappa²(g-1)`. For `g>=a`, `p'(g)=g(3g-2a)+a kappa²>0`. There is exactly one positive root. In the selected cubic regime the Gram horizontal diameter lies inside disk 0, so its attained growth is strictly below r. Hence the bisection interval `[0,r]` brackets that unique root with the correct signs. The endpoint square is positive below the plateau threshold; its square-root branch is bracketed by the same interval. Perfect rational squares are detected exactly.

All arithmetic used for these decisions and enclosures is rational. Each bisection step retains the root; the final absolute width is at most `Nr/(c 2**bits)`. Squaring nonnegative bracket endpoints preserves the enclosure. The public API intentionally accepts exact rational values rather than arbitrary real-valued numerical inputs; this input restriction does not weaken the analytic theorem for real parameters.

## 6. Degenerate inputs and finite checks

At K=0, the original factors see only signed moments through degree two, which vanish, so there is no growing root. At a=0 positivity forces both species to `N delta_0`; at a=1 it forces the same equal endpoint pair. The code handles these cases before any division by a. A returned zero is the convention for the absence of positive growth. Invalid moment ranges, negative K, unknown sector and invalid precision are rejected by the relevant input checks. Scaling and the relation `F2(k)=F1(2k)/2` are consistent with the normalization.

I ran the supplied five unittest methods by import, so their `__main__` result-file writer was not invoked. All passed, including 1,000 exact disk identities, 72 cubic witnesses, 96 endpoint witnesses, 108 plateau witnesses, 36 transition cases, and 150 random enclosure/scaling cases.

I additionally performed **19,800 exact rational geometry-to-root checks** independent of the production selector. For `a=i/100`, `i=1,...,99`, and `y=j/50`, `j=1,...,100`, I chose the width branch using the comparisons `R0²+delta²<=R1²` and `R1²+delta²<R0²`, formed rational `k=NW(y)/c` at `N=7/3` for both sectors, and checked both the reported branch and containment of the exact `gamma=ky`. All passed: 1,086 cubic, 13,814 endpoint, 4,900 plateau. This is a finite audit of branch implementation and enclosure behavior; it is not a universal proof.

The audit is retained in `independent_full_band_check.py`; run `python3 independent_full_band_check.py` from this directory. It checks the pinned algorithm and proof hashes before importing the evaluator, uses explicit failures that remain active under Python optimization, and saves `results/independent-full-band-validation.json`. The saved rerun completed all 19,800 cases successfully at 30 enclosure bits. Its checker SHA-256 is `1c5351734e973696a9b8d6822d14a086c7b2591b0114fea53b54ee36e794d75f`. Its reference branch remains geometric and does not reuse the production selector.

## 7. Integrated-paper transcription check

I also read the integrated complete-envelope section in `paper/draft.typ`, source SHA-256 `78bae00bc26cf24c3312fcced5554dd248bc8139ba325e50ea903c72a7ab9451`. Its disk formulas, all-real-frequency inequality, three widths, cubic polynomial, transition definitions, piecewise rate, sector relation, positive attainers, exact selector and zero-second-moment asymptotic match the reviewed theorem. The per-point completion quantifier and physical/smooth-density limitations are retained. This was a source transcription check, not a visual audit of the compiled 11-page PDF.

One wording precision issue was sent to the author: "Both variable terms decrease away from x=0" can imply that `R1²` decreases monotonically with `x²`. It does not: it eventually increases toward a lower asymptote. The following polynomial comparison correctly proves that it is bounded above by its x=0 value. Suggested wording: "Both variable terms are bounded above by their x=0 values." This does not invalidate the displayed proof or result, but correcting it avoids an unnecessary false monotonicity interpretation. The reviewer did not edit the paper.

## 8. Remaining issues and provenance

There is no unresolved blocking issue in the reviewed universal theorem or rational implementation. The sole integrated-source wording precision issue is described above. No scientific priority or broader practical inference was evaluated.

This verdict is tied to these SHA-256 source versions captured at review start:

| Source | SHA-256 |
|---|---|
| `finite-k-transition-audit.md` | `f076d0fb272d073810fe3aaafede69d331c998abf787d6a94790047760a223e3` |
| `finite_band_envelope.py` | `e0fd9e4dce56dbb1c115690e8d7f14f4a419e6620cac2906d61886a448a9ddb5` |
| `test_finite_band_envelope.py` | `68cad6ca913b0eeaffb3e6c3f09bfb959c25f301d620ada4171cc0b883c2f39d` |
| `paper/draft.typ` | `c8059afeb8e643e6d89aad377b4fc4e013f3df32c3b7c60eea029ab6dea2f875` |
| `independent-m2-locus.md` | `257ffd32d93462ac8ce7238b7ddf994a8a4090911498224ca40bf3734bea584c` |
