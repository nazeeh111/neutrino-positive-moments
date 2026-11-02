# What the full nonaxial cap establishes, and what it cannot predict

Status: primary analysis and exact checks, not independent review. This is a consequence of the full-frequency theorem, not a separate claim to a new instability mechanism or high-impact algorithm.

## An entire interval of indistinguishable maximum growth rates

For the supplied per-species data `(N, zero vector flux, second tensor NI/3)`, the actual maximum temporal growth over all real wave vectors can take **every value in [0,2N/3]**. Thus the theorem's upper endpoint is an exact description of the information missing from these moments, not an accurate point prediction for a particular background.

Let nu+,nu- be the five-atom axial attainer and its inversion, each with mass N and the stated full moments. They attain gamma0=2N/3 at k0=2N. Let sigma be any common positive measure with those same moments, for example the isotropic six-axis octahedron. For 0<=epsilon<=1, define

```
mu_s(epsilon) = (1-epsilon)*sigma + epsilon*nu_s.
```

Every species has exactly the same supplied moments for every epsilon, while

```
G(epsilon) = epsilon*(nu+ - nu-).
```

Only G enters the shifted fast linear operator. At k=epsilon*k0 and omega=epsilon*omega0, its polarization equals that of the original attainer. Conversely, the same signed background can be represented by the effective positive species epsilon*nu_s, with mass epsilon*N and isotropic second tensor epsilon*NI/3. Applying the proved upper bound to those effective species excludes every growth rate above epsilon*2N/3. Therefore the actual all-k maximum is **exactly** epsilon*2N/3. At epsilon=0 the signed interaction vanishes and every streaming eigenvalue is real.

This construction is within the arbitrary-positive-measure class. It does not establish exact atomic attainment under occupation or smooth-density constraints. It varies the attaining wave number and therefore does not solve a fixed-k inverse problem.

The checker `exact_growth_range.py` verifies positivity, all full moments, cancellation of the common component and the original beam-operator eigenvector for 14 rational cases, including both endpoints and two masses. The scaling argument and the existing universal theorem, rather than those finite cases alone, prove the full interval.

## A sharp lower limit on point-prediction accuracy

Any deterministic point estimator using only these identical moments must output one number c for the whole family. Its worst absolute error is at least

```
max(|c|, |c-2N/3|) >= N/3.
```

The constant midpoint c=N/3 attains that worst-case error because the full range is [0,2N/3]. Hence N/3 is the optimal worst-case absolute error over this uncertainty class. This is a mathematical accuracy limit, not a recommended physical estimator. A nonzero point estimate cannot have uniformly bounded relative error as the true growth approaches zero. Returning zero instead incurs 100% relative error on every positive member.

No decision of the form “the true maximum exceeds theta,” with 0<theta<2N/3, is determined by these moments: the family contains backgrounds on both sides. The cap can exclude impossible reported rates above its endpoint, but cannot identify the rate or certify instability below it.

## Comparison with a serious published metric

[Richers, *Evaluating Approximate Flavor Instability Metrics in Neutron Star Mergers*, arXiv2206.08444v3, Section IV.2, Eqs.21–24](https://arxiv.org/html/2206.08444v3#S4.SS2) treats a geometric-mean crossing-depth quantity as a growth-rate estimate, not a universal sharp bound. Its direct version uses angular crossing integrals; the maximum-entropy version adds a distributional closure. The discussion following Eq.24 says that, in the studied snapshot, the existence and eventual amount of conversion matter more than precise growth because the estimated rates already exceed competing rates. These are source-reported results, not a reproduction here.

For our own finite measures, define nonnegative crossing masses J+ and J- as the total positive and negative magnitudes of G, avoiding signed-negative-component notation. The attainer's two species have disjoint support, so along the mixture family J+=J-=epsilon*N. The geometric-mean magnitude proxy is then epsilon*N, whereas the exact maximum is epsilon*2N/3. Their ratio is 3/2 for epsilon>0. This is an exact constructed-family comparison, **not** an astrophysical accuracy benchmark or a refutation of the paper's deliberately approximate estimate.

Crucially, J+ and J- are not determined by the supplied moments: they vary from zero to N along this very family. It would be an invalid comparison to credit a moment-only theorem with outperforming a proxy that is given extra angular information, or to claim that an approximate proxy is a competing rigorous certificate.

## Revised usefulness verdict

The full theorem removes the axial restriction from one sharply specified uncertainty bound. It closes the earlier mathematical gap between the safe capN and the attainable2N/3 in that class. That is a substantive mathematical extension, but the capN was the project's earlier safe estimate, not an established best published operational baseline.

The demonstrated scientific content is a sharp information limit and a possible supplementary regression constraint. No operational astrophysical decision, transport accuracy improvement, runtime advantage or nonlinear-conversion prediction has been demonstrated. The old counter-flux utility test remains failed, and the old conic timing comparison remains unchanged. Replacing an upper bound by a point estimate, or equating “some completion grows” with the actual background's instability, would violate the new exact range result.

A useful external assessment should now ask whether this exact three-dimensional moment-constrained extremum is known or a routine corollary, and whether any realistic workflow needs its exclusion guarantee. Those are separate questions. Primary-source searches and exact arithmetic have not answered either with independent expert evidence. The pending expert-message authorization has not been granted; no outreach occurred.
