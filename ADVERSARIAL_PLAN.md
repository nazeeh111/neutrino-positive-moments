# Exact-boundary and scaling stress, before its run

This is an adversarial analytic check, not statistical held-out physical data. StageB504-point gate already passed. No original criterion is relaxed.

Use known sharp-envelope rational boundaries: longitudinalN=1,k=2/3,gamma=1/2; transverseN=1,k=1/3,gamma=1/4. For each, test the boundary and below/above shiftsdelta={1/100,1/10000,1/10^8,1/10^12}. Scale N,k,gamma together by{10^-50,10^-6,1,10^6,10^50}. All90 classifications must match the independently proved threshold and pass the independent exact checker. Any error fails this correctness stress. Record unresolveds, never round a boundary into an answer.

At unit scale compare32,128,512-point angular-grid LPs. Repair every claimed feasible primal with exact rational equalities as before. Where the grid claims infeasibility, attempt an exact Farkas separator: numerical dual candidate, rational coefficients, outward repair of mass multipliers, then verify ALL grid columns nonnegative and right-hand-side dot dual strictlynegative. No claim about continuous impossibility follows from such a grid certificate.

At the longitudinal boundary, the unique extremal interior angle is a=3/8; at the transverse boundary it is a=3/8 too. None of these three uniform node grids includes±3/8. Exact disk/witness construction should therefore expose a discretization limitation at these designed boundary tests. This is a deliberate adversarial example and must not be presented as a naturally sampled accuracy rate or a numerical solver failure.

SLSQP is a generic numerical proposer, not a ground-truth oracle. Record its exact-certificate outcome and unresolveds. Timing from the current development runs is exploratory; no controlled speedup claim is part of this check.
