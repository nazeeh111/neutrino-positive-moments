# Rational-frequency exclusions with a standalone verifier

This interface implements the sufficient certificate in general-moment-certificate.md for arbitrary rational real frequency, positive rational growth and nonnegative rational reduced wave number. It removes the earlier requirement to supply rational a,b,h,t. It remains an exclusion method: **unresolved does not mean possible or stable**.

The proposer and verifier are separate implementations written by the same primary agent. The verifier imports only Python's standard library and has been run in an isolated temporary directory without the proposer or repository modules. This is separately checkable computation, not independent scientific review of the theorem.

## Input and use

Input JSON has exactly five fields: plus, minus, k, wr, gamma. Each species has mass N, a three-entry vector F, and a symmetric3x3 tensor Q. All scalar values must be JSON integers or rational strings. Floats and booleans are rejected. The tensor trace must equal N and the moment Gram matrix must be positive semidefinite. The mathematical claim is conditional on the stated positive-measure model; this interface is not a general physical-data validation service.

Use the bundled nonzero-real-frequency example:

```sh
python3 -B -O nonaxial-exploration/general_frequency_certificate.py examples/nonaxial-general-frequency-input.json
python3 -B -O nonaxial-exploration/verify_general_frequency.py examples/nonaxial-general-frequency-input.json examples/nonaxial-general-frequency-result.json
```

The first command proposes and verifies an exclusion before returning it. The second replays the saved certificate independently of the proposer. A successful verifier exits0; a rejected certificate exits1. Neither command overwrites the supplied files.

The z axis is the direction of the reduced wave vector. k is its nonnegative magnitude; moments must already be expressed in that frame. Known real matter/current shifts must be applied consistently. Gamma is supplied directly, not as gamma-squared. These are certificates at **exact rational points**, not intervals around measured moments or numerical roots. Rounding an observed mode and excluding the rounded point does not exclude the actual mode.

## Exact square-root enclosures

For k>0 put x=wr/k, y=gamma/k, and v=1-x²+y². The needed square root has positive real component

```
s = sqrt(v²+4x²y²),
p = sqrt((s+v)/2),
q = -xy/p,  H2=p²+x².
```

The proposer encloses s and p with rational endpoints using integer square roots. It begins at32 bits and doubles up to the requested limit (default512, maximum4096). Rational perfect squares are recognized exactly. Cancellation that leaves a zero lower bound for p triggers refinement; it is never replaced by a floating approximation.

The moment-certificate coefficients require no further square roots:

```
alpha = p*(2-1/H2), beta=-y/p,
c_plus=(p-q)²/H2, c_minus=(p+q)²/H2,
e_plus=c_plus-1+1/H2, e_minus=c_minus-1+1/H2.
```

The known reduction matrix M0 is exact rational-complex arithmetic. Its upper base entry is `Re((1-r²)*M0_00)-(x*g0+g1z)/k`; its mixed entry is `q*Re(M0_0i)+p*Im(M0_0i)`; the lower block is `-Re(M0_ij)`. Rational interval arithmetic encloses the resulting real symmetric exclusion matrix L.

## What a certificate proves

An ordinary certificate contains a version tag, a canonical-input SHA256 and the two root intervals. It does not ask the verifier to trust a proposed matrix or solver status. The verifier checks, with exact rational comparisons,

```
0 <= s_lo <= s_hi,  s_lo² <= v²+4x²y² <= s_hi²,
0 < p_lo <= p_hi,
p_lo² <= (s_lo+v)/2,  p_hi² >= (s_hi+v)/2.
```

It then reconstructs L itself. The proposer uses interval leading principal determinants; the verifier uses interval LDL elimination and requires every pivot's lower bound to be strictly positive. Interval operations enclose the true values despite dependency between entries. The positive-matrix proof then excludes every positive completion with those moments.

The digest binds the certificate to normalized exact input values. Equivalent fractions are accepted; a certificate for changed moments or frequency is rejected. It is an input-binding mechanism, not an author signature.

Delta=0 returns unresolved without division. At k=0, the entire original polarization is fixed by the supplied moments. A separate certificate type checks its exact complex determinant directly: nonzero excludes the point, while zero returns unresolved. Failure to prove positivity before the precision limit also returns unresolved.

## Checks and limits

The retained validation checks square-root enclosures, agreement with the earlier exact-parameter formula, arbitrary rational frequency cases, common scales of10^-30 and10^30, forged/changed-input certificates, malformed moments and exact growing modes. The oscillating two-beam control has directions+x and+z, unit masses, k=6/5 and omega=3/5+4i/5; the original polarization is exactly singular and the new interface leaves it unresolved. The verifier also runs as a copied standalone file under optimized Python with an empty Python path.

No statistical accuracy rate is inferred from the finite case counts. The old frozen synthetic improvement over the species-block norm bound remains evidence of one exclusion improvement, not a general speedup or an astrophysical benchmark. Measurement-error enclosures, general wave-vector rotation in the API, physical mode data, independent proof review, novelty and operational usefulness remain open.
