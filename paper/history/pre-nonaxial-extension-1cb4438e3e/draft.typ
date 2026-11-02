#set document(title: "Exact frequency uncertainty from neutrino angular moments", author: "Research draft")
#set page(paper: "a4", margin: (x: 22mm, y: 20mm), numbering: "1")
#set text(font: "Libertinus Serif", size: 10pt)
#set par(justify: true, leading: .6em)
#set heading(numbering: "1.")
#set math.equation(numbering: "(1)")
#align(center)[
  #text(size: 20pt, weight: "bold")[Exact frequency uncertainty from neutrino angular moments]
  #v(4mm)
  #text(size: 11pt)[A four-disk diagnostic and sharp growth envelopes]
  #v(3mm)
  #text(size: 9pt)[Working research draft · 1 November 2026 · Not submitted or peer reviewed]
]
#v(4mm)
#block(inset: 8pt, fill: luma(96%))[
*Abstract.* Angular moment data do not specify a neutrino angular distribution. We study the inverse stability question: which unstable temporal frequencies can occur for some positive angular completion of given number density and flux? In the local, collisionless, axisymmetric fast limit, both longitudinal and transverse dispersion equations reduce to one complex Stieltjes value. Classical finite-interval moment theory turns each species' attainable values into two disks. A frequency is attainable exactly when four translated disks intersect. At most sixteen explicit candidates suffice, and exact algebraic witnesses or independent emptiness checks replace solver-status decisions. We derive sharp growth caps over every real frequency for equal masses and zero flux, with explicit attaining measures. Supplied second number moments also give sharp caps over all wave numbers, with essential stable edge cases. We extend the exact diagnostic to nonzero species flux and supplied second moments. For equal supplied second moments, we give the complete sharp finite-wave-number envelope: endpoint, cubic and plateau regimes, with positive attainers and a cube-root small-moment limit. A fixed published counter-flux test yields no added-moment rejection of the reproduced closure modes, so practical utility remains open. Positive smooth completions can grow despite sharing their first and second number angular moments with a stable isotropic pair. All 504 predeclared classifications and 90 boundary/scaling cases have exact independent checks. These results concern a restricted uncertainty class, not prediction of an actual supernova background. Weyl-disk theory and completion ambiguity are established precedents; priority and practical impact of the specialized envelopes remain to be assessed.
]

= Question and scope

Two-moment transport provides limited angular information. A closure chooses a completion; an uncertainty diagnostic asks what the supplied information permits. Completion ambiguity and crossing diagnostics are established topics @johns @abbar. A reduced zero-mode growth diagnostic already depends only on signed second moments @dasgupta. Quantum moment closures also support efficient stability calculations @froustey @qme. The present output is different: *some* compatible angular distribution has a specified frequency, or *none* does. A possible frequency does not imply that every completion is unstable, and different frequencies may require different completions.

Let $mu_+$ and $mu_-$ be positive measures on $u in [-1,1]$, uniform in azimuth, with known masses $N_+,N_-$ and axial number fluxes $M_+,M_-$. Their signed measure is $G=mu_+-mu_-$. Heavy-flavor cancellation is an explicit pointwise symmetry assumption, not inferred from equal moments. Absorb a common self-interaction scale into masses, frequencies and wave numbers. All moments here are number weighted; calling the second moment thermodynamic pressure would require additional spectral assumptions.

We use reduced real $k>0$ along the symmetry axis and $op("Im") omega>0$. Matter and background-current shifts are real and fixed by the supplied low moments. The model omits collisions, vacuum mixing, spatial variation and nonlinear saturation. Analytic angular integration, including its correct branch structure, is established @morinaga. No first three-dimensional solver or new generic moment theory is claimed.

= Exact reduction of the stability question

Define $g_0=N_+-N_-$, $g_1=M_+-M_-$ and
$ I_j=integral_(-1)^1 u^j/(omega-k u) dif G(u), quad J=I_0. $
The polarization tensor is $Pi=eta+integral v v^T/(omega-k u) dif G$, with the uniform azimuthal average understood, where $v=(1,bold(v))$, $abs(bold(v))=1$, and $eta=op("diag")(1,-1,-1,-1)$. Axial symmetry splits it into a longitudinal block and two equal transverse factors:
$ D_L=(1+I_0)(-1+I_2)-I_1^2, quad T=-1+(I_0-I_2)/2. $
The identities $k I_1=omega J-g_0$ and $k^2I_2=omega^2J-omega g_0-k g_1$ cancel the quadratic terms in $J$. Thus
$ k^2D_L=A_L J-B_L, quad 2k^2T=A_T J-B_T, $

#block(breakable: false)[where
$ A_L=omega^2-k^2+omega g_0-k g_1, quad B_L=k^2+omega g_0+k g_1+g_0^2, $
$ A_T=k^2-omega^2, quad B_T=2k^2-omega g_0-k g_1. $
]
With $z=omega/k$ and species transforms $S_s=integral dif mu_s/(z-u)$, the root equation becomes
$ S_+-S_-=k B_X/A_X, quad X in {L,T}. $
An upper-half-plane zero of $A_L$ cannot also zero $B_L$, since $A_L+B_L=(omega+g_0)^2$. It therefore gives no root. The transverse denominator vanishes only on the real axis. These exceptional cases are handled explicitly.

= From positive measures to four disks

Put $z=x+i y$, $y>0$, and consider one species of mass $N$ and flux $M$, with $N>=abs(M)$. Reweight its measure by $dif rho=dif mu/abs(z-u)^2$. For a proposed value $S=p+i q$, its first four moments satisfy
$ t_0=-q/y, quad t_1=x t_0-p, $
$ t_2=N+2x t_1-abs(z)^2t_0, quad t_3=M+2x t_2-abs(z)^2t_1. $
The degree-three Hausdorff conditions are
$ H_sigma=mat(t_0+sigma t_1,t_1+sigma t_2; t_1+sigma t_2,t_2+sigma t_3) >= 0, quad sigma in {-1,1}. $
These are classical moment conditions. Finite-interval Weyl loci are direct antecedents @di_dio. Heide et al. give the intersection of two matrix balls for nondegenerate fixed-point transform values @heide, Theorem 14.2. Its scalar, zero-order specialization supplies the nondegenerate mass/flux locus here, with the opposite transform sign. The generic two-disk construction is therefore a direct specialization of prior theory; the code also handles degenerate endpoint/zero-mass cases. Completing the square in each determinant gives the closed disks with
$ a_sigma=(x+sigma)^2+y^2, $
$ p_sigma=(N(2x+sigma)-M)/(2a_sigma), $
$ q_sigma=-((1+sigma x)(N+sigma M)+2N y^2)/(2a_sigma y), quad r_sigma^2=(N+sigma M)^2/(4y^2a_sigma). $
This is genuinely a disk condition, not just a determinant relaxation. The top-left entry is nonnegative throughout its disk; at its unique zero tangent, the off-diagonal vanishes and the remaining diagonal is $N+sigma M>=0$. Zero-radius cases correctly describe fixed endpoint atoms or zero mass.

*Finite diagnostic.* Translate the two minus-species disks by the target $k B_X/A_X$. Their intersection with the two plus-species disks is nonempty if and only if the frequency is attainable. A nonempty intersection of four closed disks contains either a disk center or a pairwise boundary-intersection point. There are four centers and at most twelve pair points. Coincident/concentric circles are skipped, and zero-radius centers remain candidates. If boundary arcs never switch inside the intersection, one full disk is contained in the others and its center works; otherwise an arc switch is a pair candidate. Tangential singleton intersections are also included.

Rational input gives rational centers and squared radii. Pair candidates lie in a real quadratic field. Exact sign checks decide membership, with no angular quadrature or numerical optimizer. Arithmetic-operation count is bounded independently of angular-grid size; bit complexity still depends on input size and separation from degeneracy. The implementation prefers feasible boundary candidates.

*Attainment and proof checking.* For feasible $S$, let $Delta=t_0t_2-t_1^2$. At rank one the measure is one atom. At rank two, the Gaussian quadrature polynomial $u^2+a u+b$ has
$ a=(-t_0t_3+t_1t_2)/Delta, quad b=(-t_2^2+t_1t_3)/Delta. $
Its roots and weights reconstruct $rho$; multiplying weights by $abs(z-u)^2$ reconstructs $mu$. Exact support, nonnegativity, moment and transform checks accompany the witness. Generic node formulas may introduce another square root beyond the coefficient field. A separate standard-library checker imports neither the production code nor its optimizer. It replays empty intersections and verifies witnesses; unresolved numerical proposals carry no authority.

= Sharp bounds over every real frequency

Assume $N_+=N_-=N>0$ and $M_+=M_-=0$. Every growing root $omega=Omega+i gamma$, with arbitrary real $Omega$, obeys
$ gamma_L^2 <= N^2 k/(k+2N), quad gamma_T^2 <= N^2 k/(4(k+N)). $ <envelope>
These caps are sharp over *all positive axial measures with those moments and every real frequency part*. Equality requires $Omega=0$. Every positive growth rate below either cap is attainable on that axis; this does not characterize the full two-dimensional frequency region.

*Proof of the global upper bounds.* Set $z=omega/k=x+i y$ and transform each species value to $U=(z^2-1)S-N z$. Direct substitution in the moment determinants gives
$ det H_sigma=(R_sigma^2-abs(U-C_sigma)^2)/abs(z-sigma)^2, $
$ C_sigma=-sigma N/2+i N(1-sigma x)/(2y), quad R_sigma^2=N^2((x-sigma)^2+y^2)/(4y^2). $
Positivity therefore puts both species values in the same two-disk intersection. The longitudinal root equation requires $U_+-U_-=k$, and the transverse equation requires $-2k$. In either case the lens contains a horizontal chord of length $t$, where $t=k$ or $2k$.

Write its endpoints as $m plus.minus h+i v$, $h=t/2$, and put $q=v-N/(2y)$. Both endpoints must lie in each disk, so
$ (abs(m+sigma N/2)+h)^2+(q+sigma N x/(2y))^2 <= N^2((x-sigma)^2+y^2)/(4y^2). $
Average the two sign inequalities. All $x$ terms cancel, leaving
$ ((abs(m+N/2)+h)^2+(abs(m-N/2)+h)^2)/2+q^2 <= N^2(1+y^2)/(4y^2). $
The left side is at least $(N/2+h)^2$, by the triangle inequality and the mean-square inequality. Hence $y^2<=N^2/[t(t+2N)]$, proving @envelope without assuming $Omega=0$. Equality forces $m=q=0$; the two individual disk slacks then have opposite signs proportional to $x$, so both are nonnegative only at $x=0$. No uniqueness statement about the distributions is made.

*Explicit attainers.* For $0<a<1$, take the actual species measures
$ mu_+=N/(1+a)(delta_a+a delta_(-1)), quad mu_-=N/(1+a)(delta_(-a)+a delta_1). $
Both have mass $N$ and flux zero. Their signed transform is
$ J=-2N k^3a(1-a)/((omega^2-k^2a^2)(omega^2-k^2)). $
The longitudinal root has $gamma^2=2N k a-k(k+2N)a^2$. Maximizing gives $a=N/(k+2N)$ and the first bound. Swapping species gives a transverse root with $gamma^2=N k a-k(k+N)a^2$, maximized at $a=N/(2(k+N))$. Varying $a$ below its optimizer realizes every positive sub-bound rate. For measures, equality is attained; under smooth-density restrictions, these atomic extrema are approached rather than asserted attained.

#figure(image("sharp-bands.svg", width: 100%), caption: [Sharp maximum growth over every real frequency at fixed wave number. Every shaded growth rate has a compatible completion at zero real frequency; the region is not a dispersion branch of one background.])

= Supplied second moments: sharp caps over all wave numbers

Retain equal masses $N$ and zero species fluxes, and now fix each second number moment $Q_+,Q_- in [0,N]$. Define positive weighted measures $dif lambda_s=(1-u^2)dif mu_s$ with masses $R_s=N-Q_s$. Polynomial division, using signed mass and flux zero, gives
$ (omega^2-k^2)J/k^2=-F, quad F=integral dif(lambda_+-lambda_-)/(omega-k u). $
Thus $D_L=-1-F$ and $T=-1+F/2$. For any positive measure of mass $R$, its transform at $omega=Omega+i gamma$ lies in the disk centered at $-i R/(2gamma)$ with radius $R/(2gamma)$. This follows by averaging the circle traced by $1/(omega-k u)$ over real $u$, or by Cauchy-Schwarz. The difference of the two disks has center $-i(R_+-R_-)/(2gamma)$ and radius $(R_++R_-)/(2gamma)$. Requiring the real targets $F=-1$ or $F=2$ proves
$ gamma_L^2 <= R_+R_-, quad gamma_T^2 <= R_+R_-/4. $ <second-cap>
The angular weight already appears in the axial dispersion relation @yi. The elementary transform-disk inequality is classical, not an invention. Here it supplies a sharp *moment-constrained* extremum, as the following construction proves.

Suppose $Q_+,Q_-<N$ and $Q_++Q_->0$. Put $a=Q_+/N$, $b=Q_-/N$, and take
$ mu_+=N/(1+a)(delta_a+a delta_(-1)), quad mu_-=N/(1+b)(delta_(-b)+b delta_1). $
These measures have the prescribed mass, flux and separate second moments. They attain the longitudinal cap at
$ k_L=N(2-a-b)/(a+b), quad Omega_L=N(a-b)/(a+b), quad gamma_L^2=N^2(1-a)(1-b). $
Reflect both measures under $u mapsto -u$ to attain the transverse cap at $k_T=k_L/2$, $Omega_T=-Omega_L/2$, $gamma_T^2=gamma_L^2/4$. Direct substitution in the original dispersion factors verifies both roots. The executable construction returns exact rational squared caps, wave numbers, real frequencies and atom weights, including unequal second moments.

*Essential edge cases and quantifiers.* If either $Q_s=N$, one weighted measure vanishes and no upper-half-plane root exists. If $Q_+=Q_-=0$, both species are forced to $N delta_0$ and are stable, despite the loose disk bound. In every other case @second-cap is sharp when optimizing over *all* positive $k$ and real $Omega$. It does not supply the sharp bound for each prescribed $k$ with second moments imposed. As $a+b$ tends to zero, the maximizing wave number diverges, whereas exactly $a=b=0$ is stable. A physical wave-number cutoff, smoothness or occupation constraint can change the extremum.

For isotropic second tensors, $Q_+=Q_-=N/3$, the exact worst-case rates are $2N/3$ longitudinal at $k=2N$ and $N/3$ transverse at $k=N$. A stable isotropic pair has the same supplied moments. These are information limits, not predictions of flavor conversion in a particular background.

= Complete sharp envelope in a finite wave-number band

Let each species have moments $(N,0,N a)$, $0<a<1$, and put $r=1-a$. The reduced wave number is restricted to $0<k<=K$; no physical ultraviolet cutoff is inferred. Write $c=1$ longitudinal and $c=2$ transverse. We now determine the maximum over all compatible positive measures and *every real frequency part*, for every $a,K$. The reduced zero-mode diagnostic @dasgupta is stable in this class because the signed moments through order two vanish. The result below concerns nonzero reduced wave number.

*Both moment disks matter.* Set $z=x+i y$, $y>0$, and $U=(z^2-1)S-N z$. The classical species locus lies in the intersection of disks with
$ C_0=i N r/(2y), quad R_0=N r/(2y), $
$ C_1=N/(a+x^2+y^2) [-x+i ((2-a)y^2+a(1-x^2))/(2y)], $
$ R_1^2=N^2a^2 abs(z^2-1)^2/(4y^2(a+x^2+y^2)^2). $
The original factors require $U_+-U_-=k$ or $-2k$, so two points of the same lens must have horizontal separation $c k$. An enclosing Gram disk alone can give a loose cap.

*The all-real-frequency step.* Let $L^2=abs(C_1-C_0)^2$. For every $beta in [0,1]$, the intersection is enclosed by a disk centered at $(1-beta)C_0+beta C_1$ with squared radius
$ R_beta^2=(1-beta)R_0^2+beta R_1^2-beta(1-beta)L^2 $
$ R_beta^2=(1-beta)R_0^2+beta^2R_1^2+beta(1-beta)(R_1^2-L^2). $
Both variable terms are bounded above by their $x=0$ values. For $t=y^2$, $d=x^2$, the squared-radius comparison for $R_1$ reduces to
$ (1+t)^2(a+t+d)^2-(a+t)^2[(t+d-1)^2+4t] $
$ =d [2(a+t)(1+a+(3-a)t)+(1-a)(1+a+2t)d] >= 0. $
The other term simplifies to
$ R_1^2-L^2=N^2(1-a^2)(a-x^2-y^2)/(4y^2(a+x^2+y^2)), $
which decreases with $x^2$. Nonnegative coefficients therefore give $R_beta^2(x,y)<=R_beta^2(0,y)$ for every fixed $beta$. This proves an upper bound at arbitrary real frequency, without assuming that a smaller single radius controls the lens.

At $x=0$, the unconstrained minimizer is
$ beta_0=r(y^2-a^2)(a+y^2)/(a^2+y^2)^2. $
Clamp it to $[0,1]$. If $beta_0<=0$, disk 0's horizontal extrema lie in disk 1; if $beta_0>=1$, the reverse holds. Otherwise the circle intersections have height $N r y/(a^2+y^2)$ and half-separation $N a r/(a^2+y^2)$, exactly the optimized enclosing radius. Consequently the sharp imaginary-axis chord width, and an upper bound at every real part, is
$ W(y)=cases(
 N r/y & "if" y<=a,
 (N a(1+y^2))/(y(a+y^2)) & "if" y+a/y<=r,
 (2N a r)/(a^2+y^2) & "otherwise".
) $ <width>
The first two conditions are disjoint. Each active piece decreases strictly with $y$, and they agree at their boundaries. Explicit attainers below establish sharpness without assuming an enclosing disk is wholly attainable.

#block(breakable: false)[
*Complete formula.* Put $k_*=N r/(c a)$ and $E_c(k)=sqrt(2N a r k/c-a^2k^2)$. Let $Gamma_c(k)$ be the unique positive root of
$ Gamma^3-(N a/c)Gamma^2+a k^2 Gamma-(N a/c)k^2=0. $ <cubic>
]

When $a<=3-2sqrt(2)$, define
$ y_(plus.minus)=(r plus.minus sqrt(1-6a+a^2))/2, quad k_c(y)=(N a(1+y^2))/(c y(a+y^2)), $
and $k_("lo")=k_c(y_+)$, $k_("hi")=k_c(y_-)$. The exact fixed-$k$ maximum is
$ F_c(k)=cases(
 N r/c & "if" k>=k_*,
 Gamma_c(k) & "if" k_("lo")<=k<=k_("hi"),
 E_c(k) & "otherwise".
) $ <finite-envelope>
The cubic band is absent when $a>3-2sqrt(2)$ and is a single point at equality. Otherwise $0<k_("lo")<=k_("hi")<k_*$. The endpoint piece increases to $k_*$, the cubic increases by implicit differentiation, and the transitions agree. Thus $F_c$ is nondecreasing, and the sharp maximum over $0<k<=K$ is $F_c(K)$. Reflection gives the same answer for $abs(k)<=K$. At $k=0$, $a=0$ or $a=1$, the maximum is zero. The sector relation is $F_2(k)=F_1(2k)/2$.

*Positive attainers.* In the endpoint regime use $mu_+=N/(1+a)(delta_a+a delta_(-1))$ and its reflection. In the cubic regime set $y=Gamma_c(k)/k$, $p q=a$ and $q-p=y+a/y$; the band condition ensures $p,q<=1$. Then
$ mu_+=N q/(p+q)delta_p+N p/(p+q)delta_(-q), quad mu_-=op("reflection")(mu_+). $
In the plateau regime set $y=N r/(c k)<=a$, $w=N r/(1-y^2)$, and take
$ mu_+=w delta_y+N(a-y)/(2(1-y))delta_1+N(a+y)/(2(1+y))delta_(-1), $
with its reflection. All weights are nonnegative, and every pair has the required mass, zero flux and second moment. The weighted plateau measures are $N r delta_(plus.minus y)$, so their $U$ difference is $N r/y=c k$. Each construction uses the shown species order for longitudinal; swap species for transverse. Direct substitution in the original factors verifies the roots. The plateau generally needs three atoms in this construction; no smooth exact attainment is claimed.

The endpoint family is genuinely nonoptimal in the cubic band. For $N=1$, $a=3/50$, $k=109/75$, the pair $mu_+=6/7 delta_(1/10)+1/7 delta_(-3/5)$ and its reflection attains $gamma=109/250$, above the endpoint family's squared rate $244269/1562500$. Taking the minimum of $E_c$, $Gamma_c$ and the plateau would therefore be incorrect.

The evaluator uses exact rational branch decisions. For $k<k_*$ put $t=2N a r/(c k)-a^2$. The cubic is active precisely when $(t+a)^2<r^2t$; equality may use the endpoint expression. Rational bisection encloses the positive root, with no angular or frequency grid. Its arithmetic cost grows with requested precision and input bit size; no measured solver speed advantage is asserted.

#figure(image("finite-band.svg", width: 100%), placement: top, caption: [Complete finite-band maxima for equal second moments. The prior comparator combines the mass/flux cap, Gram cubic and supplied-second-moment all-wave-number cap. The fixed endpoint family misses the interior cubic regime. Each point optimizes over a potentially different angular completion.])

*Increment over existing bounds.* For $N=1$, $Q=1/3$, combine the earlier mass/flux cap $sqrt(k/(k+2))$, second-moment cap $2/3$ and Gram cubic. The complete cap tightens their minimum by at most $1-sqrt(8/9)$, approximately $5.72%$. Indeed, below $k=8/5$ the squared ratio of the first two caps' minimum to $F_1(k)$ is $9/((k+2)(4-k))<=9/8$; between $8/5$ and $2$ it is $4/(k(4-k))<=25/24$, and beyond $2$ both plateau. Adding the cubic only reduces that ratio. Sharp attainment is mathematically stronger, but this isotropic comparison establishes only a modest numerical tightening.

At fixed finite $K,N$ as $Q=N a$ tends to zero,
$ F_c(K)/((Q K^2/c)^(1/3)) -> 1. $
Rescaling @cubic by $Gamma=(N a K^2/c)^(1/3)h$ leaves $h^3-1$ plus coefficients tending to zero. The positive root tends to one in these units, while $y=Gamma/K$ and $a/y$ tend to zero, eventually placing $K$ in the cubic band. Hence the maximum tends continuously to the stable $Q=0$ corner and the transverse/longitudinal ratio tends to $2^(-1/3)$. Neither a physical cutoff nor general unequal-second-moment finite-band extrema follow.

#block(breakable: false)[
= Exact frequency compatibility with general supplied second moments

The diagnostic also admits nonzero species flux and independently supplied $Q_s$. For one positive species, assume $N>0$, $Delta=N Q-M^2>=0$ and $Q<=N$. Set $d=Q-2M x+N abs(z)^2>0$. Its exact transform locus is the intersection of
$ abs(S-[N(N x-M)/d-i(2N^2y+Delta/y)/(2d)])^2 <= Delta^2/(4y^2d^2), $
and
$ abs((z^2-1)S-N z-M-i(N-Q)/(2y))^2 <= (N-Q)^2/(4y^2). $
]

These are classical moment-theory specializations. The first follows by the Schur complement of the Gram matrix for $1,u,1/(z-u)$; the second uses the positive measure $(1-u^2)dif mu$. Neither is a determinant-only relaxation.

For sufficiency reweight as before and add $t_4=Q+2x t_3-abs(z)^2t_2$. The degree-four Hausdorff matrices are
$ H=mat(t_0,t_1,t_2;t_1,t_2,t_3;t_2,t_3,t_4)>=0, quad
 L=mat(t_0-t_2,t_1-t_3;t_1-t_3,t_2-t_4)>=0. $
Invertible polynomial changes of basis make the first disk equivalent to $H>=0$ when $Delta>0$, and the second equivalent to $L>=0$. At $Delta=0$ positivity forces $N delta_(M/N)$; at $Q=N$ it forces the unique endpoint pair. Zero mass separately forces $M=Q=S=0$. A constructive Gauss-Radau argument supplies at most three atoms from $H,L$; singular $H$ cases reduce to one or two. Three are sometimes necessary, so the M1 two-atom guarantee does not extend unchanged.

Translation by the unchanged dispersion target again gives four disks. The code accepts $omega=w_r+i sqrt(g)$ with rational $w_r$, rational $g>0$ and positive rational $k$. Coordinates $S=p+i sqrt(g)q$ give the rational metric $op("diag")(1,g)$. Pair intersections still need only one real quadratic coefficient field. A separate standard-library checker evaluates all principal minors of $H,L$ and the original dispersion factor for positive proofs. For negative proofs it independently reconstructs ellipses from determinant evaluations. It imports no production formulas or optimizer. Returned M2 proof objects certify existence through these moment matrices; a full numerical atom listing is not required or emitted by this interface.

= Smooth ambiguity with identical low moments

The attaining family has identical per-species second scalar moment $N a$ and isotropic transverse second moments $N(1-a)/2$. At $N=k=1$, the longitudinal optimizer has $a=1/3$, so both full number second-moment tensors equal $I/3$. Two uniform species have exactly the same moments and $G=0$, hence are stable. Indeed, zero signed mass and flux make the shifted zero-mode spatial problem real symmetric; a finite-$k$ result here is not continuity from an already unstable zero mode.

Convolve the two attaining angular measures on the sphere with a normalized rotationally invariant von Mises-Fisher kernel of concentration $kappa$. Its action scales each spherical-harmonic degree without mixing degrees. The vanishing dipole and quadrupole coefficients remain exactly zero. The resulting strictly positive smooth densities therefore retain mass one, zero flux and second tensor $I/3$. Their axial density with respect to $dif u$ is
$ f_+(u)=kappa/(2sinh kappa) [3/4 e^(kappa u/3) op("I")_0(kappa sqrt(8(1-u^2)/9))+1/4 e^(-kappa u)], quad f_-(u)=f_+(-u). $
The validated integration uses an entire hypergeometric representation of the Bessel factor, so no square-root branch is hidden in the analytic callback. The real longitudinal dispersion is integrated with Arb balls. Strict opposite signs at rational positive-$gamma$ endpoints certify a growing root by continuity, then bisection gives width below $10^(-4)$.

The fixed concentrations $2,4,8,16,32,64$ were all reported. Four brackets certify rates near $0.26694$, $0.41343$, $0.49356$ and $0.53520$ for $kappa=8,16,32,64$. The first two concentrations produced no bracket in the specified search; that is not a proof of stability. Independent SciPy integration reproduced all four roots inside the saved brackets and complex determinant residuals below $4 times 10^(-16)$.

#figure(image("smooth-completions.svg", width: 100%), caption: [Left: positive smooth electron-species completions; antineutrinos are reflected. Every curve has the same first and second number moments as the stable uniform pair. Right: validated growing-root brackets approach the sharp measure envelope. Concentrations 2 and 4 had no certified bracket and are not plotted as zero growth.])

= Reproducible checks and what they establish

The predeclared screen uses three moment classes, two nonzero wave numbers, seven real frequencies, six positive imaginary frequencies and both sectors: 504 existential decisions. The classes are equal masses/zero flux, the normalized axial projection of the published fiducial counter-flux moments @froustey, and a single-positive-species control. This is not reproduction of a neutron-star-merger snapshot or its actual instability rates.

#block(breakable: false)[
#table(columns: (1.7fr,1fr,1fr), align: (left,right,right), inset: 5pt,
  table.header([Moment class],[Possible],[Impossible]),
  [Equal masses, zero flux],[34],[134],
  [Fiducial axial moments],[38],[130],
  [Single positive species],[0],[168],
  [*Total*],[*72*],[*432*],
)
]

All 504 decisions have independent exact checks. The zero-flux class contains 22 distinct attainable $(k,omega)$ pairs although it also contains the stable uniform completion. A generic numerical moment-LMI proposer with exact proof checks agrees on this screen. Angular-grid LPs at 32, 128 and 512 nodes agree at these points as well; no sampled accuracy advantage over those grids is claimed.

An additional 90 designed tests cover exact envelope boundaries, perturbations down to $10^(-12)$ and common scales from $10^(-50)$ to $10^50$. Every exact decision matches the proved threshold. At the two unit-scale boundaries the attaining interior cosine is $3/8$, absent from all three uniform grids. Exact Farkas checks certify each grid problem infeasible, while a continuous positive-measure witness exists. This is an intentional discretization stress, not a random error rate or a numerical LP failure. The numerical proposer left two unit-scale stress cases unresolved; those flags are not negative answers.

Independent tests include rational atomic moments, direct sector identities, exact radical ordering, degenerate circles, malformed certificate rejection and boundary attainers. The initial preregistration was amended before the full screen to accept algebraic disk-emptiness proofs in place of rational duals. Four development smoke points overlapped its planned grid and are disclosed. The records are task-local, not externally timestamp-registered. The first screen attempt stopped on a reference-baseline type conversion; its corrected complete run and failure log are retained.

A stronger continuous comparison uses the original Hausdorff matrices as four second-order-cone constraints, solved by Clarabel 0.11.1. It includes all 594 development cases and 128 fresh, seeded rational cases, with three interleaved timing repetitions. The exact method certifies all 722 decisions. The conic method extracts independently accepted rational certificates for 682, leaving 40 designed stress cases unresolved: 20 attainable and 20 unattainable. It certifies every original-screen and fresh case, with no certified disagreement. An ordinary floating circle method agrees on those 632 non-stress points but labels ten near-boundary impossible stress points possible under its fixed tolerance. These are numerical approximations, not contradictory exact certificates.

Generation medians are about 2.37 ms for exact geometry and 3.58 ms for the conic proposer. Independent checking is substantially cheaper for the conic certificates. On jointly certified cases, the median paired ratio of conic total generation-plus-checking cost to geometric cost is 0.737. Thus *the exact method has no demonstrated end-to-end speed superiority*. Its demonstrated benefit is resolving all specified rational boundaries and returning exact decisions, with larger checking cost. The solver's internal median of about 0.125 ms is not compared with another method's total time. Results are limited to these synthetic instances, implementation and machine. An avoidable small-scale normalization weakness in an initial conic run was corrected before aggregate inspection; both runs and the amendment are preserved.

The global-bound extension has independent formal polynomial and exact-arithmetic checks. The second-moment construction passes 1,170 exact direct substitutions in the original dispersion factors, with positive weights and all prescribed moments checked. These checks corroborate the written proofs; no finite test set proves a universal bound.

*Published counter-flux test and a failed utility gate.* The controlled linear angular input of Richers et al. @pic has normalized species moments $(1,plus.minus 1/3,1/3)$. These are exact number moments of its prescribed continuum distributions, not a merger snapshot or closure-derived pressure. The zero-species-flux sharp theorems do not apply to this input; the general M2 diagnostic does.

Before examining classifier results, we fixed physical axial $K=j/8$, $j=-8,...,8$. The current shift is reduced $k=K-2/3$; negative reduced $k$ is handled by reflecting the axis and species fluxes. Reproducing the polynomial Minerbo/FT quantum closure in @froustey gives $chi=17/45$ and exact longitudinal squared growth, using the signed reduced $k$ before axis reflection,
$ gamma^2=(k+2/3)(2/3-17k/45). $
Its independent four-by-four characteristic polynomial is checked exactly at all 17 points. Both transverse closure blocks are nilpotent and have no growing eigenvalue; this does not exclude actual angular transverse instabilities. An independent closed-form angular integral and Arb arithmetic certify 13 positive root brackets for the original linear input. The exact-entropy reconstruction is a third, different model: it has second moment about $0.37959$, whereas the input has $1/3$ and the polynomial closure uses $17/45$. Its roots remain numerical references.

The declared utility criterion required at least one reproduced growing closure mode to be impossible with supplied M2 while possible with M0/M1. All eight growing longitudinal closure modes remain exactly compatible with BOTH information classes. The gate therefore *fails*. All 26 designated rational points near the linear/entropy reference roots and all 26 linear-bracket endpoints also admit exact M2 certificates. Those point checks are not identities of rounded roots or a proof about a whole bracket; existence of the 13 actual linear-input roots rests on the separate Arb signs. No band was widened after seeing the negative utility result. Compatibility cannot validate the closure's actual growth prediction, and these results do not establish practical astrophysical benefit.

The general M2 checker passed nine independent tests, including 1,500 locus identities, 250 exact Radau completions, 250 independent ellipse comparisons and 200 exact frequency decisions. The full finite-band extension adds 1,000 exact enclosing-disk identities, 72 cubic, 96 endpoint and 108 plateau original-factor attainers, 36 exact transition cases and 150 random rational enclosures, including common scales from $10^(-50)$ to $10^50$. A fresh reviewer separately checks the universal proof and exact branch selector. Earlier partial-bound checks remain preserved. These finite checks corroborate, rather than replace, the written proofs and their independent review.

= Novelty assessment and remaining work

The underlying finite-interval transform loci are classical @di_dio @heide. In particular, the two-ball theorem specializes directly to the scalar mass/flux problem; combining two species does not make that moment theory new. Closure/completion ambiguity and low-moment crossing disagreement are also explicit predecessors @johns. Consequently, the defensible candidate contribution is the specialized exact frequency diagnostic with attaining completions, the sharp constrained growth caps over every real frequency, the second-moment extrema over all wave numbers, the complete equal-second-moment finite-band envelope, and the explicit smooth matched-second-moment instability example. Exact beam growth maxima @chakraborty and geometric-mean crossing-depth estimates @richers are also established; the square-root-product form alone is not new. The present theorem instead bounds all compatible completions from supplied species moments and proves attainment. Focused searches found no identical constrained sharp theorem, but absence from a search is not proof of priority.

The construction has survived independent mathematical and physics review within its model; those reviewers were fresh-context model agents, not external human referees. Practical importance remains unproved; the fixed counter-flux closure-exclusion test failed its declared utility gate. Needed next steps include expert priority review, realistic number-moment provenance, more physical classes and independent hardware reproduction. A justified uncertainty class must incorporate any available energy/occupancy or angular regularity information. No additional scan of this same fixed closure test is justified as a way to turn its failed gate into success. Arbitrary wave-vector directions, full three-dimensional angular uncertainty, collisions, nonlinear flavor evolution and actual astrophysical consequences remain outside this result. The two equal-mass/zero-flux extremum theorems have their stated global quantifiers; the finite classifier map does not give a general all-frequency or all-wave-number census.

Any positive angular quadrature preserving the supplied moments is itself an admissible completion. Every exact mode of that quadrature must therefore remain compatible, even when it is absent from the intended continuum background. The diagnostic cannot by itself filter such discretization-induced modes; analytic angular integration uses additional distributional information @morinaga. Nor does a worst-case finite-band envelope guarantee the fraction of growth captured for a fixed background, or establish a stable numerical time step.

This remains an unsubmitted research draft without external peer review.

#block[
#set text(size: 9pt)
#set par(leading: .45em)
#bibliography("references.bib", style: "ieee")
]
