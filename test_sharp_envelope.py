"""Independent checks of the sharp PURE-IMAGINARY equal-zero-flux envelope.

No production module, optimizer, angular grid or closure is imported. The
theorem concerns existence of one completion at omega=i*gamma and fixed k.
It makes no bound on roots with nonzero real part.
"""
import random
import unittest
from fractions import Fraction as Q
from independent_verifier import C,R,verify


def growth_squared(n,k,a,sector):
    if sector == 'longitudinal':
        return 2*n*k*a-k*(k+2*n)*a*a
    return n*k*a-k*(k+n)*a*a


def sharp_bound(n,k,sector):
    return n*n*k/(k+2*n) if sector == 'longitudinal' else n*n*k/(4*(k+n))


def maximizing_a(n,k,sector):
    return n/(k+2*n) if sector == 'longitudinal' else n/(2*(k+n))


def parameter_for_gamma(n,k,gamma,sector):
    """Smaller root of the attainer parabola, in one ordered quadratic field."""
    if sector == 'longitudinal':
        discriminant=n*n-(k+2*n)*gamma*gamma/k
        return R(n/(k+2*n),-1/(k+2*n),discriminant)
    discriminant=n*n-4*(k+n)*gamma*gamma/k
    return R(n/(2*(k+n)),-1/(2*(k+n)),discriminant)


def attaining_measures(n,a,sector):
    # Actual angular measures mu, not the auxiliary reweighted measure nu.
    mass=n/(1+a)
    plus=[(a,mass),(-Q(1),mass*a)]
    minus=[(-a,mass),(Q(1),mass*a)]
    return (plus,minus) if sector == 'longitudinal' else (minus,plus)


def direct_integrals(first,second,w,k):
    return tuple(sum((C(sign*mass*u**j)/(w-k*u)
                      for atoms,sign in ((first,1),(second,-1)) for u,mass in atoms),C())
                 for j in range(3))


def certificate_from_measures(first,second,w,k):
    transforms=[]
    witnesses=[]
    z=w/k
    for atoms in (first,second):
        transforms.append(sum((C(mass)/(z-u) for u,mass in atoms),C()))
        weights=[mass/((z.re-u)**2+z.im*z.im) for u,mass in atoms]
        t0=sum(weights,Q(0))
        t1=sum((weight*atom[0] for weight,atom in zip(weights,atoms)),Q(0))
        left,right=atoms[0][0],atoms[1][0]
        witnesses.append(dict(kind='two',p2_a=-left-right,p2_b=left*right,
                              nu_t0=t0,nu_t1=t1))
    return dict(status='possible',s1=transforms[0],s2=transforms[1],witnesses=witnesses)


class SharpImaginaryEnvelope(unittest.TestCase):
    def test_lens_support_inequality_exactly_matches_threshold(self):
        rng=random.Random(1112026)
        for _ in range(180):
            n=Q(rng.randint(1,9),rng.randint(1,7))
            k=Q(rng.randint(1,9),rng.randint(1,7))
            gamma=Q(rng.randint(1,12),rng.randint(1,12))
            y=gamma/k
            # 2(r-c)=N/(y sqrt(1+y²))-N/(1+y²).
            support=R(-n/(1+y*y),n/(y*(1+y*y)),1+y*y)
            for sector,multiplier in (('longitudinal',1),('transverse',2)):
                target_magnitude=multiplier*k/(1+y*y)
                self.assertEqual(target_magnitude <= support,
                                 gamma*gamma <= sharp_bound(n,k,sector))

    def test_every_tested_subthreshold_frequency_has_exact_attainer(self):
        rng=random.Random(2026110102)
        count=0
        for _ in range(60):
            n=Q(rng.randint(1,6),rng.randint(1,5))
            k=Q(rng.randint(1,6),rng.randint(1,5))
            gamma=Q(rng.randint(1,12),12)
            for sector in ('longitudinal','transverse'):
                data=dict(N1=n,F1=0,N2=n,F2=0,w=C(0,gamma),k=k,sector=sector)
                if gamma*gamma > sharp_bound(n,k,sector):
                    self.assertTrue(verify(data,dict(status='impossible',proof_kind='four_disks')))
                    continue
                a=parameter_for_gamma(n,k,gamma,sector)
                self.assertGreater(a,0)
                self.assertLess(a,1)
                self.assertEqual(growth_squared(n,k,a,sector),gamma*gamma)
                first,second=attaining_measures(n,a,sector)
                for atoms in (first,second):
                    self.assertEqual(sum((mass for _,mass in atoms),Q(0)),n)
                    self.assertEqual(sum((u*mass for u,mass in atoms),Q(0)),0)
                    self.assertTrue(all(-1<=u<=1 and mass>0 for u,mass in atoms))
                i0,i1,i2=direct_integrals(first,second,C(0,gamma),k)
                direct=((1+i0)*(-1+i2)-i1*i1 if sector=='longitudinal'
                        else -1+(i0-i2)/2)
                self.assertEqual(direct,C())
                self.assertTrue(verify(data,certificate_from_measures(first,second,C(0,gamma),k)))
                count+=1
        self.assertGreater(count,20)

    def test_parabola_global_maximum_on_attainer_family(self):
        for n,k in ((Q(1),Q(1)),(Q(7,3),Q(2,5)),(Q(1,7),Q(13,2))):
            for sector in ('longitudinal','transverse'):
                peak=maximizing_a(n,k,sector)
                bound=sharp_bound(n,k,sector)
                self.assertTrue(0<peak<1)
                self.assertEqual(growth_squared(n,k,peak,sector),bound)
                curvature=k*(k+2*n) if sector=='longitudinal' else k*(k+n)
                for a in (Q(0),Q(1),Q(1,11),Q(6,7),peak):
                    self.assertEqual(bound-growth_squared(n,k,a,sector),curvature*(a-peak)**2)

    def test_boundary_roots_and_direct_beam_polynomial(self):
        # Rational gamma boundary examples permit direct complex exact checks.
        for sector,k,gamma in (('longitudinal',Q(2,3),Q(1,2)),
                               ('transverse',Q(1,3),Q(1,4))):
            n=Q(1)
            a=maximizing_a(n,k,sector)
            self.assertEqual(a,Q(3,8))
            self.assertEqual(gamma*gamma,sharp_bound(n,k,sector))
            first,second=attaining_measures(n,a,sector)
            for w in (C(0,gamma),C(Q(3,7),Q(2,5))):
                i0,i1,i2=direct_integrals(first,second,w,k)
                direct=((1+i0)*(-1+i2)-i1*i1 if sector=='longitudinal'
                        else -1+(i0-i2)/2)
                factor=2 if sector=='longitudinal' else 1
                reduced=-1-factor*n*k*a*(1-a)/(w*w-k*k*a*a)
                self.assertEqual(direct,reduced)
            data=dict(N1=n,F1=0,N2=n,F2=0,w=C(0,gamma),k=k,sector=sector)
            self.assertTrue(verify(data,certificate_from_measures(first,second,C(0,gamma),k)))
            self.assertRaises(ValueError,verify,data,dict(status='impossible',proof_kind='four_disks'))


if __name__=='__main__':
    unittest.main(verbosity=2)
