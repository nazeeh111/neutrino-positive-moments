"""Independent exact checks for the global equal-M0/zero-M1 axial growth cap.

Standard library plus the separate certificate checker only. No production
geometry, numerical optimizer, or angular grid is imported. This checks
algebraic identities used in the accompanying proof and independently
rejects selected off-axis frequencies above (and on) the sharp cap.
The finite checks corroborate the written proof, rather than replace it.
"""
from fractions import Fraction as Q
import json
import random
import unittest
from pathlib import Path
from independent_verifier import C, support_matrices, verify


def det(m):
    return m[0]*m[2]-m[1]*m[1]


def transformed_circle_slack(n,x,y,u,sign):
    cx=-sign*n/2
    cy=n*(1-sign*x)/(2*y)
    r2=n*n*((x-sign)**2+y*y)/(4*y*y)
    return r2-(u.re-cx)**2-(u.im-cy)**2


def endpoint_slack(n,x,y,m,q,h,sign):
    """Least slack of the two horizontal-chord endpoints in disk sign."""
    cx=-sign*n/2
    cy=n*(1-sign*x)/(2*y)
    r2=n*n*((x-sign)**2+y*y)/(4*y*y)
    # Endpoint midpoint has imaginary part q+n/(2y).
    return r2-(abs(m-cx)+h)**2-(q+n/(2*y)-cy)**2


class GlobalCapChecks(unittest.TestCase):
    def test_transformed_disk_identity_from_moment_gram_determinants(self):
        rng=random.Random(20261101091)
        for _ in range(600):
            n=Q(rng.randint(1,20),rng.randint(1,10))
            x=Q(rng.randint(-25,25),rng.randint(1,12))
            y=Q(rng.randint(1,20),rng.randint(1,12))
            u=C(Q(rng.randint(-20,20),rng.randint(1,12)),
                Q(rng.randint(-20,20),rng.randint(1,12)))
            z=C(x,y)
            s=(u+n*z)/(z*z-1)
            for sign,gram in zip((1,-1),support_matrices(n,0,z,s)):
                self.assertEqual(det(gram),
                                 transformed_circle_slack(n,x,y,u,sign)/((x-sign)**2+y*y))

    def test_chord_endpoint_sum_cancels_real_frequency(self):
        rng=random.Random(20261101092)
        for _ in range(1000):
            n=Q(rng.randint(1,20),rng.randint(1,10))
            x=Q(rng.randint(-30,30),rng.randint(1,12))
            y=Q(rng.randint(1,20),rng.randint(1,12))
            m=Q(rng.randint(-20,20),rng.randint(1,12))
            q=Q(rng.randint(-20,20),rng.randint(1,12))
            h=Q(rng.randint(0,20),rng.randint(1,12))
            average=sum((endpoint_slack(n,x,y,m,q,h,sign) for sign in (1,-1)),Q(0))/2
            a,b=abs(m+n/2)+h,abs(m-n/2)+h
            simple=n*n*(1+y*y)/(4*y*y)-(a*a+b*b)/2-q*q
            self.assertEqual(average,simple)
            self.assertGreaterEqual((a*a+b*b)/2,(n/2+h)**2)

    def test_bound_equality_forces_zero_real_frequency(self):
        # Rational boundary examples, both sectors. Equality in the sum
        # requires midpoint m=q=0; the residual two slacks are opposite.
        for n,y,t in ((Q(1),Q(3,4),Q(2,3)),(Q(3),Q(3,4),Q(2)),
                      (Q(2),Q(4,3),Q(1,2))):
            # The last tuple is also sharp: sqrt(1+16/9)/(4/3)=5/4.
            h=t/2
            self.assertEqual((n/2+h)**2,n*n*(1+y*y)/(4*y*y))
            for x in (Q(-100),Q(-1),Q(-1,10000),Q(0),Q(1,10000),Q(1),Q(100)):
                slacks=[endpoint_slack(n,x,y,Q(0),Q(0),h,s) for s in (1,-1)]
                self.assertEqual(slacks[0],-n*n*x/(2*y*y))
                self.assertEqual(slacks[1],n*n*x/(2*y*y))
                self.assertEqual(min(slacks)>=0,x==0)

    def test_offaxis_upper_bound_with_independent_exact_emptiness(self):
        count=0
        # These are separate algebraic examples, not a general theorem by sampling.
        for sector,k,gamma in (('longitudinal',Q(2,3),Q(1,2)),
                               ('transverse',Q(1,3),Q(1,4))):
            for re in (Q(-10),Q(-1),Q(-1,10),Q(-1,10**5),Q(1,10**5),Q(1,10),Q(1),Q(10)):
                for scale in (Q(1),Q(2),Q(10**6),Q(1,10**6)):
                    for bump in (Q(0),Q(1,10**8),Q(1,100)):
                        data=dict(N1=scale,F1=0,N2=scale,F2=0,
                                  w=C(re*scale,(gamma+bump)*scale),k=k*scale,sector=sector)
                        self.assertTrue(verify(data,dict(status='impossible',proof_kind='four_disks')))
                        count+=1
        self.assertEqual(count,192)


if __name__=='__main__':
    suite=unittest.defaultTestLoader.loadTestsFromTestCase(GlobalCapChecks)
    result=unittest.TextTestRunner(verbosity=2).run(suite)
    if result.wasSuccessful():
        output=dict(passed=True,tests=result.testsRun,
                    determinant_identities=1200,chord_identities=1000,
                    exact_offaxis_emptiness_checks=192,
                    scope='Written global-cap proof corroboration; finite checks are not the proof',
                    production_or_numeric_imports=False)
        Path('results/offaxis-independent.json').write_text(json.dumps(output,indent=2)+'\n')
    raise SystemExit(not result.wasSuccessful())
