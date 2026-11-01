"""Few exact illustrative sensitivity checks; no conic runtime experiment."""
from fractions import Fraction as F
import unittest
from independent_verifier import C, verify, support_matrices
from second_moment_bounds import second_moment_bounds


def cap(n, k, transverse=False):
    return n*n*k/(4*(k+n) if transverse else k+2*n)


def det(a,b,d):
    return a*d-b*b


def quadrature(t0,t1,t2,t3):
    # Solve the two orthogonality equations for u²+a*u+b.
    denominator=t0*t2-t1*t1
    return dict(kind='two',p2_a=(t1*t2-t0*t3)/denominator,
                p2_b=(t1*t3-t2*t2)/denominator,nu_t0=t0,nu_t1=t1)


class MomentPrecision(unittest.TestCase):
    def test_four_declared_mass_shifts(self):
        for e in (F(1,100), F(1,10**4), F(1,10**8), F(1,10**12)):
            n=1+e
            t0,t1,t2,t3=F(272,225)*n,F(16,75),F(8,25)*n,-F(3,25)
            determinants=[]
            for sigma in (1,-1):
                aa,bb,dd=t0+sigma*t1,t1+sigma*t2,t2+sigma*t3
                self.assertGreater(aa,0)
                self.assertGreater(dd,0)
                determinants.append(det(aa,bb,dd))
            self.assertEqual(determinants,[F(16,225)*(n-1)*(4*n+1),
                                           F(16,225)*(n+1)*(4*n-1)])
            self.assertTrue(all(d>0 for d in determinants))
            self.assertEqual(cap(n,F(2,3))-F(1,4),e*(5+4*e)/(4*(4+3*e)))
            self.assertEqual(cap(n,F(1,3),True)-F(1,16),e*(5+4*e)/(16*(4+3*e)))
            z=C(0,F(3,4))
            sp,sm=C(-F(16,75),-F(68,75)*n),C(F(16,75),-F(68,75)*n)
            for k,transverse in ((F(2,3),False),(F(1,3),True)):
                i0=((sm-sp) if transverse else (sp-sm))/k
                i1,i2=z*i0,z*z*i0
                original=(-1+(i0-i2)/2 if transverse else (1+i0)*(-1+i2)-i1*i1)
                self.assertEqual(original,C())
                first,second=(sm,sp) if transverse else (sp,sm)
                plus=quadrature(t0,t1,t2,t3)
                minus=quadrature(t0,-t1,t2,-t3)
                witnesses=[minus,plus] if transverse else [plus,minus]
                data=dict(N1=n,F1=0,N2=n,F2=0,w=z*k,k=k,
                          sector='transverse' if transverse else 'longitudinal')
                self.assertTrue(verify(data,dict(status='possible',s1=first,s2=second,
                                                witnesses=witnesses)))
                for species_s in (first,second):
                    for aa,bb,dd in support_matrices(n,0,z,species_s):
                        self.assertGreater(aa,0)
                        self.assertGreater(dd,0)
                        self.assertGreater(aa*dd-bb*bb,0)
                lower_data={**data,'N1':1-e,'N2':1-e}
                self.assertTrue(verify(lower_data,dict(status='impossible',proof_kind='four_disks')))
                gamma0=F(1,4) if transverse else F(1,2)
                gamma=gamma0*(1+e/2)
                self.assertGreater(gamma*gamma,cap(F(1),k,transverse))
                self.assertLess(gamma*gamma,cap(n,k,transverse))

    def test_one_interior_interval_and_degenerate_corner(self):
        a,delta=F(1,3),F(1,100)
        l=a-delta
        bound=second_moment_bounds(1,l,l)
        self.assertEqual(bound.longitudinal.gamma_squared,(1-l)**2)
        self.assertEqual(bound.transverse.gamma_squared,(1-l)**2/4)
        self.assertEqual(bound.longitudinal.attainer.k,(1-l)/l)
        self.assertEqual((1-l)-(1-a),delta)
        self.assertEqual(second_moment_bounds(1,0,0).global_gamma_squared,0)
        self.assertEqual(second_moment_bounds(1,1,1).global_gamma_squared,0)
        # Positive values approach the all-k supremum; the zero pair is stable.
        q=F(1,10**8)
        small=second_moment_bounds(1,q,q)
        self.assertEqual(small.longitudinal.gamma_squared,(1-q)**2)
        self.assertEqual(small.longitudinal.attainer.k,10**8-1)


if __name__=='__main__':
    unittest.main(verbosity=2)
