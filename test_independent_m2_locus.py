"""Independent exact checks; standard library only; no production imports."""
from fractions import Fraction as Q
import random
import unittest

from independent_verifier import InvalidCertificate, R
import independent_m2_verifier as v


def transform(atoms, wr, g, k):
    z = v.CQ(wr/k,1/k,g)
    return sum((v.CQ(weight,0,g)/(z-u) for u,weight in atoms),v.CQ(0,0,g))


def species(atoms):
    return tuple(sum((weight*u**j for u,weight in atoms),Q(0)) for j in range(3))


def radau(t):
    """Implement the proof's construction, independently of any frequency code."""
    h = [[t[i+j] for j in range(3)] for i in range(3)]
    variance = t[0]*t[2]-t[1]**2
    if variance == 0:
        return [(t[1]/t[0],t[0])]
    full_rank = v.determinant(h) > 0
    a = [t[j]-t[j+1] for j in range(4)] if full_rank else list(t[:4])
    d = a[0]*a[2]-a[1]*a[1]
    b = (-a[2]*a[2]+a[1]*a[3])/d
    c = (-a[0]*a[3]+a[1]*a[2])/d
    radical = R(0,1,c*c-4*b)
    nodes = [(-c-radical)/2,(-c+radical)/2]
    weights = [(nodes[1]*a[0]-a[1])/(nodes[1]-nodes[0]),
               (a[1]-nodes[0]*a[0])/(nodes[1]-nodes[0])]
    if not full_rank:
        return list(zip(nodes,weights))
    weights = [weight/(1-node) for node,weight in zip(nodes,weights)]
    return list(zip(nodes,weights))+[(Q(1),t[0]-sum(weights,Q(0))) ]


class IndependentM2LocusTests(unittest.TestCase):
    def test_exact_disk_equivalence_and_determinants(self):
        rng = random.Random(231101)
        for case in range(1500):
            atoms = [(Q(u,4),Q(rng.randrange(1,9),7))
                     for u in rng.sample(range(-4,5),rng.randrange(1,6))]
            n,m,q = sp = species(atoms)
            x,y = Q(rng.randrange(-8,9),3),Q(rng.randrange(1,8),4)
            g = y*y
            if case % 5 == 0:
                s = transform(atoms,x,g,Q(1))
            else:
                s = v.CQ(Q(rng.randrange(-15,16),7),
                         Q(rng.randrange(-15,16),7)/y,g)
            h,l = v.support_matrices(sp,x,g,Q(1),s)
            p,im = s.p,y*s.q
            delta = n*q-m*m
            d = q-2*m*x+n*(x*x+y*y)
            fa = d*(p*p+im*im)-2*n*(n*x-m)*p+(2*n*n*y+delta/y)*im+n**3
            z = v.CQ(x,1,g)
            u = (z*z-1)*s-n*z-m
            fb = (n-q)*u.q-u.p*u.p-g*u.q*u.q
            self.assertEqual(v.determinant(h),-fa)
            self.assertEqual(v.determinant(l),fb)
            self.assertEqual(v.psd(h) and v.psd(l),fa <= 0 and fb >= 0)

    def test_constructive_exact_three_node_completion(self):
        rng = random.Random(41481)
        used_irrational = False
        for count in range(1,6):
            for trial in range(50):
                atoms = [(Q(u,7),Q(rng.randrange(1,11),9))
                         for u in rng.sample(range(-7,8),count)]
                sp = species(atoms)
                x,y = Q(rng.randrange(-5,6),4),Q(rng.randrange(1,8),5)
                s = transform(atoms,x,y*y,Q(1))
                t = v.moments(sp,x,y*y,Q(1),s)
                recovered = radau(t)
                self.assertLessEqual(len(recovered),3)
                for u,weight in recovered:
                    self.assertTrue(-1 <= u <= 1)
                    self.assertTrue(weight > 0)
                    used_irrational |= isinstance(u,R) and u.b != 0
                for j in range(5):
                    self.assertEqual(sum((weight*u**j for u,weight in recovered),Q(0)),t[j])
                mu = [(u,weight*((x-u)**2+y*y)) for u,weight in recovered]
                for observed,expected in zip(species(mu),sp):
                    self.assertEqual(observed,expected)
                back = transform(mu,x,y*y,Q(1))
                self.assertEqual(back.p,s.p)
                self.assertEqual(back.q,s.q)
        self.assertTrue(used_irrational)

    def test_three_atoms_can_be_necessary(self):
        t = (Q(2,3),Q(0),Q(1,3),Q(0),Q(1,3))
        self.assertEqual(v.determinant([[t[i+j] for j in range(3)] for i in range(3)]),Q(1,27))
        s = v.CQ(0,Q(-2,3),Q(1))
        self.assertTrue(v.feasible_s((Q(1),Q(0),Q(2,3)),Q(0),Q(1),Q(1),s))

    def test_positive_original_root_and_mutations(self):
        plus = [(Q(1,3),Q(3,4)),(Q(-1),Q(1,4))]
        minus = [(-u,weight) for u,weight in plus]
        for sector,k,wr,g in [('longitudinal',Q(2),Q(0),Q(4,9)),
                              ('transverse',Q(1),Q(0),Q(1,9))]:
            a,b = (plus,minus) if sector == 'longitudinal' else (minus,plus)
            d = dict(N1='1',F1='0',Q1='1/3',N2='1',F2='0',Q2='1/3',
                     wr=wr,gamma2=g,k=k,sector=sector)
            s1,s2 = transform(a,wr,g,k),transform(b,wr,g,k)
            proof = dict(status='possible',s1=[s1.p,s1.q],s2=[s2.p,s2.q])
            self.assertTrue(v.verify(d,proof))
            self.assertTrue(v.verify_m1_possible(d,proof))
            self.assertTrue(v.verify_m1_possible({key:value for key,value in d.items()
                                                 if key not in ('Q1','Q2')},proof))
            with self.assertRaises(InvalidCertificate):
                v.verify(d,dict(status='possible',s1=[s1.p+Q(1,100),s1.q],s2=proof['s2']))
            with self.assertRaises(InvalidCertificate):
                v.verify(d,dict(status='impossible',proof_kind='m2_four_disks'))
            with self.assertRaises(InvalidCertificate):
                v.verify(dict(d,Q1='1/4'),proof)
            self.assertTrue(v.verify_m1_possible(dict(d,Q1='1/4'),proof))

    def test_forced_atom_positive_and_zero_affine_slope(self):
        plus = [(Q(0),Q(1))]
        minus = [(Q(-1,3),Q(3,4)),(Q(1),Q(1,4))]
        d = dict(N1='1',F1='0',Q1='0',N2='1',F2='0',Q2='1/3',
                 wr='-1',gamma2='2/3',k='5',sector='longitudinal')
        s1,s2 = transform(plus,Q(-1),Q(2,3),Q(5)),transform(minus,Q(-1),Q(2,3),Q(5))
        proof = dict(status='possible',s1=[s1.p,s1.q],s2=[s2.p,s2.q])
        self.assertTrue(v.verify(d,proof))
        self.assertTrue(v.verify_m1_possible(d,proof))
        with self.assertRaises(InvalidCertificate):
            v.verify(d,dict(status='impossible',proof_kind='m2_four_disks'))
        degenerate = dict(d,F1='-1',Q1='1',F2='1',Q2='1',wr='0',gamma2='1',k='1')
        self.assertTrue(v.verify(degenerate,dict(status='impossible',proof_kind='m2_four_disks')))

    def test_independent_empty_ellipses_and_unique_cases(self):
        impossible = dict(status='impossible',proof_kind='m2_four_disks')
        for sector in ('longitudinal','transverse'):
            for wr in (Q(0),Q(2,3),Q(-7,5)):
                for g in (Q(4),Q(17,3),Q(10**30)):
                    d = dict(N1='1',F1='1/5',Q1='2/5',N2='1',F2='-1/7',Q2='3/7',
                             wr=wr,gamma2=g,k='1',sector=sector)
                    self.assertTrue(v.verify(d,impossible))
            for moments in [('0','0','0'),('1','0','0'),('1','0','1'),('1','1','1')]:
                d = dict(N1='1',F1='0',Q1='1/3',N2=moments[0],F2=moments[1],Q2=moments[2],
                         wr='0',gamma2='4',k='1',sector=sector)
                self.assertTrue(v.verify(d,impossible))

    def test_explicit_guards_and_false_controls(self):
        d = dict(N1='1',F1='0',Q1='1/3',N2='1',F2='0',Q2='1/3',
                 wr='0',gamma2='1',k='1',sector='longitudinal')
        bad_inputs = [dict(d,N1='0',Q1='-1'),dict(d,Q1='-1'),dict(d,F1='2'),
                      dict(d,gamma2='0'),dict(d,k='0'),dict(d,sector='other')]
        for case in bad_inputs:
            with self.assertRaises(InvalidCertificate):
                v.verify(case,dict(status='impossible',proof_kind='m2_four_disks'))
        for proof in [dict(status='unresolved'),dict(status='possible',s1=[0,0],s2=[0,0]),
                      dict(status='impossible',proof_kind='untrusted'),
                      dict(status='possible',s1=[float('nan'),0],s2=[0,0])]:
            with self.assertRaises(InvalidCertificate):
                v.verify(d,proof)
            with self.assertRaises(InvalidCertificate):
                v.verify_m1_possible(d,proof)


if __name__ == '__main__':
    unittest.main()
