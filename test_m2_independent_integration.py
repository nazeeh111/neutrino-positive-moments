"""Cross-check production against the independent M2 checker and old M1 route.

The checker itself is stdlib-only. This integration test additionally needs
the unchanged M1 module's existing NumPy/SciPy environment.
"""
from fractions import Fraction as Q
import random
import unittest

import independent_m2_verifier as v
import independent_verifier as old_check
import m2_frequency as production
from test_independent_m2_locus import species, transform


class IntegrationChecks(unittest.TestCase):
    def test_250_atomic_values_and_independent_ellipse_coefficients(self):
        rng = random.Random(103331)
        compared = 0
        for case in range(250):
            atoms = [(Q(u,9),Q(rng.randrange(1,10),7))
                     for u in rng.sample(range(-9,10),3)]
            sp = species(atoms)
            wr,g,k = Q(rng.randrange(-7,8),4),Q(rng.randrange(1,15),11),Q(rng.randrange(1,10),5)
            s = transform(atoms,wr,g,k)
            ellipses = production.species_ellipses(sp,wr,g,k)
            self.assertTrue(all(v.in_ellipse((s.p,s.q),e,g) for e in ellipses))
            independently = v.determinant_ellipses(sp,sp,wr,g,k,v.CQ(0,0,g))[:2]
            self.assertEqual(ellipses,independently)
            compared += 1
        self.assertEqual(compared,250)

    def test_200_classifications_and_m1_inclusion(self):
        from moments import QC, classify as old_classify
        rng = random.Random(7102026)
        counts = {'possible':0,'impossible':0,'m1_checks':0,'radical_witnesses':0}
        for case in range(200):
            if case < 100:
                sp1 = (Q(1),Q(0),Q(rng.randrange(1,10),10))
                sp2 = (Q(1),Q(0),Q(rng.randrange(1,10),10))
            else:
                atomsets = [[(Q(u,6),Q(rng.randrange(1,8),10))
                             for u in rng.sample(range(-6,7),rng.randrange(1,5))]
                            for _ in range(2)]
                sp1,sp2 = map(species,atomsets)
            wr = Q(rng.randrange(-5,6),8)
            k = Q(rng.randrange(1,20),6)
            gamma = Q(rng.randrange(1,10),15)
            g = gamma*gamma if case % 2 == 0 else Q(rng.randrange(1,16),64)
            sector = ('longitudinal','transverse')[case % 2]
            data = dict(N1=sp1[0],F1=sp1[1],Q1=sp1[2],N2=sp2[0],F2=sp2[1],Q2=sp2[2],
                        wr=wr,gamma2=g,k=k,sector=sector)
            proof = production.classify(data)
            self.assertTrue(v.verify(data,proof))
            self.assertTrue(v.verify(production.jsonable(data),production.jsonable(proof)))
            counts[proof['status']] += 1
            if proof['status'] == 'possible':
                self.assertTrue(v.verify_m1_possible(data,proof))
                counts['radical_witnesses'] += any(getattr(x,'b',0) != 0
                                                 for s in ('s1','s2') for x in proof[s])
                root = old_check.rational_sqrt(g)
                if root is not None:
                    old_proof = old_classify(sp1[0],sp1[1],sp2[0],sp2[1],QC(wr,root),k,sector)
                    self.assertEqual(old_proof['status'],'possible')
                    old_data = dict(N1=sp1[0],F1=sp1[1],N2=sp2[0],F2=sp2[1],w=[wr,root],k=k,sector=sector)
                    self.assertTrue(old_check.verify(old_data,old_proof))
                    counts['m1_checks'] += 1
        self.assertGreater(counts['possible'],0)
        self.assertGreater(counts['impossible'],0)
        self.assertGreater(counts['m1_checks'],0)
        self.assertGreater(counts['radical_witnesses'],0)
        print('Independent seeded classification counts:',counts)


if __name__ == '__main__': unittest.main()
