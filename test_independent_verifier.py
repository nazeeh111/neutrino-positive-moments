"""Independent certificate tests: no production/optimizer imports."""
import copy
import itertools
import random
import unittest
from fractions import Fraction as Q

from independent_verifier import (C, R, InvalidCertificate, check_witness,
    direct_dispersion, normalized_constraints, parse, positive_semidefinite,
    sector_coefficients, support_matrices, verify, independent_disks,
    independent_candidates, in_disk, transformed_moments)


def atomic_data(atoms, z):
    n = sum((m for u,m in atoms), Q(0))
    f = sum((m*u for u,m in atoms), Q(0))
    s = sum((m/(z-u) for u,m in atoms), C())
    t = tuple(sum((m*u**j/((z.re-u)**2+z.im**2) for u,m in atoms),Q(0)) for j in range(4))
    return n,f,s,t


def quadrature_from_central_moments(t):
    t0,t1,t2,t3 = t
    mean = t1/t0
    variance = t2/t0-mean**2
    if variance == 0:
        return {'kind':'one','u':mean,'nu_mass':t0}
    skew = t3/t0-3*mean*t2/t0+2*mean**3
    a = -2*mean-skew/variance
    b = mean**2-variance+mean*skew/variance
    return {'kind':'two','p2_a':a,'p2_b':b,'nu_t0':t0,'nu_t1':t1}


def possible_fixture(atoms1, atoms2, w, k, sector):
    z=w/k
    n1,f1,s1,t1 = atomic_data(atoms1,z)
    n2,f2,s2,t2 = atomic_data(atoms2,z)
    data=dict(N1=n1,F1=f1,N2=n2,F2=f2,w=w,k=k,sector=sector)
    proof=dict(status='possible',s1=s1,s2=s2,
               witnesses=[quadrature_from_central_moments(t1),quadrature_from_central_moments(t2)])
    return data,proof


def det(m):
    return (m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1])
            -m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0])
            +m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]))


def tiny_exact_dual(data):
    """Finite rational search for test fixtures only, with no optimizer."""
    n1,f1,n2,f2,w,k,sector=parse(data)
    a,b=sector_coefficients(n1,f1,n2,f2,w,k,sector)
    hs=normalized_constraints(n1,f1,n2,f2,w/k,k*b/a)
    pool=[]
    for j,h in enumerate(hs):
        for v in ((1,0),(0,1),(1,1),(1,-1),(1,2),(1,-2),(2,1),(2,-1)):
            row=tuple(v[0]**2*m[0]+2*v[0]*v[1]*m[1]+v[1]**2*m[2] for m in h)
            pool.append((j,v,row))
    for j,row in enumerate(((1,1,0),(1,-1,0),(1,0,1),(0,0,-1))):
        pool.append((4+j,None,tuple(map(Q,row))))
    for rows in itertools.combinations(pool,3):
        matrix=[[Q(1)]*3,[r[2][1] for r in rows],[r[2][2] for r in rows]]
        d=det(matrix)
        if d == 0:
            continue
        weights=[]
        for j in range(3):
            changed=[row[:] for row in matrix]
            for i in range(3):
                changed[i][j]=Q(i==0)
            weights.append(det(changed)/d)
        value=sum((weight*row[2][0] for weight,row in zip(weights,rows)),Q(0))
        if min(weights)>=0 and value<0:
            return dict(status='impossible',proof_kind='rational_dual',dual={
                'terms':[dict(constraint=r[0],vector=r[1],weight=weight) for r,weight in zip(rows,weights)],
                'negative_constant':value})
    raise AssertionError('no small test dual found')


class IndependentChecks(unittest.TestCase):
    def setUp(self):
        self.rng=random.Random(20261101)

    def atoms(self):
        return [(Q(self.rng.randint(-8,8),8),Q(self.rng.randint(1,9),7))
                for _ in range(self.rng.randint(1,8))]

    def test_random_atomic_roundtrips_and_direct_sectors(self):
        for _ in range(160):
            atoms1,atoms2=self.atoms(),self.atoms()
            w=C(Q(self.rng.randint(-10,10),7),Q(self.rng.randint(1,10),9))
            k=Q(self.rng.randint(1,8),5)
            n1,f1,s1,t1=atomic_data(atoms1,w/k)
            n2,f2,s2,t2=atomic_data(atoms2,w/k)
            for n,f,s,t in ((n1,f1,s1,t1),(n2,f2,s2,t2)):
                self.assertTrue(check_witness(n,f,w/k,s,quadrature_from_central_moments(t)))
                self.assertTrue(all(positive_semidefinite(h) for h in support_matrices(n,f,w/k,s)))
            longitudinal,transverse=direct_dispersion(atoms1,atoms2,w,k)
            for sector,direct,factor in (('longitudinal',longitudinal,1),('transverse',transverse,2)):
                a,b=sector_coefficients(n1,f1,n2,f2,w,k,sector)
                self.assertEqual(a*(s1-s2)-k*b,factor*k**3*direct)

    def test_longitudinal_attained_endpoint_root(self):
        data,proof=possible_fixture([(Q(1),Q(5,8))],[(-Q(1),Q(5,8))],C(0,1),Q(3,4),'longitudinal')
        self.assertEqual(direct_dispersion([(Q(1),Q(5,8))],[(-Q(1),Q(5,8))],C(0,1),Q(3,4))[0],C())
        self.assertTrue(verify(data,proof))
        bad=copy.deepcopy(proof)
        bad['witnesses'][0]['nu_mass']+=Q(1,100)
        self.assertRaises(InvalidCertificate,verify,data,bad)

    def test_transverse_root_and_sector_mutation(self):
        data,proof=possible_fixture([(-Q(1,2),Q(10,3))],[(Q(1,2),Q(10,3))],C(0,1),Q(1),'transverse')
        self.assertTrue(verify(data,proof))
        self.assertEqual(direct_dispersion([(-Q(1,2),Q(10,3))],[(Q(1,2),Q(10,3))],C(0,1),Q(1))[1],C())
        self.assertRaises(InvalidCertificate,verify,{**data,'sector':'longitudinal'},proof)
        bad=copy.deepcopy(proof)
        bad['s1']=bad['s1']+C(Q(1,10))
        self.assertRaises(InvalidCertificate,verify,data,bad)

    def test_zero_A_nonzero_B(self):
        data=dict(N1=1,F1=-1,N2=1,F2=1,w=C(0,1),k=1)
        proof=dict(status='impossible',proof_kind='zero_A_nonzero_B')
        self.assertTrue(verify(data,proof))
        self.assertRaises(InvalidCertificate,verify,{**data,'w':C(0,2)},proof)

    def test_fixed_zero_species_impossibility(self):
        data=dict(N1=0,F1=0,N2=0,F2=0,w=C(0,1),k=1)
        proof=dict(status='impossible',proof_kind='fixed_species',s1=C(),s2=C(Q(1,2)),violations=[(1,0)])
        self.assertTrue(verify(data,proof))
        bad={**proof,'violations':[(0,0)]}
        self.assertRaises(InvalidCertificate,verify,data,bad)
        self.assertRaises(InvalidCertificate,verify,{**data,'N1':1,'N2':1},proof)

    def test_rational_dual_and_false_certificates(self):
        for sector in ('longitudinal','transverse'):
            data=dict(N1=1,F1=0,N2=1,F2=0,w=C(0,4),k=1,sector=sector)
            proof=tiny_exact_dual(data)
            self.assertTrue(verify(data,proof))
            for mutation in ('constant','weight','index','vector'):
                bad=copy.deepcopy(proof)
                if mutation=='constant':
                    bad['dual']['negative_constant']-=1
                elif mutation=='weight':
                    bad['dual']['terms'][0]['weight']=-1
                elif mutation=='index':
                    bad['dual']['terms'][0]['constraint']=-1
                else:
                    term=next(t for t in bad['dual']['terms'] if t['weight']>0 and t['constraint']<4)
                    term['vector']=(0,0)
                self.assertRaises(InvalidCertificate,verify,data,bad)

    def test_bad_support_weight_and_serialized_exact_input(self):
        z=C(0,1)
        for witness in (
            dict(kind='one',u=2,nu_mass=1),
            dict(kind='one',u=0,nu_mass=-1),
            dict(kind='two',p2_a=0,p2_b=-4,nu_t0=1,nu_t1=0),
            dict(kind='two',p2_a=0,p2_b=-Q(1,4),nu_t0=1,nu_t1=1),
            dict(kind='two',p2_a=0,p2_b=1,nu_t0=1,nu_t1=0),
            dict(kind='zero')):
            self.assertRaises(InvalidCertificate,check_witness,Q(1),Q(0),z,C(0,-1),witness)
        data=dict(N1='0',F1='0',N2='0',F2='0',w={'re':'0','im':'1'},k='1')
        proof=dict(status='impossible',proof_kind='fixed_species',s1=['0','0'],s2=['1/2','0'],violations=[[1,0]])
        self.assertTrue(verify(data,proof))

    def test_invalid_inputs_and_unresolved_refuse(self):
        data=dict(N1=1,F1=0,N2=1,F2=0,w=C(0,1),k=1)
        proof=dict(status='unresolved')
        for change in ({},{'N1':-1},{'F1':2},{'N2':0,'F2':1},{'w':C(0,0)},
                       {'k':0},{'k':-1},{'sector':'wrong'},{'N1':1.0},{'N1':True},
                       {'w':{'re':'NaN','im':1}}):
            self.assertRaises(InvalidCertificate,verify,{**data,**change},proof)

    def test_quadratic_order_division_and_near_cancellation(self):
        self.assertEqual(R(0,1,8),R(0,2,2))
        self.assertEqual(R(1,2,Q(9,4)),4)
        self.assertEqual(R(3,0,0),3)
        self.assertRaises(InvalidCertificate,R,0,1,-1)
        self.assertRaises(InvalidCertificate,lambda: R(1,1,2)+R(1,1,3))
        p,q=1,0
        for n in range(1,201):
            p,q=p+2*q,p+q
            tiny=R(p,-q,2)
            self.assertEqual(tiny.sign(),1 if n%2 == 0 else -1)
            self.assertEqual(tiny*(1/tiny),1)
            self.assertEqual(tiny+R(-p,q,2),0)
        self.assertRaises(InvalidCertificate,lambda: R(1)/R())

    def test_closed_disk_degeneracies_and_pair_candidates(self):
        def feasible(ds):
            return [p for p in independent_candidates(ds) if all(in_disk(p,d) for d in ds)]
        unit=((Q(0),Q(0)),Q(1))
        self.assertTrue(feasible([unit,unit,((Q(0),Q(0)),Q(4))]))
        self.assertEqual(feasible([unit,((Q(2),Q(0)),Q(1))]),[(R(1),R(0))])
        self.assertFalse(feasible([unit,((Q(3),Q(0)),Q(1))]))
        self.assertTrue(feasible([unit,((Q(1),Q(0)),Q(0))]))
        self.assertFalse(feasible([unit,((Q(2),Q(0)),Q(0))]))
        lens=[unit,((Q(3,2),Q(0)),Q(1))]
        points=feasible(lens)
        self.assertTrue(points)
        self.assertFalse(any(all(in_disk(c,d) for d in lens) for c,_ in lens))
        # Pairwise overlaps do not imply a common intersection.
        self.assertFalse(feasible([((-Q(1),Q(0)),Q(1)),((Q(1),Q(0)),Q(1)),((Q(0),Q(3,2)),Q(1))]))

    def test_algebraic_witness_serialization_and_false_geometry(self):
        data=dict(N1=1,F1=0,N2=1,F2=0,w=C(0,Q(1,10)),k=1)
        n1,f1,n2,f2,w,k,sector=parse(data)
        a,b=sector_coefficients(n1,f1,n2,f2,w,k,sector)
        target,z=k*b/a,w/k
        disks=independent_disks(n1,f1,n2,f2,z,target)
        point=next(p for p in independent_candidates(disks)
                   if any(x.b for x in p) and all(in_disk(p,d) for d in disks))
        s1=C(*point)
        s2=s1-target
        proof=dict(status='possible',proof_kind='algebraic_four_disks',s1=s1,s2=s2,
                   witnesses=[quadrature_from_central_moments(transformed_moments(n,f,z,s))
                              for n,f,s in ((n1,f1,s1),(n2,f2,s2))])
        self.assertTrue(verify(data,proof))
        def encode(value):
            if isinstance(value,R):
                return {name:str(getattr(value,name)) for name in ('a','b','d')}
            if isinstance(value,C):
                return {'re':encode(value.re),'im':encode(value.im)}
            if isinstance(value,Q):
                return str(value)
            if isinstance(value,dict):
                return {key:encode(v) for key,v in value.items()}
            if isinstance(value,(list,tuple)):
                return [encode(v) for v in value]
            return value
        self.assertTrue(verify(encode(data),encode(proof)))
        self.assertRaises(InvalidCertificate,verify,data,dict(status='impossible',proof_kind='four_disks'))
        for sector in ('longitudinal','transverse'):
            self.assertTrue(verify({**data,'w':C(0,4),'sector':sector},
                                   dict(status='impossible',proof_kind='four_disks')))


if __name__=='__main__':
    unittest.main(verbosity=2)
