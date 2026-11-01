"""Written-proof corroboration: exact identities,original factors and enclosures.
No floating arithmetic or production geometry. Finite checks are not proofs.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json,random,unittest
from independent_m2_verifier import CQ
from finite_band_envelope import envelope

COUNTS={}
def factors(plus,minus,k,g2):
    w=CQ(0,1,g2)
    values=[sum((CQ(sign*mass*u**j,0,g2)/(w-k*u) for atoms,sign in ((plus,1),(minus,-1)) for u,mass in atoms),CQ(0,0,g2)) for j in range(3)]
    i0,i1,i2=values
    return (1+i0)*(-1+i2)-i1*i1,-1+(i0-i2)/2

def moments(test,atoms,n,a):
    test.assertEqual(sum(w for u,w in atoms),n)
    test.assertEqual(sum(u*w for u,w in atoms),0)
    test.assertEqual(sum(u*u*w for u,w in atoms),n*a)
    test.assertTrue(all(-1<=u<=1 and w>=0 for u,w in atoms))

class FullFiniteBand(unittest.TestCase):
 def test_combination_disk_identity(self):
    rng=random.Random(2026110118)
    for _ in range(1000):
        a=F(rng.randint(1,99),100);x=F(rng.randint(-30,30),11);y=F(rng.randint(1,30),13);b=F(rng.randint(0,20),20)
        t=y*y;d=x*x;D=a+t+d
        r0=(1-a)**2/(4*t);r1=a*a*((d+t-1)**2+4*t)/(4*t*D*D)
        dist=((a*a+t-d)**2+4*d*t)/(4*t*D*D)
        expr=(1-a*a)*(a-t-d)/(4*t*D)
        self.assertEqual(r1-dist,expr)
        self.assertEqual((1-b)*r0+b*r1-b*(1-b)*dist,(1-b)*r0+b*b*r1+b*(1-b)*expr)
        ra=a*a*(1+t)**2/(4*t*(a+t)**2);ea=(1-a*a)*(a-t)/(4*t*(a+t))
        self.assertGreaterEqual(ra,r1);self.assertGreaterEqual(ea,expr)
        self.assertGreaterEqual((1-b)*r0+b*b*ra+b*(1-b)*ea,(1-b)*r0+b*b*r1+b*(1-b)*expr)

 def test_endpoint_plateau_original_factors(self):
    counts={'endpoint':0,'plateau':0}
    for a in map(F,('1/1000','3/50','3/20','1/5','1/3','9/10')):
     for n in map(F,('1/3','1','7/2')):
      for c,sector in ((1,'longitudinal'),(2,'transverse')):
       threshold=n*(1-a)/(c*a)
       for ratio in map(F,('1/100','1/10','1/2','9/10','1','2','100')):
        k=threshold*ratio;result=envelope(n,n*a,k,sector);branch=result['branch']
        if branch=='cubic':continue
        if branch=='endpoint':
            plus=[(a,n/(1+a)),(-F(1),n*a/(1+a))];g2=2*n*a*(1-a)*k/c-a*a*k*k
        else:
            self.assertEqual(branch,'plateau');y=n*(1-a)/(c*k)
            plus=[(y,n*(1-a)/(1-y*y)),(F(1),n*(a-y)/(2*(1-y))),(-F(1),n*(a+y)/(2*(1+y)))];g2=(n*(1-a)/c)**2
        minus=[(-u,w) for u,w in plus]
        for atoms in (plus,minus):moments(self,atoms,n,a)
        if c==2:plus,minus=minus,plus
        dl,dt=factors(plus,minus,k,g2);self.assertTrue((dl if c==1 else dt).zero())
        self.assertEqual(F(result['exact_growth_squared']),g2);counts[branch]+=1
    self.assertGreater(counts['endpoint'],0);self.assertGreater(counts['plateau'],0);COUNTS.update(counts)

 def test_cubic_attainers(self):
    count=0
    for m in (6,7,8,10,20,100):
     p,q=F(1,m),F(6,m);a=p*q
     for n in (F(1,3),F(1),F(7,2)):
      plus=[(p,n*q/(p+q)),(-q,n*p/(p+q))];minus=[(-u,w) for u,w in plus]
      for y in (F(2,m),F(3,m)):
       for c,sector in ((1,'longitudinal'),(2,'transverse')):
        k=n*a*(1+y*y)/(c*y*(a+y*y));g=k*y;r=envelope(n,n*a,k,sector)
        self.assertIn(r['branch'],('endpoint','cubic'));self.assertLessEqual(F(r['growth_lower']),g);self.assertGreaterEqual(F(r['growth_upper']),g)
        first,second=(plus,minus) if c==1 else (minus,plus)
        dl,dt=factors(first,second,k,g*g);self.assertTrue((dl if c==1 else dt).zero())
        for atoms in (first,second):moments(self,atoms,n,a)
        count+=1
    self.assertEqual(count,72)

 def test_transitions_and_adjacent_inputs(self):
    for y in map(F,('1/10','1/5','1/3','1/2','2/3','9/10')):
     a=y*(1-y)/(1+y)
     for n in (F(1,5),F(1),F(9)):
      for c,sector in ((1,'longitudinal'),(2,'transverse')):
       k=n*a*(1+y*y)/(c*y*(a+y*y));r=envelope(n,n*a,k,sector)
       self.assertEqual(r['branch'],'endpoint');self.assertEqual(F(r['exact_growth_squared']),(k*y)**2)
       self.assertEqual(F(r['growth_lower']),k*y);self.assertEqual(F(r['growth_upper']),k*y)
       lo=envelope(n,n*a,k*(1-F(1,10**12)),sector);hi=envelope(n,n*a,k*(1+F(1,10**12)),sector)
       self.assertIn('cubic',(lo['branch'],hi['branch']));self.assertIn('endpoint',(lo['branch'],hi['branch']))

 def test_enclosures_scaling_and_corners(self):
    rng=random.Random(108113)
    for _ in range(150):
        n=F(rng.randint(1,17),rng.randint(1,13));a=F(rng.randint(1,99),100);k=F(rng.randint(1,71),17)
        sector=rng.choice(('longitudinal','transverse'));c=1 if sector=='longitudinal' else 2;r=envelope(n,n*a,k,sector,bits=50)
        lo,hi=F(r['growth_lower']),F(r['growth_upper']);self.assertLessEqual(lo,hi);self.assertLessEqual(hi-lo,n*(1-a)/(c*2**50));self.assertLessEqual(hi,n*(1-a)/c)
        if r['branch']=='cubic':
            p=lambda g:g**3-(n*a/c)*g*g+a*k*k*g-(n*a/c)*k*k
            self.assertLessEqual(p(lo),0);self.assertGreaterEqual(p(hi),0)
        else:
            square=F(r['exact_growth_squared']);self.assertLessEqual(lo*lo,square);self.assertGreaterEqual(hi*hi,square)
        for scale in (F(1,10**50),F(10**50)):
            scaled=envelope(n*scale,n*a*scale,k*scale,sector,bits=50)
            self.assertEqual(scaled['branch'],r['branch']);self.assertEqual(F(scaled['growth_lower']),lo*scale);self.assertEqual(F(scaled['growth_upper']),hi*scale)
    for n in (F(1),F(3,7)):
     for q,k in ((0,1),(n,1),(n/3,0)):
        r=envelope(n,q,k);self.assertEqual(r['branch'],'stable');self.assertEqual(F(r['growth_upper']),0)
    for args in ((0,0,1),(1,-1,1),(1,2,1),(1,0,-1)):
        with self.assertRaises(ValueError):envelope(*args)
    with self.assertRaises(ValueError):envelope(1.,'1/3',1)

if __name__=='__main__':
    run=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(FullFiniteBand));root=Path(__file__).resolve().parent
    report={'passed':run.wasSuccessful(),'tests':run.testsRun,'formal_combination_disk_checks':1000,'rational_cubic_attainers':72,'other_original_factor_attainers':COUNTS,'exact_transition_cases':36,'random_enclosures':150,'scale_extremes':['1e-50','1e50'],'scope':'Full finite-band theorem corroboration,original-factor witnesses and enclosures;not physical usefulness or priority validation','source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'algorithm_sha256':hashlib.sha256((root/'finite_band_envelope.py').read_bytes()).hexdigest()}
    (root/'results/finite-band-validation.json').write_text(json.dumps(report,indent=2)+'\n')
    if not run.wasSuccessful():raise SystemExit(1)
