"""Independent exact checks of a Gram-disk fixed-k growth cap.

These corroborate a written all-real-frequency upper bound and its sharp
low-second-moment subregime. They do NOT prove a complete fixed-k envelope.
No production or numerical module is imported.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json,random,unittest
from independent_verifier import C


def norm2(z):return z.re*z.re+z.im*z.im


def direct(atoms_plus,atoms_minus,n,k,gamma):
 w=C(0,gamma)
 vals=[sum((C(sign*mass*u**j)/(w-k*u)
   for atoms,sign in ((atoms_plus,1),(atoms_minus,-1)) for u,mass in atoms),C())
   for j in range(3)]
 i0,i1,i2=vals
 return (1+i0)*(-1+i2)-i1*i1,-1+(i0-i2)/2


class FixedKGramCap(unittest.TestCase):
 def test_transformed_radius_has_global_real_frequency_max_at_zero(self):
  rng=random.Random(2026110107)
  for _ in range(1000):
   a=F(rng.randint(1,19),20);x=F(rng.randint(-30,30),rng.randint(1,20));y=F(rng.randint(1,30),rng.randint(1,20))
   t=y*y;d=x*x
   # Multiply the radius comparison by all strictly positive denominators.
   gap=(1+t)**2*(a+t+d)**2-(a+t)**2*((t+d-1)**2+4*t)
   factored=d*(2*(a+t)*(1+a+(3-a)*t)+(1-a)*(1+a+2*t)*d)
   self.assertEqual(gap,factored)
   self.assertGreaterEqual(gap,0)
   self.assertEqual(gap==0,x==0)

 def test_gram_schur_disk_identity(self):
  rng=random.Random(117021)
  for _ in range(400):
   n=F(rng.randint(1,9),rng.randint(1,7));a=F(rng.randint(1,9),10)
   x=F(rng.randint(-10,10),7);y=F(rng.randint(1,10),9);z=C(x,y)
   s=C(F(rng.randint(-12,12),8),F(rng.randint(-12,12),8))
   D=a+x*x+y*y
   center=C(n*x/D,-n*(y+a/(2*y))/D)
   rad2=n*n*a*a/(4*y*y*D*D)
   # Necessary Gram positivity: -Im S/y - |S|²/N - |zS-N|²/(Na) >= 0.
   schur=-s.im/y-norm2(s)/n-norm2(z*s-n)/(n*a)
   self.assertEqual(schur,D/(n*a)*(rad2-norm2(s-center)))

 def test_exact_sharp_subregime_and_original_factors(self):
  count=0;improved=0
  for m in [6,7,8,10,20,100]:
   p,q=F(1,m),F(6,m);a=p*q
   for n in [F(1,3),F(1),F(7,2)]:
    plus=[(p,n*F(6,7)),(-q,n*F(1,7))];minus=[(-u,w) for u,w in plus]
    for y in [F(2,m),F(3,m)]:
     self.assertEqual(y+a/y,q-p)
     self.assertLessEqual(y+a/y,1-a)
     kl=n*a*(1+y*y)/(y*(a+y*y));gl=kl*y
     for branch,k,gamma,first,second,c in [
      ('longitudinal',kl,gl,plus,minus,F(1)),
      ('transverse',kl/2,gl/2,minus,plus,F(2))]:
      for atoms in [first,second]:
       self.assertEqual(sum(w for u,w in atoms),n)
       self.assertEqual(sum(u*w for u,w in atoms),0)
       self.assertEqual(sum(u*u*w for u,w in atoms),n*a)
       self.assertTrue(all(-1<=u<=1 and w>0 for u,w in atoms))
      dl,dt=direct(first,second,n,k,gamma)
      self.assertEqual(dl if branch=='longitudinal' else dt,C())
      # c*k <= N*a*(1+y²)/(y*(a+y²)) is saturated.
      self.assertEqual(c*k,n*a*(1+y*y)/(y*(a+y*y)))
      self.assertEqual(gamma**3-(n*a/c)*gamma**2+a*k*k*gamma-(n*a/c)*k*k,0)
      family=(2*n*(1-a)*k*a if branch=='longitudinal' else n*(1-a)*k*a)-k*k*a*a
      self.assertGreaterEqual(gamma*gamma,family)
      if gamma*gamma>family:improved+=1
      count+=1
  self.assertEqual(count,72);self.assertEqual(improved,60)

if __name__=='__main__':
 suite=unittest.defaultTestLoader.loadTestsFromTestCase(FixedKGramCap)
 run=unittest.TextTestRunner(verbosity=2).run(suite)
 data={'passed':run.wasSuccessful(),'tests':run.testsRun,
  'formal_radius_identity_rational_checks':1000,'gram_disk_identity_checks':400,
  'exact_sharp_attainers_original_dispersion_checks':72,'strict_improvements_over_endpoint_family':60,
  'scope':'All-real-frequency upper bound; exact attainment only in stated low-Q subregime. Not full fixed-k envelope.',
  'standard_library_only':True,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
 Path('results/finite-k-cap-checks.json').write_text(json.dumps(data,indent=2)+'\n')
 if not run.wasSuccessful():raise SystemExit(1)
