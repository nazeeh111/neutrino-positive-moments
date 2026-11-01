"""Exact fixed-k counterexample to the endpoint-atom sharpness conjecture.

This does not challenge the proved unrestricted-k second-moment envelope.
A rational IVT bracket certifies a purely imaginary original D_L root that
is faster than the endpoint-compensated two-atom family at the SAME k and
separate moments. Standard library and independent C arithmetic only.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json
from independent_verifier import C

n=F(1);q=F(1,10);k=F(1,2)
plus=[(F(1,8),F(32,37)),(-F(4,5),F(5,37))]
minus=[(-u,w) for u,w in plus]
lo,hi=F(37,125),F(297,1000)
family_gamma2=2*(1-q)*k*q-k*k*q*q

def require(cond,message):
 if not cond:raise ValueError(message)

def dispersion(gamma):
 omega=C(0,gamma)
 vals=[]
 for j in range(3):
  vals.append(sum((C(sign*mass*u**j)/(omega-k*u)
   for atoms,sign in ((plus,1),(minus,-1)) for u,mass in atoms),C()))
 i0,i1,i2=vals
 return (1+i0)*(-1+i2)-i1*i1

for atoms in [plus,minus]:
 require(all(-1<=u<=1 and w>0 for u,w in atoms),'invalid positive measure')
 require(sum(w for u,w in atoms)==n,'mass')
 require(sum(u*w for u,w in atoms)==0,'flux')
 require(sum(u*u*w for u,w in atoms)==q,'second moment')
dl,dh=dispersion(lo),dispersion(hi)
require(dl.im==dh.im==0,'dispersion not real')
require(dl.re>0 and dh.re<0,'no strict sign bracket')
require(lo*lo>family_gamma2,'not strictly above family')
# The direct imaginary-axis rational numerator in lambda=gamma² is:
# P(lambda)=lambda² + (617/6400)lambda -13/800.
# Its derivative is strictly positive on lambda>=0; the positive root is unique.
p=lambda v:v*v+F(617,6400)*v-F(13,800)
require(p(lo*lo)<0 and p(hi*hi)>0,'polynomial bracket')
result={'passed':True,'scope':__doc__,'N':str(n),'Q_plus':str(q),'Q_minus':str(q),'k':str(k),
 'plus_atoms':[[str(u),str(w)] for u,w in plus],
 'minus_atoms':[[str(u),str(w)] for u,w in minus],
 'gamma_bracket':[str(lo),str(hi)],'direct_D_long_endpoints':[str(dl.re),str(dh.re)],
 'endpoint_family_gamma_squared':str(family_gamma2),'lower_bracket_squared':str(lo*lo),
 'numerator_polynomial_lambda':['1',str(F(617,6400)),str(-F(13,800))],
 'meaning':'Endpoint-compensated two-atom family is not sharp at this fixed k, even on the imaginary axis.',
 'limitation':'No all-real-frequency finite-k sharp theorem is proved; original unrestricted-k cap remains valid.',
 'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
Path('results/finite-k-counterexample.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
