"""Sharp finite-band growth for positive axial species moments(N,0,Q).

See finite-k-transition-audit.md for the all-real bound and atomic attainment.
Standard library, exact rational branch decisions and rational enclosures.
No physical UV cutoff, universal instability or speed advantage is implied.
"""
from fractions import Fraction as F
from math import isqrt
import argparse,json

def rational(value):
    if isinstance(value,bool) or not isinstance(value,(int,str,F)):
        raise ValueError('Use exact integers,rational strings or Fraction values')
    return F(value)

def bracket(sign,upper,bits):
    lo,hi=F(0),upper
    while hi-lo>upper/F(2**bits):
        mid=(lo+hi)/2;value=sign(mid)
        if value==0:return mid,mid
        if value<0:lo=mid
        else:hi=mid
    return lo,hi

def envelope(mass,second,kmax,sector='longitudinal',bits=80):
    """Sharp maximum over all positive completions,all Re(omega),0<=k<=K.

    Every positive maximum is attained at k=K. Different inputs may need
    different atomic completions. This does not identify the actual profile.
    """
    n,q,k=map(rational,(mass,second,kmax))
    if n<=0 or q<0 or q>n or k<0:raise ValueError('Requires N>0,0<=Q<=N,K>=0')
    if sector not in ('longitudinal','transverse'):raise ValueError('Unknown sector')
    if isinstance(bits,bool) or not isinstance(bits,int) or not 1<=bits<=4096:
        raise ValueError('bits must be an integer in[1,4096]')
    c=1 if sector=='longitudinal' else 2
    a=q/n;r=1-a;kappa=c*k/n;square=None
    if k==0 or a in (0,1):
        branch='stable';lo=hi=F(0);square=F(0)
    elif kappa>=r/a:
        branch='plateau';lo=hi=r;square=r*r
    else:
        t=2*a*r/kappa-a*a
        if (t+a)**2<r*r*t:
            branch='cubic'
            lo,hi=bracket(lambda g:g**3-a*g*g+a*kappa*kappa*g-a*kappa*kappa,r,bits)
        else:
            branch='endpoint';square=2*a*r*kappa-a*a*kappa*kappa
            p,s=isqrt(square.numerator),isqrt(square.denominator)
            if p*p==square.numerator and s*s==square.denominator:lo=hi=F(p,s)
            else:lo,hi=bracket(lambda g:g*g-square,r,bits)
    scale=n/c
    result={'mass':str(n),'second':str(q),'kmax':str(k),'sector':sector,
            'normalized_second':str(a),'normalized_k':str(kappa),'branch':branch,
            'growth_lower':str(scale*lo),'growth_upper':str(scale*hi),
            'growth_squared_lower':str(scale*scale*lo*lo),
            'growth_squared_upper':str(scale*scale*hi*hi),'enclosure_bits':bits,
            'quantifier':'Sharp maximum over positive species completions,all real omega parts,0<=k<=kmax; not universal instability'}
    if square is not None:result['exact_growth_squared']=str(scale*scale*square)
    if branch=='cubic':result['normalized_root_polynomial']=['1',str(-a),str(a*kappa*kappa),str(-a*kappa*kappa)]
    return result

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mass',required=True);p.add_argument('--second',required=True);p.add_argument('--kmax',required=True)
    p.add_argument('--sector',choices=('longitudinal','transverse'),default='longitudinal');p.add_argument('--bits',type=int,default=80)
    a=p.parse_args();print(json.dumps(envelope(a.mass,a.second,a.kmax,a.sector,a.bits),indent=2))
if __name__=='__main__':main()
