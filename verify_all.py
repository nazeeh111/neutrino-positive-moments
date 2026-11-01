"""Replay all saved exact proofs with Python's standard library only.

No import of production moments/geometry, NumPy, SciPy or optimizer is made.
Smooth-profile Arb integrations are separate evidence, not checked here.
"""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json
from independent_verifier import verify,parse,sector_coefficients,C
ROOT=Path(__file__).resolve().parent

def atomic(data,atoms):
    n1,f1,n2,f2,w,k,sector=parse(data);z=w/k
    ss=[]
    for at,N,flux in zip(atoms,(n1,n2),(f1,f2)):
        at=[(F(u),F(m)) for u,m in at]
        if not all(-1<=u<=1 and m>=0 for u,m in at):raise ValueError('invalid grid atom')
        if sum((m for u,m in at),F(0))!=N or sum((u*m for u,m in at),F(0))!=flux:
            raise ValueError('grid moments do not match')
        ss.append(sum((m/(z-u) for u,m in at),C()))
    A,B=sector_coefficients(n1,f1,n2,f2,w,k,sector)
    if A*(ss[0]-ss[1])!=k*B:raise ValueError('grid root equation failed')

def farkas(data,n,proof):
    n1,f1,n2,f2,w,k,sector=parse(data);z=w/k
    A,B=sector_coefficients(n1,f1,n2,f2,w,k,sector);target=k*B/A
    y=list(map(F,proof['dual']));minimum=None
    for species in (0,1):
        for j in range(n):
            u=F(2*j,n-1)-1;s=1/(z-u);sign=1 if species==0 else -1
            col=(F(species==0),u if species==0 else F(0),F(species==1),u if species==1 else F(0),sign*s.re,sign*s.im)
            v=sum((a*b for a,b in zip(col,y)),F(0))
            minimum=v if minimum is None else min(minimum,v)
    obj=sum((a*b for a,b in zip((n1,f1,n2,f2,target.re,target.im),y)),F(0))
    if not minimum>=0 or not obj<0:raise ValueError('invalid grid Farkas proof')
    if minimum!=F(proof['minimum_column_product']) or obj!=F(proof['negative_rhs_product']):
        raise ValueError('stored grid proof fields changed')

def main():
    counts={'screen_exact':0,'stress_exact':0,'numeric_exact':0,'grid_atomic':0,'grid_farkas':0}
    for name,key in (('screen.json','screen_exact'),('stress.json','stress_exact')):
        result=json.loads((ROOT/'results'/name).read_text())
        for path,sha in result['source_hashes'].items():
            if hashlib.sha256((ROOT/path).read_bytes()).hexdigest()!=sha:raise ValueError('source hash mismatch '+path)
        for row in result['cases']:
            data=row['problem'];verify(data,row['result']);counts[key]+=1
            nr=row.get('numeric_result')
            if nr and nr['status']!='unresolved':verify(data,nr);counts['numeric_exact']+=1
            for g in row.get('grids',[]):
                n=g.get('nodes',g.get('nodes_per_species'))
                if g['status']=='exact_grid_witness':atomic(data,g['atoms']);counts['grid_atomic']+=1
                if g.get('farkas'):farkas(data,n,g['farkas']);counts['grid_farkas']+=1
    counts['passed']=True
    (ROOT/'results'/'independent-replay.json').write_text(json.dumps(counts,indent=2)+'\n')
    print(json.dumps(counts))
if __name__=='__main__':main()
