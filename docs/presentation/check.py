"""Replay one retained endpoint witness using exact native arithmetic; stdlib only."""
from pathlib import Path
from fractions import Fraction as F
import sys,json,hashlib,copy
P=Path(__file__).resolve().parent
K=P/'kernel' if (P/'kernel').is_dir() else P.parents[1]
sys.path.insert(0,str(K))
from finite_band_envelope import envelope
from independent_m2_verifier import CQ

def require(ok,message):
    if not ok:raise ValueError(message)

def verify(case,result):
    for name,digest in case['kernel_sha256'].items():
        require(hashlib.sha256((K/name).read_bytes()).hexdigest()==digest,'Changed native kernel: '+name)
    a=case['input'];require(envelope(a['mass'],a['second'],a['kmax'],a['sector'],a['bits'])==result,'Evaluator result changed')
    n,q,k=map(F,(a['mass'],a['second'],a['kmax']));g=F(case['gamma_squared']);require(result['branch']=='endpoint','Not the retained endpoint branch');require(g==F(result['exact_growth_squared']),'Squared growth mismatch')
    require(F(result['growth_lower'])**2<=g<=F(result['growth_upper'])**2,'Invalid rational growth enclosure')
    species=[]
    for name in ['plus','minus']:
        atoms=[(F(x['u']),F(x['weight']))for x in case['species'][name]]
        require(all(-1<=u<=1 and w>=0 for u,w in atoms),'Nonpositive or out-of-support atom')
        require([sum(w*u**j for u,w in atoms)for j in range(3)]==[n,F(0),q],'Moments do not equal (N,0,Q)')
        species.append(atoms)
    omega=CQ(0,1,g);values=[]
    for j in range(3):
        values.append(sum((CQ(sign*w*u**j,0,g)/(omega-k*u)for atoms,sign in zip(species,[1,-1])for u,w in atoms),CQ(0,0,g)))
    i0,i1,i2=values;require(((1+i0)*(-1+i2)-i1*i1).zero(),'Original longitudinal dispersion factor is nonzero')
    # Separate rational endpoint formula; does not use branch implementation.
    alpha=q/n;require(2*n*alpha*(1-alpha)*k-alpha*alpha*k*k==g,'Endpoint formula mismatch')
    return True

def main():
    if (P/'MANIFEST.json').exists():
        for name,digest in json.loads((P/'MANIFEST.json').read_text())['sha256'].items():require(hashlib.sha256((P/name).read_bytes()).hexdigest()==digest,'Bundle byte mismatch: '+name)
    case=json.loads((P/'case.json').read_text());result=json.loads((P/'evaluator.json').read_text());verify(case,result)
    bad=copy.deepcopy(case);bad['species']['plus'][0]['weight']='1/3'
    try:verify(bad,result)
    except ValueError:pass
    else:raise ValueError('Tampered weight was accepted')
    bad=copy.deepcopy(case);bad['gamma_squared']='1/2000'
    try:verify(bad,result)
    except ValueError:pass
    else:raise ValueError('Tampered squared growth was accepted')
    print(json.dumps({'checked':True,'species':2,'atoms':4,'moments_each':['1','0','1/3'],'exact_growth_squared':case['gamma_squared'],'original_longitudinal_factor':'0 exactly','rational_growth_enclosure':True,'tampered_weight_rejected':True,'tampered_growth_rejected':True,'scope':'One retained native endpoint witness; not a new theorem or physical validation'},indent=2))
if __name__=='__main__':main()
