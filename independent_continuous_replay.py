"""Standard-library-only replay of controlled continuous comparator evidence."""
from pathlib import Path
from fractions import Fraction as F
import collections
import copy
import hashlib
import json
import random
import sys
from independent_verifier import verify, InvalidCertificate

ROOT=Path(__file__).resolve().parent

def require(condition,message):
    if not condition: raise RuntimeError(message)

def digest(path):return hashlib.sha256((ROOT/path).read_bytes()).hexdigest()

def canonical_inputs():
    out=[]
    for filename,group in [('screen.json','development_screen'),('stress.json','designed_stress')]:
        rows=json.loads((ROOT/'results'/filename).read_text())['cases']
        for i,r in enumerate(rows):out.append((f'{group}:{i}',group,r['problem']))
    rng=random.Random(271828)
    masses=[F(1,2),F(1),F(3,2),F(2)]
    flux=[F(-1),F(-3,4),F(-1,3),F(0),F(1,3),F(3,4),F(1)]
    real=[F(i,3) for i in range(-9,10)]
    imag=[F(1,100),F(1,20),F(1,5),F(1,2),F(1),F(2)]
    wave=[F(1,4),F(1,2),F(1),F(2)]
    for i in range(128):
        n1,n2=rng.choice(masses),rng.choice(masses)
        p={'N1':str(n1),'F1':str(n1*rng.choice(flux)),
           'N2':str(n2),'F2':str(n2*rng.choice(flux)),
           'w':{'re':str(rng.choice(real)),'im':str(rng.choice(imag))},
           'k':str(rng.choice(wave)),'sector':rng.choice(['longitudinal','transverse'])}
        out.append((f'fresh_synthetic:{i}','fresh_synthetic',p))
    return out

def main():
    data=json.loads((ROOT/'results/continuous-baseline.json').read_text())
    require(digest('baseline_continuous_exploration.py')==data['source_sha256'],'benchmark source differs')
    known_plan_versions={digest('independent-continuous-review.plan.md'),digest('independent-continuous-review.md')}
    require(data['plan_sha256'] in known_plan_versions,'recorded plan/report version differs')
    for f,h in data['input_sha256'].items():require(digest('results/'+f)==h,'input artifact differs')
    expected=canonical_inputs();require(len(expected)==722==len(data['cases']),'wrong case count')
    counts=collections.Counter();unresolved=collections.Counter()
    for row,(name,group,problem) in zip(data['cases'],expected):
        require((row['id'],row['group'],row['problem'])==(name,group,problem),'canonical case mismatch')
        for method in ['exact','conic']:
            result=row[method]
            if result['status']=='unresolved':
                unresolved[(method,group)]+=1
                require(method=='conic','unexpected unresolved exact result')
                continue
            require(verify(problem,result) is True,'proof refused')
            counts[method]+=1
        require(row['conic']['status']=='unresolved' or row['conic']['status']==row['exact']['status'],'certified disagreement')
        for method,ts in row['generation_seconds'].items():
            require(len(ts)==3 and all(isinstance(t,(int,float)) and t>=0 for t in ts),'invalid timing samples')
    require(counts=={'exact':722,'conic':682},'unexpected certificate counts')
    require(unresolved=={('conic','designed_stress'):40},'unexpected unresolved strata')
    # Refusal controls remain active with python -O.
    positive=next(r for r in data['cases'] if r['conic']['status']=='possible')
    bad=copy.deepcopy(positive['conic']);bad['s1']['re']='999999'
    negative=next(r for r in data['cases'] if r['conic'].get('proof_kind')=='rational_dual')
    bad_dual=copy.deepcopy(negative['conic']);bad_dual['dual']['negative_constant']='0'
    controls=[(positive['problem'],bad),(negative['problem'],bad_dual),
              (positive['problem'],{'status':'impossible','proof_kind':'four_disks'})]
    refused=0
    for p,r in controls:
        try:verify(p,r)
        except InvalidCertificate:refused+=1
        else:raise RuntimeError('false certificate accepted')
    forbidden=[n for n in ['numpy','scipy','clarabel','moments','geometry'] if n in sys.modules]
    require(not forbidden,'numerical or production module imported')
    report={'passed':True,'canonical_cases':722,'exact_proofs':722,'conic_proofs':682,
        'unresolved_conic_stress_cases':40,'mutated_controls_rejected':refused,
        'standard_library_only':True,'python_optimized':not __debug__,
        'data_sha256':digest('results/continuous-baseline.json'),
        'checker_sha256':digest('independent_verifier.py'),
        'replay_sha256':digest('independent_continuous_replay.py')}
    (ROOT/'results/continuous-replay.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
if __name__=='__main__':main()
