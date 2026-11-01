"""Fixed, physically sourced extra-moment compatibility experiment.

All frequencies checked exactly are designated rational-squared-rate points.
Continuous reference roots remain numerical or independently Arb-bracketed,
not rounded exact roots. An inclusion validates compatibility only.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
from m2_frequency import (ComplexBasis, classify, contains, ellipse_candidates,
                          jsonable, sector_target)
from independent_m2_verifier import verify, verify_m1_possible

ROOT = Path(__file__).resolve().parent


def data_at(physical_k, gamma2, sector):
    reduced = physical_k-F(2,3)
    if reduced == 0:
        raise ValueError('special reduced zero mode requires separate treatment')
    axis = 1 if reduced > 0 else -1
    return dict(N1=F(1),F1=axis*F(1,3),Q1=F(1,3),
                N2=F(1),F2=-axis*F(1,3),Q2=F(1,3),
                wr=F(0),gamma2=gamma2,k=abs(reduced),sector=sector)


def m1_proposal(data):
    """Same frequency but only M0/M1 constraints; exact extended-circle path."""
    g,k,x = data['gamma2'],data['k'],data['wr']/data['k']
    sp1 = (data['N1'],data['F1'],data['Q1'])
    sp2 = (data['N2'],data['F2'],data['Q2'])
    target = sector_target(sp1,sp2,data['wr'],g,k,data['sector'])
    if target is None: return dict(status='unresolved')
    ellipses = []
    for n,m,shift in ((sp1[0],sp1[1],ComplexBasis(0,0,g)),
                      (sp2[0],sp2[1],target)):
        for sign in (1,-1):
            a = (x+sign)**2+g/(k*k)
            cx = (n*(2*x+sign)-m)/(2*a)+shift.p
            cq = -k*((1+sign*x)*(n+sign*m)+2*n*g/(k*k))/(2*a*g)+shift.q
            r2 = (n+sign*m)**2*k*k/(4*g*a)
            ellipses.append(((cx,cq),r2))
    for tag,point in ellipse_candidates(ellipses,g):
        if all(contains(point,e,g) for e in ellipses):
            s1 = ComplexBasis(*point,g)
            s2 = s1-target
            result = dict(status='possible',s1=[s1.p,s1.q],s2=[s2.p,s2.q])
            verify_m1_possible(data,result)
            return result
    return dict(status='unresolved',reason='No independently verified M1-positive proposal')


def run():
    reference_path = ROOT/'results/counterflux-modes.json'
    reference = json.loads(reference_path.read_text())
    expected = [str(F(j,8)) for j in range(-8,9)]
    if reference['physical_K_grid'] != expected:
        raise RuntimeError('reference physical grid differs from fixed plan')
    rows = []
    for mode in reference['closure_modes']:
        if not mode['growing']:
            rows.append(dict(kind='quantum_closure',physical_K=mode['physical_K'],
                             sector=mode['sector'],gamma2=mode['gamma_squared'],
                             status='no_growing_mode_from_exact_closure_formula'))
            continue
        data = data_at(F(mode['physical_K']),F(mode['gamma_squared']),mode['sector'])
        result = classify(data)
        verify(data,result)
        low = result if result['status'] == 'possible' else m1_proposal(data)
        if low['status'] == 'possible': verify_m1_possible(data,low)
        rows.append(dict(kind='quantum_closure',physical_K=mode['physical_K'],
                         sector=mode['sector'],data=data,result=result,
                         m1_result=low,independent_m2_check=True,
                         independent_m1_positive_check=low['status']=='possible'))
    # Stable closure blocks are an exact matrix statement, not proof of no
    # possible transverse instability in the original angular problem.
    transverse = reference['quantum_closure_transverse']
    controls = []
    for group in reference['continuous_modes']:
        for branch in group['branches']:
            for root in branch['roots']:
                gamma = F(str(root['gamma']))
                data = data_at(F(group['physical_K']),gamma*gamma,branch['sector'])
                result = classify(data)
                verify(data,result)
                controls.append(dict(completion=group['completion'],
                                     physical_K=group['physical_K'],sector=branch['sector'],
                                     scope='Exact designated rational point near numerical reference; not an exact identity of the reference root',
                                     data=data,result=result,independently_checked=True))
    brackets = []
    for row in reference['validated_linear_root_brackets']:
        checks=[]
        for endpoint in row['gamma_interval']:
            gamma=F(endpoint)
            data=data_at(F(row['physical_K']),gamma*gamma,row['sector'])
            result=classify(data);verify(data,result)
            checks.append(dict(data=data,result=result))
        brackets.append(dict(physical_K=row['physical_K'],sector=row['sector'],
                             interval=row['gamma_interval'],
                             endpoint_frequency_checks=checks,
                             scope='Arb certifies an actual linear-input root between endpoints; endpoint compatibility alone does not certify a whole interval'))
    growing = [r for r in rows if 'result' in r]
    exclusions = [r for r in growing if r['result']['status']=='impossible'
                  and r['m1_result']['status']=='possible']
    files=['M2_COMPATIBILITY_PLAN.md','m2_frequency.py','geometry.py',
           'independent_m2_verifier.py','independent_verifier.py',
           'counterflux_modes.py','m2_counterflux_screen.py']
    output=dict(scope='Prespecified 17 physical K values; published counterflux input; existential compatibility, not a root census or astrophysical validation',
                density_normalization='N±=1',number_moments={'plus':['1','1/3','1/3'],'minus':['1','-1/3','1/3']},
                current_shift='reduced k=physical K−2/3; reflect species fluxes when reduced k<0',
                utility_gate='At least one exact growing quantum-closure mode impossible with supplied M2 but independently possible with M0/M1',
                utility_gate_passed=bool(exclusions),extra_information_exclusions=len(exclusions),
                closure_growing_modes=len(growing),
                closure_m2_possible=sum(r['result']['status']=='possible' for r in growing),
                closure_m2_impossible=sum(r['result']['status']=='impossible' for r in growing),
                closure_modes=rows,quantum_closure_transverse=transverse,
                continuous_designated_point_controls=controls,validated_linear_bracket_endpoint_checks=brackets,
                reference_sha256=hashlib.sha256(reference_path.read_bytes()).hexdigest(),
                source_hashes={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in files},
                conclusion='No added-information closure exclusion in the fixed test; utility gate failed. Compatibility does not validate actual closure growth.')
    path=ROOT/'results/m2-counterflux-screen.json'
    path.write_text(json.dumps(jsonable(output),indent=2)+'\n')
    print(json.dumps({k:output[k] for k in ['utility_gate_passed','extra_information_exclusions','closure_growing_modes','closure_m2_possible','closure_m2_impossible','conclusion']},indent=2))
    print('Designated-point controls:',len(controls),'Arb-bracket endpoint pairs:',len(brackets))

if __name__=='__main__':run()
