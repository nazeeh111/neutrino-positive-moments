"""Independent standard-library replay of the fixed counterflux experiment."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
from independent_m2_verifier import verify, verify_m1_possible
from independent_verifier import require, InvalidCertificate

ROOT=Path(__file__).resolve().parent


def check_input(data, physical_k):
    k=F(physical_k)-F(2,3)
    require(k != 0,'special reduced zero mode is not on this grid')
    sign=1 if k>0 else -1
    require(F(data['k']) == abs(k),'wrong current shift or reflection')
    require(F(data['wr']) == 0,'reference mode is purely imaginary')
    for suffix,flux in (('1',sign*F(1,3)),('2',-sign*F(1,3))):
        require(F(data['N'+suffix])==1 and F(data['F'+suffix])==flux
                and F(data['Q'+suffix])==F(1,3),'changed physical moment class')


def run():
    path=ROOT/'results/m2-counterflux-screen.json'
    d=json.loads(path.read_text())
    require([r['physical_K'] for r in d['closure_modes']]==[str(F(j,8)) for j in range(-8,9)],'fixed grid changed')
    counts={'m2_proofs':0,'m1_positive_proofs':0,'exact_closure_modes':0,
            'designated_near_reference_points':0,'bracket_endpoint_checks':0,
            'extra_information_exclusions':0}
    for row in d['closure_modes']:
        physical=F(row['physical_K']);k=physical-F(2,3)
        g2=(k+F(2,3))*(F(2,3)-F(17,45)*k)
        if g2<=0:
            require(row['status']=='no_growing_mode_from_exact_closure_formula' and F(row['gamma2'])==g2,'incorrect stable closure flag')
            continue
        check_input(row['data'],physical)
        require(F(row['data']['gamma2'])==g2 and row['data']['sector']=='longitudinal','incorrect exact closure frequency')
        verify(row['data'],row['result']);counts['m2_proofs']+=1
        counts['exact_closure_modes']+=1
        if row['m1_result']['status']=='possible':
            verify_m1_possible(row['data'],row['m1_result']);counts['m1_positive_proofs']+=1
            counts['extra_information_exclusions']+=row['result']['status']=='impossible'
    for row in d['continuous_designated_point_controls']:
        check_input(row['data'],row['physical_K']);verify(row['data'],row['result'])
        counts['m2_proofs']+=1;counts['designated_near_reference_points']+=1
    for row in d['validated_linear_bracket_endpoint_checks']:
        require(len(row['interval'])==len(row['endpoint_frequency_checks'])==2,'missing bracket endpoints')
        require(0<F(row['interval'][0])<F(row['interval'][1]),'invalid gamma bracket')
        for endpoint,entry in zip(row['interval'],row['endpoint_frequency_checks']):
            check_input(entry['data'],row['physical_K'])
            require(F(entry['data']['gamma2'])==F(endpoint)**2,'endpoint point mismatch')
            verify(entry['data'],entry['result']);counts['m2_proofs']+=1
            counts['bracket_endpoint_checks']+=1
    matrix=[[F(x) for x in row] for row in d['quantum_closure_transverse']['matrix']]
    require(matrix==[[F(14,45),-F(14,45)],[F(14,45),-F(14,45)]],'wrong transverse block')
    require(all(sum(matrix[i][h]*matrix[h][j] for h in (0,1))==0 for i in (0,1) for j in (0,1)),'transverse block not nilpotent')
    require(d['utility_gate_passed']==(counts['extra_information_exclusions']>0),'incorrect utility-gate flag')
    require(counts['extra_information_exclusions']==d['extra_information_exclusions'],'incorrect exclusion count')
    require(len(d['validated_linear_bracket_endpoint_checks'])==13,'lost linear controls')
    false=next(r for r in d['closure_modes'] if 'result' in r)
    rejected=False
    try:verify(false['data'],{'status':'impossible','proof_kind':'m2_four_disks'})
    except InvalidCertificate:rejected=True
    require(rejected,'false impossibility control accepted')
    result=dict(passed=True,standard_library_only=True,optimized_python=not __debug__,
                counts=counts,false_certificate_rejected=rejected,
                utility_gate_passed=d['utility_gate_passed'],
                data_sha256=hashlib.sha256(path.read_bytes()).hexdigest())
    (ROOT/'results/m2-screen-replay.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':run()
