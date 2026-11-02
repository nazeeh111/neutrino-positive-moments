"""Exact algebra checks for the conditional QME quadrature obstruction.

Synthetic rational h fixtures test the general proof, not exponential integrals
or astrophysical data. Primary implementation, not an independent review.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from independent_verifier import C
from m2_frequency import classify, jsonable
from independent_m2_verifier import verify


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def moments(atoms, count):
    return [sum((w*u**n for u, w in atoms), F(0)) for n in range(count)]


def multiply(a, b):
    return [[sum(a[i][t]*b[t][j] for t in range(len(b)))
             for j in range(len(b[0]))] for i in range(len(a))]


def sector(c, a, d, mass, aa, bb):
    require(0 < a < d and abs(c)+d <= 1, 'Invalid h support')
    edge = mass*a*a/(2*d*d)
    h = [(c-d, edge), (c, mass-2*edge), (c+d, edge)]
    nodes = [c-a, c+a]
    quad = [(u, mass/2) for u in nodes]
    ii = moments(h, 4)
    require(ii == moments(quad, 4), 'Gaussian h moments differ')
    dd = ii[0]*ii[2]-ii[1]**2
    require(dd > 0, 'Degenerate h Gram matrix')
    kap = (ii[2]**2-ii[1]*ii[3])/dd
    eta = (ii[0]*ii[3]-ii[1]*ii[2])/dd
    require(all(u*u == kap+eta*u for u in nodes), 'Wrong nodes')
    source = [(u, w*(aa+bb*u)) for u, w in quad]
    g = moments([(u, w*(aa+bb*u)) for u, w in h], 3)
    require(moments(source, 3) == g, 'Source moments differ')
    require(g[2] == kap*g[0]+eta*g[1], 'Background identity fails')
    return {'h': h, 'source': source, 'g': g, 'kappa': kap, 'eta': eta}


def operators(ss, k):
    atoms = ss[0]['source']+ss[1]['source']
    signs = [1, 1, -1, -1]
    ll = [[(k*u if i == j else F(0))-w*(1-u*v)*signs[j]
           for j, (v, _) in enumerate(atoms)]
          for i, (u, w) in enumerate(atoms)]
    tt = [[F(0) for _ in range(4)] for _ in range(4)]
    mm = [[F(0) for _ in range(4)] for _ in range(4)]
    for s, data in enumerate(ss):
        n, f, q = data['g']
        for j, (u, _) in enumerate(data['source']):
            tt[2*s][2*s+j] = F(1)
            tt[2*s+1][2*s+j] = u
        mm[2*s][2*s+1] += k
        mm[2*s+1][2*s] += k*data['kappa']
        mm[2*s+1][2*s+1] += k*data['eta']
        for t, sig in enumerate((1, -1)):
            mm[2*s][2*t] -= sig*n
            mm[2*s][2*t+1] += sig*f
            mm[2*s+1][2*t] -= sig*f
            mm[2*s+1][2*t+1] += sig*q
    return ll, tt, mm


def run():
    nonzero_eta = 0
    heavy = [(F(-1, 4), F(2, 7)), (F(1, 3), F(3, 11))]
    for j in range(96):
        ss = [sector(F((j+2*s) % 5-2, 10), F(1+(j+s) % 3, 10),
                     F(1, 2), F(3+(j+3*s) % 7, 5),
                     F(2+(j+s) % 4), F((j+2*s) % 5-2, 2))
              for s in range(2)]
        require(all(w >= 0 for z in ss for _, w in z['source']),
                'Ordered-flavor fixture is not positive')
        ll, tt, mm = operators(ss, F(j % 13-6, 7))
        require(multiply(mm, tt) == multiply(tt, ll), 'Similarity fails')
        nonzero_eta += int(any(z['eta'] != 0 for z in ss))
        hm = moments(heavy, 3)
        for z in ss:
            require(moments(heavy+z['source'], 3) ==
                    [hm[i]+z['g'][i] for i in range(3)],
                    'Electron moments not preserved')
        require([x-y for x, y in zip(moments(heavy+ss[0]['source'], 3),
                                    moments(heavy+ss[1]['source'], 3))] ==
                [ss[0]['g'][i]-ss[1]['g'][i] for i in range(3)],
                'Common component failed to cancel')

    # Positive h does not ensure positive source quadrature.
    bad = sector(F(0), F(3, 5), F(4, 5), F(1), F(0), F(1))
    require(all(w > 0 for _, w in bad['h']) and
            any(w < 0 for _, w in bad['source']), 'Sign control failed')

    ss = [sector(F(0), F(3, 5), F(4, 5), F(1), F(2), F(sign))
          for sign in (1, -1)]
    k = F(16, 25)
    omega = C(0, F(102, 125))
    ii = [sum((sig*w*u**n/(omega-k*u)
               for z, sig in zip(ss, (1, -1)) for u, w in z['source']), C())
          for n in range(3)]
    require((1+ii[0])*(-1+ii[2])-ii[1]*ii[1] == C(),
            'Original growing dispersion root failed')
    data = {'N1': '2', 'F1': '9/25', 'Q1': '18/25',
            'N2': '2', 'F2': '-9/25', 'Q2': '18/25',
            'wr': '0', 'gamma2': '10404/15625', 'k': '16/25',
            'sector': 'longitudinal'}
    proof = jsonable(classify(data))
    require(proof['status'] == 'possible', 'Compatible root rejected')
    require(verify(data, proof), 'Separate M2 verifier failed')
    out = {'passed': True, 'exact_similarity_cases': 96,
           'cases_with_nonzero_eta': nonzero_eta,
           'common_component_checks': 96,
           'negative_source_weight_control': True,
           'original_growing_root_exact_zero': True,
           'root_input': data, 'root_proof': proof,
           'scope': 'Conditional longitudinal quadrature obstruction; synthetic h algebra fixtures, not physical data or a utility pass; primary checked only',
           'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    Path(__file__).with_name('qme-quadrature-gate.json').write_text(
        json.dumps(out, indent=2)+'\n')
    print(json.dumps({k: v for k, v in out.items() if k != 'root_proof'}, indent=2))


if __name__ == '__main__':
    run()
