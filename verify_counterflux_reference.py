"""Independent exact closure polynomial and Arb linear-root bracket replay.

Uses original recurrence factors and one closed-form angular integral, not
counterflux_modes.py, NumPy/SciPy or a numerical angular grid.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json
from itertools import permutations
from flint import arb,ctx,fmpq

ROOT=Path(__file__).resolve().parent


def determinant(matrix):
    value=F(0)
    for p in permutations(range(len(matrix))):
        term=F((-1)**sum(p[i]>p[j] for i in range(len(p)) for j in range(i+1,len(p))))
        for i,j in enumerate(p):term*=matrix[i][j]
        value+=term
    return value


def ball(q):
    q=F(q);return arb(fmpq(q.numerator,q.denominator))


def factor(gamma,k,sector):
    g,kb=ball(gamma),ball(k)
    # integral u/(i*gamma-k*u) du on [-1,1], real by oddness.
    j=-2/kb+2*g/(kb*kb)*(kb/g).atan()
    g1=ball(F(2,3))
    if sector=='longitudinal':
        return -1-g1/kb-j*(1+g*g/(kb*kb)+g1/kb)
    return -1+(j*(1+g*g/(kb*kb))+g1/kb)/2


def run():
    ctx.prec=192
    path=ROOT/'results/counterflux-modes.json';d=json.loads(path.read_text())
    exact=0
    for row in d['closure_modes']:
        k=F(row['physical_K'])-F(2,3);f,c=F(1,3),F(17,45)
        matrix=[[-F(1),k+f,F(1),-f],[c*k-f,c,f,-c],
                [-F(1),-f,F(1),k+f],[f,c,c*k-f,-c]]
        g2=(k+2*f)*(2*f-c*k)
        if F(row['gamma_squared'])!=g2:raise ArithmeticError('closure formula mismatch')
        # Both sides are degree-four polynomials; five exact values prove identity.
        for w in map(F,(-2,-1,0,1,2)):
            a=[[w*(i==j)-matrix[i][j] for j in range(4)] for i in range(4)]
            if determinant(a)!=(w*w-c*k*k)*(w*w+g2):raise ArithmeticError('closure characteristic polynomial mismatch')
        exact+=1
    brackets=[]
    for row in d['validated_linear_root_brackets']:
        lo,hi=map(F,row['gamma_interval']);k=F(row['reduced_k'])
        if not 0<lo<hi or k==0:raise ValueError('bad bracket domain')
        vl,vh=factor(lo,k,row['sector']),factor(hi,k,row['sector'])
        if not ((vl>0 and vh<0) or(vl<0 and vh>0)):raise ArithmeticError('independent Arb sign check failed')
        brackets.append(dict(physical_K=row['physical_K'],sector=row['sector'],gamma_interval=row['gamma_interval'],independent_endpoint_values=[str(vl),str(vh)]))
    result=dict(passed=True,precision_bits=ctx.prec,exact_closure_characteristic_identities=exact,
                validated_linear_root_brackets=len(brackets),brackets=brackets,
                source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                reference_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                scope='Exact closure polynomial identities and IVT existence for13linear-input roots; no root uniqueness, ME validation or full spectral census')
    (ROOT/'results/counterflux-reference-replay.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='brackets'},indent=2))

if __name__=='__main__':run()
