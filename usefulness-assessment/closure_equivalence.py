"""Exact conditional scalar-closure/quadrature equivalence.

This does not reproduce the full tensor closure or nonlinear advection in FLASH.
It shows why a consistent fixed scalar closure cannot be excluded by the same
positive moment class. Finite checks corroborate the matrix identity proof.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json,random,sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from independent_verifier import C
from m2_frequency import classify,jsonable
from independent_m2_verifier import verify

def require(ok,message):
    if not ok:raise RuntimeError(message)

def multiply(a,b):return [[sum(a[i][t]*b[t][j] for t in range(len(b))) for j in range(len(b[0]))] for i in range(len(a))]

def matrices(n,f,b,k):
    u=[b[0],-b[0],b[1],-b[1]]
    weights=[(n[s]+sign*f[s]/b[s])/2 for s in range(2) for sign in (1,-1)]
    signs=[1,1,-1,-1]
    A=[[(k*u[i] if i==j else F(0))-weights[i]*(1-u[i]*u[j])*signs[j] for j in range(4)] for i in range(4)]
    T=[[F(0) for _ in range(4)] for _ in range(4)]
    M=[[F(0) for _ in range(4)] for _ in range(4)]
    for s in range(2):
        T[2*s][2*s]=T[2*s][2*s+1]=F(1)
        T[2*s+1][2*s]=b[s];T[2*s+1][2*s+1]=-b[s]
        M[2*s][2*s+1]+=k;M[2*s+1][2*s]+=k*b[s]**2
        for t,sigma in enumerate((1,-1)):
            M[2*s][2*t]-=n[s]*sigma;M[2*s][2*t+1]+=f[s]*sigma
            M[2*s+1][2*t]-=f[s]*sigma;M[2*s+1][2*t+1]+=n[s]*b[s]**2*sigma
    return A,T,M,u,weights

def run():
    rng=random.Random(2026110123);boundary=0
    for _ in range(200):
        n=[F(rng.randint(1,15),7) for s in range(2)]
        b=[F(rng.randint(1,10),10) for s in range(2)]
        f=[n[s]*b[s]*F(rng.randint(-5,5),5) for s in range(2)]
        k=F(rng.randint(-10,10),9)
        A,T,M,u,w=matrices(n,f,b,k)
        require(multiply(T,A)==multiply(M,T),'Moment and quadrature operators differ')
        for s in range(2):
            at=list(zip(u[2*s:2*s+2],w[2*s:2*s+2]))
            require(all(-1<=x<=1 and z>=0 for x,z in at),'Nonpositive quadrature')
            require([sum(z*x**j for x,z in at) for j in range(3)]==[n[s],f[s],n[s]*b[s]**2],'Wrong moments')
        boundary+=int(any(z==0 for z in w))
    data={'N1':'1','F1':'1/3','Q1':'1/4','N2':'1','F2':'-1/3','Q2':'1/4','wr':'0','gamma2':'25/36','k':'1','sector':'longitudinal'}
    proof=jsonable(classify(data));require(proof['status']=='possible','Known quadrature root rejected');require(verify(data,proof),'Independent M2 checker failed')
    plus=[(F(1,2),F(5,6)),(F(-1,2),F(1,6))];minus=[(-u,m) for u,m in plus];omega=C(0,F(5,6))
    I=[sum((sign*m*u**j/(omega-u) for atoms,sign in ((plus,1),(minus,-1)) for u,m in atoms),C()) for j in range(3)]
    require((1+I[0])*(-1+I[2])-I[1]*I[1]==C(),'Original dispersion factor nonzero')
    out={'passed':True,'operator_identity_cases':200,'cases_with_zero_weight':boundary,'example':data,'example_proof':proof,'original_longitudinal_factor_exact_zero':True,'scope':'Constant scalar longitudinal closures withQ_s=chi_sN_s anddeltaP_s=chi_sdeltaN_s; not the general FLASHtensor/advectionoperator; not a new utility pass','source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (Path(__file__).resolve().parent/'closure-equivalence.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k!='example_proof'},indent=2))
if __name__=='__main__':run()
