"""Angular-grid LP inner approximation and exact zero-mode moment baseline."""
from fractions import Fraction as F
import numpy as np
from scipy.optimize import linprog
from moments import QC,root_target

def rref_solve(A,b):
    rows=[list(row)+[v] for row,v in zip(A,b)];n=len(A[0]) if A else 0
    pivots=[];rank=0
    for j in range(n):
        pivot=next((i for i in range(rank,len(rows)) if rows[i][j]),None)
        if pivot is None:continue
        rows[rank],rows[pivot]=rows[pivot],rows[rank]
        div=rows[rank][j];rows[rank]=[x/div for x in rows[rank]]
        for i in range(len(rows)):
            if i!=rank and rows[i][j]:
                c=rows[i][j];rows[i]=[x-c*y for x,y in zip(rows[i],rows[rank])]
        pivots.append(j);rank+=1
    if any(all(x==0 for x in row[:-1]) and row[-1]!=0 for row in rows):return None
    x=[F(0)]*n
    for i,j in enumerate(pivots):x[j]=rows[i][-1]
    return x

def grid_lp(N1,F1,N2,F2,w,k,sector,n):
    N1,F1,N2,F2,k=map(F,(N1,F1,N2,F2,k));target=root_target(N1,F1,N2,F2,w,k,sector)
    if target is None:return {'status':'structurally_impossible'}
    z=w/k;u=np.linspace(-1,1,n);f=1/(complex(z)-u)
    A=np.zeros((6,2*n));A[0,:n]=1;A[1,:n]=u;A[2,n:]=1;A[3,n:]=u
    A[4,:n]=f.real;A[4,n:]=-f.real;A[5,:n]=f.imag;A[5,n:]=-f.imag
    b=[N1,F1,N2,F2,target.re,target.im]
    result=linprog(np.zeros(2*n),A_eq=A,b_eq=np.array(list(map(float,b))),
                   bounds=(0,None),method='highs',options={'primal_feasibility_tolerance':1e-9})
    if not result.success:return {'status':'grid_infeasible_or_solver_failure','solver_status':int(result.status)}
    for tol in (1e-8,1e-12,0.):
        active=np.flatnonzero(result.x>tol)
        if len(active)>12:continue
        columns=[];nodes=[]
        for j in active:
            v=F(2*int(j%n),n-1)-1;s=1/(z-v);first=bool(j<n);sign=1 if first else -1
            columns.append([F(first),v if first else F(0),F(not first),v if not first else F(0),sign*s.re,sign*s.im])
            nodes.append((0 if first else 1,v))
        if not columns:continue
        exact=rref_solve(list(map(list,zip(*columns))),b)
        if exact is None or min(exact)<0:continue
        atoms=[[],[]]
        for (species,node),mass in zip(nodes,exact):
            if mass:atoms[species].append((node,mass))
        # Direct independent equations, not the LP residual/tolerance.
        if any(sum((mass for node,mass in at),F(0))!=N or
               sum((node*mass for node,mass in at),F(0))!=flux
               for at,N,flux in zip(atoms,(N1,N2),(F1,F2))):continue
        s=[sum((mass/(z-node) for node,mass in at),QC()) for at in atoms]
        if s[0]-s[1]!=target:continue
        return {'status':'exact_grid_witness','atoms':atoms,'active_atoms':sum(map(len,atoms))}
    return {'status':'grid_feasible_unvalidated','solver_status':int(result.status)}

def zero_mode(N1,F1,N2,F2):
    N1,F1,N2,F2=map(F,(N1,F1,N2,F2));g0=N1-N2;g1=F1-F2
    lo=g0+(F1*F1/N1 if N1 else 0)-N2
    hi=g0+N1-(F2*F2/N2 if N2 else 0)
    minabs=F(0) if lo<=0<=hi else min(abs(lo),abs(hi));maxabs=max(abs(lo),abs(hi))
    low=max(F(0),g1*g1-maxabs*maxabs/4);high=max(F(0),g1*g1-minabs*minabs/4)
    return {'minimum_growth_squared':low,'maximum_growth_squared':high,
            'robust_zero_mode_instability':low>0,
            'degree1_crossing_test':(g0+g1)*(g0-g1)<0}
