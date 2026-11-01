"""Bounded nonaxial tetrahedron falsifier; no universal extension claim."""
from pathlib import Path
import hashlib,json,sys
import numpy as np
ROOT=Path(__file__).resolve().parent

def rotation(q):
    w,x,y,z=q/np.linalg.norm(q)
    return np.array([[1-2*(y*y+z*z),2*(x*y-z*w),2*(x*z+y*w)],
                     [2*(x*y+z*w),1-2*(x*x+z*z),2*(y*z-x*w)],
                     [2*(x*z-y*w),2*(y*z+x*w),1-2*(x*x+y*y)]])

def run():
    rng=np.random.default_rng(2026110120)
    tetra=np.array([[1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1]])/np.sqrt(3)
    ks=np.geomspace(.001,8,45);rows=[];pairs=[];best={'ratio':0.}
    moment_error=0.;polar_error=0.;polar_count=0
    for case in range(300):
        qs=rng.integers(-9,10,(2,4));qs[~qs.any(axis=1),0]=1;pairs.append(qs.tolist())
        v=np.r_[tetra@rotation(qs[0]).T,tetra@rotation(qs[1]).T]
        for points in (v[:4],v[4:]):
            moment_error=max(moment_error,float(np.max(np.abs(points.mean(axis=0)))),float(np.max(np.abs(points.T@points/4-np.eye(3)/3))),float(np.max(np.abs(np.sum(points*points,axis=1)-1))))
        g=np.r_[np.full(4,.25),np.full(4,-.25)];coupling=-(1-v@v.T)*g[None,:]
        V=np.c_[np.ones(8),v]
        for k in ks:
            A=np.diag(k*v[:,2])+coupling
            ev=np.linalg.eigvals(A);w=ev[np.argmax(ev.imag)];gamma=float(w.imag)
            cap=np.sqrt(4*k/9-k*k/9) if k<2 else 2/3
            row={'case':case,'k':float(k),'real_frequency':float(w.real),'growth':gamma,'axial_bound':float(cap),'ratio':float(gamma/cap)}
            if gamma>1e-8:
                P=np.diag([1.,-1.,-1.,-1.]).astype(complex)
                for weight,velocity,fourvector in zip(g,v,V):P+=weight*np.outer(fourvector,fourvector)/(w-k*velocity[2])
                singular=np.linalg.svd(P,compute_uv=False);residual=float(singular[-1]/singular[0])
                row['polarization_relative_singular_value']=residual;polar_error=max(polar_error,residual);polar_count+=1
            rows.append(row)
            if row['ratio']>best['ratio']:best=dict(row,quaternions=qs.tolist())
    if moment_error>1e-12 or polar_error>1e-7:raise RuntimeError('Verification tolerance failed')
    result={'scope':'300rotatedtetrahedron pairs x45k; absenceofviolation does not prove a3D bound','cases':len(rows),'growing_roots_checked':polar_count,'maximum_moment_error':moment_error,'maximum_polarization_residual':polar_error,'best':best,'quaternions':pairs,'records':rows,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'plan_sha256':hashlib.sha256((ROOT/'PLAN.md').read_bytes()).hexdigest(),'numpy':np.__version__,'python':sys.version}
    (ROOT/'results/probe.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('quaternions','records','python')},indent=2))
if __name__=='__main__':run()
