from pathlib import Path
import numpy as np,json,hashlib
rng=np.random.default_rng(2026110120)
tet=np.array([[1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1]])/np.sqrt(3)
def rotation(q):
 w,x,y,z=q/np.linalg.norm(q)
 return np.array([[1-2*(y*y+z*z),2*(x*y-z*w),2*(x*z+y*w)],[2*(x*y+z*w),1-2*(x*x+z*z),2*(y*z-x*w)],[2*(x*z-y*w),2*(y*z+x*w),1-2*(x*x+y*y)]])
ks=np.geomspace(.001,8,45);best={'ratio':0};rows=[]
for i in range(300):
 qs=rng.integers(-9,10,(2,4));qs[~qs.any(axis=1),0]=1
 v=np.r_[tet@rotation(qs[0]).T,tet@rotation(qs[1]).T]
 g=np.r_[np.full(4,.25),np.full(4,-.25)]
 interaction=-(1-v@v.T)*g[None,:]
 for k in ks:
  matrix=np.diag(k*v[:,2])+interaction
  vals=np.linalg.eigvals(matrix);index=np.argmax(vals.imag);w=vals[index];gamma=float(w.imag)
  cap=np.sqrt(4*k/9-k*k/9) if k<2 else 2/3
  if gamma/cap>best['ratio']:best={'ratio':gamma/cap,'quaternions':qs.tolist(),'k':float(k),'omega':[float(w.real),gamma],'axial_cap':float(cap),'case':i}
print(json.dumps(best,indent=2))
Path('/private/tmp/nonaxial-probe-result.json').write_text(json.dumps({'scope':'Exploratory300 rotatedtetrahedron pairs x45positivek; eigenvalues only, not a universal proof','best':best,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2)+'\n')
