"""Exact nonaxisymmetric instability erased by azimuthal averaging.

Directions=(sqrt(2/3)*X,sqrt(2/3)*Y,Z/sqrt(3)), k=1/sqrt(3),
omega=i/3. All operator entries after common coordinate scaling are rational.
Known symmetry-breaking physics; no new instability or bound violation claim.
"""
from fractions import Fraction as F
from pathlib import Path
import sys,json,hashlib
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from independent_verifier import C

def require(ok,message):
    if not ok:raise RuntimeError(message)

def run():
    plus=[(1,0,1),(-1,0,1),(0,1,-1),(0,-1,-1)]
    minus=[(x,y,-z) for x,y,z in plus]
    for points in (plus,minus):
        require(all(2*x*x+2*y*y+z*z==3 for x,y,z in points),'unit sphere')
        require(all(sum(p[j] for p in points)==0 for j in range(3)),'zero flux')
        require(all(sum(p[i]*p[j] for p in points)==0 for i in range(3) for j in range(i)),'cross moments')
        require([F(sum(p[j]**2 for p in points),4)*s for j,s in enumerate((F(2,3),F(2,3),F(1,3)))]==[F(1,3)]*3,'isotropic second tensor')
    require(sorted(p[2] for p in plus)==sorted(p[2] for p in minus),'equal azimuthal projections')
    v=plus+minus;g=[F(1,4)]*4+[F(-1,4)]*4;omega=C(0,F(1,3))
    eigenvector=[C(p[1])/(omega-F(p[2],3)) for p in v]
    require(any(q!=C() for q in eigenvector),'nonzero eigenvector')
    for i,p in enumerate(v):
        value=C(F(p[2],3))*eigenvector[i]
        for j,q in enumerate(v):
            dot=F(2*p[0]*q[0]+2*p[1]*q[1]+p[2]*q[2],3)
            value-=(1-dot)*g[j]*eigenvector[j]
        require(value==omega*eigenvector[i],'original eight-angle eigen equation')
    # Off-diagonal polarization moments vanish before coordinate radicals.
    vectors=[(1,*p) for p in v]
    for i in range(4):
        for j in range(i):
            require(sum((g[s]*vectors[s][i]*vectors[s][j]/(omega-F(v[s][2],3)) for s in range(8)),C())==C(),'offdiagonal polarization')
    diagonals=[]
    for j,(metric,scale) in enumerate(zip((1,-1,-1,-1),(F(1),F(2,3),F(2,3),F(1,3)))):
        diagonals.append(C(metric)+sum((g[s]*scale*vectors[s][j]**2/(omega-F(v[s][2],3)) for s in range(8)),C()))
    require(diagonals==list(map(C,(1,-2,0,-1))),'original polarization zero')
    result={'passed':True,'species_mass':'1','species_flux':['0','0','0'],'species_second_tensor':'I/3','directions_encoded':{'plus':plus,'minus':minus,'coordinate_scale_squared':['2/3','2/3','1/3']},'each_atom_weight':'1/4','wave_number_squared':'1/3','frequency_real':'0','growth':'1/3','polarization_diagonal':['1','-2','0','-1'],'azimuthal_projection_signed_measure_zero':True,'original_eigen_equations_exact':8,'scope':'Explicit exact projection-loss example; not a counterexample to the axial worst-case envelope or a novel symmetry-breaking mechanism','source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (Path(__file__).resolve().parent/'results/exact-projection.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':run()
