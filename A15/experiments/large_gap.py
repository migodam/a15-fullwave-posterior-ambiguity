import sys,json
from pathlib import Path
R=Path(__file__).resolve().parents[1];sys.path.insert(0,str(R/'src'))
import numpy as np
from cylinders import CylinderModel
from physics import Family
c=json.loads((R/'configs/protocol.json').read_text());kw=c['family'].copy();kw['gap_max']=1.5;m=CylinderModel(Family(**kw),order=20);th=np.load(R/'results/discovery/calibration.npz')['theta'][:32];rows=[]
for i,t in enumerate(th):
 for gap in [float(t[0]),.4,.8,1.2]:
  tt=t.copy();tt[0]=gap;a=m.solve(tt,c['k0']);b=m.solve(tt,c['k0'],coupling=False)
  j=m.parameter_jacobian(tt,c['k0']);scale=np.array([.14,6,.079])/np.sqrt(12);A=np.vstack([j['y'][:2].real,j['y'][:2].imag])*scale/c['sigma_primary'];B=j['Q']*scale
  rows.append({'scene':i,'gap':gap,'Q':a['Q'],'Q_uncoupled':b['Q'],'Q_coupling_effect':a['Q']-b['Q'],'relative_external_coupling_effect':float(np.linalg.norm(a['y'][:2]-b['y'][:2])/np.linalg.norm(a['y'][:2])),'local_task_variance_proxy':float(B@np.linalg.solve(np.eye(3)+A.T@A,B))})
(R/'results/large_gap.json').write_text(json.dumps({'rows':rows,'full_forward_solves':m.solve_count,'scope':'physical separated-cylinder boundary check, not a new population reference. Domain is infinite Helmholtz background; translated body can lie outside the original VIE unit box.'},indent=2))
