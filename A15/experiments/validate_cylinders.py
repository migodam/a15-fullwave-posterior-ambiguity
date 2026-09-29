import sys,json,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np
from cylinders import CylinderModel
from physics import Family,pack

def main():
 rng=np.random.default_rng(15004);records=[];start=time.perf_counter()
 for idx in range(24):
  radius=[.10,.12][idx%2];eps=[8.,12.,16.][idx%3];f=Family(radius=radius,epsilon_base=eps)
  th=np.array([rng.uniform(.02,.16),rng.uniform(-.25,.25)*eps,rng.uniform(.001,.08)])
  if idx<4:th[0]=.02001
  k=rng.uniform(3,12);sols={};row={'scene':idx,'theta':th.tolist(),'k':k,'radius':radius,'epsilon':eps,'orders':[]}
  for order in [8,12,16,20,24]:sols[order]=CylinderModel(f,order=order).solve(th,k)
  ref=sols[24]
  for order in [8,12,16,20]:
   err=sols[order]['y']-ref['y'];sigma=.01*np.sqrt(np.mean(pack(ref['y'][:2])**2))
   row['orders'].append({'order':order,'whitened_rms':float(np.sqrt(np.mean(pack(err)**2))/sigma),'Q_abs':abs(sols[order]['Q']-ref['Q'])})
  mod=CylinderModel(f,order=20);row['flux']=mod.flux_diagnostics(th,k)
  # Modal coefficients are compared in a fixed family at each point. E boundary from two sides.
  sol=mod.solve(th,k,fields=True);jumps=[]
  for body in range(2):
   for ang in [.0,.7,1.9,3.1]:
    v=np.array([np.cos(ang),np.sin(ang)]);p=sol['state']['centers'][body]+radius*v;d=1e-7
    vals=mod.field(sol,np.array([p-d*v,p+d*v]));jumps.append(float(abs(vals[0]-vals[1])/max(abs(vals[0]),abs(vals[1]),1e-9)))
  row['boundary_value_jump_max']=max(jumps)
  if idx<8:
   jac=mod.parameter_jacobian(th,k);tangent=[]
   for d in range(8):
    v=rng.normal(size=3)*np.array([.02,1.,.01]);c=rng.normal(size=18)+1j*rng.normal(size=18)
    adj=abs(np.real(np.vdot(c,jac['y']@v))-v@np.real(jac['y'].conj().T@c))
    errors=[]
    for step in [1e-3,1e-4,1e-5,1e-6]:
     if th[0]-abs(step*v[0])<.02:continue
     a=mod.solve(th+step*v,k);b=mod.solve(th-step*v,k);dy=(a['y']-b['y'])/(2*step);dq=(a['Q']-b['Q'])/(2*step)
     errors.append({'step':step,'y':float(np.linalg.norm(dy-jac['y']@v)/max(np.linalg.norm(jac['y']@v),1e-12)),
                    'Q':float(abs(dq-jac['Q']@v)/max(abs(jac['Q']@v),1e-10))})
    tangent.append({'direction':d,'adjoint_abs':float(adj),'errors':errors})
   row['tangents']=tangent
  records.append(row)
 errors=[r['orders'][-1] for r in records]
 out={'scenes':records,'seconds':time.perf_counter()-start,'summary':{'order20_vs24_whitened_rms_max':max(e['whitened_rms'] for e in errors),'order20_vs24_Q_max':max(e['Q_abs'] for e in errors),'flux_balance_max':max(r['flux']['balance_relative'] for r in records),'boundary_value_jump_max':max(r['boundary_value_jump_max'] for r in records)},'scope':'physical full-wave numerical validation only'}
 (ROOT/'results/cylinder_validation.json').write_text(json.dumps(out,indent=2));print(json.dumps(out['summary'],indent=2))
if __name__=='__main__':main()
