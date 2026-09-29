import sys,time,json,platform,os
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np
from scipy.linalg import solve
from scipy.sparse.linalg import gmres
from cylinders import CylinderModel
from physics import Family
cfg=json.loads((ROOT/"configs/protocol.json").read_text());m=CylinderModel(Family(**cfg["family"]),order=20);s=m.components([.083,.7,.023],cfg["k0"]);a=s['M'];b=s['rhs'];ref=solve(a,b);out={'environment':platform.platform(),'matrix_size':len(b),'results':[]}
for name in ['scipy_dense','scipy_gmres','torch_cuda_single','torch_cuda_batch32']:
 try:
  if name.startswith('torch'):
   import torch
   aa=torch.as_tensor(a,device='cuda');bb=torch.as_tensor(b,device='cuda')
   if name.endswith('32'):aa=aa.repeat(32,1,1);bb=bb.repeat(32,1)
   torch.cuda.synchronize()
  t=time.perf_counter()
  for _ in range(100):
   if name=='scipy_dense':x=solve(a,b,check_finite=False)
   elif name=='scipy_gmres':x,info=gmres(a,b,rtol=1e-11,atol=1e-13,restart=50);assert info==0
   else:x=torch.linalg.solve(aa,bb)
  if name.startswith('torch'):torch.cuda.synchronize();x=x.cpu().numpy();x=x[0] if x.ndim==2 else x
  elapsed=(time.perf_counter()-t)/100/(32 if name.endswith('32') else 1)
  out['results'].append({'name':name,'seconds_per_solve':elapsed,'relative_error':float(np.linalg.norm(x-ref)/np.linalg.norm(ref))})
 except Exception as e:out['results'].append({'name':name,'error':str(e)})
t=time.perf_counter()
for _ in range(200):m.solve([.083,.7,.023],cfg["k0"])
out['cached_complete_seconds']=(time.perf_counter()-t)/200
out['selected']='minimum reliable wall-clock including assembly for actual table; no mandated GPU'
(ROOT/'results/benchmark20.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
