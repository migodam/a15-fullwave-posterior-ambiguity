"""Positive-weight adaptive cubature, full prior support, no failed-node pruning.
A cell's 3^3 and 5^3 rules estimate error; all current 5^3 leaf rules form a
positive reference measure. Budget-limited results NEVER silently pass.
"""
import sys,json,time,heapq,hashlib
from pathlib import Path
R=Path(__file__).resolve().parents[1];sys.path.insert(0,str(R/'src'))
import numpy as np
from scipy.special import logsumexp
from scipy.optimize import minimize
from numpy.polynomial.legendre import leggauss
from cylinders import CylinderModel
from physics import Family
c=json.loads((R/'configs/protocol.json').read_text());bounds=np.array(c['bounds']);span=np.diff(bounds).ravel();m=CylinderModel(Family(**c['family']),order=20);qd=np.load(R/'results/queries_validation_11_coupled.npz');table=R/'results/table_41_coupled';tn=np.load(table/'prior.npz')['theta'];ty=np.load(table/'y.npy')[:,1,:2]
allruns=[];t0=time.perf_counter()
for qi in [0,1]:
 obs=qd['y'][qi,1,:2]+c['sigma_primary']*qd['noise'][qi,1,:2];cache={};calls0=m.solve_count
 def ev(u):
  key=tuple(np.round(u,15))
  if key not in cache:
   s=m.solve(bounds[:,0]+span*u,c['k0']);cache[key]=(-.5*np.sum(abs(s['y'][:2]-obs)**2)/c['sigma_primary']**2,s['Q'])
  return cache[key]
 # Locate likely mass only to seed subdivision; no removed prior region.
 best=tn[np.argmin(np.sum(abs(ty-obs)**2,axis=1))];u0=(best-bounds[:,0])/span
 opt=minimize(lambda u:-ev(u)[0],u0,method='Nelder-Mead',bounds=[(0,1)]*3,options={'maxiter':100,'xatol':1e-7});peak=-opt.fun
 def cell(lo,hi):
  values=[]
  for n in [3,5]:
   x,w=leggauss(n);g=np.meshgrid(*([x]*3),indexing='ij');u=(np.stack(g,-1).reshape(-1,3)+1)/2*(hi-lo)+lo
   ww=np.prod(np.stack(np.meshgrid(w,w,w,indexing='ij'),-1).reshape(-1,3),axis=1)*np.prod(hi-lo)/8
   val=np.array([ev(v) for v in u]);weights=ww*np.exp(val[:,0]-peak)
   moment=np.array([weights.sum(),weights@val[:,1],weights@(val[:,1]**2)])
   values.append((moment,u,val,weights))
  err=np.abs(values[0][0]-values[1][0]);return dict(lo=lo,hi=hi,err=err,mom=values[1][0],u=values[1][1],val=values[1][2],w=values[1][3])
 # Put an interior boundary on each side of MAP in two most likelihood-sensitive axes.
 cuts=[]
 for i in range(3):
  if i==2:cuts.append([0.,1.])
  else:cuts.append(sorted(set([0.,float(np.clip(opt.x[i]-.06,.01,.99)),float(np.clip(opt.x[i]+.06,.01,.99)),1.])))
 cells=[]
 for i in range(len(cuts[0])-1):
  for j in range(len(cuts[1])-1):cells.append(cell(np.array([cuts[0][i],cuts[1][j],0.]),np.array([cuts[0][i+1],cuts[1][j+1],1.])))
 history=[]
 while True:
  mom=sum(x['mom'] for x in cells);err=sum(x['err'] for x in cells);mean=mom[1]/mom[0];var=max(0,mom[2]/mom[0]-mean**2)
  ww=np.concatenate([x['w'] for x in cells]);vv=np.concatenate([x['val'][:,1] for x in cells]);idx=np.argsort(vv);quant=np.interp([.025,.05,.5,.95,.975],np.cumsum(ww[idx])/ww.sum(),vv[idx])
  history.append({'solves':m.solve_count-calls0,'leaves':len(cells),'mean':float(mean),'variance':float(var),'quantiles':quant.tolist(),'relative_moment_error_estimates':(err/np.maximum(abs(mom),1e-30)).tolist()})
  if m.solve_count-calls0>=9500:break
  # Refine largest relative raw-moment error, split longest sensitivity-scaled side.
  scores=[np.max(x['err']/np.maximum(abs(mom),1e-30)) for x in cells];ci=int(np.argmax(scores));cc=cells.pop(ci)
  length=(cc['hi']-cc['lo'])*np.array([1,1,.1]);axis=int(np.argmax(length));mid=(cc['lo'][axis]+cc['hi'][axis])/2
  hi=cc['hi'].copy();hi[axis]=mid;lo=cc['lo'].copy();lo[axis]=mid;cells.extend([cell(cc['lo'],hi),cell(lo,cc['hi'])])
 np.savez_compressed(R/f'results/adaptive_q{qi}.npz',theta=bounds[:,0]+span*np.concatenate([x['u'] for x in cells]),Q=vv,weights=ww/ww.sum())
 allruns.append({'query':qi,'history':history,'solves':m.solve_count-calls0,'status':'budget-limited; local embedded-rule agreement is not an independent convergence proof','map_success':bool(opt.success)})
 (R/'results/adaptive_reference.json').write_text(json.dumps({'runs':allruns,'seconds':time.perf_counter()-t0,'solves':m.solve_count,'scope':'all prior support retained, positive tensor weights on adaptive disjoint cells; two preselected validation queries only'},indent=2));print(qi,history[-1],flush=True)
