import sys,json,time,argparse,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np
from scipy.special import logsumexp
from cylinders import CylinderModel
from physics import Family

def run():
 ap=argparse.ArgumentParser();ap.add_argument('--order',type=int,required=True);ap.add_argument('--split',choices=['validation','test'],default='validation');ap.add_argument('--coupling',type=int,default=1);a=ap.parse_args()
 cfg=json.loads((ROOT/'configs/protocol.json').read_text());ph=hashlib.sha256((ROOT/'configs/protocol.json').read_bytes()).hexdigest()
 if a.split=='test':
  gate=json.loads((ROOT/'results/reference_validation_gate.json').read_text());assert gate['passed'] and gate['protocol_hash']==ph,'test blocked by reference validation'
 folder=ROOT/'results'/f'table_{a.order}_{"coupled" if a.coupling else "uncoupled"}';meta=json.loads((folder/'status.json').read_text());assert meta['status']=='complete' and meta['protocol_hash']==ph
 y=np.load(folder/'y.npy',mmap_mode='r').reshape(-1,54);Q=np.load(folder/'Q.npy',mmap_mode='r')[:,1];prior=np.load(folder/'prior.npz');theta=prior['theta'];w0=prior['weights'];lw0=np.log(w0)
 sortedQ=np.argsort(Q);qsorted=Q[sortedQ]
 selectors={'base1':[18],'base':[18,19],'same_angles':[18,19,20,21],'base8':list(range(18,26)),
            'optimized_angles':[18,19]+[18+i for i in cfg['extra_optimized']],
            'neighbor':[18,19]+[18*cfg['neighbor_index']+i for i in [0,1]],
            'near':[18,19,34,35],'random':[18,19]+[18+i for i in cfg['extra_random']]}
 outdir=ROOT/'results'/f'{a.split}_{a.order}_{"coupled" if a.coupling else "uncoupled"}';outdir.mkdir(exist_ok=True)
 count=cfg[a.split+'_per_seed'];model=CylinderModel(Family(**cfg['family']),order=24);bounds=np.array(cfg['bounds']);allrows=[];start=time.perf_counter()
 for seed in cfg['seeds']:
  querypath=ROOT/'results'/f'queries_{a.split}_{seed}_{"coupled" if a.coupling else "uncoupled"}.npz'
  if querypath.exists():qdata=np.load(querypath);truth=qdata['theta'];clean=qdata['y'];qtrue=qdata['Q'];noise=qdata['noise']
  else:
   rng=np.random.default_rng(np.random.SeedSequence([seed,15100 if a.split=='validation' else 15200]))
   truth=rng.uniform(bounds[:,0],bounds[:,1],(count,3));noise=rng.normal(size=(count,3,18))+1j*rng.normal(size=(count,3,18));clean=np.empty((count,3,18),complex);qtrue=np.empty(count)
   for i,th in enumerate(truth):
    for j,k in enumerate(cfg['frequencies']):
     s=model.solve(th,k,coupling=bool(a.coupling));clean[i,j]=s['y']
     if j==1:qtrue[i]=s['Q']
   np.savez_compressed(querypath,theta=truth,y=clean,Q=qtrue,noise=noise)
  for noise_level in cfg['noise_levels']:
   sigma=noise_level*cfg['noise_rms'];observations=(clean+sigma*noise).reshape(count,54)
   for method,sel in selectors.items():
    pred=y[:,sel];pn=np.sum(abs(pred)**2,axis=1);rows=[]
    for i,obs in enumerate(observations[:,sel]):
     d2=pn-2*np.real(pred.conj()@obs)+np.sum(abs(obs)**2)
     # Dot expansion can create tiny roundoff negatives, not physical negative distances.
     d2=np.maximum(d2,0);lw=lw0-.5*d2/sigma**2;logZ=logsumexp(lw);w=np.exp(lw-logZ)
     mean=float(w@Q);var=float(w@((Q-mean)**2));cum=np.cumsum(w[sortedQ]);quant=np.interp([.025,.05,.5,.95,.975],cum,qsorted)
     mt=w@theta;dt=theta-mt;ct=dt.T@(w[:,None]*dt)
     py=w@pred;dc=pred-py;pc=(dc.conj().T@(w[:,None]*dc))
     row={'seed':seed,'query':i,'noise':noise_level,'method':method,'Q_true':float(qtrue[i]),'mean':mean,'variance':var,'std':float(np.sqrt(var)),
          'quantiles':quant.tolist(),'covered90':bool(quant[1]<=qtrue[i]<=quant[3]),'covered95':bool(quant[0]<=qtrue[i]<=quant[4]),
          'ambiguous':bool(np.sqrt(var)>cfg['tau_q'] and quant[3]-quant[1]>=np.log(cfg['ci90_energy_ratio_min'])),
          'theta_mean':mt.tolist(),'theta_covariance':ct.tolist(),'predictive_mean_real':py.real.tolist(),'predictive_mean_imag':py.imag.tolist(),
          'predictive_covariance_real':pc.real.tolist(),'predictive_covariance_imag':pc.imag.tolist(),
          'log_evidence_without_noise_constant':float(logZ),'integration_effective_nodes':float(1/np.sum(w*w))}
     rows.append(row)
    path=outdir/f'seed{seed}_noise{noise_level}_{method}.json';path.write_text(json.dumps(rows));allrows.extend(rows)
   print(a.split,a.order,seed,noise_level,round(time.perf_counter()-start,1),flush=True)
 summary=[]
 for level in cfg['noise_levels']:
  for method in selectors:
   rr=[r for r in allrows if r['noise']==level and r['method']==method];summary.append({'noise':level,'method':method,'count':len(rr),'risk':float(np.mean([r['variance'] for r in rr])),'prevalence':float(np.mean([r['ambiguous'] for r in rr])),'coverage90':float(np.mean([r['covered90'] for r in rr])),'mean_squared_error':float(np.mean([(r['mean']-r['Q_true'])**2 for r in rr]))})
 (outdir/'summary.json').write_text(json.dumps({'summary':summary,'seconds':time.perf_counter()-start,'protocol_hash':ph,'source_table':meta,'split':a.split,'status':'complete'},indent=2))
if __name__=='__main__':run()
