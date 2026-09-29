import sys,json,time,hashlib
from pathlib import Path
R=Path(__file__).resolve().parents[1];sys.path.insert(0,str(R/'src'))
import numpy as np
from cylinders import CylinderModel
from physics import Family
from posterior import adaptive_tempered_smc,readout_statistics
c=json.loads((R/'configs/protocol.json').read_text());m=CylinderModel(Family(**c['family']),order=20);q=np.load(R/'results/queries_validation_11_coupled.npz');start=time.perf_counter();out=[]
# Fixed first queries, independent runs. Diagnostic; convergence is not presumed.
for qi in [0,1]:
 for seed in [15171,15172]:
  obs=q['y'][qi,1,:2]+c['sigma_primary']*q['noise'][qi,1,:2]
  def ll(th):
   return np.array([-.5*np.sum(abs(m.solve(t,c['k0'])['y'][:2]-obs)**2)/c['sigma_primary']**2 for t in th])
  rr=adaptive_tempered_smc(ll,c['bounds'],n_particles=256,seed=seed,pcn_steps=3,pcn_beta=.12)
  qq=np.array([m.solve(t,c['k0'])['Q'] for t in rr['theta']]);ss=readout_statistics(qq,rr['weights'])
  np.savez_compressed(R/f'results/smc_q{qi}_seed{seed}.npz',theta=rr['theta'],weights=rr['weights'],Q=qq,log_likelihood=rr['log_likelihood'])
  row={'query':qi,'seed':seed,'n':256,'lambda':rr['lambda'],'ess':rr['ess'],'tempering_complete':rr['converged'],'diagnostics':rr['diagnostics'],'statistics':ss,'evaluations':rr['n_likelihood_evals']+256}
  # JSON conversion preserves all diagnostics and weighted values.
  row=json.loads(json.dumps(row,default=lambda x:x.tolist() if hasattr(x,'tolist') else x));out.append(row)
  (R/'results/smc_validation.json').write_text(json.dumps({'runs':out,'seconds':time.perf_counter()-start,'solves':m.solve_count,'status':'complete' if len(out)==4 else 'running','protocol_hash':hashlib.sha256((R/'configs/protocol.json').read_bytes()).hexdigest()},indent=2));print(qi,seed,ss,flush=True)
