import json,hashlib
from pathlib import Path
import numpy as np
R=Path(__file__).resolve().parents[1];c=json.loads((R/'configs/protocol.json').read_text());rng=np.random.default_rng(15191);output={'status':'PROVISIONAL VALIDATION DIAGNOSTIC — reference failed, not confirmatory','protocol_hash':hashlib.sha256((R/'configs/protocol.json').read_bytes()).hexdigest(),'methods':{},'comparisons':{}}
folder=R/'results/validation_41_coupled'
def read(method):
 return sum([json.loads((folder/f'seed{seed}_noise0.02_{method}.json').read_text()) for seed in c['seeds']],[])
b=read('base');n=len(b);ix=rng.integers(n,size=(5000,n));bv=np.array([x['variance'] for x in b]);baseboot=bv[ix].mean(1)
for method in ['base1','base','same_angles','optimized_angles','neighbor','near','random','base8']:
 rr=read(method);v=np.array([x['variance'] for x in rr]);hit=np.array([x['ambiguous'] for x in rr]);mse=np.array([(x['mean']-x['Q_true'])**2 for x in rr]);vb=v[ix].mean(1);red=1-vb/baseboot
 p=hit.mean();z=1.95996398454;den=1+z*z/n;mid=(p+z*z/2/n)/den;half=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
 output['methods'][method]={'risk':float(v.mean()),'risk_ci':np.quantile(vb,[.025,.975]).tolist(),'risk_reduction':float(1-v.mean()/bv.mean()),'risk_reduction_ci':np.quantile(red,[.025,.975]).tolist(),'prevalence':p,'prevalence_wilson':[(mid-half),(mid+half)],'seed_prevalence':[np.mean(hit[k*128:(k+1)*128]) for k in range(3)],'coverage90':np.mean([x['covered90'] for x in rr]),'mse':float(mse.mean()),'mse_minus_variance_ci':np.quantile((mse-v)[ix].mean(1),[.025,.975]).tolist()}
for target in ['neighbor','near']:
 tv=np.array([x['variance'] for x in read(target)])
 for baseline in ['same_angles','optimized_angles']:
  vv=np.array([x['variance'] for x in read(baseline)]);diff=(vv-tv);boot=diff[ix].mean(1)/baseboot
  output['comparisons'][target+'_vs_'+baseline]={'additional_risk_reduction_pp':float(diff.mean()/bv.mean()*100),'paired95_pp':(np.quantile(boot,[.025,.975])*100).tolist(),'familywise97_5_pp':(np.quantile(boot,[.0125,.9875])*100).tolist(),'seed_pp':[float(diff[k*128:(k+1)*128].mean()/bv[k*128:(k+1)*128].mean()*100) for k in range(3)],'inference':'numerical reference invalid; intervals do not include integration bias; no Holm-confirmatory PASS asserted'}
(R/'results/population_diagnostics.json').write_text(json.dumps(output,indent=2));print(json.dumps(output,indent=2))
