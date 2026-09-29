"""Strict per-query validation; diagnostic node counts never certify convergence."""
import json,hashlib
from pathlib import Path
import numpy as np
R=Path(__file__).resolve().parents[1];c=json.loads((R/'configs/protocol.json').read_text());rows=[]
for p in sorted((R/'results/validation_41_coupled').glob('seed*.json')):
 a=json.loads((R/'results/validation_33_coupled'/p.name).read_text());b=json.loads(p.read_text())
 for x,y in zip(a,b):
  assert (x['seed'],x['query'],x['noise'],x['method'])==(y['seed'],y['query'],y['noise'],y['method'])
  dm=abs(x['mean']-y['mean']);dv=abs(x['variance']-y['variance']);dq=max(abs(np.array(x['quantiles'])-y['quantiles']))
  rows.append({k:y[k] for k in ['seed','query','noise','method']}|{'mean_difference':dm,'variance_difference':dv,'variance_relative':dv/max(y['variance'],1e-30),'quantile_difference':dq,'passed':bool(dm<=c['reference_mean_tol'] and dv<=max(c['reference_variance_abs_tol'],c['reference_variance_rel_tol']*y['variance']) and dq<=c['reference_quantile_tol'])})
summary=[]
for noise in c['noise_levels']:
 for method in sorted({x['method'] for x in rows}):
  rr=[x for x in rows if x['noise']==noise and x['method']==method]
  summary.append({'noise':noise,'method':method,'passed':sum(x['passed'] for x in rr),'total':len(rr),'max_mean_difference':max(x['mean_difference'] for x in rr),'max_variance_relative':max(x['variance_relative'] for x in rr),'max_quantile_difference':max(x['quantile_difference'] for x in rr)})
out={'passed':False,'quadrature_passed':all(x['passed'] for x in rows),'protocol_hash':hashlib.sha256((R/'configs/protocol.json').read_bytes()).hexdigest(),'reason':'Independent SMC and all quadrature tolerances required. Test stays closed.','summary':summary,'branch_status':'not classified; no established physical multimodality','rows':rows}
(R/'results/reference_validation_gate.json').write_text(json.dumps(out,indent=2));print(json.dumps(summary,indent=2))
