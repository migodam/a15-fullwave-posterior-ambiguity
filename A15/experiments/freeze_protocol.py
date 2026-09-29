import sys,json,hashlib,time,itertools
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np
from cylinders import CylinderModel
from physics import Family,pack

p=ROOT/'configs/protocol.json'
if p.exists():raise RuntimeError('Protocol already frozen; refuse overwrite')
qa=json.loads((ROOT/'results/cylinder_validation.json').read_text())['summary']
assert qa['order20_vs24_whitened_rms_max']<.1 and qa['order20_vs24_Q_max']<.02 and qa['flux_balance_max']<1e-8,qa
screen=json.loads((ROOT/'results/discovery/coarse_reference_summary.json').read_text());selected=next(c for c in screen['candidates'] if c['candidate']==screen['selected'])
cfg=json.loads((ROOT/'configs/discovery.json').read_text());f=Family(radius=selected['radius'],epsilon_base=selected['epsilon_base']);m=CylinderModel(f,order=20);k0=selected['k0']
bounds=np.array([[.02,.16],[-.25*f.epsilon_base,.25*f.epsilon_base],[.001,.08]])
rng=np.random.default_rng(15005);th=rng.uniform(bounds[:,0],bounds[:,1],(1024,3));yy=np.array([m.solve(t,k0)['y'][:2] for t in th]);rms=np.sqrt(np.mean(np.r_[yy.real.ravel(),yy.imag.ravel()]**2));sigma=.02*rms
np.savez_compressed(ROOT/'results/discovery/calibration.npz',theta=th,y=yy,rms=rms)
scale=(bounds[:,1]-bounds[:,0])/np.sqrt(12);ks=np.array([.95,1,1.05])*k0
# Data-free local posterior metric is explicitly a discovery DESIGN PROXY only.
angles=list(itertools.combinations(range(2,16),2));risks={'base':[],'neighbor_minus':[],'neighbor_plus':[],'near':[]};pair_risk=np.zeros(len(angles))
for t in th[:32]:
 jac=m.parameter_jacobian(t,k0);b=jac['Q']*scale
 def realA(j):return np.vstack([j.real,j.imag])*scale/sigma
 A=realA(jac['y'][:2]);H=np.eye(3)+A.T@A
 risks['base'].append(float(b@np.linalg.solve(H,b)))
 for name,kk in [('neighbor_minus',ks[0]),('neighbor_plus',ks[2])]:
  C=realA(m.parameter_jacobian(t,kk)['y'][:2]);risks[name].append(float(b@np.linalg.solve(H+C.T@C,b)))
 C=realA(jac['y'][16:]);risks['near'].append(float(b@np.linalg.solve(H+C.T@C,b)))
 for i,pair in enumerate(angles):
  C=realA(jac['y'][list(pair)]);pair_risk[i]+=float(b@np.linalg.solve(H+C.T@C,b))/32
neighbor=0 if np.mean(risks['neighbor_minus'])<np.mean(risks['neighbor_plus']) else 2;optimal=list(angles[int(np.argmin(pair_risk))]);random_pair=np.random.default_rng(15006).choice(np.arange(2,16),2,replace=False).tolist()
plan={**cfg,'stage':'confirmatory_frozen','family':f.__dict__,'bounds':bounds.tolist(),'k0':k0,'frequencies':ks.tolist(),'backend':'converged full-wave cylinder expansion, CPU complex128','modal_order':20,'radial_order':48,'noise_rms':float(rms),'sigma_primary':float(sigma),'primary_channels':[0,1],'extra_fixed':[2,3],'extra_optimized':optimal,'extra_random':random_pair,'near_channels':[16,17],'neighbor_index':neighbor,'external_receivers':m.receivers[:16].tolist(),'near_receivers':m.receivers[16:].tolist(),'prior_box_mass_within_family':1.,'family_weight_if_uniform_six_template_prior':1/6,'selection':screen['selection_rule'],'predictions':{name:float(np.mean(v)) for name,v in risks.items()},'optimized_angle_predicted_risk':float(pair_risk.min()),'target_prediction_scope':'full-wave local parameter tangent proxy, NOT validated two-mode prediction; coupled-mode explanation evaluated separately','coupling_prediction':'deleting cross-body feedback is predicted to decrease Q posterior risk and diminish neighboring-frequency advantage; test both, retain contrary results','created_unix':time.time(),'test_opened':False}
p.write_text(json.dumps(plan,indent=2));digest=hashlib.sha256(p.read_bytes()).hexdigest();(ROOT/'configs/protocol.sha256').write_text(digest+'\n')
(ROOT/'EXPERIMENT_PROTOCOL.md').write_text(f'''# A15 frozen confirmatory protocol\nSHA256: {digest}\n\nSource of requirements: A15_prompt.md; all thresholds inherited from the approved plan and configs/discovery.json.\n\nSelected two-circle family: radius={f.radius}, epsilon_base={f.epsilon_base}, target center=({f.target_x},0), bright center=(target_x-2r-s,0). Independent uniform bounds {bounds.tolist()}. Full width retained: mass=1 within selected family; 1/6 under an explicitly uniform six-template discovery mixture. No universal Maxwell prior probability is asserted.\n\nTime exp(-i omega t), epsilon=epsilon_real+i gamma k_ref/k; unit E_z plane wave propagating +x. ROI is fixed target disk. Q=natural log of mean ROI |E|^2 relative to incident.\n\nk0={k0}; frequencies={ks.tolist()}; RMS per REAL component={rms}; primary sigma={sigma}. Shared absolute noise for all channels, 1/2/5% boundaries.\n\nBase channels 0/1; same-frequency fixed 2/3; optimized {optimal}; random {random_pair}; neighbor frequency index {neighbor}; near 16/17. Every addition has exactly 2 complex outputs. Receiver positions and design-proxy predictions in JSON. Local full-wave proxy prediction is not a successful two-mode explanation.\n\nseeds 11/29/47; validation 128 per seed; locked test 512 per seed. Truth is continuous prior, not quadrature nodes. 33^3/41^3 positive-weight quadrature, independent SMC diagnostic. No table/node failure deletion.\n\nGate A: SD(Q|y)>ln1.25 AND 90% energy interval ratio>=1.5. Each seed prevalence >=.1 and pooled Wilson lower95>=.1, conditional on valid physical/reference accuracy. Gate B also requires two-mode errors below noise/task margin, predicted targeted advantage >=10 percentage points vs fixed AND optimized same-frequency baselines, paired query bootstrap with Holm adjustment and consistent seed effects, plus coupling attribution.\n\nFrozen numerical backend adjustment: VIE retained and audited; exact-circle expansion order20 with independent order24 checks selected because VIE n<=128 did not meet the noise error budget. This does not alter physical family or gates.\n\nFlow blocked until Gate B PASS and cheap simple posterior models fail. All readouts use same physical theta. Failure/incomplete results must remain visible. Resume always rereads full A15_prompt.md before this protocol and RESUME_STATE.md.\n''')
(ROOT/'RESUME_STATE.md').write_text(f'''# A15 resume\nRead A15_prompt.md completely first.\nStage: frozen protocol, reference table implementation.\nProtocol SHA256: {digest}.\nDiscovery: complete, coarse screen only; six templates, three followups. Exact full-wave circular backend order20 passed output truncation/flux checks; VIE limitations retained.\nGate A/B NOT RUN. Test NOT OPENED. Flow BLOCKED.\nNext: generate 33^3 and41^3 tables on direct-LAN XINAN; compare validation posterior integration, no test selection.\n''')
print(json.dumps({'family':f.__dict__,'k0':k0,'sigma':sigma,'neighbor_index':neighbor,'optimized_angles':optimal,'hash':digest},indent=2))
