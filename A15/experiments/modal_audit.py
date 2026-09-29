"""Fixed-grid physical injection audit; discretization limitations kept explicit."""
import sys,json,time
from pathlib import Path
R=Path(__file__).resolve().parents[1];sys.path.insert(0,str(R/'src'))
import numpy as np
from scipy.linalg import svd,solve
from physics import Family,PhysicalModel
from cylinders import CylinderModel
c=json.loads((R/'configs/protocol.json').read_text());family=Family(**c['family']);model=CylinderModel(family,order=20);vm=PhysicalModel(64,family)
cal=np.load(R/'results/discovery/calibration.npz');print(cal.files,flush=True)
# Predeclared first 16 calibration states; not selected from test.
thetas=cal['theta'][:16];records=[];t0=time.perf_counter()
for i,theta in enumerate(thetas):
 sol=model.solve(theta,c['k0'],fields=True);M=sol['state']['M'];u,s,vh=svd(M);v=vh.conj().T
 # Resolvent truncation and response-POD projection are distinct diagnostics.
 Z=np.stack([model.solve(theta,k,fields=True)['z'] for k in np.linspace(c['k0']*.95,c['k0']*1.05,9)],axis=1)
 # All expansion z use unit-current L2 modes, but radial bases depend on frequency;
 # sample common spatial current fields for genuinely common-coordinate response POD below.
 fields=[];points=vm.points
 for k in np.linspace(c['k0']*.95,c['k0']*1.05,9):
  so=model.solve(theta,k,fields=True);ef=model.field(so,points);mat=vm.material(theta,k)['chi'];fields.append(mat*ef)
 J=np.stack(fields,axis=1)*vm.h;ur,sr,_=svd(J,full_matrices=False)
 f=vm.solve(theta,c['k0'],fields=True);mm=np.eye(len(f['chi']))-f['chi'][:,None]*f['GD'];uv,ss,vvh=svd(mm);vv=vvh.conj().T
 injection=f['E'][:,None]*vm.material_jacobian(theta,c['k0'])*np.array([.14,6,.079])[None,:]/np.sqrt(12)
 amp=(uv.conj().T@injection)/ss[:,None];vis=np.linalg.norm(f['GS'][:2]@vv,axis=0)
 L=(2*np.conj(f['E'])*vm.roi)@f['GD']/np.sum(vm.roi*abs(f['E'])**2)
 task=np.abs(L@vv);inj=np.linalg.norm(uv.conj().T@injection,axis=1)
 residual_terms=[];rom={}
 for name,V in [('resolvent',v[:,-2:])]:
  z2=V@solve(V.conj().T@M@V,V.conj().T@sol['state']['rhs'])
  yerr=float(np.linalg.norm(sol['S'][:2]@(z2-sol['z']))/2/c['sigma_primary']);en=np.sum(abs(z2[model.d:]/(sol['state']['epsilon'][1]-1))**2);q2=np.log(en/(np.pi*family.radius**2))
  rom[name]={'response_error_sigma_rms_real':yerr,'Q_abs_error':float(abs(q2-sol['Q']))}
 V=ur[:,:2]/vm.h
 jr=V@solve(V.conj().T@mm@V,V.conj().T@(f['chi']*f['incident']))
 er=f['incident']+f['GD']@jr;qr=np.log(np.sum(vm.roi*abs(er)**2)/vm.roi.sum())
 rom['response_POD']={'response_error_sigma_rms_real':float(np.linalg.norm(f['GS'][:2]@(jr-f['J']))/2/c['sigma_primary']),'Q_abs_error':float(abs(qr-f['Q'])),'response_snapshot_energy_fraction':float(np.sum(sr[:2]**2)/np.sum(sr**2))}
 rec={'scene':i,'theta':theta.tolist(),'Q':sol['Q'],'VIE_vs_circle_Q':float(f['Q']-sol['Q']),'two_mode':rom,'smallest_singular_values':ss[-8:].tolist(),'injection_last8':inj[-8:].tolist(),'measurement_visibility_last8':vis[-8:].tolist(),'task_visibility_last8':task[-8:].tolist(),'tangent_current_residual_fraction_last2':float(np.linalg.norm(amp[:-2])/np.linalg.norm(amp))}
 records.append(rec)
 if i==0:
  images=np.zeros((64*64,2),complex);images[vm.indices]=vv[:,-2:];
  np.savez_compressed(R/'results/modal_audit.npz',singular_values=ss,injection=inj,measurement_visibility=vis,task_visibility=task,current_modes=images.reshape(64,64,2),mode_Q_effect=task/ss,mode_y_effect=vis/ss,response_singular_values=sr)
 (R/'results/modal_audit.json').write_text(json.dumps({'records':records,'scope':'fixed-grid VIE physical D_E J_material injection; n64 has known forward bias so diagnostic, not quantitative certification. Two-mode response POD is optimistic in-family frequency response basis, not universal deployable model. Exact circle resolvent in mass-normalized current coordinates separately.','status':'complete' if i==15 else 'running','seconds':time.perf_counter()-t0},indent=2));print(i,rom,flush=True)
# Selected-family 81-point frequency refinement, diagnostic only; frozen k0 unchanged.
k=np.linspace(c['k0']-.45,c['k0']+.45,81);out=[model.solve(thetas[0],x) for x in k]
np.savez_compressed(R/'results/frequency_refined.npz',k=k,Q=[x['Q'] for x in out],y=[x['y'] for x in out],theta=thetas[0])
