import sys,os,json,time,platform
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np
from scipy.linalg import lu_factor,lu_solve
from physics import PhysicalModel,Family,circle_fractions,cylinder_exact,pack
from physics_backend import Solver,self_cell_green

def main():
    out={'stage':'implementation_validation','environment':{'python':sys.version,'platform':platform.platform()},'checks':{},'mesh':[]}
    started=time.perf_counter();rng=np.random.default_rng(15002)
    model=PhysicalModel(32);f=model.family;k=6.;theta=np.array([.083,.7,.023])
    # Uniform geometry integration and actual nonzero material displacement.
    area=[]
    for s in np.linspace(.02,.16,23):
        b,t=model.fractions(s);area.extend([abs(model.h**2*b.sum()-np.pi*f.radius**2),abs(model.h**2*t.sum()-np.pi*f.radius**2)])
    out['checks']['circle_area_max_abs']=max(area)
    sol=model.solve(theta,k,fields=True)
    old=Solver(n=32,k=k,n_rx=2,n_illum=1);chi=np.zeros(32**2,complex);chi[model.indices]=sol['chi'];ref=old.solve(chi)
    out['checks']['active_vs_full_E']=float(np.linalg.norm(sol['E']-ref['E'][model.indices,0])/np.linalg.norm(sol['E']))
    out['checks']['active_vs_full_J']=float(np.linalg.norm(sol['J']-ref['J'][model.indices,0])/np.linalg.norm(sol['J']))
    out['checks']['absorption_positive']=sol['absorption']>=0
    out['checks']['flux']=model.flux_diagnostics(theta,k)
    # Zero contrast and Born limits use copied kernel on original full grid.
    zero=old.solve(np.zeros(32**2));out['checks']['zero_contrast']=float(np.linalg.norm(zero['J']))
    born=old.GS@(1e-7*chi[:,None]*old.incident);weak=old.solve(1e-7*chi)['y']
    out['checks']['born_relative']=float(np.linalg.norm(weak-born)/np.linalg.norm(born))
    out['tangents']=[]
    for scene in range(8):
        th=np.array([rng.uniform(.035,.145),rng.uniform(-1.5,1.5),rng.uniform(.005,.07)])
        jac=model.parameter_jacobian(th,k)
        for direction_id in range(8):
            v=rng.normal(size=3)*np.array([.02,1.,.01]);direction=v/np.linalg.norm(v)
            d=jac['y']@direction;c=rng.normal(size=18)+1j*rng.normal(size=18)
            adj=abs(np.real(np.vdot(c,d))-direction@np.real(jac['y'].conj().T@c))/max(abs(np.real(np.vdot(c,d))),1.)
            errors=[]
            for eps in 10.**np.arange(-3,-8,-1):
                yp=model.solve(th+eps*direction,k);ym=model.solve(th-eps*direction,k)
                fd=(yp['y']-ym['y'])/(2*eps);fdq=(yp['Q']-ym['Q'])/(2*eps)
                errors.append({'eps':eps,'y':float(np.linalg.norm(fd-d)/max(np.linalg.norm(d),1e-12)),
                               'Q':float(abs(fdq-jac['Q']@direction)/max(abs(jac['Q']@direction),1e-12))})
            out['tangents'].append({'scene':scene,'direction':direction_id,'adjoint_error':adj,'errors':errors})
    # One-cylinder independent analytic check using same VIE kernel, not old solver outputs.
    for n in [32,48,64,96,128]:
        m=PhysicalModel(n);gd,gs,inc=m.operators(k);a=m.roi>0;cc=m.roi[a]*(f.epsilon_base-1+1j*.023*f.k_ref/k)
        mat=np.eye(a.sum())-cc[:,None]*gd[np.ix_(a,a)];jj=lu_solve(lu_factor(mat),cc*inc[a])
        yy=gs[:,a]@jj;exact=cylinder_exact(m.receivers,k,f.radius,f.epsilon_base+1j*.023*f.k_ref/k,(f.target_x,0))
        pair=m.solve(theta,k)
        out['mesh'].append({'n':n,'single_cylinder_relative_y':float(np.linalg.norm(yy-exact)/np.linalg.norm(exact)),
                            'single_cylinder_max_abs_y':float(abs(yy-exact).max()),'pair_y_real':pair['y'].real.tolist(),
                            'pair_y_imag':pair['y'].imag.tolist(),'pair_Q':pair['Q'],'seconds':pair['seconds'],'active_cells':pair['active_cells']})
        print('mesh',n,out['mesh'][-1]['single_cylinder_relative_y'],pair['Q'],flush=True)
    out['seconds']=time.perf_counter()-started
    path=ROOT/'results/physics_validation.json';path.write_text(json.dumps(out,indent=2));print(path)
if __name__=='__main__':main()
