import sys,json,time,itertools
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np
from scipy.stats import qmc
from scipy.linalg import svd,solve
from cylinders import CylinderModel
from physics import Family,pack


def jsonwrite(path,data):
    p=Path(path);tmp=p.with_suffix('.tmp');tmp.write_text(json.dumps(data,indent=2));tmp.replace(p)

def main():
    cfg=json.loads((ROOT/'configs/discovery.json').read_text());folder=ROOT/'results/discovery';folder.mkdir(exist_ok=True)
    klist=np.linspace(3,12,41);unit=qmc.Sobol(3,scramble=True,seed=15001).random_base2(5)[:24]
    records=[];start=time.perf_counter()
    for radius,epsilon in itertools.product(cfg['radii'],cfg['epsilon_base']):
        name=f'r{radius:.2f}_e{epsilon:.0f}';path=folder/(name+'.npz')
        theta=np.c_[.02+.14*unit[:,0],epsilon*(-.25+.5*unit[:,1]),.001+.079*unit[:,2]]
        if path.exists():
            data=np.load(path);ys=data['y'];qs=data['Q'];res=data['residual']
        else:
            model=CylinderModel(Family(radius=radius,epsilon_base=epsilon),order=12)
            ys=np.empty((24,41,18),complex);qs=np.empty((24,41));res=np.empty((24,41))
            for i,th in enumerate(theta):
                for j,k in enumerate(klist):
                    sol=model.solve(th,k);ys[i,j]=sol['y'];qs[i,j]=sol['Q'];res[i,j]=sol['residual']
            np.savez_compressed(path,theta=theta,k=klist,y=ys,Q=qs,residual=res)
        # Local scale-normalized diagnosis ONLY, never posterior truth.
        model=CylinderModel(Family(radius=radius,epsilon_base=epsilon),order=12)
        scale=np.array([.14,.5*epsilon,.079])/np.sqrt(12)
        scores=[]
        for j,k in enumerate(klist):
            sigma=.02*np.sqrt(np.mean(np.r_[ys[:,j,:2].real.ravel(),ys[:,j,:2].imag.ravel()]**2))
            local=[];ratios=[];romerr=[]
            for i in [0,8,16]:
                jac=model.parameter_jacobian(theta[i],k);A=np.vstack([jac['y'][:2].real,jac['y'][:2].imag])*scale/sigma;b=jac['Q']*scale
                cov=np.linalg.inv(np.eye(3)+A.T@A);local.append(float(b@cov@b))
                sol=model.solve(theta[i],k,fields=True);u,s,vh=svd(sol['state']['M']);V=vh.conj().T[:,-2:]
                zr=V@solve(V.conj().T@sol['state']['M']@V,V.conj().T@sol['state']['rhs'])
                romerr.append(float(np.linalg.norm(pack(sol['S']@(zr-sol['z'])))/sigma/2))
                vis=np.linalg.norm(sol['S'][:2]@V,axis=0);ratios.append(float(max(vis)/max(min(vis),1e-30)))
            scores.append({'k':float(k),'local_task_variance_proxy':float(np.median(local)),
                           'two_mode_response_error_sigma_median':float(np.median(romerr)),
                           'radiation_contrast_median':float(np.median(ratios)),
                           'Q_std_prior_probe':float(qs[:,j].std()),'sigma':float(sigma)})
        # Ranking does not exclude failure of the explanatory two-mode model.
        best=max(scores,key=lambda x:x['local_task_variance_proxy'])
        records.append({'name':name,'radius':radius,'epsilon_base':epsilon,'frequencies':scores,'best':best,
                        'max_residual':float(res.max()),'states':len(theta),'solves':24*41+41*3*2})
        jsonwrite(folder/'summary.json',{'status':'running','templates':records,'seconds':time.perf_counter()-start})
        print(name,best,flush=True)
    ranked=sorted(records,key=lambda x:x['best']['local_task_variance_proxy'],reverse=True)
    jsonwrite(folder/'summary.json',{'status':'complete','templates':records,'candidate_order':[r['name'] for r in ranked],
                                    'seconds':time.perf_counter()-start,'selection':'median local Q variance diagnostic, NOT posterior prevalence; refine top three candidates with discovery quadrature before freeze',
                                    'test_opened':False})
if __name__=='__main__':main()
