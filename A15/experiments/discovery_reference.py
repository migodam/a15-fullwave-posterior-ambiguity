import sys,json,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np
from scipy.special import logsumexp
from cylinders import CylinderModel
from physics import Family
from posterior import gauss_legendre_box_nodes

def stats(q,y,prior,observation,sigma):
    lw=np.log(prior)-.5*np.sum(abs(y-observation)**2,axis=1)/sigma**2;w=np.exp(lw-logsumexp(lw));mean=w@q
    return mean,w@(q-mean)**2

def main():
    discovery=json.loads((ROOT/'results/discovery/summary.json').read_text());ranked=sorted(discovery['templates'],key=lambda x:x['best']['local_task_variance_proxy'],reverse=True)[:3]
    records=[]
    for i,c in enumerate(ranked):
        f=Family(radius=c['radius'],epsilon_base=c['epsilon_base']);k=c['best']['k'];ks=k*np.array([.95,1,1.05]);model=CylinderModel(f)
        bounds=np.array([[.02,.16],[-.25*f.epsilon_base,.25*f.epsilon_base],[.001,.08]])
        rule=gauss_legendre_box_nodes(bounds,9);theta=rule['nodes'];path=ROOT/'results/discovery'/('coarse_'+c['name']+'.npz')
        if path.exists():data=np.load(path);ys=data['y'];qs=data['Q']
        else:
            ys=np.empty((len(theta),3,18),complex);qs=np.empty((len(theta),3))
            for a,th in enumerate(theta):
                for b,kk in enumerate(ks):s=model.solve(th,kk);ys[a,b]=s['y'];qs[a,b]=s['Q']
            np.savez_compressed(path,theta=theta,weights=rule['weights'],k=ks,y=ys,Q=qs)
        rng=np.random.default_rng(np.random.SeedSequence([15003,i]));truth=rng.uniform(bounds[:,0],bounds[:,1],size=(128,3));sigma=.02*np.sqrt(np.sum(rule['weights'][:,None]*abs(ys[:,1,:2])**2)/4)
        var=[];actual=[];cases=[]
        for th in truth:
            s=model.solve(th,k);obs=s['y'][:2]+sigma*(rng.normal(size=2)+1j*rng.normal(size=2));mu,v=stats(qs[:,1],ys[:,1,:2],rule['weights'],obs,sigma);var.append(v);actual.append((s['Q']-mu)**2)
            cases.append({'theta':th.tolist(),'Q':s['Q'],'mean_Q':mu,'variance_Q':v})
        record={'candidate':c['name'],'radius':f.radius,'epsilon_base':f.epsilon_base,'k0':k,'ks':ks.tolist(),'sigma':sigma,
                'coarse_risk':float(np.mean(var)),'coarse_frequency_sd_above':float(np.mean(np.sqrt(var)>np.log(1.25))),
                'mean_estimator_squared_error':float(np.mean(actual)),'prior_box_mass':1.,'note':'9^3 discovery screening only; NOT reference truth, coverage or population evidence',
                'cases':cases}
        records.append(record);(ROOT/'results/discovery/coarse_reference_summary.json').write_text(json.dumps({'status':'running','candidates':records},indent=2));print(record['candidate'],record['coarse_risk'],record['coarse_frequency_sd_above'],flush=True)
    # Keep broad prior; maximum discovery posterior-risk screen with local proxy tie-break.
    selected=max(records,key=lambda r:r['coarse_risk']);(ROOT/'results/discovery/coarse_reference_summary.json').write_text(json.dumps({'status':'complete','candidates':records,'selected':selected['candidate'],'selection_rule':'largest coarse task risk across three pre-ranked candidates; full discovery theta box retained; no confirmatory prevalence claim'},indent=2))
if __name__=='__main__':main()
