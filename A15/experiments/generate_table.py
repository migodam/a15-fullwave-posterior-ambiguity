import sys,json,time,hashlib,argparse,os
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np
from cylinders import CylinderModel
from physics import Family
from posterior import gauss_legendre_box_nodes

def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def atomic(path,obj):
 p=Path(path);tmp=p.with_suffix('.tmp');tmp.write_text(json.dumps(obj,indent=2));tmp.replace(p)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--order',type=int,required=True);ap.add_argument('--coupling',type=int,default=1);a=ap.parse_args()
 cfg=json.loads((ROOT/'configs/protocol.json').read_text());ph=digest(ROOT/'configs/protocol.json');assert ph==(ROOT/'configs/protocol.sha256').read_text().strip()
 codehash={p.name:digest(p) for p in (ROOT/'src').glob('*.py')};rule=gauss_legendre_box_nodes(cfg['bounds'],a.order);theta=rule['nodes'];N=len(theta)
 name=f'table_{a.order}_'+('coupled' if a.coupling else 'uncoupled');folder=ROOT/'results'/name;folder.mkdir(exist_ok=True)
 statuspath=folder/'status.json';old=json.loads(statuspath.read_text()) if statuspath.exists() else None
 if old:
  assert old['protocol_hash']==ph and old['code_hashes']==codehash,'resume provenance mismatch'
  if old['status']=='complete':print('already complete');return
 mode='r+' if old else 'w+';start_index=old['completed'] if old else 0
 y=np.lib.format.open_memmap(folder/'y.npy',mode=mode,dtype=np.complex128,shape=(N,3,18))
 Q=np.lib.format.open_memmap(folder/'Q.npy',mode=mode,dtype=np.float64,shape=(N,3))
 diagnostics=np.lib.format.open_memmap(folder/'diagnostics.npy',mode=mode,dtype=np.float64,shape=(N,3,3))
 if not old:np.savez(folder/'prior.npz',theta=theta,weights=rule['weights'])
 model=CylinderModel(Family(**cfg['family']),order=cfg['modal_order'],radial_order=cfg['radial_order']);started=time.perf_counter()
 status={'status':'running','completed':start_index,'total':N,'order':a.order,'coupling':bool(a.coupling),'protocol_hash':ph,'code_hashes':codehash,'seed':'deterministic Gauss-Legendre prior integration','host':os.environ.get('COMPUTERNAME','local'),'seconds_before_resume':old.get('seconds',0) if old else 0}
 for lo in range(start_index,N,256):
  hi=min(N,lo+256)
  try:
   for i in range(lo,hi):
    for j,k in enumerate(cfg['frequencies']):
     s=model.solve(theta[i],k,coupling=bool(a.coupling));y[i,j]=s['y'];Q[i,j]=s['Q'];diagnostics[i,j]=[s['residual'],s['absorption'],s['current_energy']]
  except Exception as exc:
   atomic(folder/f'failure_{i}.json',{'index':i,'theta':theta[i].tolist(),'error':repr(exc),'frequency_index':j});raise
  y.flush();Q.flush();diagnostics.flush();status.update(completed=hi,seconds=status['seconds_before_resume']+time.perf_counter()-started,full_wave_solves=hi*3)
  atomic(statuspath,status)
  if hi%4096==0 or hi==N:print(name,hi,N,round(status['seconds'],2),flush=True)
 status['status']='complete';status['files']={p.name:digest(p) for p in folder.glob('*.npy')};status['files']['prior.npz']=digest(folder/'prior.npz');atomic(statuspath,status)
 print('COMPLETE',name,status['seconds'],flush=True)
if __name__=='__main__':main()
