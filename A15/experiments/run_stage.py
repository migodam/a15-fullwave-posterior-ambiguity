"""Stage guard. Formal execution cannot inherit a positive gate from a stale run."""
import argparse,json,hashlib,subprocess,sys
from pathlib import Path
R=Path(__file__).resolve().parents[1]
def guard(stage):
 p=R/'configs/protocol.json';ph=hashlib.sha256(p.read_bytes()).hexdigest()
 assert ph==(R/'configs/protocol.sha256').read_text().strip(),'frozen protocol changed'
 if stage=='test':
  g=json.loads((R/'results/reference_validation_gate.json').read_text())
  assert g['passed'] and g['protocol_hash']==ph,'reference gate closed'
  for order in [33,41]:
   t=R/f'results/table_{order}_coupled';m=json.loads((t/'status.json').read_text());assert m['status']=='complete' and m['protocol_hash']==ph
   for fn,h in m['files'].items():assert hashlib.sha256((t/fn).read_bytes()).hexdigest()==h,'table hash mismatch'
 elif stage=='flow':
  g=json.loads((R/'results/gates.json').read_text())
  assert g['protocol_hash']==ph and g['Gate_B']=='PASS' and g.get('simple_methods_failed') is True,'Flow gate closed'
 return ph
if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('stage',choices=['test','flow']);x=a.parse_args();guard(x.stage)
 if x.stage=='test':subprocess.run([sys.executable,str(R/'experiments/evaluate_reference.py'),'--order','41','--split','test'],check=True)
 else:raise RuntimeError('Training implementation intentionally absent until scientific and simple-model gates pass')
