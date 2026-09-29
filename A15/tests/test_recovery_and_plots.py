import sys,unittest,tempfile,shutil,json,hashlib,importlib.util,contextlib,io
from pathlib import Path
import numpy as np
R=Path(__file__).resolve().parents[1]
def load(name):
 s=importlib.util.spec_from_file_location(name,R/'experiments'/f'{name}.py');m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
class Recovery(unittest.TestCase):
 def test_checkpoint_resume_and_hash_guard(self):
  m=load('generate_table')
  with tempfile.TemporaryDirectory() as tmp:
   root=Path(tmp);shutil.copytree(R/'src',root/'src');shutil.copytree(R/'configs',root/'configs');(root/'results').mkdir();m.ROOT=root
   old=sys.argv;sys.argv=['generate_table','--order','3']
   try:
    with contextlib.redirect_stdout(io.StringIO()):m.main()
    folder=root/'results/table_3_coupled';ref=np.load(folder/'y.npy').copy();st=json.loads((folder/'status.json').read_text());st['status']='running';st['completed']=20;(folder/'status.json').write_text(json.dumps(st))
    y=np.load(folder/'y.npy',mmap_mode='r+');y[20:]=0;y.flush();del y
    with contextlib.redirect_stdout(io.StringIO()):m.main()
    np.testing.assert_array_equal(np.load(folder/'y.npy'),ref)
    (root/'src/cylinders.py').write_text((root/'src/cylinders.py').read_text()+'\n# mutation\n')
    with self.assertRaises(AssertionError):m.main()
   finally:sys.argv=old
 def test_figure_same_query_observation(self):
  m=load('plot_results');table=m._load_table();query=m._load_queries();path,row=m._representative_query();theta,w,q,*_=m._posterior_for_query(table,query,row['query'])
  self.assertEqual(row['seed'],11);self.assertEqual(row['noise'],.02)
  self.assertAlmostEqual(float(w@q),row['mean'],places=9)
  self.assertAlmostEqual(float(w@(q-w@q)**2),row['variance'],places=9)
 def test_guard_rejects_test_and_flow(self):
  m=load('run_stage')
  for stage in ['test','flow']:
   with self.assertRaises(AssertionError):m.guard(stage)
if __name__=='__main__':unittest.main()
