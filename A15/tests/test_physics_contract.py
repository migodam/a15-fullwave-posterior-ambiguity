import sys,unittest,hashlib,json
from pathlib import Path
R=Path(__file__).resolve().parents[1];sys.path.insert(0,str(R/'src'))
import numpy as np
from scipy.linalg import solve
from scipy.integrate import quad
from scipy.special import hankel1
from cylinders import CylinderModel
from physics import Family,cylinder_exact
from physics_backend import self_cell_green
class Contract(unittest.TestCase):
 def test_frozen_hash(self):
  self.assertEqual(hashlib.sha256((R/'configs/protocol.json').read_bytes()).hexdigest(),(R/'configs/protocol.sha256').read_text().strip())
 def test_no_measurement_in_readout(self):
  import inspect,physics
  self.assertNotIn('measurement',inspect.signature(physics.physical_readout).parameters)
 def test_self_cell_disk_integral(self):
  for k in [3.,8.625,12.]:
   a=1/32/np.sqrt(np.pi);f=lambda r:2*np.pi*r*.25j*hankel1(0,k*r)
   z=quad(lambda r:f(r).real,0,a,epsabs=1e-13)[0]+1j*quad(lambda r:f(r).imag,0,a,epsabs=1e-13)[0]
   self.assertLess(abs(z-self_cell_green(k,1/32)),1e-12)
 def test_reciprocity(self):
  m=CylinderModel(Family(radius=.1));s=m.components([.08,.5,.03],8.625)
  def far(alpha,beta):
   inc=np.exp(1j*8.625*(s['centers']@np.array([np.cos(beta),np.sin(beta)])))[:,None]*(1j*np.exp(-1j*beta))**m.m
   z=solve(s['M'],(s['inject']*inc).ravel());return np.sum(z.reshape(2,m.d)*s['radiation']*np.exp(-1j*8.625*(s['centers']@np.array([np.cos(alpha),np.sin(alpha)])))[:,None]*(-1j*np.exp(1j*alpha))**m.m)
  for a,b in [(.3,1.1),(.8,-.7),(2.1,.4)]:self.assertLess(abs(far(a,b)-far(b+np.pi,a+np.pi)),1e-10)
 def test_uncoupled_independent_cylinders(self):
  m=CylinderModel(Family(radius=.1));th=[.08,.5,.03];k=8.625;s=m.solve(th,k,coupling=False,fields=True)
  y=sum(cylinder_exact(m.receivers,k,m.family.radius,eps,center) for eps,center in zip(s['state']['epsilon'],s['state']['centers']))
  self.assertLess(np.linalg.norm(y-s['y'])/np.linalg.norm(y),1e-10)
 def test_reference_blocks_test(self):
  gate=json.loads((R/'results/reference_validation_gate.json').read_text());self.assertFalse(gate['passed'])
  self.assertFalse(any((R/'results').glob('queries_test*')))
 def test_physical_validation_adjacent_steps(self):
  j=json.loads((R/'results/cylinder_validation.json').read_text())
  for scene in j['scenes']:
   for row in scene.get('tangents',[]):
    for task in ['y','Q']:
     good=[x[task]<1e-5 for x in row['errors']]
     self.assertTrue(any(a and b for a,b in zip(good,good[1:])),(scene['scene'],row['direction'],task))
if __name__=='__main__':unittest.main()
