import unittest
import numpy as np
from model import hartree,kinetic_symbol
from semilocal import Semilocal
try:
 from scf import HSEFunctional
except ImportError:HSEFunctional=None
class SCFTests(unittest.TestCase):
 def test_total_energy_directional_derivative(self):
  self.assertIsNotNone(HSEFunctional,'HSE energy/action not implemented')
  class Toy:
   shape=(4,4,4);h=.8;dv=h**3;nk=8;nlcc=np.zeros(shape);native_energies=np.zeros(7)
   weights=np.ones(8)/8
   axes=(np.arange(2)+.5-1)*2*np.pi/(2*4*h)
   k=np.stack(np.meshgrid(axes,axes,axes,indexing='ij'),axis=-1).reshape(8,3)
   def density(self,u):return 2*np.einsum('k,knxyz->xyz',self.weights,abs(u)**2)
   def core(self,u):return .15*u
   def expectation(self,u,v):return float(2*np.einsum('k,knxyz,knxyz->',self.weights,u.conj(),v).real*self.dv)
   nab=np.tile(np.array([.5,0,0,0])[:,None]/h,(1,3))
  model=Toy();rng=np.random.default_rng(88);u=(rng.normal(size=(8,1,4,4,4))+1j*rng.normal(size=(8,1,4,4,4)))*.04
  d=(rng.normal(size=u.shape)+1j*rng.normal(size=u.shape))*.03;g=np.ones((8,1,1),complex)
  with Semilocal('hse06') as xc:
   functional=HSEFunctional(model,xc);hu,e,stats=functional.evaluate(u,g,0.)
   step=1e-5;ep=functional.evaluate(u+step*d,g,0.)[1]['total'];em=functional.evaluate(u-step*d,g,0.)[1]['total']
  self.assertAlmostEqual((ep-em)/(2*step),2*model.expectation(d,hu),places=7)
 def test_trial_acceptance_rejects_uphill_and_nonfinite(self):
  import scf
  self.assertTrue(hasattr(scf,'accept_trial'),'trial acceptance guard missing')
  self.assertTrue(scf.accept_trial(-2.,.001,np.ones(2),-1.,0.))
  self.assertFalse(scf.accept_trial(-.9,.001,np.ones(2),-1.,0.))
  for energy,residual,u in [(float('nan'),.1,np.ones(2)),(-2.,float('inf'),np.ones(2)),(-2.,.1,np.array([np.nan]))]:
   with self.assertRaises(FloatingPointError):scf.accept_trial(energy,residual,u,-1.,0.)
if __name__=='__main__':unittest.main()
