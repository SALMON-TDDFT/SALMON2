import unittest
import numpy as np
try:
 from rt import rk4_step
except ImportError:rk4_step=None
class RTTests(unittest.TestCase):
 def test_order_for_hermitian_hamiltonian(self):
  self.assertIsNotNone(rk4_step,'real-time propagation missing')
  u=np.array([.6,.8j]);energy=np.array([.7,1.1]);rhs=lambda t,x:-1j*energy*x
  errors=[]
  for n in [4,8]:
   x=u.copy()
   for j in range(n):x=rk4_step(x,j/n,1/n,rhs)
   errors.append(np.linalg.norm(x-u*np.exp(-1j*energy)))
  self.assertGreater(errors[0]/errors[1],15.)
  self.assertLess(errors[0]/errors[1],17.)
 def test_nonlinear_stages_are_not_frozen(self):
  self.assertIsNotNone(rk4_step,'real-time propagation missing')
  calls=[]
  def rhs(t,x):calls.append(x.copy());return x*x
  got=rk4_step(np.array([.2]),0.,.1,rhs)
  self.assertEqual(len(calls),4);self.assertNotEqual(calls[0][0],calls[1][0])
  self.assertAlmostEqual(got[0],.2/(1-.02),places=10)
 def test_invalid_pilot_request(self):
  import rt
  self.assertTrue(hasattr(rt,'validate_request'),'pilot argument validation missing')
  for args in [(0,.08,1e-4,0.),(1,0.,1e-4,0.),(1,.08,float('nan'),0.),(1,.08,1e-4,-1.)]:
   with self.assertRaises(ValueError):rt.validate_request(*args)
if __name__=='__main__':unittest.main()
