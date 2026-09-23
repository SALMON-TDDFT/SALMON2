import unittest
import numpy as np
from ace import ACE

class ACETests(unittest.TestCase):
 def setUp(self):
  rng=np.random.default_rng(7);self.dv=.3
  self.u=rng.normal(size=(2,3,9))+1j*rng.normal(size=(2,3,9))
  a=rng.normal(size=(2,9,9))+1j*rng.normal(size=(2,9,9))
  self.h=-a@a.conj().transpose(0,2,1)-np.eye(9)
  self.w=self.u@self.h.transpose(0,2,1)
 def test_interpolation_and_arbitrary_target_hermiticity(self):
  ace=ACE(self.u,self.w,self.dv)
  np.testing.assert_allclose(ace.apply(self.u),self.w,rtol=1e-12,atol=1e-11)
  identity=np.broadcast_to(np.eye(9),(2,9,9));h=ace.apply(identity)
  np.testing.assert_allclose(h,h.conj().transpose(0,2,1),atol=1e-11)
  self.assertLess(np.linalg.eigvalsh(h).max(),1e-10)
 def test_complex_basis_covariance(self):
  v=np.linalg.qr(self.u[0,:,:3])[0]
  a=ACE(self.u,self.w,self.dv);b=ACE(v@self.u,v@self.w,self.dv)
  np.testing.assert_allclose(a.apply(self.u),b.apply(self.u),atol=1e-10)
 def test_invalid_input_rejected(self):
  for w in [self.w*1j,-self.w,np.zeros_like(self.w),self.w*np.nan]:
   with self.assertRaises(ValueError):ACE(self.u,w,self.dv)
  with self.assertRaises(ValueError):ACE(self.u,self.w,-1.)

if __name__=='__main__':unittest.main()
