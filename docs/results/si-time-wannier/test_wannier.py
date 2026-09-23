import unittest
import numpy as np
from wannier import projected_gauge,reconstruct,moments

class WannierTests(unittest.TestCase):
 def test_projected_gauge_rotation_invariance(self):
  rng=np.random.default_rng(3);u=np.linalg.qr(rng.normal(size=(40,3))+1j*rng.normal(size=(40,3)))[0]
  trial=rng.normal(size=(40,3))+1j*rng.normal(size=(40,3));v=np.linalg.qr(rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)))[0]
  g,_=projected_gauge(u,trial,1);h,_=projected_gauge(u@v,trial,1)
  np.testing.assert_allclose(u@g,u@v@h,atol=1e-13)
  np.testing.assert_allclose(g.conj().T@g,np.eye(3),atol=1e-13)
  with self.assertRaises(ValueError):projected_gauge(u,np.zeros_like(trial),1)
 def test_supercell_normalization_and_density(self):
  # Four half-shifted k points, orthonormal primitive-cell delta orbitals.
  k=np.array([[-.375,0,0],[-.125,0,0],[.125,0,0],[.375,0,0]])*2*np.pi
  r=np.array([[0.,0,0],[.5,0,0]])
  cells=np.array([[i,0,0] for i in range(4)])
  u=np.tile(np.eye(2)[None,:,:],(4,1,1)).astype(complex)
  w=reconstruct(u,k,r,cells)
  np.testing.assert_allclose(np.einsum('srn,srm->nm',w.conj(),w),np.eye(2),atol=1e-13)
  np.testing.assert_allclose(np.sum(abs(w)**2,axis=(0,2)),np.mean(np.sum(abs(u)**2,axis=2),axis=0),atol=1e-13)
 def test_periodic_center_and_variance(self):
  # A packet crossing the periodic boundary must remain centered at zero.
  pos=np.array([[.1,1,1],[9.9,1,1]]);rho=np.array([[.5],[.5]])
  m=moments(rho,pos,np.array([10.,10.,10.]),1.)
  np.testing.assert_allclose(np.minimum(m['center'][0],10-m['center'][0]),[0,1,1],atol=1e-12)
  np.testing.assert_allclose(m['spread'][0],.01,atol=1e-12)

if __name__=='__main__':unittest.main()
