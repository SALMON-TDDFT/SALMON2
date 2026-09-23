import unittest
import numpy as np
from mlwf import functional,rotate,minimize,overlap_mesh,transport_gauge

class MLWFTests(unittest.TestCase):
 def setUp(self):
  self.nb=np.array([[1,2],[2,0],[0,1]])
  self.b=np.array([[1.,0,0],[-1.,0,0]])
  self.w=np.array([.5,.5])
  rng=np.random.default_rng(2)
  self.u=np.array([np.linalg.qr(rng.normal(size=(2,2))+1j*rng.normal(size=(2,2)))[0] for _ in range(3)])
  self.m=np.tile((np.eye(2)*.8)[None,None,:,:],(3,2,1,1)).astype(complex)
 def test_analytic_gradient(self):
  f,d,_=functional(self.u,self.m,self.nb,self.b,self.w)
  eps=1e-6
  fp= functional(rotate(self.u,d,eps),self.m,self.nb,self.b,self.w)[0]
  fm= functional(rotate(self.u,d,-eps),self.m,self.nb,self.b,self.w)[0]
  self.assertAlmostEqual((fp-fm)/(2*eps),-np.sum(abs(d)**2),places=6)
 def test_spread_decomposition_and_variable_gradient(self):
  f,d,info=functional(self.u,self.m,self.nb,self.b,self.w)
  self.assertAlmostEqual(f,info['omega_invariant']+info['variable_spread'],places=12)
  rng=np.random.default_rng(17)
  z=rng.normal(size=self.u.shape)+1j*rng.normal(size=self.u.shape)
  direction=(z-z.conj().transpose(0,2,1))/2
  eps=1e-6
  fp=functional(rotate(self.u,direction,eps),self.m,self.nb,self.b,self.w)[2]['variable_spread']
  fm=functional(rotate(self.u,direction,-eps),self.m,self.nb,self.b,self.w)[2]['variable_spread']
  self.assertAlmostEqual((fp-fm)/(2*eps),-np.vdot(d,direction).real,places=7)

 def test_known_localized_minimum(self):
  u,log=minimize(self.u,self.m,self.nb,self.b,self.w,maxiter=500,tolerance=1e-7)
  f,_,_=functional(u,self.m,self.nb,self.b,self.w)
  self.assertLess(f,1.)
  self.assertLess(log['gradient_norm'],1e-5)
  np.testing.assert_allclose(u.conj().transpose(0,2,1)@u,np.tile(np.eye(2),(3,1,1)),atol=1e-12)

 def test_temporal_transport_removes_basis_rotation(self):
  rng=np.random.default_rng(25)
  previous=np.array([np.linalg.qr(rng.normal(size=(7,2))+1j*rng.normal(size=(7,2)))[0] for _ in range(3)])
  rotation=np.array([np.linalg.qr(rng.normal(size=(2,2))+1j*rng.normal(size=(2,2)))[0] for _ in range(3)])
  current=previous@rotation
  old_wannier=previous@self.u
  gauge,singular=transport_gauge(current,old_wannier,1.)
  np.testing.assert_allclose(current@gauge,old_wannier,atol=2e-15)
  np.testing.assert_allclose(singular,1,atol=2e-15)

 def test_periodic_boundary_phase(self):
  import itertools
  k=np.array(list(itertools.product([-.375,-.125,.125,.375],repeat=3)))*2*np.pi
  r=np.array([[.2,.3,.1]])
  orbitals=np.exp(-1j*k@r.T)[:,:,None]
  raw,neighbors,b,weights=overlap_mesh(orbitals,k,r,1.,1.)
  u=np.ones((64,1,1),complex)
  f,d,info=functional(u,raw,neighbors,b,weights)
  self.assertAlmostEqual(f,0.,places=12)
  np.testing.assert_allclose(info['centers'],r,atol=1e-12)
  np.testing.assert_allclose(d,0,atol=1e-12)
 def test_occupied_gauge_covariance(self):
  rng=np.random.default_rng(9)
  v=np.array([np.linalg.qr(rng.normal(size=(2,2))+1j*rng.normal(size=(2,2)))[0] for _ in range(3)])
  raw=v.conj().transpose(0,2,1)[:,None]@self.m@v[self.nb]
  uu=v.conj().transpose(0,2,1)@self.u
  f,_,info=functional(self.u,self.m,self.nb,self.b,self.w)
  g,_,other=functional(uu,raw,self.nb,self.b,self.w)
  self.assertAlmostEqual(f,g,places=12)
  self.assertAlmostEqual(info['omega_invariant'],other['omega_invariant'],places=12)

if __name__=='__main__':unittest.main()
