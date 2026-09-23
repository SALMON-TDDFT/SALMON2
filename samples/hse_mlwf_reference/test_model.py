import unittest
import numpy as np
try:
 from model import kinetic_symbol,gradient,hartree,semilocal_potential,orthonormalize
except ImportError:
 kinetic_symbol=None
class ModelTests(unittest.TestCase):
 def setUp(self):self.assertIsNotNone(kinetic_symbol,'reference Hamiltonian module missing')
 def test_kinetic_and_gradient_plane_wave(self):
  n=8;h=.7;k=np.array([.1,.2,-.3]);lap=np.tile(np.array([1.,0,0,0])[:,None]/h**2,(1,3));nab=np.tile(np.array([.5,0,0,0])[:,None]/h,(1,3))
  symbol=kinetic_symbol((n,)*3,h,k,3/h**2,lap,nab)
  self.assertAlmostEqual(symbol[0,0,0],np.dot(k,k)/2)
  q=2*np.pi/(n*h);expected=(1-np.cos(q*h))/h**2+k[0]*np.sin(q*h)/h+np.dot(k,k)/2
  self.assertAlmostEqual(symbol[1,0,0],expected)
  rho=np.sin(2*np.pi*np.arange(n)[:,None,None]/n)*np.ones((n,n,n))
  g=gradient(rho,nab)
  np.testing.assert_allclose(g[0],np.cos(2*np.pi*np.arange(n)[:,None,None]/n)*np.ones((n,n,n))*np.sin(q*h)/h,atol=1e-14)
 def test_hartree_and_gga_energy_gradient(self):
  from semilocal import Semilocal
  rng=np.random.default_rng(19);shape=(6,6,6);h=.7;rho=.04+.002*rng.random(shape);d=rng.normal(size=shape)*.001
  nab=np.tile(np.array([.5,0,0,0])[:,None]/h,(1,3));delta=1e-4
  with Semilocal('hse06') as xc:
   v,e=semilocal_potential(rho,nab,xc,h**3)
   ep=semilocal_potential(rho+delta*d,nab,xc,h**3)[1];em=semilocal_potential(rho-delta*d,nab,xc,h**3)[1]
   self.assertAlmostEqual((ep-em)/(2*delta),np.sum(v*d)*h**3,places=9)
  v,e=hartree(rho,h);ep=hartree(rho+delta*d,h)[1];em=hartree(rho-delta*d,h)[1]
  self.assertAlmostEqual((ep-em)/(2*delta),np.sum(v*d)*h**3,places=10)
 def test_orthonormality_preserves_subspace(self):
  rng=np.random.default_rng(20);u=rng.normal(size=(2,3,4,4,4))+1j*rng.normal(size=(2,3,4,4,4));dv=.5
  v=orthonormalize(u,dv)
  for x in v:np.testing.assert_allclose(x.reshape(3,-1).conj()@x.reshape(3,-1).T*dv,np.eye(3),atol=1e-13)
if __name__=='__main__':unittest.main()
