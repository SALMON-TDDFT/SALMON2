import unittest
import numpy as np
import tdelf
class Tests(unittest.TestCase):
 def test_single_orbital(self):
  rng=np.random.default_rng(2);u=rng.normal(size=(1,20,1))+1j*rng.normal(size=(1,20,1));g=rng.normal(size=(3,1,20,1))+1j*rng.normal(size=(3,1,20,1))
  z=tdelf.fields(u,g)
  np.testing.assert_allclose(z['elf'],1,atol=1e-12)
 def test_gauge_and_orbital_invariance(self):
  rng=np.random.default_rng(8);u=rng.normal(size=(2,20,3))+1j*rng.normal(size=(2,20,3));g=rng.normal(size=(3,2,20,3))+1j*rng.normal(size=(3,2,20,3))
  z=tdelf.fields(u,g)
  boost=np.array([.7,-.2,.4]);zz=tdelf.fields(u,g+1j*boost[:,None,None,None]*u)
  np.testing.assert_allclose(z['elf'],zz['elf'],atol=1e-14)
  v=np.linalg.qr(rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)))[0]
  np.testing.assert_allclose(z['elf'],tdelf.fields(u@v,g@v)['elf'],atol=1e-14)
 def test_plane_wave_derivative(self):
  n=6;L=4.;r=np.stack(np.meshgrid(*(np.arange(n)*L/n for _ in range(3)),indexing='ij'),axis=-1).reshape(-1,3,order='F')
  k=np.array([[.1,.2,-.3]]);gv=2*np.pi/L*np.array([1,-1,0]);u=np.exp(1j*r@gv)[None,:,None]
  grad=tdelf.gradients(u,k,L,n)
  np.testing.assert_allclose(grad,1j*(k[0]+gv)[:,None,None,None]*u,atol=1e-14)
 def test_real_nyquist_derivative_option(self):
  n=4;L=4.
  x=np.stack(np.meshgrid(*(np.arange(n) for _ in range(3)),indexing='ij'),axis=-1).reshape(-1,3,order='F')
  u=(-1.)**x[:,0];u=u[None,:,None].astype(complex)
  grad=tdelf.gradients(u,np.zeros((1,3)),L,n,zero_nyquist=True)
  np.testing.assert_allclose(grad,0,atol=1e-14)
 def test_uniform_gas_formula(self):
  # Six opposite momenta, chosen to give the continuum cold-gas tau at n=6.
  n=6.;p=np.sqrt(3/5*(6*np.pi**2*n)**(2/3))
  q=np.vstack([np.eye(3),-np.eye(3)])*p
  u=np.ones((6,4,1),complex)*np.sqrt(6)
  g=1j*q.T[:,:,None,None]*u
  z=tdelf.fields(u,g)
  np.testing.assert_allclose(z['elf'],.5,atol=1e-14)
if __name__=='__main__':unittest.main()
