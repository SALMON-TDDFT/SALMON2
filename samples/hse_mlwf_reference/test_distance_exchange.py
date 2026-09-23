import unittest
from types import SimpleNamespace
import numpy as np
from exchange import bloch_exchange
try:
 from distance_exchange import DistanceExchange
except ImportError:
 DistanceExchange=None


class DistanceExchangeTests(unittest.TestCase):
 def setUp(self):
  self.assertIsNotNone(DistanceExchange,'common distance exchange is not implemented')
  self.rng=np.random.default_rng(184);self.n=3;self.h=.8;self.mesh=2
  axes=(np.arange(self.mesh)+.3-self.mesh/2)*2*np.pi/(self.mesh*self.n*self.h)
  self.k=np.stack(np.meshgrid(axes,axes,axes,indexing='ij'),axis=-1).reshape(-1,3)
  self.k=self.k[self.rng.permutation(len(self.k))]
  self.u=self.random((8,2,3,3,3))*.08
 def random(self,shape):return self.rng.normal(size=shape)+1j*self.rng.normal(size=shape)
 def test_full_radius_matches_bloch(self):
  expected=bloch_exchange(self.u,self.k,self.h,.3)
  op=DistanceExchange((3,3,3),self.h,self.k,.3,block_rows=4)
  got,stats=op.apply(self.u,self.u)
  np.testing.assert_allclose(got,expected,atol=2e-12,rtol=2e-12)
  self.assertGreater(stats['workspace_bytes'],0)
 def test_fixed_source_hermiticity_linearity_and_gauge_covariance(self):
  op=DistanceExchange((3,3,3),self.h,self.k,.3,radius=1.7,block_rows=4)
  v=self.random(self.u.shape)*.1;z=self.random(self.u.shape)*.1
  kv=op.apply(self.u,v)[0];kz=op.apply(self.u,z)[0]
  np.testing.assert_allclose(np.vdot(z,kv),np.vdot(kz,v),atol=2e-12)
  np.testing.assert_allclose(op.apply(self.u,2*v+1j*z)[0],2*kv+1j*kz,atol=2e-12)
  unitary=np.array([np.linalg.qr(self.random((2,2)))[0] for _ in self.k])
  rotated=np.einsum('kab,kbxyz->kaxyz',unitary,self.u)
  np.testing.assert_allclose(op.apply(rotated,v)[0],kv,atol=2e-12)
 def test_energy_gradient(self):
  op=DistanceExchange((3,3,3),self.h,self.k,.3,radius=1.7,block_rows=4)
  d=self.random(self.u.shape);a=op.apply(self.u,self.u)[0];dv=self.h**3/len(self.k)
  def energy(u):return np.vdot(u,op.apply(u,u)[0]).real*dv
  eps=1e-6
  self.assertAlmostEqual((energy(self.u+eps*d)-energy(self.u-eps*d))/(2*eps),4*np.vdot(d,a).real*dv,places=8)
 def test_invalid_inputs(self):
  for radius in (-1,0,float('nan')):
   with self.assertRaises(ValueError):DistanceExchange((3,3,3),self.h,self.k,.3,radius=radius)
  with self.assertRaises(ValueError):DistanceExchange((3,3,3),self.h,self.k,.3,block_rows=0)
  op=DistanceExchange((3,3,3),self.h,self.k,.3)
  with self.assertRaises(ValueError):op.apply(self.u[:1],self.u)

 def test_onsite_limit_and_spatial_skipping(self):
  op=DistanceExchange((3,3,3),self.h,self.k,.3,radius=.4,block_rows=1)
  action,stats=op.apply(self.u,self.u)
  density=np.sum(abs(self.u)**2,axis=(0,1))/len(self.k)
  np.testing.assert_allclose(action,-op.real_kernel[0,0,0]*density*self.u,atol=2e-13)
  self.assertLess(stats['primitive_pair_fraction'],.1)

 def test_finite_cutoff_matches_direct_wannier_convolution(self):
  from exchange import (bloch_to_wannier,wannier_to_bloch,localized_exchange,
                        ScreenedKernel,cell_shifts)
  op=DistanceExchange((3,3,3),self.h,self.k,.3,radius=1.7)
  kernel=ScreenedKernel((6,6,6),self.h,.3);kernel.multiplier=op.multiplier
  w,twist=bloch_to_wannier(self.u,self.k,self.h)
  ref=localized_exchange(w,cell_shifts(2,3),kernel,0.)[0]
  np.testing.assert_allclose(op.apply(self.u,self.u)[0],wannier_to_bloch(ref,self.k,self.h,twist),atol=2e-12)

 def test_functional_optional_backend(self):
  from scf import HSEFunctional
  model=SimpleNamespace(shape=(3,3,3),h=self.h,nk=len(self.k))
  op=DistanceExchange(model.shape,self.h,self.k,.3)
  f=HSEFunctional(model,None,exchange_backend=op)
  try:
   np.testing.assert_allclose(f.exchange(self.u,None,0.)[0],op.apply(self.u,self.u)[0],atol=2e-13)
   with self.assertRaises(ValueError):f.exchange(self.u,None,1e-6)
  finally:f.close()

 def test_hse_rejects_indefinite_cutoff_kernel(self):
  from scf import HSEFunctional
  model=SimpleNamespace(shape=(3,3,3),h=self.h,nk=len(self.k))
  op=DistanceExchange(model.shape,self.h,self.k,.3,radius=1.7)
  self.assertLess(op.multiplier.min(),0.)
  with self.assertRaisesRegex(ValueError,'positive'):HSEFunctional(model,None,exchange_backend=op)

 def test_mesh_three_shifted_permuted(self):
  from exchange import (bloch_to_wannier,wannier_to_bloch,localized_exchange,
                        ScreenedKernel,cell_shifts)
  n=2;mesh=3
  axes=(np.arange(mesh)+.2-mesh/2)*2*np.pi/(mesh*n*self.h)
  k=np.stack(np.meshgrid(axes,axes,axes,indexing='ij'),axis=-1).reshape(-1,3)
  k=k[self.rng.permutation(len(k))];u=self.random((27,2,n,n,n))*.08
  for radius in (None,1.7):
   op=DistanceExchange((n,)*3,self.h,k,.3,radius=radius,block_rows=3)
   kernel=ScreenedKernel((n*mesh,)*3,self.h,.3);kernel.multiplier=op.multiplier
   w,twist=bloch_to_wannier(u,k,self.h)
   ref=localized_exchange(w,cell_shifts(mesh,n),kernel,0.)[0]
   np.testing.assert_allclose(op.apply(u,u)[0],wannier_to_bloch(ref,k,self.h,twist),atol=2e-12)

 def test_disjoint_row_partitions_recover_full_action(self):
  for radius in (None,1.7):
   op=DistanceExchange((3,3,3),self.h,self.k,.3,radius=radius,block_rows=4)
   reference=op.apply(self.u,self.u)[0]
   for size in (3,10):
    parts=[op.apply(self.u,self.u,row_rank=rank,row_size=size)[0] for rank in range(size)]
    np.testing.assert_allclose(sum(parts),reference,atol=0,rtol=0)
   for rank,size in ((-1,2),(2,2),(0,0),(0,1.5)):
    with self.assertRaises(ValueError):op.apply(self.u,self.u,row_rank=rank,row_size=size)

if __name__=='__main__':unittest.main()
