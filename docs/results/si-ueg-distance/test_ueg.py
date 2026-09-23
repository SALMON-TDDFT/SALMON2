import unittest
import numpy as np
import ueg

class UEGTests(unittest.TestCase):
 def test_fermi_count_and_symmetry(self):
  q=ueg.momentum_cube(12,.2)
  f=ueg.fermi_reference(q,40)
  self.assertAlmostEqual(f.sum(),40)
  np.testing.assert_allclose(ueg.mean_momentum(f,q),0,atol=1e-15)
  self.assertTrue(np.all((f>=0)&(f<=1)))
 def test_fractional_shift_matches_number_and_flow(self):
  q=ueg.momentum_cube(12,.2);f=ueg.fermi_reference(q,40)
  target=np.array([.13,-.07,.04]);g=ueg.flow_reference(f,target,.2)
  self.assertAlmostEqual(g.sum(),40)
  np.testing.assert_allclose(ueg.mean_momentum(g,q),target,atol=1e-15)
  self.assertTrue(np.all((g>=0)&(g<=1)))
  with self.assertRaises(ValueError):ueg.flow_reference(f,[10,0,0],.2)
 def test_distance_matches_explicit_matrix_and_gauge(self):
  rng=np.random.default_rng(8)
  c=np.linalg.qr(rng.normal(size=(8,3))+1j*rng.normal(size=(8,3)))[0][None]
  f=np.array([[1,1,.5,.5,0,0,0,0.]])
  p=c[0]@c[0].conj().T
  expected=np.linalg.norm(p-np.diag(f[0]))**2
  self.assertAlmostEqual(ueg.hs_distance(c,f),expected,places=12)
  v=np.linalg.qr(rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)))[0]
  self.assertAlmostEqual(ueg.hs_distance(c@v,f),expected,places=12)
  self.assertAlmostEqual(ueg.hs_distance(np.eye(8)[:,:3][None],np.array([[1,1,1,0,0,0,0,0]])),0)
 def test_fft_plane_wave_and_indexing(self):
  n=4;L=2.;dv=(L/n)**3
  r=np.stack(np.meshgrid(*(np.arange(n)*L/n for _ in range(3)),indexing='ij'),axis=-1).reshape(-1,3,order='F')
  k=np.array([[-.25,-.25,-.25],[.25,.25,.25]])*2*np.pi/L
  g=np.array([1,-1,0]);u=np.tile((np.exp(2j*np.pi*r@g/L)/L**1.5)[None,:,None],(2,1,1))
  c=ueg.fourier_coefficients(u,n,dv)
  self.assertAlmostEqual(float(np.sum(abs(c[0])**2)),1)
  idx=np.unravel_index(np.argmax(abs(c[0,:,0])),(n,n,n),order='F')
  self.assertEqual(idx,(1,3,0))
  indices,q=ueg.global_indices(k,n,L,2)
  expected=k+2*np.pi/L*g
  np.testing.assert_allclose(q.reshape(-1,3)[indices[:,np.argmax(abs(c[0,:,0]))]],expected)
 def test_ensemble_endpoint_and_rank_lower_bound(self):
  f=np.array([[1.,.7,.3,0.]])
  c=np.diag(np.sqrt(f[0]))[None]
  self.assertAlmostEqual(ueg.hs_distance(c,f),0,places=12)
  projector=np.eye(4)[:,:2][None]
  self.assertAlmostEqual(ueg.hs_distance(projector,f),ueg.rank_lower_bound(f,2),places=12)
 def test_integer_boost_invariance(self):
  rng=np.random.default_rng(4)
  c=np.linalg.qr(rng.normal(size=(10,3))+1j*rng.normal(size=(10,3)))[0][None]
  f=np.array([[0.,0,0,1,1,1,0,0,0,0]])
  self.assertAlmostEqual(ueg.hs_distance(c,f),ueg.hs_distance(np.roll(c,2,axis=1),np.roll(f,2,axis=1)),places=12)
 def test_unclipped_alpha(self):
  np.testing.assert_allclose(ueg.alpha_from_distance(np.array([0.,2.,4.]),2.),[0,.2,.4])

if __name__=='__main__':unittest.main()
