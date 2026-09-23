import unittest
import tempfile
import hashlib
from pathlib import Path
import numpy as np
import exchange as ex
try:
 from local_support import LocalSupportExchange
except ImportError:
 LocalSupportExchange = None


class LocalSupportTests(unittest.TestCase):
 def setUp(self):
  self.assertIsNotNone(LocalSupportExchange, 'variational local support is not implemented')
  self.rng=np.random.default_rng(73)
  self.kernel=ex.ScreenedKernel((8,8,8),.8,.3)
  self.w=(self.rng.normal(size=(2,8,8,8))+1j*self.rng.normal(size=(2,8,8,8)))*.03
  self.shifts=ex.cell_shifts(2,4)

 def reference(self,w,origins,width):
  mask=np.zeros(w.shape)
  for i,o in enumerate(origins):
   ix=np.ix_(*[(np.arange(width)+s)%8 for s in o]);mask[(i,)+ix]=1
  v=mask*w
  sources=np.array([np.roll(a,tuple(s),(0,1,2)) for a in v for s in self.shifts])
  return mask*ex.apply_exchange(sources,v,self.kernel)

 def test_full_and_wrapped_support(self):
  for width in (3,5,8):
   origins=np.array([[7,6,1],[1,0,7]])
   with LocalSupportExchange(self.kernel,self.shifts,origins,width) as op:
    got,stats=op.apply(self.w)
   np.testing.assert_allclose(got,self.reference(self.w,origins,width),atol=2e-13)
   self.assertLessEqual(stats['convolutions'],24)

 def test_energy_gradient_with_fixed_domains(self):
  origins=np.array([[7,6,1],[1,0,7]])
  with LocalSupportExchange(self.kernel,self.shifts,origins,5) as op:
   a,_=op.apply(self.w);d=self.rng.normal(size=self.w.shape)+1j*self.rng.normal(size=self.w.shape)
   def energy(w):return ex.closed_shell_energy(w,op.apply(w)[0],self.kernel.dv)
   eps=1e-6;fd=(energy(self.w+eps*d)-energy(self.w-eps*d))/(2*eps)
   self.assertAlmostEqual(fd,4*self.kernel.dv*np.vdot(d,a).real,places=8)

 def test_empty_pairs_skipped(self):
  origins=np.array([[0,0,0],[4,4,4]])
  with LocalSupportExchange(self.kernel,np.zeros((1,3),int),origins,2) as op:
   a,s=op.apply(self.w)
   self.assertEqual(s['convolutions'],2)
  self.assertTrue(np.isfinite(a).all())

 def test_validation(self):
  for width in (0,9,2.5):
   with self.assertRaises(ValueError):LocalSupportExchange(self.kernel,self.shifts,[[0,0,0]]*2,width)
  with self.assertRaises(ValueError):LocalSupportExchange(self.kernel,[[0,0,0],[1,0,0]],[[0,0,0]]*2,3)

 def test_fftw_parity(self):
  origins=np.array([[7,6,1],[1,0,7]])
  with LocalSupportExchange(self.kernel,self.shifts,origins,3,fftw=True) as op:
   a,_=op.apply(self.w)
  np.testing.assert_allclose(a,self.reference(self.w,origins,3),atol=2e-13)

 def test_nonself_inverse_translations_and_gradient(self):
  kernel=ex.ScreenedKernel((9,9,9),.8,.3);shifts=ex.cell_shifts(3,3)
  w=(self.rng.normal(size=(2,9,9,9))+1j*self.rng.normal(size=(2,9,9,9)))*.02
  origins=np.array([[8,7,0],[2,0,8]]);mask=np.zeros_like(w.real)
  for i,o in enumerate(origins):mask[(i,)+np.ix_(*[(np.arange(4)+s)%9 for s in o])]=1
  sources=np.array([np.roll(a,tuple(s),(0,1,2)) for a in w*mask for s in shifts])
  with LocalSupportExchange(kernel,shifts,origins,4) as op:
   action,_=op.apply(w)
   np.testing.assert_allclose(action,mask*ex.apply_exchange(sources,w*mask,kernel),atol=2e-13)
   d=self.rng.normal(size=w.shape)+1j*self.rng.normal(size=w.shape);eps=1e-6
   def energy(x):return ex.closed_shell_energy(x,op.apply(x)[0],kernel.dv)
   self.assertAlmostEqual((energy(w+eps*d)-energy(w-eps*d))/(2*eps),4*kernel.dv*np.vdot(d,action).real,places=8)

 def test_orbital_specific_support_is_not_automatically_ace_compatible(self):
  from ace import ACE
  with LocalSupportExchange(self.kernel,[[0,0,0]],[[7,6,1],[1,0,7]],3) as op:
   action,_=op.apply(self.w)
  # A real energy gradient need not have a Hermitian occupied-space metric.
  # Keep this rejection: symmetrizing the metric would hide a model change.
  with self.assertRaisesRegex(ValueError,'Hermitian'):ACE(self.w[None],action[None],self.kernel.dv)

 def test_benchmark_fingerprint_matches_loaded_archive(self):
  import benchmark_local_support as bench
  self.assertTrue(hasattr(bench,'load_snapshot'),'atomic benchmark snapshot reader missing')
  with tempfile.TemporaryDirectory() as directory:
   path=Path(directory)/'state.npz';g=np.eye(2)[None]
   np.savez(path,u=self.w,gauge=g)
   expected=hashlib.sha256(path.read_bytes()).hexdigest()
   u,gauge,digest=bench.load_snapshot(path)
   np.testing.assert_array_equal(u,self.w);np.testing.assert_array_equal(gauge,g)
   self.assertEqual(digest,expected)

if __name__=='__main__':unittest.main()
