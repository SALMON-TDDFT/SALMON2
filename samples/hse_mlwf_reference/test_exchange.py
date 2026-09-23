import unittest
import numpy as np
try:
 import exchange as ex
except ImportError:
 ex=None

class ExchangeTests(unittest.TestCase):
 def setUp(self):
  self.assertIsNotNone(ex, 'screened exchange module is not implemented')
 def test_kernel_zero_and_fourier_mode(self):
  shape=(6,6,6);h=.7;omega=.11;k=ex.ScreenedKernel(shape,h,omega)
  self.assertAlmostEqual(k.multiplier[0,0,0],np.pi/omega**2)
  rho=np.exp(2j*np.pi*np.arange(6)[:,None,None]/6)*np.ones(shape)
  q=2*np.pi/(6*h);answer=4*np.pi*(-np.expm1(-q*q/(4*omega*omega)))/q**2
  np.testing.assert_allclose(k.convolve(rho),answer*rho,atol=2e-13)
 def test_action_hermiticity_and_energy_gradient(self):
  rng=np.random.default_rng(3);shape=(4,4,4);h=.8;dv=h**3
  a=rng.normal(size=(64,3))+1j*rng.normal(size=(64,3));a=np.linalg.qr(a)[0]/np.sqrt(dv)
  w=a.T.reshape((3,)+shape);kernel=ex.ScreenedKernel(shape,h,.3)
  v=(rng.normal(size=(2,)+shape)+1j*rng.normal(size=(2,)+shape))*.1
  kv=ex.apply_exchange(w,v,kernel)
  np.testing.assert_allclose(np.vdot(v[0],kv[1]),np.vdot(kv[0],v[1]),atol=1e-12)
  kw=ex.apply_exchange(w,w,kernel);energy=ex.closed_shell_energy(w,kw,dv)
  self.assertLess(energy,0)
  d=rng.normal(size=w.shape)+1j*rng.normal(size=w.shape);step=1e-6
  def f(x):return ex.closed_shell_energy(x,ex.apply_exchange(x,x,kernel),dv)
  derivative=(f(w+step*d)-f(w-step*d))/(2*step)
  self.assertAlmostEqual(derivative,4*np.vdot(d,kw).real*dv,places=7)
 def test_occupied_unitary_covariance(self):
  rng=np.random.default_rng(4);shape=(4,4,4);h=.8
  w=(rng.normal(size=(3,)+shape)+1j*rng.normal(size=(3,)+shape))*.05
  u=np.linalg.qr(rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)))[0]
  rotated=np.einsum('ij,jxyz->ixyz',u,w);kernel=ex.ScreenedKernel(shape,h,.2)
  np.testing.assert_allclose(ex.apply_exchange(rotated,rotated,kernel),np.einsum('ij,jxyz->ixyz',u,ex.apply_exchange(w,w,kernel)),atol=1e-13)
 def test_translation_pair_selection_and_zero_threshold(self):
  rng=np.random.default_rng(5);shape=(6,6,6);h=.8
  w=(rng.normal(size=(2,)+shape)+1j*rng.normal(size=(2,)+shape))*.03
  shifts=np.array([[0,0,0],[3,0,0]])
  kernel=ex.ScreenedKernel(shape,h,.3)
  sources=np.array([np.roll(a,tuple(s),(0,1,2)) for a in w for s in shifts])
  expected=ex.apply_exchange(sources,w,kernel)
  got,stats=ex.localized_exchange(w,shifts,kernel,pair_tolerance=0.)
  np.testing.assert_allclose(got,expected,atol=1e-13)
  self.assertEqual(stats['retained_pairs'],8)
 def test_primitive_bloch_matches_wannier_supercell(self):
  # A shifted 2x2x2 mesh checks boundary twist and FFT alias conventions.
  rng=np.random.default_rng(8);n=4;mesh=2;h=.9;length=n*h
  axes=(np.arange(mesh)+.5-mesh/2)*2*np.pi/(mesh*length)
  k=np.stack(np.meshgrid(axes,axes,axes,indexing='ij'),axis=-1).reshape(-1,3)
  u=(rng.normal(size=(8,1,n,n,n))+1j*rng.normal(size=(8,1,n,n,n)))*.05
  w,twist=ex.bloch_to_wannier(u,k,h)
  shifts=ex.cell_shifts(mesh,n);kernel=ex.ScreenedKernel((n*mesh,)*3,h,.3)
  action,_=ex.localized_exchange(w,shifts,kernel,0.)
  back=ex.wannier_to_bloch(action,k,h,twist)
  expected=ex.bloch_exchange(u,k,h,.3)
  np.testing.assert_allclose(back,expected,atol=1e-11,rtol=1e-11)
 def test_absolute_closed_shell_normalization(self):
  shape=(4,4,4);h=.8;omega=.3;volume=np.prod(shape)*h**3
  w=np.ones((1,)+shape)/np.sqrt(volume);kernel=ex.ScreenedKernel(shape,h,omega)
  energy=ex.closed_shell_energy(w,ex.apply_exchange(w,w,kernel),h**3)
  self.assertAlmostEqual(energy,-np.pi/(omega**2*volume),places=13)
 def test_finite_threshold_translation_symmetry(self):
  rng=np.random.default_rng(33);w=rng.random((2,6,6,6))**8;shifts=ex.cell_shifts(2,3)
  kernel=ex.ScreenedKernel((6,6,6),.8,.3);scores=ex.pair_overlaps(w,shifts,.8**3);cut=.05
  got,stats=ex.localized_exchange(w,shifts,kernel,cut);manual=np.zeros_like(w,dtype=complex)
  for i in range(2):
   for j in range(2):
    for ir,shift in enumerate(shifts):
     inverse=np.where(np.all(shifts==(-shift)%6,axis=1))[0][0]
     self.assertEqual(scores[i,j,ir]>=cut,scores[j,i,inverse]>=cut)
     if scores[i,j,ir]>=cut:
      source=np.roll(w[i],tuple(shift),(0,1,2));manual[j]-=source*kernel.convolve(source.conj()*w[j])
  np.testing.assert_allclose(got,manual,atol=1e-13)
 def test_restricted_convolution_matches_periodic_kernel(self):
  self.assertTrue(hasattr(ex,'LocalConvolution'),'local convolution not implemented')
  rng=np.random.default_rng(32);kernel=ex.ScreenedKernel((12,12,12),.8,.2)
  rho=rng.normal(size=(4,4,4))+1j*rng.normal(size=(4,4,4));full=np.zeros(kernel.shape,complex);full[:4,:4,:4]=rho
  local=ex.LocalConvolution(kernel,4)
  np.testing.assert_allclose(local.convolve(rho),kernel.convolve(full)[:4,:4,:4],atol=1e-13)
 def test_local_source_action_and_hermiticity(self):
  self.assertTrue(hasattr(ex,'source_box_exchange'),'source box exchange not implemented')
  rng=np.random.default_rng(34);w=(rng.normal(size=(2,8,8,8))+1j*rng.normal(size=(2,8,8,8)))*.02
  kernel=ex.ScreenedKernel((8,8,8),.8,.2);origins=np.array([[6,7,5],[1,2,0]]);width=3
  sources=np.zeros_like(w)
  for i,origin in enumerate(origins):
   ix=np.ix_(*[(np.arange(width)+s)%8 for s in origin]);sources[(i,)+ix]=w[(i,)+ix]
  targets=(rng.normal(size=(2,8,8,8))+1j*rng.normal(size=(2,8,8,8)))*.02
  got=ex.source_box_exchange(w,targets,np.zeros((1,3),int),origins,width,kernel)
  np.testing.assert_allclose(got,ex.apply_exchange(sources,targets,kernel),atol=1e-13)
  np.testing.assert_allclose(np.vdot(targets[0],got[1]),np.vdot(got[0],targets[1]),atol=1e-13)
 def test_reciprocal_pair_reuse(self):
  self.assertTrue(hasattr(ex,'symmetric_exchange'),'symmetric pair exchange missing')
  rng=np.random.default_rng(42);w=(rng.normal(size=(2,8,8,8))+1j*rng.normal(size=(2,8,8,8)))*.03
  shifts=ex.cell_shifts(2,4);kernel=ex.ScreenedKernel((8,8,8),.8,.2)
  for threshold in [0.,.05]:
   expected,_=ex.localized_exchange(w,shifts,kernel,threshold)
   got,s=ex.symmetric_exchange(w,shifts,kernel,threshold)
   np.testing.assert_allclose(got,expected,atol=1e-13)
if __name__=='__main__':unittest.main()
