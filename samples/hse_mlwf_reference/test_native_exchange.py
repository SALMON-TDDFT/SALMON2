"""Compiled Fortran vs independent Python exchange; optional local toolchain paths."""
import os,subprocess,tempfile,unittest,shutil
from pathlib import Path
import numpy as np
from distance_exchange import DistanceExchange
from ace import ACE
ROOT=Path(__file__).resolve().parents[2]

def packed(u):
 return u.transpose(2,3,4,1,0).reshape((-1,u.shape[1],u.shape[0]),order='F')

def write_fixture(path,u,t,k,h,omega=.11):
 with Path(path).open('wb') as f:
  f.write(np.array([u.shape[2],round(len(k)**(1/3)),u.shape[1],t.shape[1]],np.int32).tobytes())
  f.write(np.array([h,omega],np.float64).tobytes())
  f.write(np.asarray(k.T,order='F').tobytes(order='F'))
  for a in (u,t):f.write(packed(a).tobytes(order='F'))

class NativeExchangeTest(unittest.TestCase):
 @classmethod
 def setUpClass(cls):
  cls.tmp=tempfile.TemporaryDirectory();cls.directory=Path(cls.tmp.name);cls.exe=cls.directory/'probe'
  source=ROOT/'src/xc/hse_exchange.f90'
  if not source.exists():raise AssertionError('Native HSE exchange module has not been implemented')
  if not (ROOT/'src/xc/hse_ace.f90').exists():raise AssertionError('Native ACE has not been implemented')
  fftw=Path(os.environ.get('FFTW_ROOT','/opt/homebrew/opt/fftw'))
  if not shutil.which(os.environ.get('FC','gfortran')) or not (fftw/'include/fftw3.f03').exists():
   raise unittest.SkipTest('Native kernel tests require gfortran and FFTW_ROOT')
  blas=Path(os.environ.get('OPENBLAS_ROOT','/opt/homebrew/opt/openblas'))
  cmd=[os.environ.get('FC','gfortran'),'-O2','-fcheck=all','-I'+str(fftw/'include'),str(source),str(ROOT/'src/xc/hse_ace.f90'),str(ROOT/'samples/hse_mlwf_reference/native_exchange_probe.f90'),'-L'+str(fftw/'lib'),'-lfftw3','-L'+str(blas/'lib'),'-lopenblas','-o',str(cls.exe)]
  result=subprocess.run(cmd,cwd=cls.directory,capture_output=True,text=True)
  if result.returncode:raise AssertionError(result.stderr)
 @classmethod
 def tearDownClass(cls):cls.tmp.cleanup()
 def test_ptcn_matches_reference(self):
  from ptcn import ptcn_step
  source=ROOT/'src/rt/hse_ptcn_core.f90'
  self.assertTrue(source.exists(),'Native PT-CN module has not been implemented')
  exe=self.directory/'ptcn'
  blas=Path(os.environ.get('OPENBLAS_ROOT','/opt/homebrew/opt/openblas'))
  result=subprocess.run([os.environ.get('FC','gfortran'),str(source),str(ROOT/'samples/hse_mlwf_reference/native_ptcn_probe.f90'),'-L'+str(blas/'lib'),'-lopenblas','-o',str(exe)],cwd=self.directory,capture_output=True,text=True)
  self.assertEqual(result.returncode,0,result.stderr)
  values=np.loadtxt(subprocess.run([str(exe)],check=True,capture_output=True,text=True).stdout.splitlines())
  u=np.array([1,.2+.1j,-.1+.3j]).reshape(1,1,3);u/=np.linalg.norm(u)
  matrix=np.array([[.2,.1,0],[.1,.7,0],[0,0,1.1]])
  def local(x):return x@matrix+.05*abs(x)**2*x
  def exchange(x):return (lambda a:np.zeros_like(a)),np.zeros_like(x)
  expected,_=ptcn_step(u,.1,local,exchange,1.)
  np.testing.assert_allclose(values[:,0]+1j*values[:,1],expected.ravel(),rtol=1e-13,atol=1e-13)
 def test_semilocal_matches_libxc_reference(self):
  from semilocal import Semilocal
  source=ROOT/'src/xc/hse_semilocal.f90'
  self.assertTrue(source.exists(),'Native HSE semilocal module has not been implemented')
  lib=Path(os.environ.get('LIBXC_ROOT','/opt/homebrew/opt/libxc'))
  exe=self.directory/'semilocal'
  subprocess.run([os.environ.get('FC','gfortran'),str(source),str(ROOT/'samples/hse_mlwf_reference/native_semilocal_probe.f90'),'-L'+str(lib/'lib'),'-lxc','-o',str(exe)],cwd=self.directory,check=True,capture_output=True)
  values=np.loadtxt(subprocess.run([str(exe)],check=True,capture_output=True,text=True).stdout.splitlines())
  with Semilocal() as xc:ref=np.array(xc.evaluate([0,1e-4,.1,2],[0,.001,.2,3])).T
  np.testing.assert_allclose(values,ref,rtol=2e-14,atol=1e-14)
 def test_shifted_shuffled_mesh_and_row_partition(self):
  rng=np.random.default_rng(792);n=3;mesh=2;h=.7
  k=np.array(list(np.ndindex((mesh,)*3)))*2*np.pi/(mesh*n*h)+[.031,-.02,.017]
  k=k[rng.permutation(len(k))]
  u=rng.normal(size=(len(k),2,n,n,n))+1j*rng.normal(size=(len(k),2,n,n,n))
  t=rng.normal(size=(len(k),3,n,n,n))+1j*rng.normal(size=(len(k),3,n,n,n))
  ref=packed(DistanceExchange((n,)*3,h,k).apply(u,t)[0]);path=self.directory/'input.bin';write_fixture(path,u,t,k,h)
  for size in (1,3,9):
   total=np.zeros_like(ref)
   for rank in range(size):
    out=self.directory/'output.bin'
    subprocess.run([str(self.exe),str(path),str(out),str(rank),str(size)],check=True,capture_output=True)
    total+=np.fromfile(out,np.complex128).reshape(ref.shape,order='F')
   np.testing.assert_allclose(total,ref,rtol=2e-12,atol=2e-12)
   if size==1:
    w=DistanceExchange((n,)*3,h,k).apply(u,u)[0]
    ace_ref=packed(ACE(u,w,h**3).apply(t))
    ace_out=np.fromfile(str(out)+'.ace',np.complex128).reshape(ref.shape,order='F')
    np.testing.assert_allclose(ace_out,ace_ref,rtol=3e-12,atol=3e-12)
    interpolation=np.fromfile(str(out)+'.occupied',np.complex128).reshape(packed(u).shape,order='F')
    np.testing.assert_allclose(interpolation,packed(w),rtol=3e-12,atol=3e-12)

if __name__=='__main__':unittest.main()
