"""Opt-in real-MPI kernel parity, including uneven k and idle row owners."""
import os, subprocess, tempfile, unittest
from pathlib import Path
import numpy as np
from distance_exchange import DistanceExchange
from test_native_exchange import packed, write_fixture, ROOT

@unittest.skipUnless(os.environ.get('SALMON_TEST_MPI') == '1', 'set SALMON_TEST_MPI=1 for real MPI tests')
class DistributedExchangeTest(unittest.TestCase):
 def test_distributed_arbitrary_targets(self):
  with tempfile.TemporaryDirectory() as tmp:
   d=Path(tmp); exe=d/'probe'
   fftw=Path(os.environ.get('FFTW_ROOT','/opt/homebrew/opt/fftw'))
   blas=Path(os.environ.get('OPENBLAS_ROOT','/opt/homebrew/opt/openblas'))
   cmd=['mpifort','-O2','-fopenmp','-fcheck=all','-I'+str(fftw/'include'),str(ROOT/'src/xc/hse_exchange.f90'),str(ROOT/'samples/hse_mlwf_reference/native_distributed_exchange_probe.f90'),'-L'+str(fftw/'lib'),'-lfftw3','-L'+str(blas/'lib'),'-lopenblas','-o',str(exe)]
   result=subprocess.run(cmd,cwd=d,capture_output=True,text=True)
   self.assertEqual(result.returncode,0,result.stderr)
   rng=np.random.default_rng(804)
   for n,m in ((2,2),(3,2),(2,4)):
    h=.7;k=np.array(list(np.ndindex((m,)*3)))*2*np.pi/(m*n*h)+[.03,-.02,.01]
    k=k[rng.permutation(len(k))]
    u=rng.normal(size=(m**3,2,n,n,n))+1j*rng.normal(size=(m**3,2,n,n,n))
    t=rng.normal(size=(m**3,3,n,n,n))+1j*rng.normal(size=(m**3,3,n,n,n))
    ref=packed(DistanceExchange((n,)*3,h,k).apply(u,t)[0]);write_fixture(d/'in',u,t,k,h)
    for layout in (None,'auto','strided','contiguous'):
     for np_,threads,block in ((1,1,2),(3,1,2),(8,1,2),(3,2,2),(3,4,2),(1,12,2),(3,1,5),(8,1,7)):
      env=dict(os.environ,OMP_NUM_THREADS=str(threads),OMP_DYNAMIC='FALSE',OPENBLAS_NUM_THREADS=str(threads),SALMON_HSE_BLOCK_ROWS=str(block),SALMON_HSE_PROFILE='1',SALMON_HSE_FFT_LAYOUT=layout or 'auto')
      if layout is None:env.pop('SALMON_HSE_FFT_LAYOUT',None)
      r=subprocess.run(['mpiexec','-n',str(np_),str(exe),str(d/'in'),str(d/'out')],env=env,capture_output=True,text=True,timeout=90)
      self.assertEqual(r.returncode,0,r.stderr)
      out=np.fromfile(d/'out',np.complex128).reshape(ref.shape,order='F')
      np.testing.assert_allclose(out,ref,rtol=3e-12,atol=3e-12)
   for invalid in ('nan','shape','block'):
    guard_env=dict(env)
    guard_env.pop('SALMON_HSE_BLOCK_ROWS',None)
    r=subprocess.run(['mpiexec','-n','3',str(exe),str(d/'in'),str(d/'out'),invalid],env=guard_env,capture_output=True,text=True,timeout=30)
    self.assertEqual(r.returncode,0,r.stderr)

   for invalid_block in ('0','-1','2oops','999999999999999999999999','9999'):
    bad_env=dict(env,SALMON_HSE_BLOCK_ROWS=invalid_block)
    r=subprocess.run(['mpiexec','-n','3',str(exe),str(d/'in'),str(d/'out')],env=bad_env,capture_output=True,text=True,timeout=30)
    self.assertNotEqual(r.returncode,0,invalid_block)

   bad_env=dict(env,SALMON_HSE_FFT_LAYOUT='invalid')
   r=subprocess.run(['mpiexec','-n','3',str(exe),str(d/'in'),str(d/'out')],env=bad_env,capture_output=True,text=True,timeout=30)
   self.assertNotEqual(r.returncode,0)
