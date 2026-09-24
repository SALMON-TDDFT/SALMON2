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
   cmd=['mpifort','-O2','-fcheck=all','-I'+str(fftw/'include'),str(ROOT/'src/xc/hse_exchange.f90'),str(ROOT/'samples/hse_mlwf_reference/native_distributed_exchange_probe.f90'),'-L'+str(fftw/'lib'),'-lfftw3','-L'+str(blas/'lib'),'-lopenblas','-o',str(exe)]
   result=subprocess.run(cmd,cwd=d,capture_output=True,text=True)
   self.assertEqual(result.returncode,0,result.stderr)
   rng=np.random.default_rng(804)
   for n in (2,3):
    m=2;h=.7;k=np.array(list(np.ndindex((m,)*3)))*2*np.pi/(m*n*h)+[.03,-.02,.01]
    k=k[rng.permutation(len(k))]
    u=rng.normal(size=(8,2,n,n,n))+1j*rng.normal(size=(8,2,n,n,n))
    t=rng.normal(size=(8,3,n,n,n))+1j*rng.normal(size=(8,3,n,n,n))
    ref=packed(DistanceExchange((n,)*3,h,k).apply(u,t)[0]);write_fixture(d/'in',u,t,k,h)
    for np_ in (1,3,8):
     env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
     r=subprocess.run(['mpiexec','-n',str(np_),str(exe),str(d/'in'),str(d/'out')],env=env,capture_output=True,text=True,timeout=90)
     self.assertEqual(r.returncode,0,r.stderr)
     out=np.fromfile(d/'out',np.complex128).reshape(ref.shape,order='F')
     np.testing.assert_allclose(out,ref,rtol=3e-12,atol=3e-12)

   for invalid in ('nan','shape'):
    r=subprocess.run(['mpiexec','-n','3',str(exe),str(d/'in'),str(d/'out'),invalid],env=env,capture_output=True,text=True,timeout=30)
    self.assertEqual(r.returncode,0,r.stderr)
