"""Finite-range overlap-save versus direct periodic complex-kernel convolution."""
import argparse,math,os,subprocess,tempfile
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--build',type=Path,required=True);a=p.parse_args()
b=a.build.resolve();obj=b/'src/CMakeFiles/salmon.dir'
with tempfile.TemporaryDirectory() as tmp:
 exe=Path(tmp)/'probe'
 subprocess.run(['mpifort','-fopenmp','-fcheck=all','-I'+str(b),str(Path(__file__).with_name('probe.f90')),
  *[str(obj/('xc/'+x+'.f90.o')) for x in ['exx_spatial_local','exx_local_fft','exx_batch_backend','fftw_blocks','fftw_pencils']],
  str(obj/'parallel/communication.f90.o'),str(obj/'misc/nvtx_wrapper.f90.o'),
  '-L/opt/homebrew/lib','-lfftw3','-o',str(exe)],check=True)
 for dims in [(1,1,1),(2,1,1),(3,1,1),(2,2,1),(2,2,2)]:
  for omp in (1,2):
   subprocess.run(['mpiexec','-n',str(math.prod(dims)),str(exe),' '.join(map(str,dims))],
    env=dict(os.environ,OMP_NUM_THREADS=str(omp)),check=True,timeout=120)
