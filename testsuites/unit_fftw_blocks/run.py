"""Native Cartesian FFT and localized HSE/Coulomb action versus legacy reference."""
import argparse,os,subprocess,tempfile,shlex
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--build',type=Path,required=True);p.add_argument('--action',action='store_true')
p.add_argument('--slab',action='store_true')
p.add_argument('--mixed',action='store_true')
a=p.parse_args();b=a.build.resolve();obj=b/'src/CMakeFiles/salmon.dir'
names=['fftw_blocks','fftw_pencils']
libs=[]
if a.action:
 names+=['exx_ace','exx_blas_threads','exx_distributed_gauge','exx_sparse_orbitals','exx_distributed_metric','exx_orbitals',
  'exx_pair_candidates','exx_spatial','exx_batch_backend','exx_spatial_local','exx_wannier_gauge','exx_local_fft']
 libs=['-L/opt/homebrew/opt/openblas/lib','-lopenblas']
 if 'USE_SCALAPACK:BOOL=ON' in (b/'CMakeCache.txt').read_text():
  libs+=shlex.split(os.environ.get('SCALAPACK_LIBS','-L/opt/homebrew/opt/scalapack/lib -lscalapack'))
with tempfile.TemporaryDirectory(prefix='block-fft-') as tmp:
 exe=Path(tmp)/'probe'
 subprocess.run(['mpifort','-fopenmp','-fcheck=all','-I'+str(b),str(Path(__file__).with_name('action_probe.f90' if a.action else 'probe.f90')),
  *[str(obj/('xc/'+n+'.f90.o')) for n in names],str(obj/'parallel/communication.f90.o'),str(obj/'misc/nvtx_wrapper.f90.o'),
  '-L/opt/homebrew/lib','-lfftw3',*libs,'-o',str(exe)],check=True)
 for dims in [(1,1,1),(2,1,1),(2,2,1),(8,1,1),(2,2,2),(1,2,4)]:
  for threads in (1,2,4):
   subprocess.run(['mpiexec','-n',str(dims[0]*dims[1]*dims[2]),str(exe),' '.join(map(str,dims)),*(['mixed'] if a.mixed else ['slab'] if a.slab else [])],
    check=True,timeout=120,env=dict(os.environ,OMP_NUM_THREADS=str(threads),OPENBLAS_NUM_THREADS='1'))

 if not a.action:
  subprocess.run(['mpiexec','-n','8',str(exe),'8 1 1','24 10 6'],check=True,timeout=120,
   env=dict(os.environ,OMP_NUM_THREADS='2',OPENBLAS_NUM_THREADS='1'))
  subprocess.run(['mpiexec','-n','6',str(exe),'3 2 1','15 8 6'],check=True,timeout=120,
   env=dict(os.environ,OMP_NUM_THREADS='2',OPENBLAS_NUM_THREADS='1'))

  # One-field cache probe has fewer axis lines than peers (no full slot).
  subprocess.run(['mpiexec','-n','8',str(exe),'8 1 1','16 2 2'],check=True,timeout=120,
   env=dict(os.environ,OMP_NUM_THREADS='4',OPENBLAS_NUM_THREADS='1'))
