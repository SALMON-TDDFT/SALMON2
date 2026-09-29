"""Check threaded FFT chunks and the same module compiled without OpenMP."""
import argparse,os,shlex,subprocess,tempfile
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--build',type=Path,required=True)
p.add_argument('--no-openmp',action='store_true')
p.add_argument('--action',action='store_true')
p.add_argument('--benchmark',action='store_true')
p.add_argument('--fftw-prefix',type=Path,default=Path('/opt/homebrew'))
p.add_argument('--mpifort',default='mpifort');p.add_argument('--mpiexec',default='mpiexec')
a=p.parse_args()
if a.benchmark and (a.action or a.no_openmp):p.error('--benchmark uses the configured FFT object without other modes')
build=a.build.resolve();root=Path(__file__).resolve().parents[2]
obj=build/'src/CMakeFiles/salmon.dir'
with tempfile.TemporaryDirectory(prefix='fftw-omp-') as t:
 t=Path(t)
 cmd=[a.mpifort,'-O2','-fcheck=all','-I'+str(build),'-I'+str(a.fftw_prefix/'include'),'-J'+str(t)]
 if not a.no_openmp:cmd+=['-fopenmp']
 subprocess.run(cmd+['-c',str(root/'src/xc/fftw_pencils.f90'),'-o',str(t/'fft.o')],check=True,cwd=t)
 extra=[];libraries=[]
 if a.action:
  names=['exx_ace','exx_blas_threads','exx_distributed_gauge','exx_sparse_orbitals','exx_distributed_metric','exx_orbitals',
         'exx_pair_candidates','exx_spatial','exx_batch_backend','exx_spatial_local','exx_wannier_gauge','exx_local_fft']
  extra=[str(obj/('xc/'+n+'.f90.o')) for n in names]
  libraries=['-L/opt/homebrew/opt/openblas/lib','-lopenblas']
  if 'USE_SCALAPACK:BOOL=ON' in (build/'CMakeCache.txt').read_text():
   libraries+=shlex.split(os.environ.get('SCALAPACK_LIBS','-L/opt/homebrew/opt/scalapack/lib -lscalapack'))
 subprocess.run([a.mpifort,'-fopenmp','-fcheck=all','-I'+str(t),'-I'+str(build),
  str(Path(__file__).with_name('benchmark.f90' if a.benchmark else 'action_probe.f90' if a.action else 'probe.f90')),str(obj/'xc/fftw_pencils.f90.o','xc/fftw_blocks.f90.o' if a.benchmark else t/'fft.o'),*extra,
  str(obj/'parallel/communication.f90.o'),str(obj/'misc/nvtx_wrapper.f90.o'),
  '-L'+str(a.fftw_prefix/'lib'),'-lfftw3',*libraries,'-o',str(t/'probe')],check=True)
 for ranks in ((4,) if a.action or a.benchmark else (1,2,4)):
  for scale in ((1,) if a.action or a.benchmark else (1,5)):
   subprocess.run([a.mpiexec,'-n',str(ranks),str(t/'probe'),'serial' if a.no_openmp else 'omp',str(scale)],
   env=dict(os.environ,OMP_NUM_THREADS='4',OPENBLAS_NUM_THREADS='1'),check=True,timeout=120)
