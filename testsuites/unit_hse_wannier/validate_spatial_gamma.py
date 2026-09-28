"""Convergence and unitarity of mixed localized subspaces across MPI layouts."""
import argparse,os,subprocess,tempfile
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--build',type=Path,required=True);a=p.parse_args()
b=a.build.resolve();obj=b/'src/CMakeFiles/salmon.dir'
names=['xc/exx_orbitals.f90.o','xc/hse_spatial.f90.o','xc/exx_spatial_local.f90.o','xc/hse_wannier_gauge.f90.o','xc/exx_local_fft.f90.o','xc/fftw_pencils.f90.o','parallel/communication.f90.o','misc/nvtx_wrapper.f90.o']
with tempfile.TemporaryDirectory() as tmp:
 exe=Path(tmp)/'probe'
 subprocess.run(['mpifort','-fopenmp','-fcheck=all','-ffree-line-length-none','-I'+str(b),str(Path(__file__).with_name('spatial_gamma_probe.f90')),*[str(obj/n) for n in names],'-L/opt/homebrew/opt/fftw/lib','-lfftw3','-L/opt/homebrew/opt/openblas/lib','-lopenblas','-o',str(exe)],check=True)
 for ranks,orbitals in [(1,1),(2,1),(4,1),(2,2),(4,2)]:
  subprocess.run(['mpiexec','-n',str(ranks),str(exe),str(orbitals)],check=True,env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'))
