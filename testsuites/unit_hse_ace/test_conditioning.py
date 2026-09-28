import argparse,os,subprocess,tempfile
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--build',type=Path,required=True);a=p.parse_args();b=a.build.resolve();o=b/'src/CMakeFiles/salmon.dir'
with tempfile.TemporaryDirectory() as temp:
 exe=Path(temp)/'conditioning'
 subprocess.run(['mpifort','-fopenmp','-fcheck=all','-ffree-line-length-none','-I'+str(b),str(Path(__file__).with_name('conditioning_driver.f90')),*[str(o/n) for n in ('xc/hse_ace.f90.o','xc/exx_orbitals.f90.o','parallel/communication.f90.o','misc/nvtx_wrapper.f90.o')],'-L/opt/homebrew/opt/openblas/lib','-lopenblas','-o',str(exe)],check=True)
 subprocess.run(['mpiexec','-n','1',str(exe)],check=True,env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1'))
