"""Compile the distributed solver with bounds checks; compare to LAPACK."""
import argparse,os,subprocess,tempfile
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--build',type=Path,required=True);p.add_argument('--quick',action='store_true');a=p.parse_args()
b=a.build.resolve();root=Path(__file__).resolve().parents[3]
with tempfile.TemporaryDirectory() as tmp:
    exe=Path(tmp)/'probe'
    subprocess.run(['mpifort','-fopenmp','-fcheck=all','-ffree-line-length-none','-J'+tmp,'-I'+str(b),
        str(root/'src/gs/dc/lcfo_scalapack.f90'),str(Path(__file__).with_name('driver.f90')),
        str(b/'src/CMakeFiles/salmon.dir/parallel/communication.f90.o'),
        str(b/'src/CMakeFiles/salmon.dir/misc/nvtx_wrapper.f90.o'),
        '-L/opt/homebrew/opt/scalapack/lib','-lscalapack','-L/opt/homebrew/opt/openblas/lib','-lopenblas','-o',str(exe)],check=True)
    for ranks in ((1,4) if a.quick else (1,2,4,6)):
        subprocess.run(['mpiexec','-n',str(ranks),str(exe)],check=True,timeout=30,
                       env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'))
