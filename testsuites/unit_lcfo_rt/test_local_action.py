from pathlib import Path
import os,subprocess,tempfile
root=Path(__file__).resolve().parents[2];here=Path(__file__).parent
with tempfile.TemporaryDirectory() as folder:
 exe=Path(folder)/'probe'
 subprocess.run(['mpifort','-O2','-fexternal-blas','-fno-tree-loop-vectorize','-fcheck=all',
   str(here/'local_action_stubs.f90'),str(root/'src/xc/hse_ace.f90'),
   str(root/'src/xc/lcfo_ace_local.f90'),str(here/'local_action_probe.f90'),
   '-L/opt/homebrew/opt/openblas/lib','-lopenblas','-o',str(exe)],cwd=folder,check=True)
 subprocess.run(['mpirun','-np','2',str(exe)],cwd=folder,check=True,
   env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1'))
