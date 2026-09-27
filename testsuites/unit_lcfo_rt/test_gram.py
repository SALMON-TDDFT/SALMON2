from pathlib import Path
import os,subprocess,tempfile
here=Path(__file__).parent;root=here.resolve().parents[1]
with tempfile.TemporaryDirectory() as folder:
 exe=Path(folder)/'probe'
 subprocess.run(['mpifort','-O2','-fexternal-blas','-fno-tree-loop-vectorize','-fcheck=all',str(here.resolve()/'local_action_stubs.f90'),str(root/'src/xc/lcfo_gram.f90'),str(here.resolve()/'gram_probe.f90'),'-L/opt/homebrew/opt/openblas/lib','-lopenblas','-o',str(exe)],cwd=folder,check=True)
 for n in (2,4):
  subprocess.run(['mpirun','-np',str(n),str(exe)],cwd=folder,check=True,env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'))
