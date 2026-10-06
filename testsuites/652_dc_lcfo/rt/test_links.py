from pathlib import Path
import os,subprocess,tempfile
here=Path(__file__).resolve().parent;root=here.parents[2]
with tempfile.TemporaryDirectory() as folder:
 exe=Path(folder)/'probe'
 subprocess.run(['mpifort','-O2','-fexternal-blas','-fno-tree-loop-vectorize','-fcheck=all',str(here/'local_action_stubs.f90'),str(root/'src/xc/lcfo_mlwf_links.f90'),str(here/'link_probe.f90'),'-L/opt/homebrew/opt/openblas/lib','-lopenblas','-o',str(exe)],cwd=folder,check=True)
 for n in (2,4):
  subprocess.run(['mpirun','-np',str(n),str(exe)],cwd=folder,check=True,env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'))
