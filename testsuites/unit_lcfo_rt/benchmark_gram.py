"""Isolated Gram timing; one MPI job at a time, no timing pass/fail threshold."""
from pathlib import Path
import os,subprocess,tempfile
here=Path(__file__).resolve().parent;root=here.parents[1]
with tempfile.TemporaryDirectory() as folder:
 exe=Path(folder)/'probe'
 subprocess.run(['mpifort','-O3','-fexternal-blas','-fno-tree-loop-vectorize',str(here/'local_action_stubs.f90'),str(root/'src/xc/lcfo_gram.f90'),str(here/'gram_benchmark.f90'),'-L/opt/homebrew/opt/openblas/lib','-lopenblas','-o',str(exe)],cwd=folder,check=True)
 for n in (4,8,16):
  subprocess.run(['mpirun','--bind-to','none','--map-by','slot','--nooversubscribe','-np',str(n),str(exe)],cwd=folder,check=True,env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'))
