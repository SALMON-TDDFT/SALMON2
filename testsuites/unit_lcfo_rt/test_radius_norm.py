from pathlib import Path
import os,subprocess,tempfile
here=Path(__file__).resolve().parent;repo=here.parents[1]
with tempfile.TemporaryDirectory() as d:
 exe=Path(d)/'probe'
 subprocess.run(['gfortran','-O2','-fcheck=all',str(repo/'src/xc/lcfo_wf_support.f90'),str(here/'radius_norm_probe.f90'),'-L/opt/homebrew/opt/openblas/lib','-lopenblas','-o',str(exe)],cwd=d,check=True)
 subprocess.run([str(exe)],cwd=d,check=True,env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1'))
