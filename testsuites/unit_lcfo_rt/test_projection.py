from pathlib import Path
import subprocess,tempfile,os
r=Path(__file__).resolve().parents[2]
with tempfile.TemporaryDirectory() as d:
 exe=Path(d)/'probe'
 subprocess.run(['gfortran-15','-O2','-fcheck=all','-fexternal-blas','-fno-tree-loop-vectorize',str(r/'src/xc/lcfo_projection.f90'),str(Path(__file__).with_name('projection_probe.f90')),'-L/opt/homebrew/opt/openblas/lib','-lopenblas','-o',str(exe)],cwd=d,check=True)
 subprocess.run([str(exe)],cwd=d,env=dict(os.environ,OPENBLAS_NUM_THREADS='1'),check=True)
