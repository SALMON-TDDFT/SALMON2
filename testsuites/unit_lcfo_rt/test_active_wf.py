from pathlib import Path
import os,subprocess,tempfile,argparse
p=argparse.ArgumentParser();p.add_argument('--compile-only',action='store_true');args=p.parse_args()
r=Path(__file__).resolve().parents[2];h=Path(__file__).parent
with tempfile.TemporaryDirectory() as d:
 exe=Path(d)/'probe'
 subprocess.run(['gfortran-15','-O2','-fcheck=all','-fexternal-blas','-fno-tree-loop-vectorize',str(r/'src/xc/lcfo_wf_support.f90'),str(h/'active_wf_probe.f90'),'-L/opt/homebrew/opt/openblas/lib','-lopenblas','-o',str(exe)],cwd=d,check=True)
 if not args.compile_only:subprocess.run([str(exe)],cwd=d,check=True,env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'))
