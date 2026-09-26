from pathlib import Path
import subprocess,tempfile,os
r=Path(__file__).resolve().parents[2]
with tempfile.TemporaryDirectory() as d:
 exe=Path(d)/'probe'
 subprocess.run(['gfortran-15','-O2','-fcheck=all','-fopenmp','-fno-tree-loop-vectorize','-I/opt/homebrew/opt/fftw/include',str(r/'src/xc/hse_wannier_gauge.f90'),str(r/'src/xc/hse_wannier.f90'),str(Path(__file__).with_name('exact_pair_probe.f90')),'-L/opt/homebrew/opt/fftw/lib','-lfftw3','-L/opt/homebrew/opt/openblas/lib','-lopenblas','-o',str(exe)],cwd=d,check=True)
 for n in (1,2,4):subprocess.run([str(exe)],cwd=d,env=dict(os.environ,OMP_NUM_THREADS=str(n),OPENBLAS_NUM_THREADS='1'),check=True)
