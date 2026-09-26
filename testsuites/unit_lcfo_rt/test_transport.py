from pathlib import Path
import os,subprocess,tempfile
root=Path(__file__).resolve().parents[2];here=Path(__file__).parent
with tempfile.TemporaryDirectory() as d:
    exe=Path(d)/'probe'
    cmd=[os.environ.get('FC','gfortran'),'-fopenmp','-O2','-fexternal-blas','-fno-tree-loop-vectorize','-fcheck=all',
         str(here/'transport_stubs.f90'),str(root/'src/xc/hse_wannier_gauge.f90'),
         str(root/'src/xc/lcfo_rt_wannier.f90'),str(here/'transport_probe.f90'),
         '-L/opt/homebrew/opt/openblas/lib','-lopenblas','-o',str(exe)]
    subprocess.run(cmd,cwd=d,check=True)
    for threads in ['1','2','4']:
        subprocess.run([str(exe)],cwd=d,check=True,env=dict(os.environ,OPENBLAS_NUM_THREADS='1',
                       OMP_NUM_THREADS=threads,SALMON_LCFO_RT_MLWF='1',SALMON_LCFO_RT_RADIUS='1'))
