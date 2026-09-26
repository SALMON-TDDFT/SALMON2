from pathlib import Path
import os, subprocess, tempfile
root=Path(__file__).resolve().parents[2]
with tempfile.TemporaryDirectory() as d:
    cmd=[os.environ.get('FC','gfortran'),'-O2','-fexternal-blas','-fno-tree-loop-vectorize','-fcheck=all',
         str(root/'src/xc/hse_wannier_gauge.f90'),str(Path(__file__).with_name('gamma_probe.f90')),
         '-L/opt/homebrew/opt/openblas/lib','-lopenblas','-o',str(Path(d)/'probe')]
    subprocess.run(cmd,cwd=d,check=True)
    subprocess.run([str(Path(d)/'probe')],cwd=d,check=True,
                   env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1'))
