"""Compile and exercise the standalone projected RT core (one BLAS thread)."""
from pathlib import Path
import os, subprocess, tempfile
root=Path(__file__).resolve().parents[2]
source=root/'src/rt/lcfo_rt_core.f90'
assert source.exists(), 'missing LCFO real-time propagation kernel'
with tempfile.TemporaryDirectory(prefix='lcfo-rt-') as work:
    exe=Path(work)/'probe'
    subprocess.run(['gfortran','-O2','-fcheck=all','-fbacktrace','-fno-tree-loop-vectorize',
                    str(source),str(Path(__file__).with_name('probe.f90')),
                    '-L/opt/homebrew/opt/openblas/lib','-lopenblas','-o',str(exe)],cwd=work,check=True)
    subprocess.run([str(exe)],env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1'),check=True)
