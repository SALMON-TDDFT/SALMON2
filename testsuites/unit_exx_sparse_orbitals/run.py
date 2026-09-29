"""Exact sparse support storage; no MPI or external libraries required."""
import os
from pathlib import Path
import subprocess
import tempfile
root = Path(__file__).resolve().parents[2]
with tempfile.TemporaryDirectory(prefix='exx-sparse-') as tmp:
    exe = Path(tmp) / 'probe'
    subprocess.run([os.environ.get('FC', 'gfortran'), '-O2', '-fcheck=all',
                    str(root / 'src/xc/exx_sparse_orbitals.f90'),
                    str(Path(__file__).with_name('probe.f90')), '-o', str(exe)], cwd=tmp, check=True)
    subprocess.run([str(exe)], check=True)
