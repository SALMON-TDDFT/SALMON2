"""Distributed screened/global exchange versus the retained full-support oracle."""
import os
from pathlib import Path
import subprocess
import tempfile

root = Path(__file__).resolve().parents[3]
fftw = Path(os.environ.get('FFTW_ROOT', '/opt/homebrew/opt/fftw'))
blas = Path(os.environ.get('OPENBLAS_ROOT', '/opt/homebrew/opt/openblas'))
modules = ['exx_k_backend', 'exx_k_exchange', 'exx_local_fft', 'lcfo_wf_support',
           'exx_wannier_gauge', 'exx_wannier']
with tempfile.TemporaryDirectory() as folder:
    exe = Path(folder) / 'probe'
    subprocess.run([os.environ.get('MPIFC', 'mpifort'), '-O2', '-fcheck=all', '-fopenmp',
                    '-I'+str(fftw/'include'),
                    *[str(root/'src/xc'/f'{m}.f90') for m in modules],
                    str(Path(__file__).with_name('probe.f90')), '-L'+str(fftw/'lib'), '-lfftw3',
                    '-L'+str(blas/'lib'), '-lopenblas', '-o', str(exe)], cwd=folder, check=True)
    for ranks in (1, 2, 3, 4):
        subprocess.run([os.environ.get('MPIEXEC', 'mpiexec'), '-n', str(ranks), str(exe)],
                       env=dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1'),
                       check=True, timeout=120)
