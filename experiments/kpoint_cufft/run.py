"""MPI convolution-boundary parity; --gpu requires NVHPC and a real GPU per rank."""
import argparse
import os
from pathlib import Path
import shlex
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]
p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--gpu', action='store_true')
p.add_argument('--ranks', type=int, nargs='+', default=[1, 2, 3, 4])
p.add_argument('--threads', type=int, default=1)
p.add_argument('--mpi-args', default='', help='site-specific launcher mapping options')
a = p.parse_args()
if any(n < 1 for n in a.ranks):
    p.error('positive rank counts required')
fftw = Path(os.environ.get('FFTW_ROOT', '/opt/homebrew/opt/fftw'))
blas = Path(os.environ.get('OPENBLAS_ROOT', '/opt/homebrew/opt/openblas'))
includes = shlex.split(os.environ.get('FFTW_FFLAGS', '-I' + str(fftw / 'include')))
libraries = shlex.split(os.environ.get('FFTW_LIBS', '-L' + str(fftw / 'lib') + ' -lfftw3'))
libraries += shlex.split(os.environ.get('BLAS_LIBS', '-L' + str(blas / 'lib') + ' -lopenblas'))
flags = ['-O2', '-cpp', '-fopenmp', '-fcheck=all']
sources = ['src/xc/exx_k_backend.f90']
if a.gpu:
    flags = ['-O2', '-cpp', '-mp', '-acc', '-cudalib=cufft', '-Mbounds', '-DUSE_EXX_CUFFT']
    flags += shlex.split(os.environ.get('CUFFT_TEST_FFLAGS', '-gpu=cc90'))
    sources += ['src/xc/exx_batch_backend.f90', 'src/xc/exx_cufft.f90', 'src/xc/exx_k_cufft.f90']
sources += ['src/xc/exx_k_exchange.f90', 'experiments/kpoint_cufft/probe.f90']
with tempfile.TemporaryDirectory(prefix='salmon-k-cufft-') as tmp:
    folder = Path(tmp)
    (folder / 'config.h').write_text('#define USE_MPI\n')
    exe = folder / 'probe'
    subprocess.run([os.environ.get('MPIFC', 'mpifort'), *flags, '-I' + tmp, *includes,
                    *[str(ROOT / f) for f in sources], *libraries, '-o', str(exe)],
                   cwd=tmp, check=True, timeout=180)
    for ranks in a.ranks:
        subprocess.run([os.environ.get('MPIEXEC', 'mpiexec'), *shlex.split(a.mpi_args), '-n', str(ranks), str(exe)],
                       check=True, timeout=120,
                       env=dict(os.environ, OMP_NUM_THREADS=str(a.threads), OPENBLAS_NUM_THREADS='1'))
