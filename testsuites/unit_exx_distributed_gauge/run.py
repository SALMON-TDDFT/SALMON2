"""Compare replicated and tiled Gamma MLWF localization and transport."""
import argparse
import os
from pathlib import Path
import shlex
import subprocess
import tempfile

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--build', type=Path, required=True)
p.add_argument('--no-scalapack', action='store_true')
p.add_argument('--empty-only', action='store_true')
a = p.parse_args()
b = a.build.resolve()
root = Path(__file__).resolve().parents[2]
objects = b / 'src/CMakeFiles/salmon.dir'
env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')
with tempfile.TemporaryDirectory(prefix='gauge-tiles-') as directory:
    tmp = Path(directory)
    (tmp / 'config.h').write_text('' if a.no_scalapack else '#define USE_SCALAPACK\n')
    libs = shlex.split(os.environ.get('BLAS_LIBS', '-L/opt/homebrew/opt/openblas/lib -lopenblas'))
    libs += shlex.split(os.environ.get('FFTW_LIBS', '-L/opt/homebrew/opt/fftw/lib -lfftw3'))
    if not a.no_scalapack:
        libs += shlex.split(os.environ.get('SCALAPACK_LIBS', '-L/opt/homebrew/opt/scalapack/lib -lscalapack'))
    sources = ['src/xc/exx_sparse_orbitals.f90', 'src/xc/exx_ace.f90', 'src/xc/exx_distributed_metric.f90', 'src/xc/exx_orbitals.f90',
               'src/xc/exx_distributed_gauge.f90', 'src/xc/exx_spatial.f90',
               'testsuites/unit_exx_distributed_gauge/driver.f90']
    deps = ['xc/exx_pair_candidates.f90.o', 'xc/exx_batch_backend.f90.o', 'xc/exx_spatial_local.f90.o',
            'xc/exx_wannier_gauge.f90.o', 'xc/exx_local_fft.f90.o', 'xc/fftw_pencils.f90.o','xc/fftw_blocks.f90.o',
            'parallel/communication.f90.o', 'misc/nvtx_wrapper.f90.o']
    subprocess.run([os.environ.get('MPIFC', 'mpifort'), '-cpp', '-fopenmp', '-fcheck=all', '-g',
                    '-ffree-line-length-none', '-I'+str(tmp), '-I'+str(b)] +
                   [str(root / f) for f in sources] + [str(objects / f) for f in deps] +
                   libs + ['-o', str(tmp / 'driver')], cwd=tmp, check=True)
    cases = [(1, 1, 32), (2, 1, 32), (2, 2, 32), (4, 1, 32),
                                  (4, 2, 32), (4, 4, 32), (8, 1, 2), (8, 8, 2)]
    if a.empty_only:
        cases = [(8, 8, 2)]
    for ranks, orbital, size in cases:
        subprocess.run([os.environ.get('MPIEXEC', 'mpiexec'), '-n', str(ranks), str(tmp / 'driver'), str(orbital), str(size)],
                       cwd=tmp, env=env, timeout=120, check=True)
