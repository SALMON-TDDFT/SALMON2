"""Manual MPI regression; pass --no-scalapack to exercise the fallback."""
import argparse
import os
from pathlib import Path
import shlex
import subprocess
import tempfile

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--build', type=Path, required=True)
p.add_argument('--no-scalapack', action='store_true')
p.add_argument('--sizes', nargs='+', type=int, default=[2, 7, 65])
p.add_argument('--ranks', nargs='+', type=int, default=[1, 2, 3, 4, 8])
a = p.parse_args()
b = a.build.resolve()
root = Path(__file__).resolve().parents[2]
objects = b / 'src/CMakeFiles/salmon.dir'
env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')
with tempfile.TemporaryDirectory(prefix='ace-metric-') as directory:
    tmp = Path(directory)
    (tmp / 'config.h').write_text('' if a.no_scalapack else '#define USE_SCALAPACK\n')
    libs = shlex.split(os.environ.get('BLAS_LIBS', '-L/opt/homebrew/opt/openblas/lib -lopenblas'))
    if not a.no_scalapack:
        libs = shlex.split(os.environ.get('SCALAPACK_LIBS', '-L/opt/homebrew/opt/scalapack/lib -lscalapack')) + libs
    cmd = [os.environ.get('MPIFC', 'mpifort'), '-cpp', '-fopenmp', '-fcheck=all', '-g',
           '-I'+str(tmp), '-I'+str(b)]
    sources = ['src/xc/exx_sparse_orbitals.f90', 'src/xc/exx_ace.f90', 'src/xc/exx_distributed_metric.f90', 'src/xc/exx_orbitals.f90',
               'testsuites/unit_exx_distributed_metric/driver.f90']
    subprocess.run(cmd + [str(root / f) for f in sources] +
                   [str(objects / f) for f in ['parallel/communication.f90.o', 'misc/nvtx_wrapper.f90.o']] +
                   libs + ['-o', str(tmp / 'driver')], cwd=tmp, check=True)
    for ranks in a.ranks:
        for size in a.sizes:
            subprocess.run([os.environ.get('MPIEXEC', 'mpiexec'), '-n', str(ranks), str(tmp / 'driver'), str(size)],
                           cwd=tmp, env=env, timeout=120, check=True)
