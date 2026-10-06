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
p.add_argument('--threads', nargs='+', type=int, default=[1, 2])
p.add_argument('--ranks', nargs='+', type=int, default=[1, 2, 3, 4, 8])
a = p.parse_args()
b = a.build.resolve()
root = Path(__file__).resolve().parents[3]
objects = b / 'src/CMakeFiles/salmon.dir'
env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')
with tempfile.TemporaryDirectory(prefix='ace-metric-') as directory:
    tmp = Path(directory)
    (tmp / 'config.h').write_text(('' if a.no_scalapack else '#define USE_SCALAPACK\n') +
        ('#define HAVE_EXX_OPENBLAS_THREADS\n' if '#define HAVE_EXX_OPENBLAS_THREADS' in (b/'config.h').read_text() else ''))
    backend=(root/'src/xc/exx_blas_threads.f90').read_text().replace('module exx_blas_threads','module gauge_test_backend')
    (tmp/'backend.f90').write_text(backend)
    libs = shlex.split(os.environ.get('BLAS_LIBS', '-L/opt/homebrew/opt/openblas/lib -lopenblas'))
    if not a.no_scalapack:
        libs = shlex.split(os.environ.get('SCALAPACK_LIBS', '-L/opt/homebrew/opt/scalapack/lib -lscalapack')) + libs
    cmd = [os.environ.get('MPIFC', 'mpifort'), '-cpp', '-fopenmp', '-fcheck=all', '-g',
           '-I'+str(tmp), '-I'+str(b)]
    sources = ['testsuites/651_hybrid_exchange/gauge/thread_probe.f90', 'src/xc/exx_sparse_orbitals.f90', 'src/xc/exx_ace.f90', 'src/xc/exx_distributed_metric.f90', 'src/xc/exx_orbitals.f90',
               'testsuites/651_hybrid_exchange/metric/driver.f90']
    subprocess.run(cmd + [str(tmp/'backend.f90')] + [str(root / f) for f in sources] +
                   [str(objects / f) for f in ['parallel/communication.f90.o', 'misc/nvtx_wrapper.f90.o']] +
                   libs + ['-o', str(tmp / 'driver')], cwd=tmp, check=True)
    for ranks in a.ranks:
        for size in a.sizes:
          for threads in a.threads:
            env['OMP_NUM_THREADS']=str(threads)
            subprocess.run([os.environ.get('MPIEXEC', 'mpiexec'), '-n', str(ranks), str(tmp / 'driver'), str(size)],
                           cwd=tmp, env=env, timeout=120, check=True)
