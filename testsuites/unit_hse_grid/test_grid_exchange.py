"""Sequential MPI partition/parity, high-mode target and Hermiticity checks."""
import argparse
import os
from pathlib import Path
import subprocess
import tempfile

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--ranks', type=int, nargs='+', default=[1, 2, 4])
parser.add_argument('--compile-only', action='store_true')
args = parser.parse_args()
root = Path(__file__).resolve().parents[2]
with tempfile.TemporaryDirectory(prefix='hse-grid-') as directory:
    work = Path(directory)
    (work / 'config.h').write_text('#define USE_MPI\n')
    fftw = Path(os.environ.get('FFTW_ROOT', '/opt/homebrew/opt/fftw'))
    blas = Path(os.environ.get('OPENBLAS_ROOT', '/opt/homebrew/opt/openblas'))
    executable = work / 'probe'
    command = [os.environ.get('MPIFC', 'mpifort'), '-cpp', '-fopenmp', '-O2',
               '-fcheck=all', '-fno-tree-loop-vectorize', '-I'+str(work), '-I'+str(fftw/'include')]
    command += [str(root/'src/xc'/name) for name in
                ['hse_wannier_gauge.f90', 'hse_wannier.f90', 'hse_grid_exchange.f90']]
    command += [str(Path(__file__).with_name('probe.f90')), '-L'+str(fftw/'lib'),
                '-lfftw3', '-L'+str(blas/'lib'), '-lopenblas', '-o', str(executable)]
    subprocess.run(command, cwd=work, check=True)
    for ranks in ([] if args.compile_only else args.ranks):
        subprocess.run(['mpiexec', '-n', str(ranks), str(executable)], cwd=work,
                       env=dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1'),
                       check=True, timeout=60)
