"""Regress Gamma Jacobi saddle escape, projected seeds and deferred reseeding.

Run with --build PATH or SALMON_BUILD=PATH after building SALMON with HSE/MPI.
The probe compiles the current spatial module against the build's dependencies.
"""
import shlex
import argparse
import os
from pathlib import Path
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build', type=Path, default=os.environ.get('SALMON_BUILD'))
    args = parser.parse_args()
    if args.build is None:
        parser.error('provide --build PATH or set SALMON_BUILD')
    build = args.build.resolve()
    objects = build / 'src/CMakeFiles/salmon.dir'
    dependencies = [
        'xc/exx_distributed_gauge.f90.o','xc/exx_sparse_orbitals.f90.o','xc/exx_distributed_metric.f90.o','xc/exx_orbitals.f90.o', 'xc/exx_ace.f90.o',
        'xc/exx_batch_backend.f90.o','xc/exx_spatial_local.f90.o', 'xc/exx_wannier_gauge.f90.o',
        'xc/exx_local_fft.f90.o', 'xc/fftw_pencils.f90.o','xc/fftw_blocks.f90.o',
        'parallel/communication.f90.o', 'misc/nvtx_wrapper.f90.o',
    ]
    fftw = Path(os.environ.get('FFTW_ROOT', '/opt/homebrew/opt/fftw'))
    blas = Path(os.environ.get('OPENBLAS_ROOT', '/opt/homebrew/opt/openblas'))
    with tempfile.TemporaryDirectory() as tmp:
        executable = Path(tmp) / 'probe'
        subprocess.run([
            os.environ.get('MPIFC', 'mpifort'), '-fopenmp', '-fcheck=all',
            '-ffree-line-length-none', '-I' + str(build),
            str(ROOT / 'src/xc/exx_pair_candidates.f90'),
            str(ROOT / 'src/xc/exx_spatial.f90'),
            str(Path(__file__).with_name('projected_seed_probe.f90')),
            *[str(objects / name) for name in dependencies],
            '-L' + str(fftw / 'lib'), '-lfftw3',
            '-L' + str(blas / 'lib'), *(shlex.split(os.environ.get('SCALAPACK_LIBS', '-L/opt/homebrew/opt/scalapack/lib -lscalapack')) if 'USE_SCALAPACK:BOOL=ON' in (build / 'CMakeCache.txt').read_text() else []),'-lopenblas', '-o', str(executable),
        ], cwd=tmp, check=True)
        for ranks in (1, 2, 4):
            subprocess.run([
                os.environ.get('MPIEXEC', 'mpiexec'), '-n', str(ranks), str(executable),
            ], cwd=tmp, check=True, timeout=30,
                env=dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1'))


if __name__ == '__main__':
    main()
