"""Opt-in direct resident k-mesh cuFFT versus FFTW parity; no MPI required.

Use --gpu or SALMON_TEST_CUFFT_GPU=1 to require NVHPC and a working NVIDIA GPU.
Without that request GPU parity is explicitly skipped, never reported passed.
"""
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
GPU_REQUESTED = os.environ.get('SALMON_TEST_CUFFT_GPU', '').lower() in ('1', 'yes', 'true')
if '--gpu' in sys.argv:
    GPU_REQUESTED = True
    sys.argv.remove('--gpu')


def fftw_flags():
    if 'FFTW_FFLAGS' in os.environ or 'FFTW_LIBS' in os.environ:
        if not all(name in os.environ for name in ('FFTW_FFLAGS', 'FFTW_LIBS')):
            raise RuntimeError('Set both FFTW_FFLAGS and FFTW_LIBS')
        return shlex.split(os.environ['FFTW_FFLAGS']), shlex.split(os.environ['FFTW_LIBS'])
    if os.environ.get('FFTW_ROOT'):
        root = Path(os.environ['FFTW_ROOT']).expanduser().resolve()
        library = root / 'lib'
        if not library.is_dir() and (root / 'lib64').is_dir():
            library = root / 'lib64'
        return ['-I' + str(root / 'include')], ['-L' + str(library), '-lfftw3']
    if shutil.which('pkg-config'):
        cflags = subprocess.run(['pkg-config', '--cflags', 'fftw3'], text=True, capture_output=True)
        libs = subprocess.run(['pkg-config', '--libs', 'fftw3'], text=True, capture_output=True)
        if cflags.returncode == libs.returncode == 0:
            return shlex.split(cflags.stdout), shlex.split(libs.stdout)
    raise RuntimeError('FFTW not configured: set FFTW_ROOT or FFTW_FFLAGS and FFTW_LIBS')


class KpointCufftGPU(unittest.TestCase):
    @unittest.skipUnless(GPU_REQUESTED, 'GPU parity not requested; use --gpu or SALMON_TEST_CUFFT_GPU=1')
    def test_resident_tiles_match_fftw(self):
        compiler = os.environ.get('NVFC', 'nvfortran')
        self.assertTrue(shutil.which(compiler), 'GPU test requested but NVHPC compiler missing: ' + compiler)
        includes, libraries = fftw_flags()
        with tempfile.TemporaryDirectory(prefix='salmon-k-cufft-gpu-') as tmp:
            folder = Path(tmp)
            (folder / 'config.h').write_text('')
            exe = folder / 'probe'
            flags = ['-cpp', '-DUSE_EXX_CUFFT', '-acc', '-cudalib=cufft', '-O2', '-Mbounds',
                     *shlex.split(os.environ.get('CUFFT_TEST_FFLAGS', ''))]
            sources = [ROOT / 'src/xc' / name for name in
                       ('exx_batch_backend.f90', 'exx_cufft.f90', 'exx_k_backend.f90', 'exx_k_cufft.f90')]
            command = [compiler, '-I' + tmp, *flags, *includes, *map(str, sources),
                       str(HERE / 'gpu_probe.f90'), *libraries, '-o', str(exe)]
            build = subprocess.run(command, cwd=tmp, text=True, capture_output=True, timeout=180)
            self.assertEqual(build.returncode, 0, ' '.join(command) + '\n' + build.stdout + build.stderr)
            run = subprocess.run([str(exe)], cwd=tmp, text=True, capture_output=True, timeout=120,
                                 env=dict(os.environ, OMP_NUM_THREADS='1'))
            self.assertEqual(run.returncode, 0, run.stdout + run.stderr)
            self.assertIn('PASS k-cuFFT/FFTW parity', run.stdout)
            print(run.stdout.strip())


if __name__ == '__main__':
    unittest.main()
