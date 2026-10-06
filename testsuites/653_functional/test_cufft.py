"""CPU stub checks; opt-in real cuFFT/FFTW parity with --gpu or SALMON_TEST_CUFFT_GPU=1.

An explicitly requested GPU test fails if its compiler, FFTW, or CUDA device
is unavailable. The default reports the GPU test as skipped, never as passed.
FFTW: set FFTW_ROOT, or both FFTW_FFLAGS and FFTW_LIBS, or install pkg-config fftw3.
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


class CufftBackend(unittest.TestCase):
    def compile_and_run(self, compiler, flags, sources, libraries, success):
        with tempfile.TemporaryDirectory(prefix='salmon-cufft-test-') as tmp:
            exe = Path(tmp) / 'probe'
            (Path(tmp) / 'config.h').write_text('')
            command = [compiler, '-I' + tmp, *flags, *map(str, sources), *libraries, '-o', str(exe)]
            build = subprocess.run(command, cwd=tmp, text=True, capture_output=True, timeout=120)
            self.assertEqual(build.returncode, 0, ' '.join(command) + '\n' + build.stdout + build.stderr)
            run = subprocess.run([str(exe)], cwd=tmp, text=True, capture_output=True, timeout=120,
                                 env=dict(os.environ, OMP_NUM_THREADS='1'))
            self.assertEqual(run.returncode, 0, run.stdout + run.stderr)
            self.assertIn(success, run.stdout)
            print(run.stdout.strip())

    def test_cpu_disabled_stub(self):
        compiler = os.environ.get('GNU_FC', 'gfortran')
        if not shutil.which(compiler):
            self.skipTest('GNU Fortran unavailable; set GNU_FC to test the CPU stub')
        self.compile_and_run(compiler, ['-cpp', '-O0', '-g', '-fcheck=all'],
                             [ROOT / 'src/xc/exx_batch_backend.f90', ROOT / 'src/xc/exx_cufft.f90', HERE / 'cufft_stub_probe.f90'], [],
                             'PASS CPU cuFFT stub')

    @unittest.skipUnless(os.environ.get('SALMON_TEST_MPIEXEC'), 'set SALMON_TEST_MPIEXEC for abort test')
    def test_fatal_callback_aborts_mpi_job(self):
        compiler = os.environ.get('MPIFC', 'mpifort')
        self.assertTrue(shutil.which(compiler), 'MPI Fortran compiler missing')
        with tempfile.TemporaryDirectory(prefix='salmon-cufft-abort-') as tmp:
            folder = Path(tmp)
            (folder / 'config.h').write_text('#define USE_MPI\n')
            driver = folder / 'abort.f90'
            driver.write_text("""program abort_probe
  use mpi
  use iso_c_binding, only: c_ptr,c_null_ptr
  implicit none
  interface
    subroutine fatal(message) bind(C,name='salmon_exx_cufft_abort')
      import c_ptr
      implicit none
      type(c_ptr),value :: message
    end subroutine
  end interface
  integer :: rank,ierr
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  if(rank==0)call fatal(c_null_ptr)
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  print *, 'UNEXPECTED SUCCESS'
  call MPI_Finalize(ierr)
end program
""")
            exe = folder / 'abort'
            build = subprocess.run([compiler, '-cpp', '-I' + tmp, str(ROOT / 'src/xc/exx_batch_backend.f90'),
                                    str(ROOT / 'src/xc/exx_cufft.f90'),
                                    str(driver), '-o', str(exe)], cwd=tmp, text=True, capture_output=True)
            self.assertEqual(build.returncode, 0, build.stdout + build.stderr)
            run = subprocess.run([os.environ['SALMON_TEST_MPIEXEC'], '-n', '2', str(exe)],
                                 cwd=tmp, text=True, capture_output=True, timeout=30)
            self.assertNotEqual(run.returncode, 0)
            self.assertIn('EXX cuFFT: fatal OpenACC runtime error', run.stdout + run.stderr)
            self.assertNotIn('UNEXPECTED SUCCESS', run.stdout + run.stderr)

    @unittest.skipUnless(GPU_REQUESTED, 'GPU parity not requested; use --gpu or SALMON_TEST_CUFFT_GPU=1')
    def test_gpu_matches_fftw(self):
        compiler = os.environ.get('NVFC', 'nvfortran')
        self.assertTrue(shutil.which(compiler), 'GPU test requested but NVHPC compiler missing: ' + compiler)
        include, libraries = fftw_flags()
        flags = ['-cpp', '-DUSE_EXX_CUFFT', '-acc', '-cudalib=cufft', '-O2', '-Mbounds',
                 *shlex.split(os.environ.get('CUFFT_TEST_FFLAGS', '')), *include]
        self.compile_and_run(compiler, flags,
                             [ROOT / 'src/xc/exx_batch_backend.f90', ROOT / 'src/xc/exx_cufft.f90', ROOT / 'src/xc/exx_local_fft.f90',
                              HERE / 'cufft_gpu_probe.f90'], libraries, 'PASS cuFFT/FFTW complex128 parity')


if __name__ == '__main__':
    unittest.main()
