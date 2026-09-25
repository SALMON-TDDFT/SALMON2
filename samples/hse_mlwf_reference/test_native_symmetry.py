"""Small compiled star/phase/projector tests; no production run or MPI required."""
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]

class NativeSymmetryTest(unittest.TestCase):
    def test_star_projectors_and_bloch_phases(self):
        compiler = os.environ.get('FC', shutil.which('gfortran') or shutil.which('gfortran-15'))
        if not compiler:
            self.skipTest('Fortran compiler unavailable')
        with tempfile.TemporaryDirectory() as tmp:
            exe = Path(tmp) / 'probe'
            run = subprocess.run([compiler, '-O0', '-fcheck=all', '-Wall',
                str(ROOT / 'src/xc/hse_symmetry.f90'),
                str(ROOT / 'samples/hse_mlwf_reference/native_symmetry_probe.f90'),
                '-o', str(exe)], cwd=tmp, capture_output=True, text=True)
            self.assertEqual(run.returncode, 0, run.stderr)
            run = subprocess.run([str(exe)], capture_output=True, text=True, timeout=20)
            self.assertEqual(run.returncode, 0, run.stdout + run.stderr)
            self.assertIn('symmetry probe passed', run.stdout)

    def test_density_callback_matches_expanded_exchange(self):
        compiler = os.environ.get('FC', shutil.which('gfortran') or shutil.which('gfortran-15'))
        fftw = Path(os.environ.get('FFTW_ROOT', '/opt/homebrew/opt/fftw'))
        blas = Path(os.environ.get('OPENBLAS_ROOT', '/opt/homebrew/opt/openblas'))
        if not compiler or not (fftw / 'include/fftw3.f03').exists():
            self.skipTest('Fortran compiler and FFTW required')
        with tempfile.TemporaryDirectory() as tmp:
            exe = Path(tmp) / 'probe'
            run = subprocess.run([compiler, '-O0', '-fcheck=all', '-Wall',
                '-I' + str(fftw / 'include'), str(ROOT / 'src/xc/hse_symmetry.f90'),
                str(ROOT / 'src/xc/hse_exchange.f90'),
                str(ROOT / 'samples/hse_mlwf_reference/native_symmetry_exchange_probe.f90'),
                '-L' + str(fftw / 'lib'), '-lfftw3', '-L' + str(blas / 'lib'),
                '-lopenblas', '-o', str(exe)], cwd=tmp, capture_output=True, text=True)
            self.assertEqual(run.returncode, 0, run.stderr)
            run = subprocess.run([str(exe)], capture_output=True, text=True, timeout=30,
                                 env=dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1'))
            self.assertEqual(run.returncode, 0, run.stdout + run.stderr)
            self.assertIn('symmetry exchange callback passed', run.stdout)

if __name__ == '__main__':
    unittest.main()
