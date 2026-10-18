"""Compiled fixed-core variational and full-spinor exchange regressions.

Set RUN_MPI_TESTS=1 to include two-rank k distribution (MPI sockets required).
Toolchain paths follow test_native_exchange.py.
"""
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]

class NativeSpinNLCCTest(unittest.TestCase):
    def run_probe(self, name, mpi=False):
        compiler = os.environ.get('MPIFC' if mpi else 'FC', 'mpifort' if mpi else 'gfortran')
        fftw = Path(os.environ.get('FFTW_ROOT', '/opt/homebrew/opt/fftw'))
        blas = Path(os.environ.get('OPENBLAS_ROOT', '/opt/homebrew/opt/openblas'))
        libxc = Path(os.environ.get('LIBXC_ROOT', '/opt/homebrew/opt/libxc'))
        if not shutil.which(compiler):
            self.skipTest(f'{compiler} not installed')
        spinor = 'spinor' in name
        sources = ['exx_k_backend', 'exx_k_exchange', 'exx_spinor'] if spinor else ['hybrid_semilocal']
        if spinor and not mpi:
            sources.append('exx_ace')
        with tempfile.TemporaryDirectory() as directory:
            exe = Path(directory) / 'probe'
            cmd = [compiler, '-O1', '-g', '-fcheck=all', '-fopenmp']
            if spinor:
                cmd += ['-I' + str(fftw / 'include')]
            cmd += [str(ROOT / 'src/xc' / (s + '.f90')) for s in sources]
            cmd += [str(ROOT / 'samples/hse_mlwf_reference' / (name + '.f90'))]
            if spinor:
                cmd += ['-L' + str(fftw / 'lib'), '-lfftw3', '-L' + str(blas / 'lib'), '-lopenblas']
            else:
                cmd += ['-L' + str(libxc / 'lib'), '-lxc']
            cmd += ['-o', str(exe)]
            result = subprocess.run(cmd, cwd=directory, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            env = dict(os.environ, OMP_NUM_THREADS='2', OPENBLAS_NUM_THREADS='1')
            cmd = [str(exe)]
            if mpi:
                cmd = [os.environ.get('MPIEXEC', 'mpiexec'), '-n', '2', str(exe)]
            for suffix in ([], ['nan']) if mpi else ([],):
                result = subprocess.run(cmd + suffix, cwd=directory, env=env, capture_output=True,
                                        text=True, timeout=120)
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                self.assertIn('PASS', result.stdout)

    def test_fixed_core_libxc_and_variational_derivatives(self):
        self.run_probe('native_nlcc_probe')

    def test_spin_core_variational_derivatives(self):
        self.run_probe('native_spin_nlcc_probe')

    def test_spinor_offdiagonal_exchange_rotation_and_ace(self):
        self.run_probe('native_spinor_probe')

    @unittest.skipUnless(os.environ.get('RUN_MPI_TESTS') == '1', 'Set RUN_MPI_TESTS=1 for MPI validation')
    def test_spinor_k_distribution(self):
        self.run_probe('native_spinor_mpi_probe', mpi=True)

if __name__ == '__main__':
    unittest.main()
