"""Validate the production k-mesh cuFFT stub and owning lifecycle with GNU Fortran."""
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent


class KpointCufftStub(unittest.TestCase):
    def test_validation_and_lifecycle(self):
        compiler = os.environ.get('GNU_FC', 'gfortran')
        self.assertTrue(shutil.which(compiler), 'GNU Fortran required; set GNU_FC')
        with tempfile.TemporaryDirectory(prefix='salmon-k-cufft-stub-') as tmp:
            folder = Path(tmp)
            (folder / 'config.h').write_text('')
            exe = folder / 'probe'
            build = subprocess.run([compiler, '-cpp', '-O0', '-g', '-fcheck=all', '-I' + tmp,
                                    str(ROOT / 'src/xc/exx_k_backend.f90'),
                                    str(ROOT / 'src/xc/exx_k_cufft.f90'), str(HERE / 'stub_probe.f90'),
                                    '-o', str(exe)], cwd=tmp, capture_output=True, text=True, timeout=120)
            self.assertEqual(build.returncode, 0, build.stdout + build.stderr)
            run = subprocess.run([str(exe)], cwd=tmp, capture_output=True, text=True, timeout=30)
            self.assertEqual(run.returncode, 0, run.stdout + run.stderr)
            self.assertIn('k-cuFFT CPU stub validation passed', run.stdout)
            print(run.stdout.strip())


if __name__ == '__main__':
    unittest.main()
