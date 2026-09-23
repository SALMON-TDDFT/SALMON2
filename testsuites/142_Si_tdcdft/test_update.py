"""Compile and exercise the production xc-field update against analytic solutions."""
from pathlib import Path
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]

class FieldUpdateTest(unittest.TestCase):
    def test_analytic_solutions(self):
        source = ROOT / 'src/rt/tdcdft_lrc.f90'
        self.assertTrue(source.exists(), 'TDCDFT field update is not implemented')
        with tempfile.TemporaryDirectory() as work:
            subprocess.run(['gfortran', '-std=f2008', '-Wall', '-Wextra', '-fcheck=all',
                            str(source), str(Path(__file__).with_name('test_update.f90')),
                            '-o', 'test_update'], cwd=work, check=True)
            subprocess.run([str(Path(work) / 'test_update')], cwd=work, check=True)

if __name__ == '__main__':
    unittest.main()
