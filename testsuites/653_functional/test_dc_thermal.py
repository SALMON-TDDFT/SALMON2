"""Fixed-charge thermal occupation response, independent of the SCF implementation."""
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]

@unittest.skipUnless(shutil.which('gfortran'), 'gfortran required')
class DCThermalTest(unittest.TestCase):
    def test_weighted_thermal_response(self):
        source = ROOT/'src/gs/dc/dc_thermal.f90'
        self.assertTrue(source.exists(), 'weighted thermal response module missing')
        with tempfile.TemporaryDirectory() as tmp:
            exe = Path(tmp)/'probe'
            build = subprocess.run(['gfortran', '-O2', '-fcheck=all', '-ffpe-trap=invalid,zero,overflow',
                str(source), str(Path(__file__).with_name('dc_thermal_probe.f90')), '-o', str(exe)],
                cwd=tmp, text=True, capture_output=True)
            self.assertEqual(build.returncode, 0, build.stdout+build.stderr)
            run = subprocess.run([str(exe)], cwd=tmp, text=True, capture_output=True)
            self.assertEqual(run.returncode, 0, run.stdout+run.stderr)

if __name__ == '__main__': unittest.main()
