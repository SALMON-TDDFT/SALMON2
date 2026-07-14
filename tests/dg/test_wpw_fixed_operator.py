from pathlib import Path
import shutil, subprocess, tempfile, unittest

ROOT = Path(__file__).resolve().parents[2]

class FixedOperatorTest(unittest.TestCase):
    def test_driver(self):
        fc = shutil.which("gfortran")
        self.assertIsNotNone(fc)
        source = ROOT / "src/common/dg_wpw_fixed_operator.f90"
        self.assertTrue(source.exists(), "missing shared fixed-operator context")
        with tempfile.TemporaryDirectory(prefix="wpw-fixed-") as td:
            td = Path(td)
            exe = td / "check"
            result = subprocess.run([
                fc, str(source), str(ROOT / "tests/dg/fixtures/wpw_fixed_operator_driver.F90"),
                "-o", str(exe)
            ], cwd=td, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            self.assertEqual(result.returncode, 0, result.stdout)
            result = subprocess.run([str(exe)], text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            self.assertEqual(result.returncode, 0, result.stdout)
            self.assertIn("PASS WPW fixed operator", result.stdout)

if __name__ == "__main__":
    unittest.main()
