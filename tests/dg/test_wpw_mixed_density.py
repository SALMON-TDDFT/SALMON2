from pathlib import Path
import shutil,subprocess,tempfile,unittest
ROOT=Path(__file__).resolve().parents[2]
class MixedDensityTest(unittest.TestCase):
 def test_driver(self):
  fc=shutil.which('gfortran');self.assertIsNotNone(fc)
  with tempfile.TemporaryDirectory(prefix='wpw-density-') as td:
   td=Path(td);(td/'CMakeLists.txt').write_text(f'''cmake_minimum_required(VERSION 3.16)
project(wpw_density LANGUAGES Fortran)
add_executable(check "{ROOT/'src/common/dg_wpw_density.f90'}" "{ROOT/'tests/dg/fixtures/wpw_density_driver.F90'}")
''')
   for cmd in (["cmake","-S",str(td),"-B",str(td/'b'),f"-DCMAKE_Fortran_COMPILER={fc}"],["cmake","--build",str(td/'b'),"-j","2"]):
    r=subprocess.run(cmd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True);self.assertEqual(r.returncode,0,r.stdout)
   r=subprocess.run([str(td/'b/check')],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True);self.assertEqual(r.returncode,0,r.stdout);self.assertIn('PASS mixed density',r.stdout)
if __name__=='__main__':unittest.main()
