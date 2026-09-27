import os
from pathlib import Path
import subprocess
import tempfile
import unittest
ROOT=Path(__file__).resolve().parents[2]
class SpatialACE(unittest.TestCase):
    def test_mpi_partitions(self):
        with tempfile.TemporaryDirectory() as tmp:
            exe=Path(tmp)/'probe'
            subprocess.run([os.environ.get('MPIFC','mpifort'),'-O0','-g','-fcheck=all','-ffree-line-length-none',str(ROOT/'src/xc/hse_ace.f90'),str(Path(__file__).with_name('spatial_driver.f90')),'-L'+os.environ.get('OPENBLAS_ROOT','/opt/homebrew/opt/openblas')+'/lib','-lopenblas','-o',str(exe)],cwd=tmp,check=True)
            for ranks in (1,2,4):
                with self.subTest(ranks=ranks):
                    run=subprocess.run([os.environ.get('MPIEXEC','mpiexec'),'-n',str(ranks),str(exe)],capture_output=True,text=True,timeout=30,env=dict(os.environ,OPENBLAS_NUM_THREADS='1'))
                    self.assertEqual(run.returncode,0,run.stdout+run.stderr)
                    self.assertIn('PASS spatial ACE',run.stdout)
                    print(run.stdout.strip())
if __name__=='__main__':unittest.main()
