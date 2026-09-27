"""Integration rejection tests; set SALMON_TEST_EXE to an HSE-enabled executable."""
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest
ROOT=Path(__file__).resolve().parents[2]

@unittest.skipUnless(os.environ.get('SALMON_TEST_EXE'),'SALMON_TEST_EXE is required')
class InputTest(unittest.TestCase):
    def run_input(self, inp, expected):
        with tempfile.TemporaryDirectory() as tmp:
            for atom in ('H','O'):shutil.copy(ROOT/f'testsuites/pseudo/{atom}_rps.dat',tmp)
            result=subprocess.run([os.environ['SALMON_TEST_EXE']],input=inp,cwd=tmp,text=True,capture_output=True,
                timeout=90,env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',OMPI_MCA_btl='self,vader'))
            self.assertNotEqual(result.returncode,0,result.stdout[-500:])
            self.assertIn(expected,result.stdout+result.stderr)
    def base(self):
        return (ROOT/'samples/pbeh40_rvv10/water.inp').read_text().replace('nscf = 300','nscf = 1')
    def test_unconverged_bomd_is_rejected(self):
        inp=self.base().replace("theory='dft'","theory='dft_md'")
        inp+="\n&tgrid\n dt=.01\n nt=1\n/\n&md\n ensemble='NVE'\n yn_set_ini_velocity='y'\n temperature0_ion_k=300\n/\n"
        self.run_input(inp,'SCF not converged; ionic step rejected')
    def test_legacy_snapshot_is_rejected(self):
        inp=self.base().replace('rvv10_nq=32',"rvv10_nq=32\n yn_hse_wannier_snapshot='y'")
        self.run_input(inp,'PBEh40: legacy HSE Wannier snapshot cannot encode Coulomb cutoff')
    def test_dc_md_is_rejected(self):
        self.run_input(self.base().replace("theory='dft'","theory='dft_md'\n yn_dc='y'"),'DC MD is not yet supported')
    def test_restart_is_rejected(self):
        self.run_input(self.base().replace("sysname = 'H2O'","sysname = 'H2O'\n yn_restart='y'"),'checkpoint parameter validation')
    def test_invalid_rvv10_parameters_are_rejected(self):
        self.run_input(self.base().replace('rvv10_nq=32','rvv10_nq=3'),'invalid rVV10 parameters')
    def test_hse_md_stays_rejected(self):
        self.run_input(self.base().replace("theory='dft'","theory='dft_md'").replace("xc ='pbeh40_rvv10'","xc ='hse06'"),'unsupported calculation theory')
if __name__=='__main__':unittest.main()
