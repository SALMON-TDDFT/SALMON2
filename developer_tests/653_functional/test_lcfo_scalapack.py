"""Distributed full LCFO diagonalization and native mesh reconstruction parity."""
import os,re,unittest
import numpy as np
import test_ehrenfest as helpers

@unittest.skipUnless(os.environ.get('SALMON_TEST_SCALAPACK')=='1', 'ScaLAPACK build required')
class LCFOScalapack(unittest.TestCase):
    functional='hse06'
    setUpClass=classmethod(helpers.RealspaceEhrenfest.setUpClass.__func__)
    execute=classmethod(helpers.RealspaceEhrenfest.execute.__func__)
    rt_input=helpers.RealspaceEhrenfest.rt_input

    def test_distributed_diagonalization(self):
        reference=np.loadtxt(next((self.root/'gs/data_dcdft/total').glob('*_eigen.data')))[:,3]
        rt=self.rt_input(nt=20,moving=False).replace('nproc_rgrid=1,1,1','nproc_rgrid=1,2,1')
        ref_folder,ref_run=self.execute('reference_rt',rt,ranks=2,rt=True)
        self.assertEqual(ref_run.returncode,0,ref_run.stdout[-3000:]+ref_run.stderr)
        ref_energy=np.loadtxt(next(ref_folder.glob('*_rt_energy.data')))
        ref_current=np.loadtxt(next(ref_folder.glob('*_rt.data')))
        for ranks,layout,orbitals,total in [(2,'1,1,1',1,'2,1,1'),(8,'1,2,1',2,'4,2,1')]:
            inp=self.base.replace("yn_dc_lcfo='y'","yn_dc_lcfo='y'\n lcfo_eigensolver='scalapack'")
            inp=inp.replace('nproc_rgrid=1,1,1','nproc_rgrid='+layout).replace('nproc_ob=1','nproc_ob='+str(orbitals))
            inp=inp.replace('nproc_rgrid_tot=2,1,1','nproc_rgrid_tot='+total)
            folder,run=self.execute('scalapack'+str(ranks),inp,ranks=ranks)
            self.assertEqual(run.returncode,0,run.stdout[-4000:]+run.stderr)
            self.assertIn('end DC-LCFO complex',run.stdout)
            self.assertIn('LCFO_SCALAPACK',run.stdout)
            eigen=np.loadtxt(next((folder/'data_dcdft/total').glob('*_eigen.data')))[:,3]
            np.testing.assert_allclose(eigen,reference,atol=1e-7,rtol=0)
            tiles=re.findall(r'LCFO_DENSE rank/local_rows/local_cols/global:\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)',run.stdout)
            self.assertEqual(len(tiles),ranks)
            n=int(tiles[0][3])
            self.assertEqual(sum(int(r[1])*int(r[2]) for r in tiles),n*n)
            self.assertTrue(all(int(r[1])*int(r[2])<n*n for r in tiles))
            out,runrt=self.execute('from_scalapack'+str(ranks),rt,ranks=2,rt='scalapack'+str(ranks))
            self.assertEqual(runrt.returncode,0,runrt.stdout[-3000:]+runrt.stderr)
            energy=np.loadtxt(next(out.glob('*_rt_energy.data')))
            current=np.loadtxt(next(out.glob('*_rt.data')))
            np.testing.assert_allclose(energy,ref_energy,atol=1e-7,rtol=0)
            np.testing.assert_allclose(current,ref_current,atol=1e-7,rtol=0)
            print('LCFO ScaLAPACK functional/ranks/eigen/RT error:',self.functional,ranks,
                  np.max(abs(eigen-reference)),np.max(abs(energy-ref_energy)))

class PBEhLCFOScalapack(LCFOScalapack):
    functional='pbeh40_rvv10'
