"""Fixed sphere takes precedence over norm diagnostics in distributed SCF."""
import os
import unittest
import test_adaptive_scf as helper

@unittest.skipUnless(os.environ.get('SALMON_TEST_EXE') and os.environ.get('SALMON_TEST_MPIEXEC'),'MPI executable required')
class RadiusPriority(unittest.TestCase):
    run_case=helper.AdaptiveSCF.run_case
    def test_spatial_radius_priority(self):
        for functional in ('hse06','pbeh40'):
            a=self.run_case(0,ranks=2,functional=functional,radius=5,pre_scf=1e-4)
            b=self.run_case(.999,ranks=2,functional=functional,radius=5,pre_scf=1e-4)
            c=self.run_case(1,ranks=4,functional=functional,radius=5,pre_scf=1e-4)
            self.assertLess(abs(a-b),2e-6)
            self.assertLess(abs(b-c),2e-6)

    def test_full_cell_radius_needs_no_converged_gauge(self):
        self.run_case(0,ranks=2,functional='hse06',radius=100,pre_scf=1e-4,
                      mlwf_maxiter=1,mlwf_tolerance=1e-30)

if __name__=='__main__': unittest.main()
