import sys
import unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'samples/dc_hse'))
from locality import eigen_diagnostics, axial_tail, support_sweep

class LocalityTest(unittest.TestCase):
    def test_eigen_residual_and_overlap(self):
        psi=np.eye(3,dtype=complex)/np.sqrt(.5)
        hpsi=psi@np.diag([1.,2.,4.])
        r=eigen_diagnostics(psi,hpsi,.5,np.array([2.,2.,0.]))
        self.assertLess(r['max_residual_Ha'],1e-14)
        self.assertLess(r['max_orthogonality_error'],1e-14)
        hpsi[2,0]+=.1/np.sqrt(.5)
        r=eigen_diagnostics(psi,hpsi,.5,np.array([2.,2.,0.]))
        self.assertAlmostEqual(r['max_residual_Ha'],.1)
    def test_periodic_tail(self):
        q=np.zeros((1,8,2,2),complex)
        q[0,0,0,0]=np.sqrt(.75);q[0,7,0,0]=.5
        result=axial_tail(q,[1.,1.,1.],[0.,1.,4.])
        self.assertEqual(result['tail_fraction'].shape,(3,1))
        self.assertAlmostEqual(result['tail_fraction'][1,0],0.)
        self.assertAlmostEqual(result['tail_fraction'][2,0],0.)

    def test_support_sweep_full_limit_and_hermiticity(self):
        rng=np.random.default_rng(21)
        q=rng.normal(size=(2,6,2,2))+1j*rng.normal(size=(2,6,2,2))
        q/=np.sqrt(np.sum(abs(q)**2,axis=(1,2,3)))[:,None,None,None]
        reports=support_sweep(q,[.8,.7,.6],.11,[2.4,.8],2)
        self.assertAlmostEqual(reports[0]['exchange_error_Ha'],0.,places=12)
        self.assertAlmostEqual(reports[0]['relative_action_error'],0.,places=12)
        self.assertGreater(reports[1]['discarded_norm_fraction'],0.)
        self.assertLess(reports[1]['target_metric_antihermitian_relative'],1e-12)
        self.assertFalse(reports[1]['certified_for_scf'])

if __name__=='__main__':unittest.main()
