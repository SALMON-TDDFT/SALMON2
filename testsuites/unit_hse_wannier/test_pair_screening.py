import sys
from pathlib import Path
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'samples/dc_hse'))
from pair_screening import analyze

class PairScreening(unittest.TestCase):
    def fixture(self):
        rng=np.random.default_rng(621)
        q=rng.normal(size=(3,4,5,6))+1j*rng.normal(size=(3,4,5,6))
        q/=np.sqrt(np.sum(abs(q)**2,axis=(1,2,3)))[:,None,None,None]
        return q*np.sqrt(np.array([1.,.37,.02]))[:,None,None,None]

    def test_zero_budget_is_exact(self):
        r=analyze(self.fixture(),np.array([.6,.8,.9]),.11,0)
        self.assertEqual(r['kept_ordered_pairs'],9)
        self.assertEqual(r['exchange_error_Ha'],0)
        self.assertEqual(r['relative_action_error'],0)
        self.assertLess(r['screened_metric_antihermitian_relative'],1e-12)
        self.assertGreater(r['screened_metric_min_eigenvalue'],0)

    def test_positive_budget_bound(self):
        r=analyze(self.fixture(),np.ones(3),.11,1e6)
        self.assertEqual(r['kept_ordered_pairs'],3)
        self.assertGreaterEqual(r['exchange_error_Ha'],-1e-12)
        self.assertLessEqual(r['exchange_error_Ha'],r['discarded_energy_bound_Ha']+1e-12)
        self.assertLessEqual(r['discarded_energy_bound_Ha'],1e6)
        self.assertGreater(r['relative_action_error'],0)

    def test_disjoint_support_can_skip_exactly(self):
        q=np.zeros((2,4,4,4),complex);q[0,0,0,0]=1;q[1,2,2,2]=1j
        r=analyze(q,np.ones(3),.11,0)
        self.assertEqual(r['kept_ordered_pairs'],2)
        self.assertEqual(r['exchange_error_Ha'],0)
        self.assertEqual(r['relative_action_error'],0)

    def test_exact_against_dense_real_space_operator(self):
        q=self.fixture();spacing=np.array([.6,.8,.9]);omega=.11;shape=q.shape[1:]
        axes=[2*np.pi*np.fft.fftfreq(s,d=h) for s,h in zip(shape,spacing)]
        g2=sum(x*x for x in np.meshgrid(*axes,indexing='ij'))
        kernel=np.full(shape,np.pi/omega**2)
        np.divide(4*np.pi*(-np.expm1(-g2/(4*omega**2))),g2,out=kernel,where=g2>0)
        spatial=np.fft.ifftn(kernel)
        points=np.indices(shape).reshape(3,-1).T
        displacement=(points[:,None,:]-points[None,:,:])%np.array(shape)
        dense_kernel=spatial[tuple(displacement.transpose(2,0,1))]
        flat=q.reshape(len(q),-1).T
        density_matrix=flat@flat.conj().T
        action=-(density_matrix*dense_kernel)@flat
        reference=.25*np.prod(spacing)*np.vdot(flat,action).real
        result=analyze(q,spacing,omega,0)
        self.assertAlmostEqual(result['full_exchange_Ha'],reference,places=12)
        rng=np.random.default_rng(2)
        u=np.linalg.qr(rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)))[0]
        rotated=(flat@u).T.reshape(q.shape)
        self.assertAlmostEqual(analyze(rotated,spacing,omega,0)['full_exchange_Ha'],reference,places=12)

    def test_invalid_parameters(self):
        for omega,budget in [(0,0),(.11,-1),(float('nan'),0),(.11,float('nan'))]:
            with self.assertRaises(ValueError): analyze(self.fixture(),np.ones(3),omega,budget)

if __name__=='__main__': unittest.main()
