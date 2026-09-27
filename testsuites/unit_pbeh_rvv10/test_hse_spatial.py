"""HSE DC initialization -> fixed-ion mesh RT, serial/distributed parity."""
import os
import unittest
import numpy as np
import test_ehrenfest as helpers

@unittest.skipUnless(os.environ.get('SALMON_TEST_EXE') and os.environ.get('SALMON_TEST_MPIEXEC'),
                     'SALMON_TEST_EXE and SALMON_TEST_MPIEXEC required')
class HSESpatial(unittest.TestCase):
    functional='hse06'
    setUpClass=classmethod(helpers.RealspaceEhrenfest.setUpClass.__func__)
    execute=classmethod(helpers.RealspaceEhrenfest.execute.__func__)
    rt_input=helpers.RealspaceEhrenfest.rt_input
    pulse_input=helpers.RealspaceEhrenfest.pulse_input

    def test_serial_spatial_parity(self):
        for field in ('impulse','pulse'):
            results=[]
            for ranks,layout in [(1,'1,1,1'),(2,'1,2,1'),(4,'1,2,2')]:
                inp=(self.rt_input(nt=20,moving=False) if field=='impulse' else
                     self.pulse_input(dt=.08,nt=120,moving=False))
                if ranks==1:
                    inp=inp.replace("xc='hse06'","xc='hse06'\n yn_hse_wannier='y'")
                inp=inp.replace('nproc_rgrid=1,1,1','nproc_rgrid='+layout)
                folder,run=self.execute(field+str(ranks),inp,ranks=ranks,rt=True)
                self.assertEqual(run.returncode,0,run.stdout[-3000:]+run.stderr)
                self.assertIn('end SALMON',run.stdout)
                self.assertIn('EXX_SPATIAL' if ranks>1 else 'HSE_WANNIER',run.stdout)
                self.assertNotIn('Native LCFO RT active',run.stdout)
                energy=np.loadtxt(next(folder.glob('*_rt_energy.data')))
                data=np.loadtxt(next(folder.glob('*_rt.data')))
                self.assertTrue(np.isfinite(energy).all() and np.isfinite(data).all())
                results.append((energy,data))
                for actual,reference in zip(results[-1],results[0]):
                    self.assertLess(np.max(abs(actual-reference)),2e-8)
                print('HSE parity field/ranks/energy/current:',field,ranks,
                      np.max(abs(energy-results[0][0])),np.max(abs(data-results[0][1])))

    def test_hse_md_still_rejected(self):
        _,run=self.execute('md_rejected',self.rt_input(moving=True),rt=True)
        self.assertNotEqual(run.returncode,0)
        self.assertIn('HSE06 requires fixed-ion',run.stdout+run.stderr)

    def test_legacy_full_propagator_admission(self):
        inp=self.rt_input(nt=2,moving=False)+"\n&propagation\n propagator='hse_taylor4_full'\n/\n"
        _,run=self.execute('legacy_full',inp,rt=True)
        # This noncubic fixture reaches the legacy kernel's existing grid guard.
        # It must not be rejected earlier by the new Wannier propagator admission.
        self.assertIn('end complex DC-LCFO wavefunction reconstruction',run.stdout)
        self.assertIn('HSE06: cubic grid and full cubic k mesh required',run.stdout+run.stderr)
        self.assertNotIn('HSE Wannier RT:',run.stdout+run.stderr)
