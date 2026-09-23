import unittest
import numpy as np
from analyze import response, difference, peak_metrics, HARTREE_EV

class AnalysisTests(unittest.TestCase):
    def test_damped_oscillator_and_amplitude(self):
        t=np.arange(1,20001)*0.05
        j=np.exp(-0.02*t)*np.cos(0.15*t)
        energy=np.linspace(2,6,401)
        eps=response(t,0.001*j,0.001,energy)
        np.testing.assert_allclose(eps,response(t,0.002*j,0.002,energy))
        metrics=peak_metrics(energy,eps.imag,3.5,4.5)
        self.assertLess(abs(metrics['peak_eV']-0.15*HARTREE_EV),0.08)
        self.assertGreater(metrics['area_eV'],0)
    def test_pump_cancellation_and_delay(self):
        t=np.arange(1,1001)*0.1
        pump=np.sin(t); probe=np.where(t>20,np.exp(-(t-20)/10),0)
        np.testing.assert_allclose(difference(t,pump+probe,t,pump),probe,atol=1e-15)
        energy=np.array([2.,3.])
        np.testing.assert_allclose(response(t,probe,1,energy,20),
                                   response(t[t>20]-20,probe[t>20],1,energy))
    def test_invalid_grids(self):
        with self.assertRaises(ValueError): difference(np.arange(3),np.ones(3),np.arange(3)+0.1,np.ones(3))
        with self.assertRaises(ValueError): response(np.array([1.,2.,4.]),np.ones(3),1,np.array([2.]))
        with self.assertRaises(ValueError): response(np.arange(3.),np.ones(3),0,np.array([2.]))
        with self.assertRaises(ValueError): peak_metrics(np.array([1.,2.]),np.ones(2),3,4)

if __name__=='__main__': unittest.main()
