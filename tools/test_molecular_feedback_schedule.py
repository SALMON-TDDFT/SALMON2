"""Fast independent checks for the development oscillator scheduling fixture."""
import unittest
import numpy as np
from salmonpy.dipole_source import _midpoint_oscillator


class OscillatorChecks(unittest.TestCase):
    def test_free_analytic_order(self):
        errors=[]
        for dt in (.04,.02,.01,.005):
            p=np.array([1.,0.,0.]);v=np.zeros(3)
            for _ in range(round(1/dt)):
                p,v,_,stats=_midpoint_oscillator(p,v,np.zeros(3),np.zeros(3),dt=dt,omega=2.)
                balance=stats[5]
                self.assertLess(abs(balance),1e-12)
            errors.append(abs(p[0]-np.cos(2.)))
        self.assertTrue(np.all(np.log2(np.array(errors[:-1])/errors[1:])>1.99))

    def test_work_identity_and_sign_negative_control(self):
        p=np.array([.2,.1,0.]);v=np.array([.3,0.,0.]);E=np.array([.5,.2,0.])
        new_p,new_v,_,stats=_midpoint_oscillator(p,v,E,E,dt=.01,omega=2.)
        balance=stats[5]
        energy_change=.5*np.sum(new_v**2+4*new_p**2-v**2-4*p**2)
        work=np.dot(E,new_p-p)
        self.assertLess(abs(energy_change-work),1e-12)
        self.assertGreater(abs(energy_change+work),1e-5)
        self.assertLess(abs(balance),1e-12)


if __name__=='__main__':unittest.main()
