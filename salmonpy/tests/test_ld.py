import unittest
import numpy as np
from salmonpy import LorentzDrude

class LDTests(unittest.TestCase):
    def material(self):
        return LorentzDrude(plasma_frequency=1., strength=[1., .2], resonance=[0., 2.], damping=[.1, .2])

    def test_vacuum(self):
        m = LorentzDrude(plasma_frequency=0., strength=[1.], resonance=[0.], damping=[0.])
        r = m.respond(np.ones((5,3)), dt=.1)
        np.testing.assert_array_equal(r.current, 0.)
        np.testing.assert_array_equal(r.polarization, 0.)

    def test_drude_constant_field(self):
        m = LorentzDrude(plasma_frequency=2., strength=[1.], resonance=[0.], damping=[0.])
        r = m.respond(np.ones((11,3)), dt=.01)
        np.testing.assert_allclose(r.current[:,0], r.time / np.pi, atol=1e-15)
        np.testing.assert_allclose(r.polarization[:,0], r.time**2 / (2*np.pi), atol=1e-15)

    def test_batch_step_and_independence(self):
        m=self.material(); e=np.arange(33).reshape(11,3)*.001
        r=m.respond(e, dt=.01)
        s=m.new_state(dt=.01, initial_field=e[0]); other=m.new_state(dt=.01)
        for i in range(1,len(e)):
            p,j=s.step(e[i])
            np.testing.assert_array_equal(j,r.current[i])
            np.testing.assert_array_equal(p,r.polarization[i])
        np.testing.assert_array_equal(other.current,0.)
        np.testing.assert_array_equal(m.respond(e,dt=.01).current,r.current)

    def test_dielectric(self):
        m=LorentzDrude(plasma_frequency=2.,strength=[.5],resonance=[3.],damping=[.2])
        np.testing.assert_allclose(m.permittivity([1.,2.]),1+2/(9-np.array([1.,2.])**2-1j*.2*np.array([1.,2.])))

    def test_damped_drude_convergence(self):
        m=LorentzDrude(plasma_frequency=1.,strength=[1.],resonance=[0.],damping=[.7])
        errors=[]
        for dt in [.1,.05]:
            r=m.respond(np.ones((round(2/dt)+1,3)),dt=dt)
            exact=(1-np.exp(-.7*r.time[-1]))/(4*np.pi*.7)
            errors.append(abs(r.current[-1,0]-exact))
        self.assertGreater(errors[0]/errors[1],3.9)
        self.assertLess(errors[0]/errors[1],4.1)

    def test_complex_rejected(self):
        with self.assertRaises(ValueError):
            self.material().respond(np.ones((2,3),dtype=complex)*(1+1j),dt=.1)
        with self.assertRaises(ValueError):
            LorentzDrude(plasma_frequency=1.,strength=np.array([1+1j]),resonance=[0.],damping=[0.])

    def test_invalid(self):
        for x in [float('nan'), float('inf'), -1.]:
            with self.assertRaises(ValueError):
                LorentzDrude(plasma_frequency=x,strength=[1.],resonance=[0.],damping=[0.])
        m=self.material()
        for dt in [0.,-1.,float('nan')]:
            with self.assertRaises(ValueError): m.respond(np.zeros((2,3)),dt=dt)
        with self.assertRaises(ValueError): m.respond([[0,0,float('inf')]],dt=.1)
        with self.assertRaises(ValueError): m.respond([1,2,3],dt=.1)

if __name__=='__main__': unittest.main()
