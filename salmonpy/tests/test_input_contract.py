import unittest
from dataclasses import FrozenInstanceError, asdict
from unittest.mock import patch
import numpy as np
from salmonpy import LorentzDrude


class InputContractTests(unittest.TestCase):
    def model(self):
        return LorentzDrude(plasma_frequency=1., strength=[1.], resonance=[0.], damping=[.1])

    def test_resolved_config(self):
        m=self.model(); s=m.new_state(dt=.1)
        self.assertEqual(m.config.strength, (1.,))
        self.assertEqual(m.config.units, 'atomic_gaussian')
        self.assertEqual(s.config.initial_field, (0.,0.,0.))
        self.assertEqual(s.config.model, m.config)
        self.assertEqual(s.config.dt, .1)
        self.assertEqual(s.config.initial_time, 0.)
        self.assertEqual(asdict(s.config)['model']['plasma_frequency'],1.)
        with self.assertRaises(FrozenInstanceError): m.config.plasma_frequency=2.
        with self.assertRaises(FrozenInstanceError): s.config.dt=.2

    def test_public_mutation_rejected(self):
        m=self.model(); s=m.new_state(dt=.1)
        for obj,name,value in [(m,'plasma_frequency',2.),(m,'strength',[2.]),
                               (m,'config',None),(s,'dt',.2),(s,'config',None),
                               (s,'material',self.model())]:
            with self.subTest(name=name), self.assertRaises(AttributeError):
                setattr(obj,name,value)
        with self.assertRaises(TypeError): s.step([0,0,0],dt=.2)

    def test_inputs_copied(self):
        f=np.array([1.]); e=np.array([1.,0.,0.])
        m=LorentzDrude(plasma_frequency=1.,strength=f,resonance=[0.],damping=[.1])
        s=m.new_state(dt=.1,initial_field=e)
        f[0]=9.;e[0]=9.
        self.assertEqual(m.config.strength,(1.,))
        self.assertEqual(s.config.initial_field,(1.,0.,0.))
        ref=self.model().new_state(dt=.1,initial_field=[1,0,0])
        np.testing.assert_array_equal(s.step([1,0,0]),ref.step([1,0,0]))

    def test_invalid_step_is_transactional(self):
        s=self.model().new_state(dt=.1); ref=self.model().new_state(dt=.1)
        s.step([1,0,0]);ref.step([1,0,0])
        before=(s.time,s.current.copy(),s.polarization.copy(),s.config)
        for field in ([1,2], [np.nan,0,0], [np.inf,0,0], [1j,0,0]):
            with self.assertRaises(ValueError): s.step(field)
            self.assertEqual(s.time,before[0])
            np.testing.assert_array_equal(s.current,before[1])
            np.testing.assert_array_equal(s.polarization,before[2])
            self.assertEqual(s.config,before[3])
        with patch('salmonpy.ld._native',side_effect=RuntimeError('unavailable')):
            with self.assertRaises(RuntimeError): s.step([2,0,0])
        np.testing.assert_array_equal(s.step([2,0,0]),ref.step([2,0,0]))

    def test_scalar_validation(self):
        for dt in (True,[.1],np.array([.1]),'0.1'):
            with self.subTest(dt=dt), self.assertRaises(ValueError): self.model().new_state(dt=dt)

    def test_nonfinite_native_output_preserves_history(self):
        s=self.model().new_state(dt=.1); ref=self.model().new_state(dt=.1)
        s.step([1,0,0]); ref.step([1,0,0])
        before=(s.time,s.current.copy(),s.polarization.copy())
        def bad_native(*args):
            args[-4][0]=99.  # modified temporary per-pole P
            args[-2][0]=float('nan')  # invalid output P
        with patch('salmonpy.ld._native',return_value=bad_native):
            with self.assertRaises(FloatingPointError): s.step([2,0,0])
        self.assertEqual(s.time,before[0])
        np.testing.assert_array_equal(s.current,before[1])
        np.testing.assert_array_equal(s.polarization,before[2])
        np.testing.assert_array_equal(s.step([2,0,0]),ref.step([2,0,0]))
