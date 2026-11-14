import unittest
from types import SimpleNamespace
from unittest.mock import patch
import numpy as np
from salmonpy.coupled import MultiscaleField
from salmonpy.material_response import _Measurements, _StepAdapter
from salmonpy.maxwell import Multiscale1D


class MultiscaleCouplingTests(unittest.TestCase):
    def context(self, size=1):
        return SimpleNamespace(group=SimpleNamespace(role='field'),
            solver_comm=SimpleNamespace(size=size), plan=SimpleNamespace(point_ids=(19, 4)))

    def test_field_preserves_update_and_zero_clock(self):
        model = Multiscale1D(2, 5., vacuum=(3, 4))
        state = model.new_state(dt=.02, steps=3)
        reference = model.new_state(dt=.02, steps=3)
        field = MultiscaleField(self.context(), state)
        for j in (np.ones((2, 3)), np.arange(6).reshape(2, 3), np.zeros((2, 3))):
            np.testing.assert_array_equal(field.advance(j), reference.step(j))
            np.testing.assert_array_equal(field.E, reference.grid_E[3:5])
        self.assertEqual(field.electronic_clock_offset, '0')
        self.assertEqual(field.point_ids, (19, 4))

    def test_material_E_indices_use_integers(self):
        for spacing in (.1, .3, 1.1):
            state = Multiscale1D(2, spacing, vacuum=(300, 50)).new_state(dt=.0001, steps=2)
            field = MultiscaleField(self.context(), state)
            field.advance(np.zeros((2, 3)))
            field.advance(np.array([[1., 2., 3.], [4., 5., 6.]]))
            np.testing.assert_array_equal(field.E, state.grid_E[300:302])
        context = self.context()
        context.plan.point_ids = (19,)
        state = Multiscale1D(1, 1.).new_state(dt=.0001, steps=1)
        field = MultiscaleField(context, state)
        np.testing.assert_array_equal(field.E, state.grid_E)

    def test_nonzero_A_is_rejected_without_mutation(self):
        model = Multiscale1D(2, 5., incident=lambda t: np.broadcast_to([0., 0., 1.], (len(t), 3)))
        state = model.new_state(dt=.02, steps=2)
        before = state.grid_A.copy()
        with self.assertRaisesRegex(ValueError, 'point 19.*vacuum.*delay'):
            MultiscaleField(self.context(), state)
        np.testing.assert_array_equal(state.grid_A, before)
        self.assertEqual(state.step_index, 0)
        with self.assertRaisesRegex(ValueError, 'one field rank'):
            MultiscaleField(self.context(2), state)

    def test_legacy_adapter_only_commits_successful_endpoint(self):
        seen = []
        class State:
            def advance(self, old, new):
                seen.append(old.copy())
                if new[0] < 0: raise ValueError('failed')
                return new*2
        state = _StepAdapter(State())
        np.testing.assert_array_equal(state.step(np.ones(3)), np.ones(3)*2)
        with self.assertRaises(ValueError): state.step(-np.ones(3))
        state.step(np.ones(3)*3)
        np.testing.assert_array_equal(seen[-1], np.ones(3))

    def test_measurement_failure_does_not_fail_work(self):
        measurement = _Measurements()
        with patch('salmonpy.material_response.time.perf_counter', side_effect=OSError('clock unavailable')):
            measurement.stop('A_distribution', measurement.start())
        self.assertIsNone(measurement.seconds['A_distribution'])
        self.assertEqual(measurement.calls['A_distribution'], 1)
        self.assertTrue(measurement.errors)
