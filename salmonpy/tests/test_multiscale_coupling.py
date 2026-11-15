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

    def test_packet_session_rejected_before_field_advances(self):
        from salmonpy.coupled import CoupledStepper
        from salmonpy.errors import InputError
        context=self.context()
        state=Multiscale1D(2,5.).new_state(dt=.02,steps=2)
        field=MultiscaleField(context,state)
        def agree(error):
            if error: raise InputError(error)
        response=SimpleNamespace(context=context,closed=False,failed=False,
                                 step_index=0,dt=.02,_field_contracts={},agree=agree)
        with self.assertRaisesRegex(InputError,'explicit sampling'):
            CoupledStepper(response,field)
        self.assertEqual(state.step_index,0)

    def test_measurement_failure_does_not_fail_work(self):
        measurement = _Measurements()
        with patch('salmonpy.material_response.time.perf_counter', side_effect=OSError('clock unavailable')):
            measurement.stop('A_distribution', measurement.start())
        self.assertIsNone(measurement.seconds['A_distribution'])
        self.assertEqual(measurement.calls['A_distribution'], 1)
        self.assertTrue(measurement.errors)


class FDTDCommonClockTests(unittest.TestCase):
    """Clock admission only: no field/native solver or MPI process is started."""
    def check(self, records):
        from salmonpy.coupled import CoupledStepper, NativeEHField
        class LocalComm:
            def bcast(self, value, root=0): return value
        context = SimpleNamespace(group=SimpleNamespace(role='field'),
            plan=SimpleNamespace(point_ids=(19, 4)), is_leader=True, field_group=0,
            leaders=LocalComm(), solver_comm=LocalComm())
        def agree(error):
            if error: raise ValueError(error)
        response = SimpleNamespace(context=context, timing_records=lambda: records, agree=agree)
        field = object.__new__(NativeEHField)
        field.dt = .02
        return CoupledStepper._check_timing(response, field)

    def records(self, second_dt=.02, second_clock=0):
        from salmonpy.timing import TimeSpec, InputPort, OutputPort
        return {19: TimeSpec(.02, (InputPort('A'),), (OutputPort('j'),),
                    integrator='material one').to_dict(),
                4: TimeSpec(second_dt, (InputPort('A'),), (OutputPort('j'),),
                    clock_offset=second_clock, integrator='material two').to_dict()}

    def test_distinct_materials_share_fdtd_clock_and_half_step(self):
        record = self.check(self.records())
        self.assertEqual(record['electronic_clock_offset'], '1/2')
        self.assertEqual(record['field']['inputs']['j']['offset'], '1/2')
        self.assertEqual(record['field']['inputs']['j']['lag'], 1)
        self.assertEqual(set(record['matter']), {19, 4})
        self.assertNotEqual(record['matter'][19]['integrator'], record['matter'][4]['integrator'])

    def test_integer_subcycle_ratios_are_rejected_without_an_adapter(self):
        for dt in (.01, .04, .015):
            with self.subTest(dt=dt), self.assertRaisesRegex(ValueError, 'point 4:.*dt differ'):
                self.check(self.records(second_dt=dt))

    def test_equal_dt_does_not_hide_wrong_half_step(self):
        with self.assertRaisesRegex(ValueError, 'point 4:.*offset'):
            self.check(self.records(second_clock='1/2'))
