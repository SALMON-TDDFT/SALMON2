import json
import unittest
from fractions import Fraction
import numpy as np
from salmonpy.timing import TimeSpec, InputPort, OutputPort, crystal_default, check_loop
from salmonpy.maxwell import Longitudinal, Transverse, couple

FIELD = TimeSpec(.1, (InputPort('j', lag=1, uses_t0=False),), (OutputPort('A'),))


class _State:
    def __init__(self, owner, dt):
        self.owner = owner; self.dt = dt; self.k = 0
    @property
    def time(self): return self.k*self.dt
    @property
    def current(self): return np.array([0., 0., 1e-3])
    def step(self, A):
        self.k += 1; return np.array([0., 0., 1e-3])
    def close(self): self.owner.closed += 1


class Declared:
    """Matter test double with an explicit timing declaration built from ``ports``."""
    backend = 'declared-test'
    def __init__(self, inputs=(InputPort('A'),), outputs=(OutputPort('J_matter'),), clock_offset=0):
        self.inputs = inputs; self.outputs = outputs; self.clock_offset = clock_offset; self.closed = 0
    def new_state(self, *, dt, time0, initial_A):
        state = _State(self, dt)
        state.timing = TimeSpec(dt, self.inputs, self.outputs, clock_offset=self.clock_offset)
        return state


class TimeSpecTests(unittest.TestCase):
    def test_validation(self):
        self.assertEqual(InputPort('A', offset='1/2').offset, Fraction(1, 2))
        self.assertEqual(OutputPort('A', offset=.5).offset, Fraction(1, 2))
        for bad in (.3, 1, -0.5, True, 'x'):
            with self.assertRaises(ValueError): InputPort('A', offset=bad)
        with self.assertRaises(ValueError): InputPort('B')
        with self.assertRaises(ValueError): InputPort('A', sampling='linear')
        with self.assertRaises(ValueError): InputPort('A', lag=-1)
        with self.assertRaises(ValueError): OutputPort('A', at_t0='maybe')
        with self.assertRaises(ValueError): TimeSpec(.1, (InputPort('A'), InputPort('A')), ())
        with self.assertRaises(ValueError): TimeSpec(0., (), ())
        with self.assertRaises(TypeError): TimeSpec(.1, (OutputPort('A'),), ())

    def test_record_is_json(self):
        record = TimeSpec(.1, (InputPort('A', offset='1/2'),), (OutputPort('j'),), clock_offset=.5).to_dict()
        self.assertEqual(json.loads(json.dumps(record))['inputs']['A']['offset'], '1/2')
        self.assertEqual(record['schema'], 'salmonpy.time/1')


class CheckLoopTests(unittest.TestCase):
    def test_crystal_default_couples_with_field(self):
        self.assertEqual(check_loop(FIELD, crystal_default(.1)), ())
        legacy = TimeSpec(.1, (InputPort('A'),), (OutputPort('J_matter'),))
        self.assertEqual(check_loop(FIELD, legacy), ('j = -J_matter (electronic charge -1)',))

    def test_rejections(self):
        cases = {
            'offset': TimeSpec(.1, (InputPort('A'),), (OutputPort('J_matter'),), clock_offset='1/2'),
            'sampling': TimeSpec(.1, (InputPort('A', sampling='interval_average'),), (OutputPort('J_matter'),)),
            'adapter': TimeSpec(.1, (InputPort('A'),), (OutputPort('E'),)),
            'dt': crystal_default(.2),
            't0': crystal_default(.1, t0=1.),
        }
        for name, spec in cases.items():
            with self.subTest(name), self.assertRaises(ValueError):
                check_loop(FIELD, spec, matter_out=spec.outputs[0].quantity)

    def test_zero_total_lag_is_rejected(self):
        implicit = TimeSpec(.1, (InputPort('j', lag=0, uses_t0=False),), (OutputPort('A'),))
        with self.assertRaisesRegex(ValueError, 'lag 0'): check_loop(implicit, crystal_default(.1))

    def test_placeholder_at_t0_is_rejected_when_used(self):
        field = TimeSpec(.1, (InputPort('j', lag=1, uses_t0=True),), (OutputPort('A'),))
        matter = TimeSpec(.1, (InputPort('A'),), (OutputPort('J_matter', at_t0='zero'),))
        with self.assertRaisesRegex(ValueError, 'placeholder'): check_loop(field, matter)
        check_loop(FIELD, matter)  # the field does not use t=0: accepted


class CoupleTimingTests(unittest.TestCase):
    def test_mismatch_rejected_before_stepping_and_state_closed(self):
        matter = Declared(clock_offset='1/2')
        with self.assertRaises(ValueError): couple(matter, Longitudinal(), dt=.1, steps=3)
        self.assertEqual(matter.closed, 1)

    def test_charge_current_output_is_not_negated(self):
        matter = Declared(outputs=(OutputPort('j'),))
        result = couple(matter, Transverse(), dt=.1, steps=2)
        np.testing.assert_array_equal(result.current[:, 0, 2], 1e-3)
        self.assertNotIn('J_matter', result.timing['matter'][0]['outputs'])

    def test_matter_current_is_converted_and_recorded(self):
        result = couple(Declared(), Transverse(), dt=.1, steps=2)
        np.testing.assert_array_equal(result.current[:, 0, 2], -1e-3)
        self.assertEqual(result.timing['field']['inputs']['j']['lag'], 1)
        self.assertFalse(result.timing['field']['inputs']['j']['uses_t0'])
        self.assertEqual(result.timing['matter'][0]['outputs']['J_matter']['at_t0'], 'physical')

    def test_backend_without_declaration_uses_crystal_default(self):
        class Plain(Declared):
            def new_state(self, **kw):
                return _State(self, kw['dt'])
        with self.assertWarnsRegex(UserWarning, 'no timing declaration'):
            result = couple(Plain(), Transverse(), dt=.1, steps=1)
        np.testing.assert_array_equal(result.current[:, 0, 2], 1e-3)  # crystal contract: already j
        self.assertIn('crystal protocol default', result.timing['matter'][0]['integrator'])


if __name__ == '__main__':
    unittest.main()

class FieldSubstepTests(unittest.TestCase):
    def make(self,dt=.01,quantum_dt=.02):
        from types import SimpleNamespace as NS
        from unittest.mock import patch
        from salmonpy.coupled import NativeEHField
        from salmonpy.material_response import _response_timing
        field=object.__new__(NativeEHField)
        comm=NS(allgather=lambda v:[v],bcast=lambda v,root=0:v)
        context=NS(group=NS(role='field'),plan=NS(point_ids=(1,)),is_leader=True,
                   solver_comm=comm,leaders=comm)
        field.context=context;field.dt=dt;field.step_index=0;field.closed=False;field.failed=False
        field.point_ids=(1,);field.A=np.zeros((1,3));field.E=field.A.copy()
        field.capacity=100;field.current_rows_by_rank=None;field.manifest={};field._write_manifest=lambda:None
        response=NS(context=context,dt=quantum_dt,step_index=0,closed=False,failed=False,
                    _field_contracts=None,timing_records=lambda:{1:_response_timing(quantum_dt,'test').to_dict()})
        def agree(error):
            if error:raise ValueError(error)
        response.agree=agree;events=[];held=[]
        def advance(j):
            events.append('field');held.append(j.copy());field.step_index+=1;field.A=field.A+.01
            return field.A.copy()
        field.advance_owned=advance;field.advance=advance
        def material(r,A,error,**kw):
            agree(error);events.append('electrons');r.step_index+=1;return np.full((1,3),r.step_index)
        self.addCleanup(patch.stopall)
        patch('salmonpy.coupled._select_response_admission',return_value=False).start()
        patch('salmonpy.coupled._response_from_field',side_effect=material).start()
        return response,field,events,held

    def test_field_updates_hold_current_and_use_one_electronic_endpoint(self):
        from salmonpy.coupled import FieldSubcycledStepper
        for count in (1,2,4,10):
            with self.subTest(field_steps=count):
                r,f,events,held=self.make(.02/count,.02)
                stepper=FieldSubcycledStepper(r,f,field_steps=count)
                for n in range(3):
                    frame=stepper.step()
                    self.assertEqual(frame['field_step'],(n+1)*count)
                    self.assertEqual(frame['electronic_step'],n+1)
                    self.assertAlmostEqual(frame['electronic_field_clock'],(n+1)*.02+.01/count)
                self.assertEqual(events,(['field']*count+['electrons'])*3)
                for i,j in enumerate(held):np.testing.assert_array_equal(j,float(i//count))

    def test_invalid_ratio_or_partial_failure_is_not_silently_reused(self):
        from salmonpy.coupled import FieldSubcycledStepper
        for count in (0,True,1.5,3):
            r,f,_,_=self.make()
            with self.subTest(count=count),self.assertRaises(ValueError):FieldSubcycledStepper(r,f,field_steps=count)
            self.assertEqual(f.step_index,0)
        r,f,events,_=self.make();stepper=FieldSubcycledStepper(r,f,field_steps=2)
        def fail(field):raise ValueError('observation failed')
        with self.assertRaisesRegex(ValueError,'observation failed'):stepper.step(after_field_step=fail)
        self.assertTrue(r.failed);self.assertEqual(events,['field'])
        with self.assertRaises(ValueError):stepper.step()
