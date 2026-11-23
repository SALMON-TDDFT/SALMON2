import unittest
from types import SimpleNamespace
import numpy as np
from salmonpy import material_response as mr
from salmonpy.errors import InputError
from salmonpy.field_input import FieldPacket, FieldValue
from salmonpy.timing import TimeSpec, InputPort, OutputPort


class State:
    def __init__(self, dt):
        self.timing = TimeSpec(dt, (InputPort('A'),), (OutputPort('j'),))
        self.time = 0.
        self.A = np.zeros(3)
        self.closed = 0
    def advance(self, old, new):
        np.testing.assert_array_equal(old, self.A)
        self.A = new.copy()
        self.time += self.timing.dt
        return -2*new
    def close(self): self.closed += 1


def binding(factory=State):
    backend = SimpleNamespace(identity='test/1', execution='leader_only', reduction='leader',
        create_state=lambda initial, point, dt, context: factory(dt))
    return mr.MaterialBinding('test', object(), backend)


def packet(n, phase=0., value=(1., 0., 0.)):
    return FieldPacket(A=FieldValue(value, (n+phase)*.1, [0,0,0], 'A/c:au'))


class PacketBinding(unittest.TestCase):
    def test_both_phases_independent_history_and_sign(self):
        self.assertTrue(hasattr(mr, 'packetize_A_binding'))
        for phase, offset in [('endpoint',0.),('yee_half_step',.5)]:
            original = binding()
            wrapped = mr.packetize_A_binding(original, field_phase=phase)
            self.assertIs(wrapped.provider, original.provider)
            a = wrapped.backend.create_state(None, 1, .1, None)
            b = wrapped.backend.create_state(None, 2, .1, None)
            np.testing.assert_array_equal(a.step_fields(packet(1, offset)), [-2,0,0])
            self.assertEqual(b.time, 0.)
            np.testing.assert_array_equal(b.step_fields(packet(1, offset)), [-2,0,0])
            self.assertEqual(a.timing.outputs[0].offset, offset)
            self.assertEqual(a.time, .1)
            a.close(); a.close(); b.close()

    def test_bad_input_before_update_and_clock_failure(self):
        a = mr.packetize_A_binding(binding(), field_phase='endpoint').backend.create_state(None,1,.1,None)
        with self.assertRaises(InputError): a.step_fields(packet(1,.5))
        self.assertEqual(a.time,0.)
        np.testing.assert_array_equal(a.step_fields(packet(1)),[-2,0,0])
        a.close()
        with self.assertRaises(InputError): a.step_fields(packet(2))
        with self.assertRaises(InputError): mr.packetize_A_binding(binding(),field_phase='implicit')

    def test_invalid_created_state_closed(self):
        raw=State(.1)
        raw.A[0]=1
        with self.assertRaises(InputError):
            mr.packetize_A_binding(binding(lambda dt:raw),field_phase='endpoint').backend.create_state(None,1,.1,None)
        self.assertEqual(raw.closed,1)

    def test_callback_failure_invalidates(self):
        raw=State(.1)
        def fail(old,new): raise RuntimeError('backend failure')
        raw.advance=fail
        a=mr.packetize_A_binding(binding(lambda dt:raw),field_phase='endpoint').backend.create_state(None,1,.1,None)
        with self.assertRaises(RuntimeError): a.step_fields(packet(1))
        with self.assertRaises(InputError): a.step_fields(packet(1))
        a.close()

    def test_native_sbe_matches_direct(self):
        from salmonpy.sbe import SBEInitialState
        from salmonpy.materials import bind_material
        p=np.zeros((2,2,3,1),complex)
        p[0,1,0,0]=p[1,0,0,0]=.3
        model=SBEInitialState(eigen=[[-.5],[.4]],momentum=p,
            nonlocal_velocity=p*.1,weights=[1.],lattice=np.eye(3)*4,electrons=2)
        for phase, offset in [('endpoint',0.),('yee_half_step',.5)]:
            base=bind_material(model,key='sbe',identity='test-sbe/1')
            wrapped=mr.packetize_A_binding(base,field_phase=phase)
            prepared=wrapped.provider.prepare(None)
            a=wrapped.backend.create_state(prepared,1,.1,None)
            b=wrapped.backend.create_state(prepared,2,.1,None)
            direct=model.new_state(dt=.1)
            try:
                for n in range(1,11):
                    values=[n*.001,0.,0.]
                    expected=direct.step(values)
                    np.testing.assert_array_equal(a.step_fields(packet(n,offset,values)),expected)
                    np.testing.assert_array_equal(b.step_fields(packet(n,offset,values)),expected)
                self.assertEqual(a.clock_mapping['original']['outputs']['j']['offset'],'0')
            finally:
                a.close();b.close();direct.close()

    def test_reused_backend_state_rejected_without_closing_owner(self):
        raw=State(.1)
        backend=mr.packetize_A_binding(binding(lambda dt:raw),field_phase='endpoint').backend
        first=backend.create_state(None,1,.1,None)
        with self.assertRaises(InputError): backend.create_state(None,2,.1,None)
        self.assertEqual(raw.closed,0)
        first.close()
        self.assertEqual(raw.closed,1)

    def test_collective_backend_refused_before_callbacks(self):
        original=binding()
        original.backend.execution='collective'
        with self.assertRaises(InputError):
            mr.packetize_A_binding(original,field_phase='endpoint')

    def test_cross_binding_alias_refused_without_double_close(self):
        raw=State(.1)
        first_backend=mr.packetize_A_binding(binding(lambda dt:raw),field_phase='endpoint').backend
        second_backend=mr.packetize_A_binding(binding(lambda dt:raw),field_phase='endpoint').backend
        first=first_backend.create_state(None,1,.1,None)
        second=second_backend.create_state(None,2,.1,None)
        with self.assertRaises(InputError): mr._admit_response_state(second,{1:first})
        second.close()
        self.assertEqual(raw.closed,0)
        first.close()
        self.assertEqual(raw.closed,1)
