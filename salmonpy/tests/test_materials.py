"""One serial material contract, conformance checks and point-MPI plumbing."""
from dataclasses import dataclass
from types import SimpleNamespace
import importlib
import unittest
import numpy as np
from salmonpy.fields import FieldSamples
from salmonpy.errors import InputError
from salmonpy.timing import TimeSpec,InputPort,OutputPort


@dataclass(frozen=True)
class Model:
    scale: float = 2.
    backend: str = 'test-memory/1'
    output: str = 'j'

    def new_state(self, *, dt, time0=0., initial_A=(0.,0.,0.)):
        scale=self.scale;output=self.output
        class State:
            def __init__(self):
                self.time=time0;self.current=np.array([.2,.3,.4]);self.last=np.array(initial_A,copy=True)
                self.timing=TimeSpec(dt,(InputPort('A'),),(OutputPort(output),),t0=time0)
                self.closed=False
            def step(self,A):
                self.current= .8*self.current+scale*(A-self.last)
                self.last=A.copy();self.time+=dt
                return self.current.copy()
            def close(self):self.closed=True
        return State()


class Materials(unittest.TestCase):
    def setUp(self):
        self.api=importlib.import_module('salmonpy.materials')
        self.drive=FieldSamples(np.arange(6)*.1,np.arange(18).reshape(6,3)*.01,'A')

    def test_checker_and_binding_match_serial(self):
        from salmonpy.crystal import respond
        initial=Model();report=self.api.check_material(initial,self.drive)
        self.assertEqual(report.steps,5);self.assertEqual(report.maximum_repeat_difference,0.)
        with self.assertRaises(AttributeError):report.steps=1
        binding=self.api.bind_material(initial,key='custom',identity='scale2/algorithm1')
        self.assertEqual(binding.backend.execution,'leader_only');self.assertEqual(binding.backend.reduction,'leader')
        state=binding.backend.create_state(binding.provider.prepare(None),7,.1,None)
        zero_drive=FieldSamples(self.drive.time,self.drive.values-self.drive.values[0],'A')
        reference=respond(initial,zero_drive)
        np.testing.assert_array_equal(state.current,reference.current[0])
        for A,J in zip(zero_drive.values[1:],reference.current[1:]):np.testing.assert_array_equal(state.step(A),J)
        state.close();state.close()
        with self.assertRaises(InputError):state.step(np.zeros(3))

    def test_declared_matter_current_is_converted_once(self):
        from salmonpy.crystal import respond
        model=Model(output='J_matter');binding=self.api.bind_material(model,key='x',identity='matter/1')
        state=binding.backend.create_state(binding.provider.prepare(None),1,.1,None)
        self.addCleanup(state.close)
        reference=respond(model,FieldSamples(self.drive.time,self.drive.values-self.drive.values[0],'A'))
        np.testing.assert_array_equal(state.current,reference.current[0])
        self.assertIsNotNone(state.timing.output('j'));self.assertIsNone(state.timing.output('J_matter'))
        for A,J in zip((self.drive.values-self.drive.values[0])[1:],reference.current[1:]):np.testing.assert_array_equal(state.step(A),J)

    def test_binding_preserves_input_lag_and_initial_usage(self):
        class Delayed:
            backend='declared-delay'
            def new_state(self,**kw):
                state=Model().new_state(**kw)
                state.timing=TimeSpec(kw['dt'],(InputPort('A',lag=1,uses_t0=False),),(OutputPort('j'),))
                return state
        binding=self.api.bind_material(Delayed(),key='lag',identity='lag1/1')
        state=binding.backend.create_state(binding.provider.prepare(None),1,.1,None)
        self.addCleanup(state.close)
        self.assertEqual(state.timing.input('A'),InputPort('A',lag=1,uses_t0=False))

    def test_explicit_identity_and_single_rank_scope(self):
        for kwargs in ({'key':'x','identity':''},{'key':'','identity':'v1'},{'key':'x','identity':None}):
            with self.assertRaises(InputError):self.api.bind_material(Model(),**kwargs)
        binding=self.api.bind_material(Model(),key='x',identity='v1')
        context=SimpleNamespace(solver_comm=SimpleNamespace(size=2))
        with self.assertRaises(InputError):binding.backend.create_state(binding.provider.prepare(None),1,.1,context)

    def test_bad_drive_rejected_before_creation(self):
        class NoCreate:
            backend='none'
            def new_state(self,**kw):raise AssertionError('must reject before creation')
        for drive in (FieldSamples([0,1,3],np.zeros((3,3)),'A'),FieldSamples([0,1],np.zeros((2,3)),'E')):
            with self.assertRaises(InputError):self.api.check_material(NoCreate(),drive)

    def test_cumulative_drive_clock_mismatch_matches_serial_rejection(self):
        from salmonpy.crystal import respond
        times=np.r_[0.,np.cumsum(np.r_[.1,np.full(8,.1+5e-10)])]
        drive=FieldSamples(times,np.zeros((10,3)),'A')
        with self.assertRaises(InputError):respond(Model(),drive)
        with self.assertRaises(InputError):self.api.check_material(Model(),drive)

    def test_shared_state_and_shared_arrays_rejected(self):
        for same_object in (True,False):
            raw=Model().new_state(dt=.1)
            class Shared:
                backend='bad-shared'
                def new_state(self,**kw):
                    if same_object:return raw
                    state=Model().new_state(**kw);state.current=raw.current
                    step=state.step
                    def mutate(A):
                        result=step(A);raw.current[:]=result;return result
                    state.step=mutate;return state
            with self.assertRaises(InputError):self.api.check_material(Shared(),self.drive)

    def test_wrong_clock_timing_nonfinite_and_mutation_rejected(self):
        for mode in ('clock','timing','nan','shape','mutation','cached-current'):
            class Bad:
                backend='bad-'+mode
                def new_state(self,**kw):
                    state=Model().new_state(**kw);step=state.step
                    if mode=='timing':state.timing=TimeSpec(kw['dt'],(InputPort('A'),),(OutputPort('j',offset='1/2'),))
                    def wrong(A):
                        result=step(A)
                        if mode=='clock':state.time+=.01
                        if mode=='nan':return np.full(3,np.nan)
                        if mode=='shape':return np.zeros(2)
                        if mode=='mutation':A[:]=999
                        if mode=='cached-current':state.current=np.zeros(3)
                        return result
                    state.step=wrong;return state
            with self.subTest(mode=mode),self.assertRaises(InputError):self.api.check_material(Bad(),self.drive)

    def test_creating_another_state_must_not_change_first(self):
        class Intrusive:
            backend='constructor-side-effect';previous=None
            def new_state(self,**kw):
                state=Model().new_state(**kw)
                if self.previous is not None:
                    state.current+=1;self.previous.current[:]=state.current
                self.previous=state;return state
        with self.assertRaises(InputError):self.api.check_material(Intrusive(),self.drive)

    def test_cleanup_on_failure_and_partial_factory_failure(self):
        states=[]
        class Failing:
            backend='factory-failure'
            def new_state(self,**kw):
                if states:raise RuntimeError('cannot create second state')
                state=Model().new_state(**kw);states.append(state);return state
        with self.assertRaises(RuntimeError):self.api.check_material(Failing(),self.drive)
        self.assertTrue(states[0].closed)

    def test_initial_input_mutation_rejected_and_closed(self):
        states=[]
        class Mutating:
            backend='mutating-initial'
            def new_state(self,**kwargs):
                state=Model().new_state(**kwargs);states.append(state)
                kwargs['initial_A'][:]=999
                return state
        with self.assertRaises(InputError):self.api.check_material(Mutating(),self.drive)
        self.assertTrue(states[0].closed)

    def test_cleanup_failure_does_not_hide_invalid_contract(self):
        class Invalid:
            backend='invalid-timing'
            def new_state(self,**kwargs):
                state=Model().new_state(**kwargs);state.timing=None
                def close():raise RuntimeError('cleanup also failed')
                state.close=close;return state
        with self.assertRaises(InputError):self.api.check_material(Invalid(),self.drive)

    def test_invalid_tolerances_rejected_before_creation(self):
        for value in (None,'wrong',complex(1),-1,float('nan'),True):
            with self.subTest(value=value),self.assertRaises(InputError):
                self.api.check_material(Model(),self.drive,rtol=value)

    def test_nondeterministic_initial_values_rejected(self):
        class Unstable:
            backend='unstable';counter=0
            def new_state(self,**kw):
                self.counter+=1;state=Model().new_state(**kw);state.current+=self.counter
                return state
        with self.assertRaises(InputError):self.api.check_material(Unstable(),self.drive)

    def test_module_import_does_not_initialize_mpi(self):
        import subprocess,sys
        result=subprocess.run([sys.executable,'-c',
            "import sys;import salmonpy.materials;assert 'mpi4py.MPI' not in sys.modules"],capture_output=True,text=True)
        self.assertEqual(result.returncode,0,result.stderr)

    def test_duplicate_live_state_is_not_hidden_by_wrapper(self):
        state=Model().new_state(dt=.1)
        class Shared:
            backend='shared'
            def new_state(self,**kw):return state
        binding=self.api.bind_material(Shared(),key='shared',identity='bad')
        prepared=binding.provider.prepare(None)
        one=binding.backend.create_state(prepared,1,.1,None);self.addCleanup(one.close)
        with self.assertRaises(InputError):binding.backend.create_state(prepared,2,.1,None)
        self.assertFalse(state.closed)

if __name__=='__main__':unittest.main()
