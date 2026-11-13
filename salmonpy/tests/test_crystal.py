import unittest
import numpy as np
from salmonpy.fields import FieldSamples
from salmonpy.crystal import respond

class TestInitialState:
    """Interface test double only: no TDDFT/SBE physics is simulated."""
    backend='test'
    def new_state(self,*,dt,time0,initial_A):
        class State:
            def __init__(self): self.time=time0; self.current=np.array([1.,2.,3.]); self.last=initial_A.copy()
            def step(self,A):
                self.current=self.current+(A-self.last); self.last=A.copy(); self.time+=dt
                return self.current.copy()
        return State()

class CrystalTests(unittest.TestCase):
    def test_A_only_and_reusable_initial_state(self):
        t=np.arange(5)*.1+2.;a=np.arange(15).reshape(5,3)*.01
        drive=FieldSamples(t,a,'A');initial=TestInitialState()
        r=respond(initial,drive)
        np.testing.assert_allclose(r.current,[1,2,3]+a-a[0])
        np.testing.assert_array_equal(respond(initial,drive).current,r.current)
        with self.assertRaises(ValueError): respond(initial,FieldSamples(t,a,'E'))

    def test_nonuniform_rejected_before_initialization(self):
        class Initial:
            def new_state(self,**kwargs): raise AssertionError('must not initialize')
        with self.assertRaises(ValueError): respond(Initial(),FieldSamples([0,1,3],np.zeros((3,3)),'A'))

    def test_accumulated_time_roundoff(self):
        t=1e6+np.arange(1025)*.1
        result=respond(TestInitialState(),FieldSamples(t,np.zeros((len(t),3)),'A'))
        self.assertEqual(len(result.time),len(t))

    def test_wrong_backend_timestep_rejected(self):
        class Wrong(TestInitialState):
            def new_state(self,**kwargs):
                kwargs['dt']*=1.01
                return super().new_state(**kwargs)
        with self.assertRaisesRegex(ValueError,'backend time'):
            respond(Wrong(),FieldSamples(np.arange(4)*.1,np.zeros((4,3)),'A'))

    def test_batch_closes_owned_backend_on_success_and_failure(self):
        closed=[]
        class Initial:
            backend='lifecycle-test'
            def new_state(self,*,dt,time0,initial_A):
                class State:
                    time=time0
                    current=np.zeros(3)
                    def step(self,A):
                        if A[0]<0: raise ValueError('test failure')
                        self.time+=dt
                        return self.current
                    def close(self): closed.append(True)
                return State()
        respond(Initial(),FieldSamples([0.,.1],np.zeros((2,3)),'A'))
        self.assertEqual(len(closed),1)
        with self.assertRaisesRegex(ValueError,'test failure'):
            respond(Initial(),FieldSamples([0.,.1],[[0,0,0],[-1,0,0]],'A'))
        self.assertEqual(len(closed),2)
