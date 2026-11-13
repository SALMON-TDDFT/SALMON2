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
