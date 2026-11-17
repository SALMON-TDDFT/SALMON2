import unittest
from types import SimpleNamespace
import numpy as np
from salmonpy.coupled import NativeEHField
from salmonpy.errors import NumericalError,InputError


class Function:
    def __init__(self,fn):self.fn=fn
    def __call__(self,*args):return self.fn(*args)


class FieldRecoveryTests(unittest.TestCase):
    def test_declared_failure_preserves_field_and_permits_retry(self):
        comm=SimpleNamespace(allgather=lambda x:[x])
        f=NativeEHField.__new__(NativeEHField)
        f.context=SimpleNamespace(solver_comm=comm)
        f.closed=f.failed=False;f.step_index=1;f.capacity=3;f.dt=.01
        f.A=np.zeros((1,3));f.E=np.zeros((1,3));f.point_ids=(1,);f.manifest={}
        status=[20]
        def step(n,J,A,E):E[:]=0 if status[0] else 2;return status[0]
        f.lib=SimpleNamespace(sp_eh_step_checked=Function(step),sp_eh_recovery_protocol=Function(lambda:1))
        f._call=lambda fn,*args:fn(*args)
        with self.assertRaises(NumericalError):f.advance(np.ones((1,3)))
        self.assertFalse(f.failed);self.assertEqual(f.step_index,1)
        np.testing.assert_array_equal(f.A,0);np.testing.assert_array_equal(f.E,0)
        status[0]=0
        np.testing.assert_array_equal(f.advance(np.ones((1,3))),np.full((1,3),-.02))
        self.assertEqual(f.step_index,2)

    def test_unknown_protocol_error_remains_conservative(self):
        comm=SimpleNamespace(allgather=lambda x:[x])
        f=NativeEHField.__new__(NativeEHField)
        f.context=SimpleNamespace(solver_comm=comm)
        f.closed=f.failed=False;f.step_index=1;f.capacity=3;f.dt=.01
        f.A=np.zeros((1,3));f.E=np.zeros((1,3));f.point_ids=(1,);f.manifest={}
        def step(n,J,A,E):E[:]=0;return 99
        f.lib=SimpleNamespace(sp_eh_step_checked=Function(step),sp_eh_recovery_protocol=Function(lambda:1))
        f._call=lambda fn,*args:fn(*args)
        with self.assertRaises(InputError):f.advance(np.ones((1,3)))
        self.assertTrue(f.failed);self.assertEqual(f.step_index,1)
