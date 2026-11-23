import unittest
from types import SimpleNamespace
import numpy as np
from salmonpy.coupled import NativeEHField
from salmonpy.errors import NumericalError,InputError,RuntimeFailure


class Function:
    def __init__(self,fn):self.fn=fn
    def __call__(self,*args):return self.fn(*args)


class FieldRecoveryTests(unittest.TestCase):
    def test_declared_failure_preserves_field_and_permits_retry(self):
        comm=SimpleNamespace(allgather=lambda x:[x])
        f=NativeEHField.__new__(NativeEHField);f._execution_mode="checked"
        f.context=SimpleNamespace(solver_comm=comm)
        f.closed=f.failed=False;f.step_index=1;f.capacity=3;f.dt=.01
        f.A=np.zeros((1,3));f.E=np.zeros((1,3));f.point_ids=(1,);f.manifest={}
        status=[20]
        def step(n,J,A,E):E[:]=0 if status[0] else 2;return status[0]
        f.lib=SimpleNamespace(sp_eh_step_checked=Function(step),sp_eh_recovery_protocol=Function(lambda:1))
        f._call=lambda fn,*args:fn(*args)
        with self.assertRaises(NumericalError):f.advance(np.ones((1,3)))
        status[0]=24
        with self.assertRaisesRegex(RuntimeFailure,'allocation'):f.advance(np.ones((1,3)))
        self.assertFalse(f.failed);self.assertEqual(f.step_index,1)
        np.testing.assert_array_equal(f.A,0);np.testing.assert_array_equal(f.E,0)
        status[0]=0
        np.testing.assert_array_equal(f.advance(np.ones((1,3))),np.full((1,3),-.02))
        self.assertEqual(f.step_index,2)

    def test_runtime_and_step_protocol_admission(self):
        from salmonpy.coupled import _require_field_record
        _require_field_record({'schema':'salmonpy.embedded-tddft/1'})
        _require_field_record({'schema':'salmonpy.fast-hse-mpi/1',
            'field_abi':{'protocol':1,'scheme':NativeEHField.scheme}})
        for declaration in (None,{}, {'protocol':True,'scheme':NativeEHField.scheme},
                            {'protocol':2,'scheme':NativeEHField.scheme}, {'protocol':1,'scheme':'other'}):
            with self.subTest(declaration=declaration),self.assertRaises(InputError):
                _require_field_record({'schema':'salmonpy.fast-hse-mpi/1','field_abi':declaration})
        with self.assertRaises(InputError):_require_field_record({'schema':'unknown'})
        comm=SimpleNamespace(allgather=lambda x:[x])
        f=NativeEHField.__new__(NativeEHField);f._execution_mode="checked"
        f.context=SimpleNamespace(solver_comm=comm)
        f.closed=f.failed=False;f.step_index=1;f.capacity=3;f.dt=.01
        f.A=np.zeros((1,3));f.E=np.zeros((1,3));f.point_ids=(1,);f.manifest={}
        def step(n,J,A,E):E[:]=0;return 99
        f.lib=SimpleNamespace(sp_eh_step_checked=Function(step),sp_eh_recovery_protocol=Function(lambda:1))
        f._call=lambda fn,*args:fn(*args)
        from pathlib import Path
        from unittest.mock import patch
        events=[]
        f.native=SimpleNamespace(_live=lambda:None);f.directory=Path.cwd()
        f.lib.sp_io_directory_protocol=Function(lambda:3)
        f.lib.sp_io_directory_enter=Function(lambda path:0)
        f.lib.sp_io_directory_leave=Function(lambda:events.append('leave'))
        with patch('os.chdir',side_effect=AssertionError('field host chdir forbidden')):
            self.assertEqual(NativeEHField._call(f,lambda:9),9)
        self.assertEqual(events,['leave'])
        f.lib.sp_io_directory_protocol=Function(lambda:2)
        with patch('os.chdir',side_effect=AssertionError('host chdir forbidden')):
            with self.assertRaises(InputError):NativeEHField._call(f,lambda:9)
        self.assertEqual(events,['leave'])
        with self.assertRaises(InputError):f.advance(np.ones((1,3)))
        self.assertTrue(f.failed);self.assertEqual(f.step_index,1)
