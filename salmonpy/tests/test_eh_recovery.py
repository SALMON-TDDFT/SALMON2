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


class ComponentCurrentIndexTests(unittest.TestCase):
    def test_local_lists_and_operation_contract(self):
        from types import SimpleNamespace
        from salmonpy.coupled import NativeEHField
        from salmonpy.errors import InputError
        from salmonpy.api import _field_operations
        field=object.__new__(NativeEHField)
        field.calculation=SimpleNamespace(geometry=SimpleNamespace(grid=SimpleNamespace(counts=(40,40,40))))
        empty=field._edge_indices(np.empty((0,3),dtype=int))
        self.assertEqual(empty.shape,(0,3))
        np.testing.assert_array_equal(field._edge_indices([[1,2,3]]),[[1,2,3]])
        for bad in ([[1.,2.,3.]],[[40,0,0]],[[-1,0,0]],[[True,False,True]]):
            with self.assertRaises(InputError):field._edge_indices(bad)
        operations=_field_operations()
        self.assertIn('rank-local',operations['sp_eh_step_edges']['routing'])
        self.assertEqual(operations['sp_eh_edges_read']['collective'],'local')


class ComponentEHQueryTests(unittest.TestCase):
    def test_geometry_and_registered_read_only_operation(self):
        from salmonpy.api import _field_operations
        self.assertEqual(_field_operations()['sp_eh_fields_read']['collective'],'local')
        self.assertEqual(_field_operations()['sp_eh_fields_read']['arguments'][-1],'array')
        f=NativeEHField.__new__(NativeEHField)
        f.electric_edge_geometry=lambda:dict(counts=(40,40,40),origin=np.full((3,3),-39.)+np.eye(3),spacing=np.full(3,2.))
        geometry=f.component_field_geometry()
        np.testing.assert_array_equal(geometry['origin'][0],np.full((3,3),-39.)+np.eye(3))
        np.testing.assert_array_equal(geometry['origin'][1],np.full((3,3),-38.)-np.eye(3))

        f.calculation=SimpleNamespace(geometry=SimpleNamespace(grid=SimpleNamespace(counts=(40,40,40))))
        f.context=SimpleNamespace(solver_comm=SimpleNamespace(allgather=lambda value:[value],Allreduce=lambda local,out:np.copyto(out,local)))
        f.closed=f.failed=False
        def read(n,indices,values):
            self.assertEqual(values.shape,(2,n,3));values[0]=2.;values[1]=3.;return 0
        f.lib=SimpleNamespace(sp_eh_fields_read=Function(read));f._call=lambda fn,*args:fn(*args)
        result=f.component_fields([[2,3,4]])
        np.testing.assert_array_equal(result,np.array([[[2.,2.,2.]],[[3.,3.,3.]]]))
