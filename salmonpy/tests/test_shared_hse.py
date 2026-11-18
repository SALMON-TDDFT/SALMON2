import tempfile
from pathlib import Path
from types import SimpleNamespace
import unittest

from salmonpy.native_steps import _AppendJournal


class SharedHSETests(unittest.TestCase):
    def test_shared_model_entry_exists(self):
        from salmonpy.fast_hse_mpi import NativeFastHSEBackend, NativeFastHSEModel
        self.assertTrue(callable(NativeFastHSEBackend.open_model))
        self.assertTrue(callable(NativeFastHSEModel.new_state))

    def test_journal_does_not_truncate_other_state(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder)
            (root/'a.data').write_bytes(b'a')
            (root/'b.data').write_bytes(b'b')
            journal=_AppendJournal.capture(root,files=('a.data',))
            (root/'a.data').write_bytes(b'a-tail')
            (root/'b.data').write_bytes(b'b-tail')
            (root/'new.data').write_bytes(b'new')
            journal.restore()
            self.assertEqual((root/'a.data').read_bytes(),b'a')
            self.assertEqual((root/'b.data').read_bytes(),b'b-tail')
            self.assertEqual((root/'new.data').read_bytes(),b'new')

    def test_model_close_keeps_live_leases(self):
        from salmonpy._shared_hse import NativeFastHSEModel,_Lease
        from salmonpy.errors import InputError
        calls=[]
        comm=SimpleNamespace(allgather=lambda x:[x])
        native=SimpleNamespace(closed=False,context=SimpleNamespace(solver_comm=comm),_live=lambda:None)
        def detach():
            calls.append('detach');native.closed=True
        native.close=detach
        model=NativeFastHSEModel(native,None,Path('.'),.01,8,'hse_taylor4','identity',())
        model._attempted=True
        a,b=object(),object();model._states={a,b}
        first,second=_Lease(model),_Lease(model);first.state=a;second.state=b
        model.close()
        self.assertEqual(calls,[])
        first._live();first.close()
        self.assertEqual(calls,[])
        second._live()
        with self.assertRaises(InputError):model.new_state()
        second.close();second.close();model.close()
        self.assertEqual(calls,['detach'])

    def test_handle_binding_does_not_overwrite_shared_signature(self):
        from salmonpy._shared_hse import _HandleLibrary
        calls=[]
        def fn(*args):calls.append(args);return 0
        fn.argtypes=['actual ABI']
        lib=SimpleNamespace(**{'sp_fast_hse_handle_'+key:fn for key in
                              ('step','checkpoint','info','close','history','transaction')})
        a,b=_HandleLibrary(lib,21),_HandleLibrary(lib,22)
        a.sp_fast_hse_step.argtypes=['compatibility signature']
        a.sp_fast_hse_step(1);b.sp_fast_hse_step(2)
        self.assertEqual(calls,[(21,1),(22,2)])
        self.assertEqual(fn.argtypes,['actual ABI'])

    def test_empty_model_release_does_not_reset_a_later_owner(self):
        from unittest.mock import patch
        from salmonpy import fast_hse_mpi
        from salmonpy._shared_hse import NativeFastHSEModel
        native=SimpleNamespace(closed=False,_live=lambda:None)
        def detach():native.closed=True
        native.close=detach
        model=NativeFastHSEModel(native,None,Path('.'),.01,8,'hse_taylor4','identity',())
        with patch.object(fast_hse_mpi,'_used',True),patch.object(fast_hse_mpi.embedded,'_used',True):
            model.close()
            self.assertFalse(fast_hse_mpi._used)
            fast_hse_mpi._used=fast_hse_mpi.embedded._used=True
            model.close()
            self.assertTrue(fast_hse_mpi._used)
            self.assertTrue(fast_hse_mpi.embedded._used)

    def test_session_factory_attempts_all_returned_state_cleanup(self):
        from salmonpy._shared_hse import _SessionGroup
        from salmonpy.errors import RuntimeFailure
        events=[]
        def fail():events.append('first');raise RuntimeError('close failed')
        group=_SessionGroup(SimpleNamespace(close=lambda:events.append('model')))
        group.states=[SimpleNamespace(close=fail),SimpleNamespace(close=lambda:events.append('unadmitted'))]
        with self.assertRaisesRegex(RuntimeFailure,'close failed'):group.close()
        self.assertEqual(events,['first','unadmitted','model'])

    def test_shared_symmetry_requires_new_capability_but_full_k_stays_compatible(self):
        from salmonpy._shared_hse import _configure_handles
        from salmonpy.errors import InputError
        class Fn:
            def __init__(self,value=0):self.value=value
            def __call__(self,*args):return self.value
        names=('init','step','info','history','transaction','checkpoint','close','label')
        lib=SimpleNamespace(**{'sp_fast_hse_handle_'+name:Fn() for name in names},
            sp_fast_hse_handle_protocol=Fn(1),sp_fast_hse_model_protocol=Fn(1),sp_fast_hse_model_new=Fn())
        _configure_handles(lib)
        with self.assertRaisesRegex(InputError,'shared symmetry HSE'):_configure_handles(lib,symmetry=True)
        lib.sp_fast_hse_shared_symmetry_protocol=Fn(0)
        with self.assertRaisesRegex(InputError,'unsupported'):_configure_handles(lib,symmetry=True)
        lib.sp_fast_hse_shared_symmetry_protocol=Fn(1)
        _configure_handles(lib,symmetry=True)
