"""Native adapter admission/state tests with no MPI initialization or compilation."""
import ctypes as c
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np

from salmonpy.errors import InputError,RuntimeFailure
from salmonpy.fast_hse_mpi import NativeFastHSEBackend,_State


class FakeComm:
    size=1
    rank=0
    def allgather(self,value):return [value]


class NativeFastHSETests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)
        self.library=self.root/'native';self.library.write_bytes(b'not-loaded-by-constructor')
        self.record=self.root/'build.json'
        self.data={'schema':'salmonpy.fast-hse-mpi/1','status':'built-and-audited',
            'library':str(self.library),'sha256':hashlib.sha256(self.library.read_bytes()).hexdigest(),
            'origin_sha256':'a'*64}
        self.record.write_text(json.dumps(self.data))

    def backend(self):return NativeFastHSEBackend(self.record,self.root)

    def state(self):
        comm=FakeComm()
        native=SimpleNamespace(context=SimpleNamespace(solver_comm=comm),_live=lambda:None,close=lambda:None)
        calls=[]
        def step(old,new,current,energy):
            calls.append((old.copy(),new.copy()));current[:]=[1,2,3]
            c.cast(energy,c.POINTER(c.c_double)).contents.value=7
            return 0
        def checkpoint(path,restore):return 0
        def info(step,A,current):
            c.cast(step,c.POINTER(c.c_int)).contents.value=2
            A[:]=[.2,.3,.4];current[:]=[4,5,6];return 0
        lib=SimpleNamespace(sp_fast_hse_step=step,sp_fast_hse_checkpoint=checkpoint,
                            sp_fast_hse_info=info,sp_fast_hse_close=lambda:0)
        state=_State(native,lib,self.root,.01,3,'hse_taylor4')
        state.initialized=True
        return state,calls

    def test_symmetry_opt_in_is_lazy_bounded_and_identity_bound(self):
        with patch('salmonpy.fast_hse_mpi.c.CDLL', side_effect=AssertionError('lazy')):
            ordinary = self.backend()
            symmetric = NativeFastHSEBackend(self.record, self.root, allow_symmetry=True)
        self.assertEqual(ordinary.supported_group_sizes, (1, 2))
        self.assertEqual(symmetric.supported_group_sizes, (1,2))
        self.assertNotEqual(ordinary.identity, symmetric.identity)
        for value in (1, 'yes', None):
            with self.assertRaises(InputError):
                NativeFastHSEBackend(self.record, self.root, allow_symmetry=value)

    def test_field_rejection_requires_specific_protocol_and_preserves_state(self):
        for protocol in (None, 0, 1):
            state, _ = self.state()
            code = [34]
            def reject(old, new, current, energy):
                current[:] = 0
                c.cast(energy, c.POINTER(c.c_double)).contents.value = 0
                return code[0]
            state.lib.sp_fast_hse_step = reject
            if protocol is not None:
                state.lib.sp_fast_hse_symmetry_protocol = lambda: protocol
            with self.subTest(protocol=protocol), self.assertRaises(InputError if protocol == 1 else RuntimeFailure):
                state.step([1e-250, 0, 0])
            self.assertEqual(state.failed, protocol != 1)
            self.assertEqual(state.steps, 0)
            np.testing.assert_array_equal(state.A, np.zeros(3))
            if protocol == 1:
                code[0] = 0
                state.step([0, 0, .1])
                self.assertEqual(state.steps, 1)
            state.close()

    def test_symmetry_entry_requires_capability_and_preserves_legacy_selection(self):
        from salmonpy.fast_hse_mpi import _initializer
        old = lambda *a: 0
        old.argtypes = []
        new = lambda *a: 0
        legacy = SimpleNamespace(sp_fast_hse_init=old)
        self.assertIs(_initializer(legacy,False),old)
        with self.assertRaises(InputError):
            _initializer(legacy,True)
        for protocol in (0,1,2):
            lib = SimpleNamespace(sp_fast_hse_init=old,sp_fast_hse_init_symmetry=new,
                                  sp_fast_hse_symmetry_protocol=lambda:protocol)
            if protocol == 1:
                self.assertIs(_initializer(lib,True),new)
            else:
                with self.assertRaises(InputError):_initializer(lib,True)

    def test_symmetry_kmpi_requires_separate_native_capability(self):
        from salmonpy.fast_hse_mpi import _initializer
        old = lambda *a: 0
        old.argtypes = []
        new = lambda *a: 0
        lib = SimpleNamespace(sp_fast_hse_init=old,sp_fast_hse_init_symmetry=new,
                              sp_fast_hse_symmetry_protocol=lambda:1)
        with self.assertRaisesRegex(InputError,'k-MPI'):
            _initializer(lib,True,2)
        for value in (0,1,2):
            lib.sp_fast_hse_symmetry_kmpi_protocol = lambda:value
            if value == 1:
                self.assertIs(_initializer(lib,True,2),new)
            else:
                with self.assertRaisesRegex(InputError,'k-MPI'):_initializer(lib,True,2)

    def test_missing_symmetry_asset_rejected_before_native_lifetime(self):
        from salmonpy.fast_hse import FastHSEMaterialProvider
        bundle = self.root/'bundle';bundle.mkdir()
        restart = bundle/'data_for_restart';restart.mkdir()
        for name in ('info.bin','occupation.bin','wfn.bin'):
            (restart/name).write_bytes(b'synthetic')
        (bundle/'rt.inp').write_text("&calculation theory='tddft_response' /\n"
            "&system yn_symmetry='yyn' /\n&functional xc='hse06',exx_factor_mode='off' /\n"
            "&propagation propagator='hse_taylor4' /\n")
        provider = FastHSEMaterialProvider(bundle,max_steps=3)
        backend = NativeFastHSEBackend(self.record,self.root,allow_symmetry=True)
        ctx = SimpleNamespace(solver_comm=FakeComm(),point_ids=(1,))
        with patch('salmonpy.fast_hse_mpi.c.CDLL',side_effect=AssertionError('no native load')):
            with self.assertRaisesRegex(InputError,'sym.dat'):
                backend.create_state(provider.prepare(ctx),1,.01,ctx)

    def test_malformed_symmetry_rejected_without_native_or_work_copy(self):
        from salmonpy.fast_hse import FastHSEMaterialProvider
        bundle = self.root/'bundle';bundle.mkdir()
        restart = bundle/'data_for_restart';restart.mkdir()
        for name in ('info.bin','occupation.bin','wfn.bin'):
            (restart/name).write_bytes(b'synthetic')
        (bundle/'rt.inp').write_text("&calculation theory='tddft_response' /\n"
            "&system yn_symmetry='yyn' /\n&functional xc='hse06',exx_factor_mode='off' /\n"
            "&propagation propagator='hse_taylor4' /\n")
        (bundle/'sym.dat').write_text('NaN 0 0 0\n0 1 0 0\n0 0 1 0\n')
        provider = FastHSEMaterialProvider(bundle,max_steps=3)
        backend = NativeFastHSEBackend(self.record,self.root,allow_symmetry=True)
        ctx = SimpleNamespace(solver_comm=FakeComm(),point_ids=(1,))
        before=set(self.root.iterdir())
        with patch('salmonpy.fast_hse_mpi.c.CDLL',side_effect=AssertionError('no native load')):
            with self.assertRaisesRegex(InputError,'sym.dat'):
                backend.create_state(provider.prepare(ctx),1,.01,ctx)
        self.assertEqual(set(self.root.iterdir()),before)

        # The generic scalar reader lowercases values; flags must instead
        # retain their spelling so unsupported native forms reject unchanged.
        inp=bundle/'rt.inp';inp.write_text(inp.read_text().replace("'yyn'","'YYN'"))
        provider=FastHSEMaterialProvider(bundle,max_steps=3)
        with patch('salmonpy.fast_hse_mpi.c.CDLL',side_effect=AssertionError('no native load')):
            with self.assertRaisesRegex(InputError,'yn_symmetry'):
                backend.create_state(provider.prepare(ctx),1,.01,ctx)
        self.assertEqual(set(self.root.iterdir()),before)

    def test_constructor_is_lazy_and_hash_bound(self):
        with patch('salmonpy.fast_hse_mpi.c.CDLL',side_effect=AssertionError('must not load')):
            backend=self.backend()
        self.assertEqual((backend.execution,backend.reduction),('collective','leader'))
        self.library.write_bytes(b'changed')
        with self.assertRaises(InputError):self.backend()

    def test_schema_origin_and_work_root_admission(self):
        for key,value in [('schema','salmonpy.fast-hse-worker/1'),('status','failed'),('origin_sha256','x')]:
            record={**self.data,key:value};self.record.write_text(json.dumps(record))
            with self.assertRaises(InputError):self.backend()
        self.record.write_text(json.dumps(self.data))
        with self.assertRaises(InputError):NativeFastHSEBackend(self.record,self.root/'missing')

    def test_k_only_group_and_one_lifetime_guard_before_loading(self):
        backend=self.backend();comm=FakeComm();comm.size=3
        ctx=SimpleNamespace(solver_comm=comm,point_ids=(1,))
        with patch('salmonpy.fast_hse_mpi.c.CDLL',side_effect=AssertionError('must not load')):
            with self.assertRaisesRegex(InputError,'one or two ranks'):
                backend.create_state(object(),1,.01,ctx)
        comm.size=1
        with patch('salmonpy.fast_hse_mpi._used',True):
            with self.assertRaisesRegex(InputError,'one native lifetime'):
                backend.create_state(object(),1,.01,ctx)

    def test_endpoint_sign_clock_and_capacity(self):
        state,calls=self.state();before=Path.cwd()
        old=np.zeros(3);new=np.array([.1,.2,.3])
        np.testing.assert_array_equal(state.advance(old,new),[1,2,3])
        self.assertEqual(state.time,.01);self.assertEqual(state.energy,7)
        new[:]=9;np.testing.assert_array_equal(state.A,[.1,.2,.3])
        with self.assertRaises(InputError):state.advance(old,[.2,.2,.2])
        self.assertEqual(len(calls),1)
        state.step([0,0,0]);state.step([0,0,0])
        with self.assertRaises(InputError):state.step([0,0,0])
        self.assertEqual(len(calls),3);self.assertEqual(Path.cwd(),before)

    def test_checkpoint_info_restores_A_current_and_clock_without_steps(self):
        state,calls=self.state();path=self.root/'saved.bin';path.write_bytes(b'caller-verified')
        state.checkpoint(path,restore=True)
        self.assertEqual(state.steps,2);self.assertEqual(state.time,.02)
        np.testing.assert_array_equal(state.A,[.2,.3,.4]);np.testing.assert_array_equal(state.current,[4,5,6])
        self.assertEqual(calls,[])
        with self.assertRaises(InputError):state.checkpoint(path,restore=True)
        state.close();state.close()

    def test_native_failure_invalidates_and_cwd_is_restored(self):
        state,calls=self.state();state.lib.sp_fast_hse_step=lambda *a:1
        before=Path.cwd()
        with self.assertRaises(RuntimeFailure):state.step([.1,.2,.3])
        self.assertTrue(state.failed);self.assertEqual(Path.cwd(),before)
        with self.assertRaises(InputError):state.step([.1,.2,.3])
        state.close()

    def test_endpoint_disagreement_rejected_before_native_call(self):
        state,calls=self.state()
        def allgather(value):
            if isinstance(value,tuple):return [value,([9,0,0],[1,2,3],0)]
            return [value,value]
        state._native.context.solver_comm.allgather=allgather
        with self.assertRaisesRegex(InputError,'endpoints differ'):state.step([.1,.2,.3])
        self.assertEqual(calls,[])

    def test_complex_endpoints_rejected_without_mutation(self):
        state,calls=self.state()
        for old,new in ((np.zeros(3,dtype=complex),np.ones(3)),
                        (np.zeros(3),np.array([1j,0,0]))):
            with self.assertRaises(InputError):state.advance(old,new)
        self.assertEqual(calls,[])
        self.assertEqual(state.steps,0)
        np.testing.assert_array_equal(state.A,np.zeros(3))

    def test_checkpoint_validation_does_not_call_native(self):
        state,calls=self.state()
        def forbidden(*args):raise AssertionError('invalid checkpoint reached native')
        state.lib.sp_fast_hse_checkpoint=forbidden
        existing=self.root/'existing';existing.write_bytes(b'saved')
        for path,restore in ((existing,False),(self.root/'missing',True),
                             (self.root/'new',1),(self.root/'bad\nname',False)):
            with self.assertRaises(InputError):state.checkpoint(path,restore)
        self.assertFalse(state.failed)
        self.assertEqual(state.steps,0)

    def test_restore_invalid_metadata_invalidates_state(self):
        state,calls=self.state();path=self.root/'saved.bin';path.write_bytes(b'saved')
        def bad_info(step,A,current):
            c.cast(step,c.POINTER(c.c_int)).contents.value=4
            A[:]=0;current[:]=0;return 0
        state.lib.sp_fast_hse_info=bad_info
        with self.assertRaisesRegex(RuntimeFailure,'metadata invalid'):
            state.checkpoint(path,True)
        self.assertTrue(state.failed)
        self.assertEqual(calls,[])
        state.close()

    def test_checkpoint_native_failure_invalidates_state(self):
        state,calls=self.state()
        state.lib.sp_fast_hse_checkpoint=lambda *args:1
        with self.assertRaisesRegex(RuntimeFailure,'checkpoint failed'):
            state.checkpoint(self.root/'new',False)
        self.assertTrue(state.failed)
        state.close()

    def test_declared_deep_recovery_preserves_state_and_allows_retry(self):
        from salmonpy.errors import NumericalError
        for protocol in (0,1):
            state,_=self.state()
            code=[30]
            def reject(old,new,current,energy):
                current[:]=0
                c.cast(energy,c.POINTER(c.c_double)).contents.value=0
                return code[0]
            state.lib.sp_fast_hse_step=reject
            state.lib.sp_fast_hse_recovery_protocol=lambda:protocol
            with self.subTest(protocol=protocol),self.assertRaises(NumericalError if protocol else RuntimeFailure):
                state.step([.1,0,0])
            np.testing.assert_array_equal(state.A,np.zeros(3))
            self.assertEqual(state.steps,0)
            self.assertEqual(state.failed,not bool(protocol))
            if protocol:
                code[0]=0
                state.step([.1,0,0])
                self.assertEqual(state.steps,1)
