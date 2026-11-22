"""Standard TDDFT/HSE shared-model ownership; numerical parity lives in MPI probes."""
import tempfile
from pathlib import Path
from types import SimpleNamespace
import unittest
from unittest.mock import patch
from salmonpy.errors import InputError
from salmonpy.native_steps import _AppendJournal


class SharedTDDFTTests(unittest.TestCase):
    def test_public_entry_and_legacy_factory_opt_in(self):
        from salmonpy.embedded_tddft import NativeTDDFTBackend,NativeTDDFTModel
        self.assertTrue(callable(NativeTDDFTModel.new_state))
        backend=object.__new__(NativeTDDFTBackend);backend.shared_model=False
        self.assertIsNone(backend.open_group(None,None,None))
        backend.shared_model=True
        with self.assertRaisesRegex(InputError,'open_group'):
            backend.create_state(None,None,None,None)

    def model(self):
        from salmonpy._shared_tddft import NativeTDDFTModel
        events=[]
        native=SimpleNamespace(context=SimpleNamespace(solver_comm=SimpleNamespace(allgather=lambda x:[x])),
                               closed=False,_live=lambda:None)
        def close():events.append('detach');native.closed=True
        native.close=close
        model=NativeTDDFTModel(native,None,Path('.'),SimpleNamespace(max_steps=32,manifest='manifest'),
                              .003,(1,1,1,1,1),'id',Path('.'))
        return model,events

    def test_model_close_keeps_siblings_live(self):
        from salmonpy._shared_tddft import _Lease
        model,events=self.model();model._attempted=True
        a,b=object(),object();model._states={a,b}
        first,second=_Lease(model),_Lease(model);first.state=a;second.state=b
        model.close();first._live();first.close();self.assertEqual(events,[])
        second._live()
        with self.assertRaises(InputError):model.new_state()
        second.close();second.close();model.close();self.assertEqual(events,['detach'])

    def test_common_parameters_are_read_only_and_empty_close_is_idempotent(self):
        from salmonpy import embedded_tddft
        model,events=self.model()
        for name,value in (('dt',1.),('capacity',5),('decomposition',(1,2,1,1,1)),('initial',None)):
            with self.assertRaises(AttributeError):setattr(model,name,value)
        with patch.object(embedded_tddft,'_used',True):
            model.close();self.assertFalse(embedded_tddft._used)
            embedded_tddft._used=True;model.close();self.assertTrue(embedded_tddft._used)
        self.assertEqual(events,['detach'])

    def test_bound_handle_cannot_overwrite_abi_signature(self):
        from salmonpy.api import native_contract, _bind_state_operations
        from salmonpy.errors import InputError
        for component in ('tddft','fast_hse'):
            for route in ('shared','registry'):
                operations=native_contract(component,route=route)['state_operations']
                self.assertTrue(any(name.endswith('handle_step') for name in operations))
                step_spec=next(v for k,v in operations.items() if k.endswith('handle_step'))
                self.assertEqual([p['name'] for p in step_spec['parameters']][:4],['state','old_A','new_A','j'])
                self.assertEqual(step_spec['parameters'][3]['direction'],'out')
                self.assertTrue(any(name.endswith('handle_checkpoint') for name in operations))
                class Fn:
                    def __call__(self,*args):raise AssertionError('binding must not execute state operations')
                functions={name:Fn() for name in operations}
                for missing in functions:
                    partial=SimpleNamespace(**{k:v for k,v in functions.items() if k!=missing})
                    with self.subTest(component=component,route=route,missing=missing),self.assertRaisesRegex(InputError,missing):
                        _bind_state_operations(partial,component,route=route)
                    self.assertFalse(any(hasattr(f,'argtypes') for f in functions.values()))
                _bind_state_operations(SimpleNamespace(**functions),component,route=route)
                import ctypes as c
                for name,fn in functions.items():self.assertIs(fn.restype,c.c_int,name)
                checkpoint=next(v for k,v in operations.items() if k.endswith('handle_checkpoint'))
                self.assertEqual(checkpoint['modes'],{0:'save',1:'fresh load',2:'stage',3:'commit',4:'cancel'})
                checkpoint['modes'][0]='edited'
                fresh=native_contract(component,route=route)['state_operations']
                self.assertEqual(next(v for k,v in fresh.items() if k.endswith('handle_checkpoint'))['modes'][0],'save')
        from salmonpy._shared_tddft import _Library
        calls=[]
        def fn(*args):calls.append(args);return 0
        fn.argtypes=['native ABI']
        lib=SimpleNamespace(**{'sp_tddft_handle_'+key:fn for key in ('step','close','checkpoint','transaction')})
        a,b=_Library(lib,31),_Library(lib,32)
        a.sp_tddft_step.argtypes=['compatibility ABI'];a.sp_tddft_step(1);b.sp_tddft_step(2)
        self.assertEqual(calls,[(31,1),(32,2)]);self.assertEqual(fn.argtypes,['native ABI'])

    def test_state_owned_timer_query_and_old_runtime(self):
        from salmonpy._shared_tddft import _SharedState
        from salmonpy.errors import RuntimeFailure
        import ctypes as c
        def probe():return 1
        def timers(token, address):
            self.assertEqual(token,17)
            values=c.cast(address,c.POINTER(c.c_double))
            for i in range(4):values[i]=float(i+1)
            return 0
        lib=SimpleNamespace(sp_tddft_timer_protocol=probe,sp_tddft_handle_timers=timers)
        state=object.__new__(_SharedState);state.closed=False;state.failed=False;state.handle=17
        state.native=SimpleNamespace(model=SimpleNamespace(lib=lib),
            context=SimpleNamespace(solver_comm=SimpleNamespace(allgather=lambda x:[x])))
        state._call=lambda fn,*args:fn(*args)
        self.assertEqual(state.timing_seconds,dict(total=1.,initialization=2.,density=3.,hartree=4.))
        del lib.sp_tddft_timer_protocol
        with self.assertRaisesRegex(RuntimeFailure,'timer'): _=state.timing_seconds
        from salmonpy.fast_hse_mpi import _State
        hse=object.__new__(_State);hse.closed=False;hse.failed=False;hse.handle=17
        hse._native=SimpleNamespace(context=state.native.context)
        hse.lib=SimpleNamespace(sp_fast_hse_timer_protocol=probe,sp_fast_hse_handle_timers=timers)
        hse._call=state._call
        self.assertEqual(hse.timing_seconds,dict(total=1.,initialization=2.,density=3.,hartree=4.))
        del hse.lib.sp_fast_hse_timer_protocol
        with self.assertRaisesRegex(RuntimeFailure,'timer'): _=hse.timing_seconds
        from salmonpy.embedded_tddft import _NativeState
        def legacy_timers(address):return timers(17,address)
        legacy=object.__new__(_NativeState);legacy.closed=False;legacy.failed=False
        legacy.native=state.native;legacy._call=state._call
        legacy.lib=SimpleNamespace(sp_tddft_legacy_timer_protocol=probe,sp_tddft_legacy_timers=legacy_timers)
        self.assertEqual(legacy.timing_seconds,dict(total=1.,initialization=2.,density=3.,hartree=4.))
        del legacy.lib.sp_tddft_legacy_timer_protocol
        with self.assertRaisesRegex(RuntimeFailure,'timer'): _=legacy.timing_seconds
        hse_legacy=object.__new__(_State);hse_legacy.closed=False;hse_legacy.failed=False
        hse_legacy._native=hse._native;hse_legacy._call=state._call
        hse_legacy.lib=SimpleNamespace(sp_fast_hse_legacy_timer_protocol=probe,sp_fast_hse_legacy_timers=legacy_timers)
        self.assertEqual(hse_legacy.timing_seconds,dict(total=1.,initialization=2.,density=3.,hartree=4.))
        del hse_legacy.lib.sp_fast_hse_legacy_timer_protocol
        with self.assertRaisesRegex(RuntimeFailure,'timer'): _=hse_legacy.timing_seconds

    def test_unsupported_output_and_unsafe_paths_reject_before_attach(self):
        from salmonpy._shared_tddft import _profile
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'rt.inp';initial=SimpleNamespace(directory=Path(d),input_file='rt.inp')
            for extra in ("yn_out_dns_rt='y'", "base_directory='../escape'", "trans_longi='lo'", "xc='hse'"):
                p.write_text("&calculation\nyn_periodic='y', "+extra+"\n/\n")
                with self.assertRaises(InputError):_profile(initial)
            p.write_text("&calculation\nyn_periodic='y', xc='libxc_pbe'\n/\n")
            self.assertEqual(_profile(initial),Path('.'))

    def test_restored_epoch_disagreement_is_not_published(self):
        from salmonpy._shared_tddft import _SharedState
        from salmonpy.errors import RuntimeFailure
        import numpy as np
        calls=[]
        def gather(value):
            calls.append(value)
            return [value,value] if value is None else [value,(9,(0.,0.,0.),(0.,0.,0.))]
        def info(token,step,A,j):step._obj.value=8;return 0
        model=SimpleNamespace(lib=SimpleNamespace(sp_tddft_handle_info=info))
        native=SimpleNamespace(model=model,context=SimpleNamespace(solver_comm=SimpleNamespace(allgather=gather)))
        state=object.__new__(_SharedState);state.native=native;state.handle=5;state.capacity=32
        state.failed=False;state.steps=0;state.A=np.zeros(3);state.current=np.zeros(3)
        state._call=lambda fn,*args:fn(*args)
        with self.assertRaisesRegex(RuntimeFailure,'endpoints differ'):state._refresh()
        self.assertTrue(state.failed);self.assertEqual(state.steps,0)

    def test_profile_ignores_coordinate_blocks_after_last_namelist(self):
        from salmonpy._shared_tddft import _profile
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);(root/'rt.inp').write_text("&system\nyn_periodic='y'\n/\n&emfield\ntrans_longi='tr'\n/\n&atomic_red_coor\n'Si' .3 .4 .45 1\n/\n")
            self.assertEqual(_profile(SimpleNamespace(directory=root,input_file='rt.inp')),Path('.'))

    def test_session_registry_shared_across_bindings_and_owned_by_context(self):
        from salmonpy.embedded_tddft import NativeTDDFTBackend
        events=[]
        registry=SimpleNamespace(backend=None,open_model=lambda initial,dt,**kw:SimpleNamespace(
            new_state=lambda p:events.append(('state',initial,p)),close=lambda:events.append(('model-close',initial))))
        backend=object.__new__(NativeTDDFTBackend)
        backend.registry_model=True;backend.shared_model=False
        backend.sha256='sha';backend.work_root=Path('/tmp/work');backend.decomposition=(1,1,1,1,1)
        registry.backend=backend;ctx=SimpleNamespace(group_resources={})
        with patch.object(backend,'open_registry',return_value=registry) as opener:
            first=backend.open_group('a',.01,ctx);second=backend.open_group('b',.01,ctx)
            first.new_state(1);second.new_state(2)
            self.assertEqual(opener.call_count,1)
            self.assertIs(ctx.group_resources['salmonpy.standard-tddft.registry/1'],registry)
            self.assertEqual(events,[('state','a',1),('state','b',2)])


class SharedHSETests(unittest.TestCase):

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

        # A never-initialized registered model must not release its live registry.
        from salmonpy._hse_registry import NativeFastHSERegistry,_ModelLease,_ModelLibrary,_RegisteredModel
        comm=SimpleNamespace(allgather=lambda value:[value])
        native=SimpleNamespace(closed=False,_live=lambda:None,context=SimpleNamespace(solver_comm=comm))
        def retire():native.closed=True
        native.close=retire
        registry=NativeFastHSERegistry(None,native,None)
        lease=_ModelLease(registry)
        child=_RegisteredModel(lease,_ModelLibrary(registry),Path('.'),.01,8,'hse_taylor4','id',())
        child._registry_id=0;lease.model=child;registry._models.add(child)
        with patch.object(fast_hse_mpi,'_used',True),patch.object(fast_hse_mpi.embedded,'_used',True):
            child.close()
            self.assertTrue(fast_hse_mpi._used);self.assertTrue(fast_hse_mpi.embedded._used)
            self.assertFalse(native.closed)
            registry.close();self.assertTrue(native.closed)

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
