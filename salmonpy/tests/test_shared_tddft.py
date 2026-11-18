"""Shared native ownership contracts; numerical parity is covered by the MPI probe."""
import tempfile
from pathlib import Path
from types import SimpleNamespace
import unittest
from unittest.mock import patch
from salmonpy.errors import InputError

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
        from salmonpy._shared_tddft import _Library
        calls=[]
        def fn(*args):calls.append(args);return 0
        fn.argtypes=['native ABI']
        lib=SimpleNamespace(**{'sp_tddft_handle_'+key:fn for key in ('step','close','checkpoint','transaction')})
        a,b=_Library(lib,31),_Library(lib,32)
        a.sp_tddft_step.argtypes=['compatibility ABI'];a.sp_tddft_step(1);b.sp_tddft_step(2)
        self.assertEqual(calls,[(31,1),(32,2)]);self.assertEqual(fn.argtypes,['native ABI'])

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
