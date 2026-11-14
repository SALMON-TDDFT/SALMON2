"""Persistent LD ownership, reference arithmetic and transactional state."""
import gc
import unittest
from unittest.mock import patch
import numpy as np
from salmonpy import LorentzDrude
from salmonpy.errors import RuntimeFailure,SnapshotError,NumericalError


def model(n=8):
    return LorentzDrude(plasma_frequency=1.,strength=np.ones(n)/n,
                       resonance=np.linspace(0,2,n),damping=np.full(n,.1))


class LDHandles(unittest.TestCase):
    def test_explicit_lifetime_and_no_pole_export_on_step(self):
        s=model().new_state(dt=.01)
        self.assertTrue(hasattr(s,'close'),'LD needs an owned persistent state')
        self.addCleanup(s.close)
        from salmonpy import _ld_handles as handles
        with patch.object(handles.NativeState,'export',side_effect=AssertionError('pole export on hot path')):
            s.step([.01,.02,.03]);s.current;s.polarization
        native=s._native_state;s.close();s.close();self.assertFalse(native.alive)
        with self.assertRaises(RuntimeFailure):s.step([0,0,0])
        with self.assertRaises(SnapshotError):s.snapshot()
        with self.assertRaises(SnapshotError):s.restore(None)
        other=model().new_state(dt=.01);finalizer=other._native_state.finalizer
        del other;gc.collect();self.assertFalse(finalizer.alive)

    def test_stateless_parity_and_snapshots(self):
        import ctypes
        from salmonpy.ld import _native
        rng=np.random.default_rng(188)
        for n in (1,8,64):
            m=model(n);fields=rng.normal(size=(31,3))*.02
            pp=np.zeros((n,3));jj=np.zeros_like(pp);po=np.empty_like(fields);jo=np.empty_like(fields)
            ptr=ctypes.POINTER(ctypes.c_double)
            _native()(n,len(fields),.013,m._wp,
                *(a.ctypes.data_as(ptr) for a in (m._f,m._w,m._g,fields,pp,jj,po,jo)))
            with m.new_state(dt=.013,initial_field=fields[0]) as s:
                for i,e in enumerate(fields[1:],1):
                    p,j=s.step(e)
                    self.assertEqual(p.tobytes(),po[i].tobytes())
                    self.assertEqual(j.tobytes(),jo[i].tobytes())
                    if i==10:token=s.snapshot()
                final=s.snapshot();s.restore(token)
                for e in fields[11:]:s.step(e)
                self.assertEqual(s.snapshot()._data,final._data)
                arrays=s._snapshot_arrays()
                self.assertEqual(arrays['_p'].tobytes(),pp.tobytes())
                self.assertEqual(arrays['_j'].tobytes(),jj.tobytes())
            response=m.respond(fields,dt=.013)
            self.assertEqual(response.polarization.tobytes(),po.tobytes())
            self.assertEqual(response.current.tobytes(),jo.tobytes())

    def test_failed_batch_and_import_preserve_state(self):
        with model().new_state(dt=.01) as s:
            s.step([.01,0,0]);token=s.snapshot()
            with self.assertRaises(ValueError):s.step([np.nan,0,0])
            with self.assertRaises(NumericalError):s._advance(np.array([s._field,[1e308]*3,[1e308]*3]))
            self.assertEqual(s.snapshot()._data,token._data);self.assertEqual(s.time,token.time)
            from dataclasses import replace
            with self.assertRaises(SnapshotError):s.restore(replace(token,time=token.time+1))
            self.assertEqual(s.snapshot()._data,token._data)
            s.step([0,0,0])

    def test_outputs_do_not_alias_state_and_all_import_elements_checked(self):
        with model().new_state(dt=.01) as s:
            p,j=s.step([.01,.02,.03]);accepted=s.snapshot()
            p[:]=99;j[:]=99
            self.assertEqual(s.snapshot()._data,accepted._data)
            fields=np.array([s._field,[.02,.03,.04]])
            p,j=s._advance(fields);accepted=s.snapshot();p[:]=99;j[:]=99
            self.assertEqual(s.snapshot()._data,accepted._data)
            raw=s._snapshot_arrays()
            for name in ('_p','_j'):
                for value in (np.nan,np.inf,-np.inf):
                    for index in np.ndindex(raw[name].shape):
                        a=raw['_p'].copy();b=raw['_j'].copy()
                        (a if name=='_p' else b)[index]=value
                        with self.assertRaises(ValueError):s._native_state.replace(raw['_field'],a,b)
                        self.assertEqual(s.snapshot()._data,accepted._data)

    def test_response_backend_closes_owned_native_state(self):
        from salmonpy.material_response import LDResponseBackend
        material=model();raw=material.new_state(dt=.01)
        self.addCleanup(raw.close)
        with patch.object(LorentzDrude,'new_state',return_value=raw):
            wrapped=LDResponseBackend().create_state(material,0,.01,None)
        wrapped.close()
        self.assertFalse(raw._native_state.alive)

    def test_finalizer_registration_failure_cleanup(self):
        from salmonpy import _ld_handles as handles
        from salmonpy.ld import _handles_api
        api=_handles_api()
        with patch.object(api,'call',wraps=api.call) as calls:
            with patch.object(handles.weakref,'finalize',side_effect=MemoryError('registry full')):
                with self.assertRaises(MemoryError):model().new_state(dt=.01)
            self.assertTrue(any(c.args[0]=='sp_ld_destroy' for c in calls.call_args_list))
        with model().new_state(dt=.01) as s:s.step([0,0,.01])

if __name__=='__main__':unittest.main()
