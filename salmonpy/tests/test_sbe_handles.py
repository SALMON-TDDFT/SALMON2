"""Persistent model/state ownership and no hot-path density export."""
import gc
import unittest
from unittest.mock import patch
import numpy as np
from salmonpy.errors import RuntimeFailure,SnapshotError
from salmonpy.sbe import SBEInitialState


def model():
    p=np.zeros((2,2,3,1),complex);p[0,1,2,0]=p[1,0,2,0]=.3
    return SBEInitialState([[-.5],[.4]],p,p*.1,[1.],np.eye(3)*4,2)


class SBEHandles(unittest.TestCase):
    def test_model_shared_states_independent_and_density_on_demand(self):
        initial=model();a=initial.new_state(dt=.02);b=initial.new_state(dt=.02)
        self.addCleanup(a.close);self.addCleanup(b.close)
        self.assertIs(a._native_state.model,b._native_state.model)
        from salmonpy import _sbe_handles as handles
        with patch.object(handles.NativeState,'density',side_effect=AssertionError('step must not export rho')):
            a.step([.01,0,.02]);a.step([.02,0,.03])
        np.testing.assert_array_equal(b.density,np.diag([2.,0.])[:,:,None])
        self.assertFalse(a.density.flags.writeable)

    def test_exact_stateless_parity_all_corrections(self):
        rng=np.random.default_rng(38219);nb=3;nk=7
        p=rng.normal(size=(nb,nb,3,nk))+1j*rng.normal(size=(nb,nb,3,nk))
        p=(p+p.swapaxes(0,1).conj())*.05
        eigen=np.linspace(-.5,.7,nb)[:,None]+np.arange(nk)[None,:]*.01
        weights=rng.uniform(.1,1.,nk);values=rng.normal(size=(15,3))*.02
        for order in range(4):
            for vnl in (False,True):
                initial=SBEInitialState(eigen,p,p*.13,weights,np.diag([4.,5.,6.]),2,order,vnl)
                with initial.new_state(dt=.013,initial_A=values[0]) as state:
                    api=state._native_state.model.api;rho=np.zeros((nb,nb,nk),complex,order='F');rho[0,0,:]=2
                    for index,A in enumerate(values):
                        old=np.ascontiguousarray(values[max(index-1,0)]);new=np.ascontiguousarray(A);j=np.empty(3)
                        arrays=(initial.eigen,initial.weights,initial.momentum,initial.nonlocal_velocity,rho,old,new)
                        api.legacy(nb,nk,2,order,int(vnl),initial.volume,.013,
                            *(a.ctypes.data for a in arrays),int(index>0),j.ctypes.data)
                        actual=state.matter_current if index==0 else -state.step(A)
                        with self.subTest(order=order,vnl=vnl,step=index):
                            self.assertEqual(actual.tobytes(),j.tobytes())
                            self.assertEqual(state.density.tobytes(order='F'),rho.tobytes(order='F'))

    def test_used_initial_model_remains_data_only_and_serializable(self):
        import copy,pickle
        from dataclasses import asdict
        initial=model()
        with initial.new_state(dt=.02) as state:state.step([.01,0,0])
        data=asdict(initial)
        self.assertNotIn('_native_model',data)
        for cloned in (pickle.loads(pickle.dumps(initial)),copy.deepcopy(initial)):
            self.assertFalse(cloned.momentum.flags.writeable)
            with cloned.new_state(dt=.02) as a,initial.new_state(dt=.02) as b:
                self.assertEqual(a.step([0,0,.01]).tobytes(),b.step([0,0,.01]).tobytes())
                self.assertIsNot(a._native_state.model,b._native_state.model)

    def test_finalizer_registration_failure_releases_native_resources(self):
        from salmonpy import _sbe_handles as handles
        from salmonpy.sbe import _handles_api
        api=_handles_api();initial=model()
        with patch.object(api,'call',wraps=api.call) as calls:
            with patch.object(handles.weakref,'finalize',side_effect=MemoryError('registry full')):
                with self.assertRaises(MemoryError):initial.new_state(dt=.02)
            self.assertTrue(any(c.args[0]=='sp_sbe_model_destroy' for c in calls.call_args_list))
        with initial.new_state(dt=.02):pass
        with patch.object(api,'call',wraps=api.call) as calls:
            with patch.object(handles.weakref,'finalize',side_effect=MemoryError('registry full')):
                with self.assertRaises(MemoryError):initial.new_state(dt=.02)
            self.assertTrue(any(c.args[0]=='sp_sbe_state_destroy' for c in calls.call_args_list))
        with initial.new_state(dt=.02) as state:state.step([0,0,.01])

    def test_weak_cache_does_not_keep_released_model_or_state_alive(self):
        import weakref
        initial=model();state=initial.new_state(dt=.02)
        owner=weakref.ref(initial);model_finalizer=state._native_state.model.finalizer
        state_finalizer=state._native_state.finalizer
        del initial,state;gc.collect()
        self.assertIsNone(owner());self.assertFalse(model_finalizer.alive);self.assertFalse(state_finalizer.alive)

    def test_close_and_finalization_release_state_handles(self):
        initial=model();state=initial.new_state(dt=.02);native=state._native_state
        state.close();state.close();self.assertFalse(native.alive)
        with self.assertRaises(RuntimeFailure):state.step([0,0,0])
        with self.assertRaises(SnapshotError):state.snapshot()
        with self.assertRaises(SnapshotError):state.restore(None)
        other=initial.new_state(dt=.02);finalizer=other._native_state.finalizer
        del other;gc.collect();self.assertFalse(finalizer.alive)

    def test_nonfinite_or_overflow_rejection_preserves_native_state(self):
        initial=model()
        with initial.new_state(dt=.02) as state:
            state.step([.01,.02,.03]);rho=state.density.copy();current=state.current.copy();time=state.time
            with self.assertRaises(ValueError):state.step([np.nan,0,0])
            with self.assertRaises((ValueError,FloatingPointError)):state.step([1e308,1e308,1e308])
            self.assertEqual(state.time,time);np.testing.assert_array_equal(state.current,current)
            np.testing.assert_array_equal(state.density,rho)
            state.step([0,0,0])

if __name__=='__main__':unittest.main()
