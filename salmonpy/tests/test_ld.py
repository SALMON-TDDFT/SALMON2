"""LD: analytic response, owned state, immutable inputs and failure behavior."""
import gc
import unittest
from dataclasses import FrozenInstanceError, asdict
from unittest.mock import patch
import numpy as np
from salmonpy import LorentzDrude
from salmonpy.errors import RuntimeFailure, SnapshotError, NumericalError

class LDTests(unittest.TestCase):
    def material(self):
        return LorentzDrude(plasma_frequency=1., strength=[1., .2], resonance=[0., 2.], damping=[.1, .2])

    def test_vacuum(self):
        m = LorentzDrude(plasma_frequency=0., strength=[1.], resonance=[0.], damping=[0.])
        r = m.respond(np.ones((5,3)), dt=.1)
        np.testing.assert_array_equal(r.current, 0.)
        np.testing.assert_array_equal(r.polarization, 0.)

    def test_drude_constant_field(self):
        m = LorentzDrude(plasma_frequency=2., strength=[1.], resonance=[0.], damping=[0.])
        r = m.respond(np.ones((11,3)), dt=.01)
        np.testing.assert_allclose(r.current[:,0], r.time / np.pi, atol=1e-15)
        np.testing.assert_allclose(r.polarization[:,0], r.time**2 / (2*np.pi), atol=1e-15)


    def test_dielectric(self):
        m=LorentzDrude(plasma_frequency=2.,strength=[.5],resonance=[3.],damping=[.2])
        np.testing.assert_allclose(m.permittivity([1.,2.]),1+2/(9-np.array([1.,2.])**2-1j*.2*np.array([1.,2.])))
        bg=LorentzDrude(plasma_frequency=2.,strength=[.5],resonance=[3.],damping=[.2],epsilon_infinity=2.25)
        self.assertEqual(bg.config.epsilon_infinity,2.25)
        self.assertEqual(bg.config.background_permittivity,2.25)
        np.testing.assert_allclose(bg.permittivity([1.,2.]),m.permittivity([1.,2.])+1.25)
        fields=np.array([[.01,.02,0.],[.02,0.,.03],[0.,.01,.02]])
        ordinary=m.respond(fields,dt=.01);with_bg=bg.respond(fields,dt=.01)
        # P/J retain the dispersive-pole meaning: instantaneous background is in D.
        np.testing.assert_array_equal(ordinary.polarization,with_bg.polarization)
        np.testing.assert_array_equal(ordinary.current,with_bg.current)
        np.testing.assert_allclose(bg.displacement(fields,with_bg.polarization),2.25*fields+4*np.pi*with_bg.polarization)
        np.testing.assert_allclose(bg.background_energy(fields,np.zeros_like(fields)),2.25*np.sum(fields**2,axis=-1)/(8*np.pi))
        with bg.new_state(dt=.01) as state, bg.new_state(dt=.01) as branch:
            state.step(fields[0]);saved=state.snapshot()
            expected=state.step(fields[1]);branch.restore(saved)
            np.testing.assert_array_equal(branch.step(fields[1]),expected)
            state.restore(saved);np.testing.assert_array_equal(state.step(fields[1]),expected)
        for bad in (0.,-1.,float('nan'),float('inf'),True):
            with self.subTest(epsilon_infinity=bad),self.assertRaises(ValueError):
                LorentzDrude(plasma_frequency=1.,strength=[1.],resonance=[1.],damping=[0.],epsilon_infinity=bad)


    def test_damped_drude_convergence(self):
        m=LorentzDrude(plasma_frequency=1.,strength=[1.],resonance=[0.],damping=[.7])
        errors=[]
        for dt in [.1,.05]:
            r=m.respond(np.ones((round(2/dt)+1,3)),dt=dt)
            exact=(1-np.exp(-.7*r.time[-1]))/(4*np.pi*.7)
            errors.append(abs(r.current[-1,0]-exact))
        self.assertGreater(errors[0]/errors[1],3.9)
        self.assertLess(errors[0]/errors[1],4.1)

    def test_complex_rejected(self):
        with self.assertRaises(ValueError):
            self.material().respond(np.ones((2,3),dtype=complex)*(1+1j),dt=.1)
        with self.assertRaises(ValueError):
            LorentzDrude(plasma_frequency=1.,strength=np.array([1+1j]),resonance=[0.],damping=[0.])

    def test_invalid(self):
        for x in [float('nan'), float('inf'), -1.]:
            with self.assertRaises(ValueError):
                LorentzDrude(plasma_frequency=x,strength=[1.],resonance=[0.],damping=[0.])
        m=self.material()
        for dt in [0.,-1.,float('nan')]:
            with self.assertRaises(ValueError): m.respond(np.zeros((2,3)),dt=dt)
        with self.assertRaises(ValueError): m.respond([[0,0,float('inf')]],dt=.1)
        with self.assertRaises(ValueError): m.respond([1,2,3],dt=.1)


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


class InputContractTests(unittest.TestCase):
    def model(self):
        return LorentzDrude(plasma_frequency=1., strength=[1.], resonance=[0.], damping=[.1])

    def test_resolved_config(self):
        m=self.model(); s=m.new_state(dt=.1)
        self.assertEqual(m.config.strength, (1.,))
        self.assertEqual(m.config.units, 'atomic_gaussian')
        self.assertEqual(s.config.initial_field, (0.,0.,0.))
        self.assertEqual(s.config.model, m.config)
        self.assertEqual(s.config.dt, .1)
        self.assertEqual(s.config.initial_time, 0.)
        self.assertEqual(asdict(s.config)['model']['plasma_frequency'],1.)
        with self.assertRaises(FrozenInstanceError): m.config.plasma_frequency=2.
        with self.assertRaises(FrozenInstanceError): s.config.dt=.2

    def test_public_mutation_rejected(self):
        m=self.model(); s=m.new_state(dt=.1)
        for obj,name,value in [(m,'plasma_frequency',2.),(m,'strength',[2.]),
                               (m,'config',None),(s,'dt',.2),(s,'config',None),
                               (s,'material',self.model())]:
            with self.subTest(name=name), self.assertRaises(AttributeError):
                setattr(obj,name,value)
        with self.assertRaises(TypeError): s.step([0,0,0],dt=.2)

    def test_inputs_copied(self):
        f=np.array([1.]); e=np.array([1.,0.,0.])
        m=LorentzDrude(plasma_frequency=1.,strength=f,resonance=[0.],damping=[.1])
        s=m.new_state(dt=.1,initial_field=e)
        f[0]=9.;e[0]=9.
        self.assertEqual(m.config.strength,(1.,))
        self.assertEqual(s.config.initial_field,(1.,0.,0.))
        ref=self.model().new_state(dt=.1,initial_field=[1,0,0])
        np.testing.assert_array_equal(s.step([1,0,0]),ref.step([1,0,0]))

    def test_invalid_step_is_transactional(self):
        s=self.model().new_state(dt=.1); ref=self.model().new_state(dt=.1)
        s.step([1,0,0]);ref.step([1,0,0])
        before=(s.time,s.current.copy(),s.polarization.copy(),s.config)
        for field in ([1,2], [np.nan,0,0], [np.inf,0,0], [1j,0,0]):
            with self.assertRaises(ValueError): s.step(field)
            self.assertEqual(s.time,before[0])
            np.testing.assert_array_equal(s.current,before[1])
            np.testing.assert_array_equal(s.polarization,before[2])
            self.assertEqual(s.config,before[3])
        with patch.object(s._native_state,'advance',side_effect=RuntimeError('unavailable')):
            with self.assertRaises(RuntimeError): s.step([2,0,0])
        np.testing.assert_array_equal(s.step([2,0,0]),ref.step([2,0,0]))

    def test_scalar_validation(self):
        for dt in (True,[.1],np.array([.1]),'0.1'):
            with self.subTest(dt=dt), self.assertRaises(ValueError): self.model().new_state(dt=dt)

if __name__=='__main__':unittest.main()
