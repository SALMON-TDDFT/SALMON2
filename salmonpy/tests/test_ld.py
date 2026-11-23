"""LD: analytic response, owned state, immutable inputs and failure behavior."""
import gc
import unittest
from dataclasses import FrozenInstanceError, asdict
from unittest.mock import patch
import numpy as np
from salmonpy import LorentzDrude
from salmonpy.errors import RuntimeFailure, SnapshotError, NumericalError

class LDTests(unittest.TestCase):
    def test_performance_modes_match_batch_steps_and_explicit_restore(self):
        from salmonpy.ld import TensorLorentzDrude
        materials=(self.material(),TensorLorentzDrude(plasma_frequency=.7,
            strength=[[[1.,.2,.1],[.2,.6,.05],[.1,.05,.4]]],resonance=[.3],damping=[.1]))
        time=np.arange(65)*.02;drive=np.column_stack((np.sin(time),np.cos(time)-1,np.sin(.7*time)))*.01
        for material in materials:
            with self.subTest(material=type(material).__name__):
                normal=material.respond(drive,dt=.02)
                fast=material.respond(drive,dt=.02,execution_mode='performance')
                np.testing.assert_array_equal(fast.current,normal.current)
                np.testing.assert_array_equal(fast.polarization,normal.polarization)
                a=material.new_state(dt=.02);b=material.new_state(dt=.02,execution_mode='performance')
                other=material.new_state(dt=.02,execution_mode='performance')
                try:
                    self.assertEqual(b.execution_mode,'performance')
                    with self.assertRaises(AttributeError):b.execution_mode='checked'
                    for n,e in enumerate(drive[1:],1):
                        p,j=a.step(e);q,k=b.step(e);other.step(-e)
                        np.testing.assert_array_equal(q,p);np.testing.assert_array_equal(k,j)
                        if n==24:saved=b.snapshot()
                    b.restore(saved)
                    for e in drive[25:]:b.step(e)
                    np.testing.assert_array_equal(a.current,b.current)
                    b.step([np.nan,0.,0.]);self.assertTrue(np.isnan(b.current).any())
                    np.testing.assert_array_equal(a.current,normal.current[-1])
                    b.restore(saved);self.assertEqual(b.execution_mode,'performance')
                finally:a.close();b.close();other.close()



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


class LDFittingTests(unittest.TestCase):
    def test_passive_fit_background_band_and_provenance(self):
        from salmonpy.ld import fit_lorentz_drude
        import tempfile,hashlib
        from pathlib import Path
        omega=np.linspace(.05,1.2,180)
        target=LorentzDrude(plasma_frequency=1.,strength=[.03,.4],resonance=[0.,.8],damping=[.04,.07],epsilon_infinity=2.5)
        epsilon=target.permittivity(omega);original=epsilon.copy()
        with tempfile.TemporaryDirectory() as tmp:
            source=Path(tmp)/'synthetic.txt';source.write_text('synthetic reference only\n')
            cfg,report=fit_lorentz_drude(omega,epsilon,poles=2,drude=True,units='atomic_gaussian',source_file=source,reference='synthetic Drude + Lorentz',rtol=.005)
            self.assertTrue(report.accepted,report.as_dict())
            self.assertLess(report.relative_rms_error,.005)
            self.assertEqual(report.source_sha256,hashlib.sha256(source.read_bytes()).hexdigest())
            self.assertEqual(report.reference,'synthetic Drude + Lorentz')
            self.assertEqual(cfg.valid_frequency_range,(.05,1.2))
            self.assertTrue(all(v>=0 for v in (*cfg.strength,*cfg.resonance,*cfg.damping)))
            self.assertGreater(cfg.epsilon_infinity,0.)
            fitted=LorentzDrude.from_config(cfg)
            error=np.linalg.norm(fitted.permittivity(omega)-epsilon)/np.linalg.norm(epsilon)
            self.assertLess(error,.005)
            np.testing.assert_array_equal(epsilon,original)
            with self.assertRaises(FrozenInstanceError):report.accepted=False
            # n+i k -> epsilon conversion uses the same fitter and provenance.
            refractive=np.sqrt(epsilon)
            other,second=fit_lorentz_drude(omega,nk=np.column_stack((refractive.real,refractive.imag)),poles=2,drude=True,units='atomic_gaussian',initial=cfg,starts=1,rtol=.005)
            self.assertTrue(second.accepted)
            self.assertEqual(second.data_kind,'nk')
            from salmonpy.analysis import OpticalSpectrum
            mask=np.ones(len(omega),dtype=bool);mask[0]=False
            observable=OpticalSpectrum(omega,epsilon,mask,{'observable':'dielectric_function',
                'units':'atomic_gaussian','sign':'exp(+i omega t)','window':'none'})
            _,optical_report=fit_lorentz_drude(observable,poles=2,drude=True,units='atomic_gaussian',initial=cfg,starts=1,rtol=.005)
            self.assertTrue(optical_report.accepted)
            self.assertEqual(optical_report.data_kind,'optical_spectrum')
            self.assertEqual(optical_report.as_dict()['provenance']['excluded_masked_samples'],1)
            # Finite DC data are valid for Lorentz-only fits with a fixed background.
            dc_omega=np.linspace(0,1.2,50)
            dc_model=LorentzDrude(plasma_frequency=1.,strength=[.4],resonance=[.8],damping=[.07],epsilon_infinity=2.5)
            dc_config,dc_report=fit_lorentz_drude(dc_omega,dc_model.permittivity(dc_omega),poles=1,
                units='atomic_gaussian',drude=False,epsilon_infinity=2.5,initial=dc_model,starts=1)
            self.assertTrue(dc_report.accepted);self.assertFalse(dc_report.fitted_background)
            self.assertEqual(dc_config.epsilon_infinity,2.5)


            self.assertEqual(other.valid_frequency_range,cfg.valid_frequency_range)
            from salmonpy.analysis import spectrum
            time=np.arange(512)*.1
            drive=spectrum(time,np.cos(2*np.pi*3*np.arange(512)/512),window='none')
            self.assertTrue(fitted.check_excitation(drive)['within_band'])
            impulse=np.zeros(512);impulse[0]=1
            with self.assertWarnsRegex(RuntimeWarning,'calibrated'):
                self.assertFalse(fitted.check_excitation(spectrum(time,impulse,window='none'))['within_band'])


    def test_fit_failure_is_reported_without_active_coefficients(self):
        from salmonpy.ld import fit_lorentz_drude
        omega=np.linspace(.1,1.,50);epsilon=np.full(50,2.-1j)
        with self.assertWarns(RuntimeWarning):
            cfg,report=fit_lorentz_drude(omega,epsilon,poles=1,units='atomic_gaussian',max_iterations=1,starts=1,rtol=1e-8)
        self.assertFalse(report.accepted)
        if cfg is not None:self.assertTrue(all(v>=0 for v in (*cfg.strength,*cfg.resonance,*cfg.damping)))
        with self.assertRaises(ValueError):fit_lorentz_drude(omega,epsilon,poles=1,units='SI')
        with self.assertRaisesRegex(ValueError,'workspace'):
            fit_lorentz_drude(omega,epsilon,poles=700_000_000,units='atomic_gaussian')
        bad=epsilon.copy();bad[3]=np.nan
        with self.assertRaises(ValueError):fit_lorentz_drude(omega,bad,poles=1,units='atomic_gaussian')


class ScipyFitTests(unittest.TestCase):
    def test_explicit_fit_damping_floor(self):
        from salmonpy.ld import fit_lorentz_drude
        import warnings
        omega=np.linspace(.03,.4,90)
        initial=LorentzDrude(plasma_frequency=1.,strength=[.12,.08],
                            resonance=[0.,.22],damping=[.025,.04],epsilon_infinity=2.)
        target=initial.permittivity(omega)
        optimizers=['native']
        try:
            import scipy
            optimizers.append('scipy')
        except ImportError:
            pass
        for optimizer in optimizers:
            with self.subTest(optimizer=optimizer),warnings.catch_warnings():
                warnings.simplefilter('ignore')
                default,old=fit_lorentz_drude(omega,target,poles=2,units='atomic_gaussian',
                    initial=initial,optimizer=optimizer,starts=1)
                zero,explicit=fit_lorentz_drude(omega,target,poles=2,units='atomic_gaussian',
                    initial=initial,optimizer=optimizer,starts=1,minimum_damping=0.)
                self.assertEqual(default,zero)
                self.assertEqual(old.data_sha256,explicit.data_sha256)
                constrained,report=fit_lorentz_drude(omega,target,poles=2,units='atomic_gaussian',
                    initial=initial,optimizer=optimizer,starts=1,minimum_damping=.03)
                self.assertTrue(all(value>=.03-1e-14 for value in constrained.damping))
                self.assertEqual(report.as_dict()['minimum_damping_au'],.03)
        for invalid in (-1.,float('nan'),float('inf')):
            with self.assertRaises(ValueError):
                fit_lorentz_drude(omega,target,poles=2,units='atomic_gaussian',minimum_damping=invalid)

    def test_scipy_optional_native_objective(self):
        try:
            import scipy
        except ImportError:
            self.skipTest("optional scipy is not installed")
        from salmonpy.ld import fit_lorentz_drude
        import warnings
        omega=np.linspace(.03,.4,90)
        original=LorentzDrude(plasma_frequency=1.,strength=[.12,.08],resonance=[0.,.22],damping=[.025,.04],epsilon_infinity=2.)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            cfg,report=fit_lorentz_drude(omega,original.permittivity(omega),poles=2,units="atomic_gaussian",optimizer="scipy",starts=2,max_iterations=500)
        self.assertTrue(report.accepted)
        self.assertTrue(report.optimizer_converged)
        self.assertLess(report.relative_rms_error,1e-6)
        self.assertIn("scipy",report.as_dict()["algorithm"])
        self.assertEqual(report.as_dict()["provenance"]["optimizer"]["version"],scipy.__version__)
        np.testing.assert_allclose(LorentzDrude.from_config(cfg).permittivity(omega),original.permittivity(omega),rtol=1e-5,atol=1e-6)
        with self.assertRaises(ValueError):
            fit_lorentz_drude(omega,original.permittivity(omega),poles=2,units="atomic_gaussian",optimizer="unknown")
        from unittest.mock import patch
        from salmonpy.errors import RuntimeFailure
        with patch.dict("sys.modules", {"scipy":None}):
            with self.assertRaisesRegex(RuntimeFailure,"optional"):
                fit_lorentz_drude(omega,original.permittivity(omega),poles=2,units="atomic_gaussian",optimizer="scipy")


class AnisotropicLDTests(unittest.TestCase):
    def test_rotated_response_and_owned_restore(self):
        from salmonpy.ld import TensorLorentzDrude
        rotation=np.array([[.6,-.8,0.],[.8,.6,0.],[0.,0.,1.]])
        strength=rotation@np.diag([.1,.2,.3])@rotation.T
        background=rotation@np.diag([2.,3.,4.])@rotation.T
        tensor=TensorLorentzDrude(plasma_frequency=1.,strength=[strength],resonance=[.4],damping=[.03],epsilon_infinity=background)
        w=np.linspace(.01,.8,40)
        actual=tensor.permittivity_tensor(w)
        expected=background[None,:,:]+strength[None,:,:]/(.4**2-w[:,None,None]**2-1j*.03*w[:,None,None])
        np.testing.assert_allclose(actual,expected,rtol=2e-13,atol=1e-13)
        self.assertLess(actual[0,0,1].real,0.)
        with self.assertRaises(ValueError):actual.setflags(write=True)
        time=np.arange(80)*.01
        electric=np.column_stack((np.sin(time),np.cos(time),np.sin(2*time)))
        local=electric@rotation
        pp=np.empty_like(electric);jj=np.empty_like(electric)
        for axis,f in enumerate([.1,.2,.3]):
            scalar=LorentzDrude(plasma_frequency=1.,strength=[f],resonance=[.4],damping=[.03])
            result=scalar.respond(local,dt=.01);pp[:,axis]=result.polarization[:,axis];jj[:,axis]=result.current[:,axis]
        result=tensor.respond(electric,dt=.01)
        np.testing.assert_allclose(result.polarization,pp@rotation.T,rtol=1e-12,atol=1e-14)
        np.testing.assert_allclose(result.current,jj@rotation.T,rtol=1e-12,atol=1e-14)
        np.testing.assert_allclose(tensor.displacement(electric,result.polarization),electric@background.T+4*np.pi*result.polarization,rtol=2e-14,atol=1e-14)
        np.testing.assert_allclose(tensor.background_energy(electric,np.zeros_like(electric)),np.sum(electric*(electric@background.T),axis=1)/(8*np.pi),rtol=2e-14)
        self.assertEqual(tensor.config.representation,'reciprocal_tensor')
        self.assertEqual(LorentzDrude.from_config(tensor.config).config,tensor.config)
        with tensor.new_state(dt=.01,initial_field=electric[0]) as state, tensor.new_state(dt=.01,initial_field=electric[0]) as independent:
            for e in electric[1:30]:state.step(e)
            saved=state.snapshot()
            first=[state.step(e) for e in electric[30:]]
            state.restore(saved)
            second=[state.step(e) for e in electric[30:]]
            np.testing.assert_array_equal(first,second)
            self.assertEqual(independent.time,0.)
            np.testing.assert_array_equal(independent.current,np.zeros(3))
            before=state.snapshot()
            with self.assertRaises(ValueError):state.step([float('inf'),0.,0.])
            np.testing.assert_array_equal(state.current,second[-1][1])
        with tensor.new_state(dt=.01) as state:
            state.step([1e308,1e308,1e308])
            previous=state.snapshot();previous_current=state.current.copy();previous_time=state.time
            with self.assertRaises(NumericalError):state.step([1e308,1e308,1e308])
            np.testing.assert_array_equal(state.current,previous_current)
            self.assertEqual(state.time,previous_time)
            retry=state.step([0.,0.,0.]);state.restore(previous)
            np.testing.assert_array_equal(retry,state.step([0.,0.,0.]))
        from salmonpy.material_response import LDMaterialProvider
        self.assertEqual(LDMaterialProvider(tensor.config).prepare(None).config,tensor.config)

    def test_anisotropic_input_rejection(self):
        from salmonpy.ld import TensorLorentzDrude
        options=dict(plasma_frequency=1.,strength=[np.diag([.1,.2,.3])],resonance=[.4],damping=[.03])
        for bad in [np.diag([1.,-1.,1.]),np.array([[1.,.2,0.],[0.,1.,0.],[0.,0.,1.]]),np.full((3,3),float('nan'))]:
            with self.assertRaises(ValueError):TensorLorentzDrude(**options,epsilon_infinity=bad)
        with self.assertRaises(ValueError):TensorLorentzDrude(**dict(options,strength=[np.diag([.1,-.2,.3])]))
        with self.assertRaises(ValueError):TensorLorentzDrude(**dict(options,damping=[-.1]))


class TensorFitTests(unittest.TestCase):
    def test_passive_non_diagonal_roundtrip(self):
        try:import scipy
        except ImportError:self.skipTest('optional SciPy unavailable')
        from salmonpy.ld import TensorLorentzDrude,fit_lorentz_drude_tensor
        f=np.array([[[.5,-.12,0],[-.12,.3,.02],[0,.02,.2]]])
        background=np.diag([1.2,1.3,1.4])
        model=TensorLorentzDrude(plasma_frequency=1.,strength=f,resonance=[.7],damping=[.08],epsilon_infinity=background)
        initial=TensorLorentzDrude(plasma_frequency=1.,strength=f*.9,resonance=[.65],damping=[.1],epsilon_infinity=background)
        w=np.linspace(.1,1.2,100);target=model.permittivity_tensor(w)
        candidate,report=fit_lorentz_drude_tensor(w,target,initial=initial,units='atomic_gaussian',rtol=1e-6)
        np.testing.assert_allclose(candidate.permittivity_tensor(w),target,rtol=1e-6,atol=1e-6)
        self.assertTrue(report['accepted'])
        self.assertTrue(report['optimizer_converged'])
        self.assertGreaterEqual(np.linalg.eigvalsh(candidate.config.strength[0]).min(),-1e-14)
        np.testing.assert_array_equal(initial.config.epsilon_infinity,background)
