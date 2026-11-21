import tempfile
import unittest
from pathlib import Path
import numpy as np
from salmonpy.fields import AcosPulse, sample_pulses
from salmonpy.geometry import Grid, Shape, Geometry
from salmonpy.analysis import read_table, spectrum, optical_conductivity, dielectric_function, optical_linearity, optical_sum_rule

class UtilityTests(unittest.TestCase):
    def test_fdtd_point_work_ports_and_energy_budget(self):
        from salmonpy.analysis import absorbed_energy_density,reflection_transmission
        from salmonpy.fdtd import energy_balance
        t=np.arange(128)*.1
        f=np.sin(2*np.pi*7*np.arange(128)/128)
        E=np.zeros((len(t),3));H=E.copy();E[:,0]=f;H[:,1]=f
        R,T=.25,.5
        ports=reflection_transmission(t,E,H,E*np.sqrt(R),-H*np.sqrt(R),E*np.sqrt(T),H*np.sqrt(T),
                                      normal=(0,0,1),units='atomic_gaussian')
        np.testing.assert_allclose(ports.R[ports.valid],R,atol=1e-14)
        np.testing.assert_allclose(ports.T[ports.valid],T,atol=1e-14)
        self.assertEqual(ports.valid.dtype,np.dtype(bool))
        impedance=reflection_transmission(t,E,H,E*np.sqrt(R),-H*np.sqrt(R),
            E*np.sqrt(T/2),H*np.sqrt(T*2),normal=(0,0,1),units='atomic_gaussian')
        np.testing.assert_allclose(impedance.T[impedance.valid],T,atol=1e-14)
        volume=np.array([2.,3.]);area=2.
        local_E=np.zeros((len(t),2,3));local_E[:,:,0]=1
        j=np.zeros_like(local_E)
        j[:,:,0]=f[:,None]**2*(137.035999084/(4*np.pi)*(1-R-T)*area)*np.array([.2,.8])/volume
        old=j.copy()
        work=absorbed_energy_density(t,local_E,j,units='atomic_gaussian',current_convention='charge',field_kind='total')
        expected=137.035999084/(4*np.pi)*np.trapezoid(f*f,t)*(1-R-T)*area*np.array([.2,.8])/volume
        np.testing.assert_allclose(work,expected,rtol=1e-14)
        np.testing.assert_array_equal(j,old)
        left=absorbed_energy_density(t[:65],local_E[:65],j[:65],units='atomic_gaussian',current_convention='charge',field_kind='total')
        right=absorbed_energy_density(t[64:],local_E[64:],j[64:],units='atomic_gaussian',current_convention='charge',field_kind='total')
        np.testing.assert_allclose(left+right,work,rtol=1e-14)
        report=energy_balance(work,volume,ports,material_ids=('film','film'),transverse_area=area,
                              field_energy_change=0.,units='atomic_gaussian',rtol=1e-12)
        self.assertTrue(report.within_tolerance)
        self.assertAlmostEqual(report.relative_residual,0.,places=14)
        for a in (work,ports.R,ports.T,ports.valid,report.absorbed_energy_density):
            with self.assertRaises(ValueError):a.setflags(write=True)
        # Stored fields and other outgoing flux are separate from material work.
        missing=report.incident_energy*.1
        with self.assertWarnsRegex(UserWarning,'energy imbalance'):
            bad=energy_balance(work,volume,ports,material_ids=('a','b'),transverse_area=area,
                               field_energy_change=missing,units='atomic_gaussian')
        self.assertFalse(bad.within_tolerance)
        self.assertAlmostEqual(bad.relative_residual,-.1,places=14)
        good=energy_balance(work*.5,volume,ports,material_ids=('a','b'),transverse_area=area,
                            field_energy_change=report.material_work*.25,
                            other_outgoing_energy=report.material_work*.25,units='atomic_gaussian',rtol=1e-12)
        self.assertTrue(good.within_tolerance)
        vacuum=reflection_transmission(t,E,H,E*0,H*0,E,H,normal=(0,0,2),units='atomic_gaussian',window='hann')
        self.assertEqual(vacuum.incident_fluence,vacuum.transmitted_fluence)
        np.testing.assert_allclose(vacuum.R[vacuum.valid],0)
        np.testing.assert_allclose(vacuum.T[vacuum.valid],1)

    def test_energy_diagnostics_reject_bad_contracts_and_native_nonfinite(self):
        import ctypes as c
        from salmonpy.analysis import absorbed_energy_density,reflection_transmission,_energy_native
        t=np.arange(4,dtype=float);a=np.ones((4,1,3))
        args=dict(units='atomic_gaussian',current_convention='charge',field_kind='total')
        for value in (float('nan'),float('inf')):
            bad=a.copy();bad[1,0,1]=value
            with self.assertRaises(ValueError):absorbed_energy_density(t,bad,a,**args)
            # Direct native path also checks every element using scalar IEEE.
            out=np.empty(1)
            fn=_energy_native('sp_absorbed_work',[c.c_int,c.c_int]+[c.c_void_p]*4)
            self.assertNotEqual(fn(4,1,t.ctypes.data,bad.ctypes.data,a.ctypes.data,out.ctypes.data),0)
        with self.assertRaisesRegex(ValueError,'charge'):absorbed_energy_density(t,a,a,**dict(args,current_convention='matter'))
        with self.assertRaises(ValueError):absorbed_energy_density(t,a[:,:,0],a,**args)
        with self.assertRaises(ValueError):absorbed_energy_density(t[::-1],a,a,**args)
        fields=np.ones((4,3))
        with self.assertRaises(ValueError):reflection_transmission(t,*([fields]*6),normal=(0,0,0),units='atomic_gaussian')
        with self.assertWarnsRegex(UserWarning,'unresolved'):
            empty=reflection_transmission(t,*([fields*0]*6),normal=(0,0,1),units='atomic_gaussian')
        self.assertFalse(np.any(empty.valid))
        self.assertTrue(np.isnan(empty.R).all())

    def test_pulse_and_composition(self):
        pulse=AcosPulse(amplitude=.1,omega=2.,duration=4.,start=1.,polarization=(1,0,0))
        t=np.linspace(0,6,601);a=pulse.sample(t)
        np.testing.assert_array_equal(a.values[[0,-1]],0.)
        np.testing.assert_allclose(sample_pulses([pulse,pulse],t).values,2*a.values)
        e=a.electric_field()
        np.testing.assert_allclose(e.values,-np.diff(a.values,axis=0)/.01,atol=1e-13)
        self.assertEqual(a.kind,'A');self.assertEqual(e.kind,'E')
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'pulse.data';a.write_salmon(p)
            np.testing.assert_allclose(np.loadtxt(p)[:,1:],a.values)

    def test_grid_and_shapes(self):
        grid=Grid(counts=(4,4,4),spacing=(1,1,1))
        np.testing.assert_array_equal(grid.coordinates[0],[-1.5,-.5,.5,1.5])
        geom=Geometry(grid,shapes=(Shape('box',(1,1,1),material=1),))
        self.assertEqual(np.count_nonzero(geom.material_ids()),8)
        self.assertEqual(geom.material_ids().shape,(4,4,4))
        with self.assertRaises(ValueError): Grid(counts=(4.2,4,4),spacing=(1,1,1))
        with self.assertRaises(ValueError): Shape('box',(-1,1,1))

    def test_read_fortran_and_units(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'obs.data'
            p.write_text('# 1:Time[a.u.] 2:E_x[a.u.]\n0d0 1D0\n1d0 2d0\n')
            table=read_table(p)
            self.assertEqual(table.columns,('Time','E_x'))
            self.assertEqual(table.units,('a.u.','a.u.'))
            np.testing.assert_array_equal(table.column('E_x'),[1,2])

    def test_fourier_convention(self):
        n=64;dt=.1;t=np.arange(n)*dt+1.2;omega=2*np.pi*3/(n*dt)
        x=np.cos(omega*t)
        s=spectrum(t,x,window='none')
        self.assertAlmostEqual(s.omega[3],omega)
        np.testing.assert_allclose(s.values[3],n*dt/2,atol=1e-12)
        self.assertEqual(s.config['sign'],'exp(+i omega t)')
        with self.assertRaises(ValueError): spectrum([0,1,3],[0,1,0])
        with self.assertRaises(ValueError): spectrum([0,1],[0,1],window='unknown')

    def test_optical_observable_contract(self):
        t=np.arange(64)*.1+1.2
        w=2*np.pi*3/(64*.1)
        e=np.cos(w*t);j=2*e+.3*np.sin(w*t)
        args=dict(current_convention='charge',field_kind='total',units='atomic_gaussian',window='none')
        original=j.copy()
        result=optical_conductivity(t,j,e,**args)
        eps=dielectric_function(result)
        np.testing.assert_allclose(result.values[3],2+.3j,atol=1e-13)
        np.testing.assert_allclose(eps.values[3],1+4*np.pi*1j*(2+.3j)/w,atol=1e-12)
        self.assertTrue(result.valid[3]);self.assertFalse(eps.valid[0])
        self.assertFalse(result.valid[4]);self.assertTrue(np.isnan(result.values[4]))
        self.assertEqual(eps.config['field_kind'],'total')
        np.testing.assert_array_equal(j,original)
        with self.assertRaises(ValueError): result.values.setflags(write=True)
        with self.assertRaises(ValueError): result.valid.setflags(write=True)
        zero=optical_conductivity(t,j,np.zeros(64),**args)
        self.assertFalse(zero.valid.any())
        for key,value in [('current_convention','matter'),('field_kind','external'),('units','SI'),('drive_rtol',1)]:
            with self.subTest(key=key),self.assertRaises(ValueError):
                optical_conductivity(t,j,e,**dict(args,**{key:value}))
        with self.assertRaises(ValueError): optical_conductivity(t,j[:,None],e,**args)
        bad=e.copy();bad[10]=np.nan
        with self.assertRaises(ValueError): optical_conductivity(t,j,bad,**args)
        # A broadband impulse has no weak-drive bins; sigma remains constant.
        impulse=np.zeros(64);impulse[0]=1
        broad=optical_conductivity(t,2*impulse,impulse,**args)
        np.testing.assert_allclose(broad.values,2,atol=1e-14)
        self.assertTrue(broad.valid.all())

        doubled=optical_conductivity(t,4*impulse,2*impulse,**args)
        check=optical_linearity(broad,doubled,rtol=1e-12)
        self.assertTrue(check['accepted'])
        altered=optical_conductivity(t,6*impulse,2*impulse,**args)
        self.assertFalse(optical_linearity(broad,altered,rtol=.01)['accepted'])
        rule=optical_sum_rule(broad,expected_weight=2*broad.omega[-1],normalization='synthetic constant sigma=2')
        self.assertAlmostEqual(rule['captured_fraction'],1.)
        self.assertFalse(rule['certified'])


class TensorOpticalTests(unittest.TestCase):
    def test_mixed_drives_tensor_and_rank_deficiency(self):
        from salmonpy.analysis import dielectric_tensor
        t=np.arange(128)*.1
        drive=np.array([[1,.2,.1],[.1,1,.2],[.3,.1,1.]])
        sigma=np.array([[.3,.05,0],[.05,.4,.02],[0,.02,.5]])
        e=np.zeros((3,len(t),3));e[:,0,:]=drive.T
        j=np.einsum('ij,rtj->rti',sigma,e)
        args=dict(current_convention='charge',field_kind='total',units='atomic_gaussian',window='none')
        out=dielectric_tensor(t,j,e,**args)
        expected=np.eye(3)[None]+4j*np.pi*sigma[None]/out.omega[1:,None,None]
        np.testing.assert_allclose(out.values[1:],expected,rtol=1e-13,atol=1e-13)
        self.assertFalse(out.valid[0].any())
        self.assertFalse(out.values.flags.writeable)
        bad=e.copy();bad[2]=bad[1]
        result=dielectric_tensor(t,j,bad,**args)
        self.assertFalse(result.valid.any())
        self.assertTrue(np.isnan(result.values).all())
        for invalid in (e[:2], np.full_like(e,np.nan)):
            with self.assertRaises(ValueError):dielectric_tensor(t,j,invalid,**args)


class WaveplateRotationTests(unittest.TestCase):
    def test_rotated_analytic_phasors(self):
        from salmonpy.analysis import waveplate_diagnostic
        omega=1.;length=432.;indices=(1.5,2.)
        transmission=[]
        for n in indices:
            phase=np.exp(1j*n*omega*length/137.035999084)
            transmission.append(4*n*phase/((n+1)**2-(n-1)**2*phase**2))
        t=np.linspace(0,20*np.pi,4001)
        principal=np.zeros((len(t),3))
        principal[:,:2]=np.real(np.exp(-1j*t[:,None])*np.array(transmission)[None,:])
        reference=None
        for angle in (0.,np.pi/6):
            cs,sn=np.cos(angle),np.sin(angle)
            rotation=np.array([[cs,-sn,0],[sn,cs,0],[0,0,1.]])
            result=waveplate_diagnostic(t,principal@rotation.T,omega=omega,indices=indices,
                        thickness=length,units='atomic_gaussian',angle=angle)
            self.assertLess(result['stokes_error_norm'],1e-13)
            if reference is None:reference=result
            else:self.assertAlmostEqual(result['stokes_normalized'][2],reference['stokes_normalized'][2],places=13)


class LongWavelengthTests(unittest.TestCase):
    def test_a4_scalar_correction_and_explicit_absorptive_extension(self):
        from salmonpy.analysis import OpticalSpectrum, long_wavelength_correction
        w=np.linspace(0.,1.,101)
        true_eps=5.+2*w*w+1j*.1*w
        true_sigma=-1j*w*(true_eps-1)/(4*np.pi)
        data=true_sigma.copy();bad=(w>0)&(w<.2)
        data[bad]+=.01+1j*.02/(w[bad]+.01)
        cfg=dict(observable='optical_conductivity',units='atomic_gaussian',sign='exp(+i omega t)',current_convention='charge',field_kind='total')
        source=OpticalSpectrum(w,data,np.ones(w.shape,dtype=bool),cfg)
        args=dict(fit_range=(.3,.6),cutoff=.2,system='gapped_insulator',units='atomic_gaussian')
        result=long_wavelength_correction(source,**args)
        np.testing.assert_allclose(result.dielectric.values[bad].real,true_eps[bad].real,atol=1e-12)
        np.testing.assert_array_equal(result.conductivity.values.real,data.real)
        np.testing.assert_array_equal(result.conductivity.values[w>=.2],data[w>=.2])
        self.assertAlmostEqual(result.report['components'][0]['static_dielectric_real'],5.,places=12)
        self.assertFalse(result.dielectric.valid[0])
        full=long_wavelength_correction(source,absorption='analytic',**args)
        np.testing.assert_allclose(full.dielectric.values[w<.2],true_eps[w<.2],atol=1e-12)
        self.assertTrue(full.dielectric.valid[0])
        self.assertEqual(full.conductivity.values[0],0j)
        self.assertFalse(full.dielectric.values.flags.writeable)
        np.testing.assert_array_equal(source.values,data)
        with self.assertRaises(ValueError):long_wavelength_correction(source,**dict(args,system='metal'))

    def test_tensor_identity_masks_and_insufficient_fit(self):
        from salmonpy.analysis import OpticalSpectrum,long_wavelength_correction
        w=np.linspace(0.,1.,101)
        data=np.zeros((len(w),3,3),complex)
        for i in range(3):data[:,i,i]=i+2.+w*w+.02j*w
        data[:,0,1]=data[:,1,0]=.2+.3*w*w
        data[:,0,2]=complex(float('nan'),float('nan'))
        valid=np.ones(data.shape,bool);valid[0]=False;valid[:,0,2]=False
        original=data.copy()
        cfg=dict(observable='dielectric_tensor',units='atomic_gaussian',sign='exp(+i omega t)')
        source=OpticalSpectrum(w,data,valid,cfg)
        with self.assertWarns(RuntimeWarning):
            out=long_wavelength_correction(source,fit_range=(.3,.6),cutoff=.2,system='gapped_insulator',units='atomic_gaussian')
        np.testing.assert_array_equal(out.dielectric.values[w>=.2],original[w>=.2])
        self.assertFalse(out.dielectric.valid[:,0,2].any())
        k=np.flatnonzero(w==.1)[0]
        np.testing.assert_allclose(out.conductivity.values[k,0,1],-1j*.1*(.2+.3*.1**2)/(4*np.pi),atol=1e-14)
        np.testing.assert_array_equal(source.values,original)
        with self.assertWarns(RuntimeWarning):
            unchanged=long_wavelength_correction(source,fit_range=(.301,.302),system='gapped_insulator',units='atomic_gaussian')
        np.testing.assert_array_equal(unchanged.dielectric.values,source.values)


class TrajectoryDiagnosticsTests(unittest.TestCase):
    def test_prefix_windows_and_shared_frequency_transform(self):
        from salmonpy.analysis import optical_sensitivity
        t=.7+np.arange(64)*.1
        e=np.zeros(64);e[0]=1.;j=2*e
        args=dict(frequencies=[0.,.7,1.3],sample_counts=[32,64],windows=['none','cubic'],
                  current_convention='charge',field_kind='total',units='atomic_gaussian')
        out=optical_sensitivity(t,j,e,**args)
        self.assertEqual(len(out.spectra),4)
        for item in out.spectra:
            np.testing.assert_allclose(item.values,2.,atol=1e-13)
            self.assertFalse(item.values.flags.writeable)
        self.assertEqual(out.reference_index,2)
        self.assertFalse(out.report['certified'])
        from salmonpy.analysis import optical_linearity
        self.assertTrue(optical_linearity(out.spectra[0],out.spectra[0],rtol=0.)['accepted'])
        tail=np.exp(-.3*(t-t[0]))
        changed=optical_sensitivity(t,tail,e,**args)
        for index,item in enumerate(changed.spectra):
            count=args['sample_counts'][index//2]
            u=np.linspace(0.,1.,count)
            win=np.ones(count) if index%2==0 else 1.-3*u*u+2*u*u*u
            expected=np.array([np.sum(win*tail[:count]*np.exp(1j*w*t[:count]))/np.exp(1j*w*t[0]) for w in args['frequencies']])
            np.testing.assert_allclose(item.values,expected,rtol=3e-13,atol=3e-13)
        self.assertGreater(changed.report['cases'][0]['relative_rms_difference'],.01)
        aliased=optical_sensitivity(t,j,e,**dict(args,frequencies=[100.]))
        self.assertFalse(aliased.spectra[0].valid.any())
        with self.assertRaises(ValueError):optical_sensitivity(t[::-1],j,e,**args)
        silent=optical_sensitivity(t,j,e,**dict(args,windows=['hann']))
        self.assertFalse(silent.spectra[0].valid.any())
        for bad in (float('nan'),float('inf')):
            broken=j.copy();broken[5]=bad
            with self.assertRaises(ValueError):optical_sensitivity(t,broken,e,**args)
        with self.assertRaises(ValueError):optical_sensitivity(t,j,e,**dict(args,sample_counts=[65]))

    def test_weighted_k_reduction_sign_cancellation_and_native_spectra(self):
        from salmonpy.outputs import DecomposedCurrent
        from salmonpy.analysis import k_current_contributions
        t=np.arange(16)*.2
        raw=np.zeros((16,3,3));raw[:,0,0]=2.;raw[:,1,0]=1.;raw[:,2,0]=-2.
        source=DecomposedCurrent(np.arange(16),np.array([[3,1,1],[3,2,1],[7,1,1]]),raw,('a.u.',)*3,())
        out=k_current_contributions(source,dt=.2,current_convention='charge',units='atomic_gaussian',
                                    frequencies=[0.,1.],window='none',total_current=raw.sum(axis=1))
        np.testing.assert_array_equal(out.k_ids,[3,7]);np.testing.assert_allclose(out.current[:,0,0],3.)
        np.testing.assert_allclose(out.current[:,1,0],-2.);np.testing.assert_allclose(out.total[:,0],1.)
        np.testing.assert_allclose(out.time,t);self.assertAlmostEqual(out.report['components'][0]['relative_rms_to_total'],3.)
        self.assertEqual(out.report['total_current_comparison']['relative_rms_error'],0.)
        expected=.2*np.exp(1j*t).sum()*3
        self.assertAlmostEqual(out.spectrum.values[1,0,0].real,expected.real,places=12)
        neg=k_current_contributions(source,dt=.2,current_convention='matter',units='atomic_gaussian')
        np.testing.assert_array_equal(neg.total,-out.total)
        self.assertFalse(out.current.flags.writeable);self.assertFalse(out.time.flags.writeable)
        np.testing.assert_array_equal(source.current,raw)
        for invalid in (float('nan'),float('inf')):
            bad=raw.copy();bad[2,1,1]=invalid
            bad_source=DecomposedCurrent(source.step,source.indices,bad,source.units,())
            with self.assertRaises(ValueError):k_current_contributions(bad_source,dt=.2,current_convention='charge',units='atomic_gaussian')
        huge_steps=DecomposedCurrent(np.arange(16,dtype=np.int64)+2**60,source.indices,raw,source.units,())
        with self.assertRaises(ValueError):k_current_contributions(huge_steps,dt=.2,current_convention='charge',units='atomic_gaussian')
        duplicate=DecomposedCurrent(source.step,np.array([[3,1,1],[3,1,1],[7,1,1]]),raw,source.units,())
        with self.assertRaises(ValueError):k_current_contributions(duplicate,dt=.2,current_convention='charge',units='atomic_gaussian')
