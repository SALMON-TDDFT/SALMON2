import tempfile
from pathlib import Path
import unittest
import numpy as np
from salmonpy.field_sampling import SampledField


class SolidFieldTests(unittest.TestCase):
    def test_uniform_quadrupole_and_curl_residual(self):
        from salmonpy.near_field import ElectricPotentialSeries
        points=np.array([[x,y,z] for x in (-1.,0.,1.) for y in (-1.,0.,1.) for z in (-1.,0.,1.)])
        gradient=np.array([[.02,.003,0.],[.003,-.01,.004],[0.,.004,-.01]])
        E=points@gradient+np.array([.01,-.02,.03])
        samples=SampledField('E',[0.,.02,.04],points,np.array([E,2*E,3*E]),'E:au')
        result=ElectricPotentialSeries.from_samples(samples,center=(0,0,0),radius=2.,lmax=2)
        probes=np.array([[.13,-.25,.41],[0.,0.,0.]])
        np.testing.assert_allclose(result.electric_field(probes,step_index=0),probes@gradient+[.01,-.02,.03],atol=1.e-14)
        expected=probes@np.array([.01,-.02,.03])+.5*np.sum((probes@gradient)*probes,axis=1)
        np.testing.assert_allclose(result.potential(probes,step_index=0),expected,atol=1.e-14)
        self.assertLess(result.relative_residual.max(),1.e-12)
        difference=ElectricPotentialSeries.from_samples(samples,center=(0,0,0),radius=2.,lmax=2,reference=samples)
        np.testing.assert_array_equal(difference.coefficients,0.)
        self.assertEqual(difference.field_component,'difference')
        with tempfile.TemporaryDirectory() as tmp:
            output=Path(tmp)/'salmonpy-ylm.dat';result.write(output)
            with self.assertRaises(FileExistsError):result.write(output)
        curl=np.stack((-points[:,1],points[:,0],np.zeros(len(points))),axis=1)
        with self.assertWarnsRegex(UserWarning,'residual'):
            bad=ElectricPotentialSeries.from_samples(SampledField('E',[0.,.02],points,np.array([curl,curl]),'E:au'),center=(0,0,0),radius=2.,lmax=2)
        self.assertGreater(bad.relative_residual[0],.9)

    def test_octupole_normalization_and_high_order_roundtrip(self):
        from salmonpy.near_field import ElectricPotentialSeries
        points=np.random.default_rng(73).uniform(-.7,.7,(120,3))
        coefficients=np.zeros((2,15));coefficients[:,11]=1.
        series=ElectricPotentialSeries((0,0,0),1.,3,[0.,.02],coefficients,[0.,0.])
        x,y,z=points.T;norm=np.sqrt(7/(4*np.pi))
        expected=.5*norm*(5*z**3-3*z*(x*x+y*y+z*z))
        np.testing.assert_allclose(series.potential(points,step_index=0),expected,atol=4.e-16)
        coefficients=np.random.default_rng(74).normal(size=(2,80))*.01
        known=ElectricPotentialSeries((0,0,0),1.,8,[0.,.02],coefficients,[0.,0.])
        fields=np.array([known.electric_field(points,step_index=i) for i in range(2)])
        fitted=ElectricPotentialSeries.from_samples(SampledField('E',[0.,.02],points,fields,'E:au'),center=(0,0,0),radius=1.,lmax=8)
        np.testing.assert_allclose(fitted.coefficients,coefficients,atol=1.e-14,rtol=1.e-12)

    def test_eh_fit_center_and_multipolar_potentials(self):
        from salmonpy.near_field import NearFieldSeries
        center=np.array([2.,-1.,.5])
        offsets=np.array([[x,y,z] for x in (-1.,0.,1.) for y in (-1.,0.,1.) for z in (-1.,0.,1.)])
        points=offsets+center
        gradient=np.array([[.01,.02,0.],[-.03,-.005,.004],[0.,.001,-.005]])
        electric=offsets@gradient.T+[.01,.02,-.03]
        magnetic=np.broadcast_to([0.,0.,.04],electric.shape)
        E=SampledField('E',[0.,.02],points,np.array([electric,electric]),'E:au')
        H=SampledField('H',[0.,.02],points,np.array([magnetic,magnetic]),'H:au')
        result=NearFieldSeries.from_samples(E,H,lmax=1,radius=2.)
        np.testing.assert_allclose(result.center,center,rtol=0,atol=1.e-15)
        self.assertLess(result.relative_residual.max(),1.e-13)
        probe=center+np.array([[.2,-.3,.4],[-.1,.5,.2]])
        scalar,vector=result.potentials(probe,step_index=0)
        r=probe-center
        expected=r@np.array([.01,.02,-.03])+.5*np.sum(r*(r@gradient.T),axis=1)
        np.testing.assert_allclose(scalar,expected,atol=1.e-15)
        np.testing.assert_allclose(vector,-np.cross(r,[0.,0.,.04])/(2*137.035999084),atol=1.e-17)
        left=SampledField('H',[-.01,.01],points,np.zeros((2,len(points),3)),'H:au')
        right=SampledField('H',[.01,.03],points,np.array([2*magnetic,2*magnetic]),'H:au')
        aligned=NearFieldSeries.from_yee_samples(E,left,right,center=center,radius=2.,lmax=1)
        np.testing.assert_allclose(aligned.H_coefficients,result.H_coefficients,atol=1.e-15)
        sampled_E,sampled_H=result.fields(probe,step_index=0)
        np.testing.assert_allclose(sampled_E,r@gradient.T+[.01,.02,-.03],atol=1.e-15)
        np.testing.assert_allclose(sampled_H,np.broadcast_to([0.,0.,.04],sampled_H.shape),atol=1.e-15)
        # The molecular center is deliberately not one of the FDTD sample positions.
        off_grid = center + np.array([.17, -.23, .31])
        shifted = NearFieldSeries.from_samples(E, H, center=off_grid, radius=2., lmax=1)
        self.assertFalse(np.any(np.all(points == off_grid, axis=1)))
        np.testing.assert_array_equal(shifted.center, off_grid)
        shifted_E, shifted_H = shifted.fields(probe, step_index=0)
        np.testing.assert_allclose(shifted_E, sampled_E, atol=1.e-15)
        np.testing.assert_allclose(shifted_H, sampled_H, atol=1.e-15)
        single=NearFieldSeries.from_samples(E,H,center=center,radius=2.,lmax=1)
        np.testing.assert_allclose(single.frame(0).E_coefficients,result.frame(0).E_coefficients,atol=1.e-15)
        from salmonpy.near_field import NearFieldFrame,_kernel,_ptr
        frame=NearFieldFrame.from_samples(SampledField('E',[0.],points,[electric],'E:au'),
            SampledField('H',[0.],points,[magnetic],'H:au'),step_index=0,center=center,radius=2.,lmax=1)
        self.assertLess(frame.relative_residual.max(),1.e-13)
        for bad in (float('nan'),float('inf')):
            field=np.stack((electric,magnetic));field[1,0,0]=bad
            geometry=np.array([*center,2.,0.,0.,0.]);ce=np.full((4,3),7.);ch=np.full((4,3),8.);stats=np.full(2,9.)
            status=_kernel('sp_near_fit')(len(points),1,_ptr(points),_ptr(field),_ptr(geometry),0,0,_ptr(ce),_ptr(ch),_ptr(stats))
            self.assertNotEqual(status,0);np.testing.assert_array_equal(ce,7.);np.testing.assert_array_equal(ch,8.)
        with tempfile.TemporaryDirectory() as tmp:
            frame.write_live_input(Path(tmp)/'live.dat',max_steps=10,dt=.02)
            path=result.write(Path(tmp)/'salmonpy-near-eh.dat')
            self.assertTrue(path.read_text().startswith('SALMONPY_NEAR_EH_V3'))

    def test_radial_terms_at_off_grid_positions(self):
        from salmonpy.near_field import NearFieldSeries
        rng=np.random.default_rng(952)
        center=np.array([.17,-.23,.31])
        points=rng.uniform(-2.,2.,(120,3))+center
        probes=rng.uniform(-1.5,1.5,(40,3))+center
        def fields(p):
            r=p-center
            radial=np.sum(r*r,axis=1)
            return radial[:,None]*np.array([.01,.02,-.03]),radial[:,None]*np.array([.003,-.004,.002])
        E,H=fields(points)
        fitted=NearFieldSeries.from_samples(
            SampledField('E',[0.,.02],points,[E,E],'E:au'),
            SampledField('H',[0.,.02],points,[H,H],'H:au'),center=center,radius=3.,lmax=2)
        self.assertEqual(fitted.E_coefficients.shape,(2,10,3))
        pe,ph=fitted.fields(probes,step_index=0)
        expected_E,expected_H=fields(probes)
        np.testing.assert_allclose(pe,expected_E,atol=2.e-15)
        np.testing.assert_allclose(ph,expected_H,atol=2.e-15)
        scalar,vector=fitted.potentials(probes,step_index=0)
        r=probes-center
        np.testing.assert_allclose(scalar,np.sum(r*expected_E,axis=1)/3,atol=2.e-15)
        np.testing.assert_allclose(vector,-np.cross(r,expected_H)/(4*137.035999084),atol=2.e-17)
        # Independent off-sample checks through total degree8, including r^4/r^6/r^8.
        points=rng.uniform(-1.,1.,(400,3))+center
        probes=rng.uniform(-.8,.8,(50,3))+center
        for degree in (4,6,8):
            r=points-center;v=np.sum(r*r,axis=1)**(degree//2)
            values=v[:,None]*np.array([.01,.02,-.03])
            zero=np.zeros_like(values)
            fit=NearFieldSeries.from_samples(
                SampledField('E',[0.,.02],points,[values,values],'E:au'),
                SampledField('H',[0.,.02],points,[zero,zero],'H:au'),
                center=center,radius=2.,lmax=degree)
            r=probes-center;expected=np.sum(r*r,axis=1)[:,None]**(degree//2)*[.01,.02,-.03]
            np.testing.assert_allclose(fit.fields(probes,step_index=0)[0],expected,atol=2.e-14)
            np.testing.assert_allclose(fit.potentials(probes,step_index=0)[0],
                np.sum(r*expected,axis=1)/(degree+1),atol=2.e-14)

    def test_exponential_envelope_at_independent_positions(self):
        from salmonpy.near_field import NearFieldSeries
        rng=np.random.default_rng(1953)
        center=np.array([.17,-.23,.31]);decay=np.array([.1,-.2,.7])
        points=rng.uniform(-1.,1.,(100,3))+center
        probes=rng.uniform(-.8,.8,(40,3))+center
        E=np.exp(-(points-center)@decay)[:,None]*[.01,.02,-.03]
        H=np.exp(-(points-center)@decay)[:,None]*[.003,-.004,.002]
        fitted=NearFieldSeries.from_samples(
            SampledField('E',[0.,.02],points,[E,E],'E:au'),
            SampledField('H',[0.,.02],points,[H,H],'H:au'),
            center=center,radius=2.,lmax=0,decay_vector=decay)
        expected=np.exp(-(probes-center)@decay)[:,None]
        e,h=fitted.fields(probes,step_index=0)
        np.testing.assert_allclose(e,expected*[.01,.02,-.03],atol=1.e-15)
        np.testing.assert_allclose(h,expected*[.003,-.004,.002],atol=1.e-15)
        plain_a=NearFieldSeries.from_samples(
            SampledField('E',[0.,.02],points,[E,E],'E:au'),
            SampledField('H',[0.,.02],points,[H,H],'H:au'),
            center=center,radius=2.,lmax=2,residual_warning=1.)
        plain_b=NearFieldSeries.from_samples(
            SampledField('E',[0.,.02],points,[E,E],'E:au'),
            SampledField('H',[0.,.02],points,[H,H],'H:au'),
            center=center,radius=2.,lmax=2,residual_warning=1.,decay_vector=(0.,0.,0.))
        np.testing.assert_array_equal(plain_a.E_coefficients,plain_b.E_coefficients)
        for left,right in zip(plain_a.potentials(probes,step_index=0),plain_b.potentials(probes,step_index=0)):
            np.testing.assert_array_equal(left,right)
        np.testing.assert_array_equal(fitted.frame(0).decay_vector,decay)
        self.assertFalse(fitted.decay_vector.flags.writeable)
        from salmonpy.near_field import couple_path_samples
        r=np.array([.2,-.3,.4])
        E=np.broadcast_to([.01,.02,-.03],(3,3));H=np.broadcast_to([.003,-.004,.002],(3,3))
        scalar,vector=couple_path_samples(r,[0.,.5,1.],[1/6,4/6,1/6],E,H)
        self.assertAlmostEqual(scalar,float(r@E[0]),places=15)
        np.testing.assert_allclose(vector,-np.cross(r,H[0])/(2*137.035999084),atol=1.e-18)
        self.assertFalse(vector.flags.writeable)
        from salmonpy.errors import InputError
        with self.assertRaisesRegex(InputError,'parameters'):
            couple_path_samples(r,[-.1,.5,1.],[1/6,4/6,1/6],E,H)


        with tempfile.TemporaryDirectory() as tmp:
            path=fitted.write(Path(tmp)/'field.dat')
            self.assertEqual(path.read_text().splitlines()[0],'SALMONPY_NEAR_EH_V3')

    def test_direct_grid_off_grid_translation_and_paths(self):
        from salmonpy.near_field import GridFieldSeries
        times=[0.,.02]
        origin=np.array([-3.,-3.,-3.]);spacing=np.array([1.,1.,1.])
        # C layout: time,z,y,x,component. Independent analytic affine field.
        positions=np.array([[[[x,y,z] for x in range(7)] for y in range(7)] for z in range(7)],dtype=float)+origin
        E=positions*.01+[.01,.02,-.03];H=np.broadcast_to([.003,-.004,.002],E.shape)
        grid=GridFieldSeries(origin,spacing,times,[E,E],[H,H],molecule_position=(.17,-.23,.31))
        points=np.array([[.2,-.3,.4],[-.1,.5,.2]])
        e,h=grid.fields(points,step_index=0)
        np.testing.assert_allclose(e,(points+[.17,-.23,.31])*.01+[.01,.02,-.03],atol=1e-15)
        scalar,vector=grid.potentials(points,step_index=0)
        expected=np.sum(points*(np.array([.17,-.23,.31])*.01+[.01,.02,-.03]),axis=1)+.005*np.sum(points*points,axis=1)
        np.testing.assert_allclose(scalar,expected,atol=1e-15)
        np.testing.assert_allclose(vector,-np.cross(points,[.003,-.004,.002])/(2*137.035999084),atol=1e-17)
        with self.assertRaises(ValueError):grid.fields([[100,0,0]],step_index=0)
        order=np.random.default_rng(1954).permutation(len(positions.reshape(-1,3)))
        points=positions.reshape(-1,3)[order]
        e=E.reshape(-1,3)[order];h=H.reshape(-1,3)[order]
        packed=GridFieldSeries.from_samples(
            SampledField('E',times,points,[e,e],'E:au'),SampledField('H',times,points,[h,h],'H:au'),
            counts=(7,7,7),origin=origin,spacing=spacing,molecule_position=(.17,-.23,.31))
        np.testing.assert_array_equal(packed.E,grid.E)
        half_left=SampledField('H',[-.01,.01],points,np.zeros((2,len(points),3)),'H:au')
        half_right=SampledField('H',[.01,.03],points,[2*h,2*h],'H:au')
        aligned=GridFieldSeries.from_yee_samples(SampledField('E',times,points,[e,e],'E:au'),half_left,half_right,
            counts=(7,7,7),origin=origin,spacing=spacing,molecule_position=(.17,-.23,.31))
        np.testing.assert_array_equal(aligned.H,grid.H)
        structured=GridFieldSeries.from_yee_arrays(origin,spacing,times,[E,E],np.zeros_like(grid.H),2*grid.H,
            H_left_times=[-.01,.01],H_right_times=[.01,.03],molecule_position=(.17,-.23,.31))
        np.testing.assert_array_equal(structured.H,grid.H)

        with tempfile.TemporaryDirectory() as tmp:
            p=grid.write(Path(tmp)/'grid.dat')
            self.assertTrue(p.read_text().startswith('SALMONPY_NEAR_GRID_V1'))


class ElectricPathResponseTests(unittest.TestCase):
    def test_spatial_work_and_integrated_dipole(self):
        from salmonpy.near_field import GridFieldSeries
        rng=np.random.default_rng(1955)
        E=rng.normal(size=(2,9,9,9,3))*.01;H=np.zeros_like(E)
        origins=np.full((2,3,3),-2.)
        origins[0,0,0]+=.03;origins[0,1,1]-=.04;origins[0,2,2]+=.02
        field=GridFieldSeries(origins,(.5,)*3,[0.,.02],E,H,molecule_position=(.17,-.23,.31),center=(.2,-.3,.1))
        points=rng.uniform(-1.,1.,(37,3));volumes=rng.uniform(.001,.003,37)
        old=-rng.uniform(.1,1.,37);new=old+rng.uniform(-.01,.01,37)
        J=field.charge_current(points,old,new,volume_elements=volumes,dt=.02)
        V,_=field.potentials(points,step_index=0)
        work=np.sum(E[0]*J)*.125*.02
        np.testing.assert_allclose(work,np.sum((new-old)*volumes*V),rtol=2e-13,atol=1e-18)
        P=field.charge_polarization(points,new,volume_elements=volumes)
        np.testing.assert_allclose(np.sum(P,axis=(0,1,2))*.125,np.sum(new[:,None]*volumes[:,None]*(points-field.center),axis=0),atol=2e-17)
        Pold=field.charge_polarization(points,old,volume_elements=volumes)
        np.testing.assert_allclose((P-Pold)/.02,J,atol=2e-16)
        self.assertFalse(J.flags.writeable)
        with self.assertRaises(ValueError):field.charge_current([[100.,0.,0.]],[1.],[2.],volume_elements=[1.],dt=.02)


class GridReferenceTests(unittest.TestCase):
    def test_fortran_reference_subtraction_and_geometry(self):
        from salmonpy.near_field import GridFieldSeries,_kernel
        from salmonpy.errors import InputError
        E=np.ones((2,2,2,2,3));H=E*.2
        total=GridFieldSeries((-40.,)*3,(80.,)*3,[0.,.01],E,H,time_origin=.03)
        own=GridFieldSeries((-40.,)*3,(80.,)*3,[0.,.01],E*.3,H*.5,time_origin=.03)
        result=total.difference(own)
        np.testing.assert_array_equal(result.E,E-E*.3);np.testing.assert_array_equal(result.H,H-H*.5)
        self.assertEqual(result.time_origin,.03)
        with self.assertRaises(InputError):total.difference(GridFieldSeries((-41.,)*3,(80.,)*3,[0.,.01],E,H,time_origin=.03))
        a=np.array([np.finfo(float).max]);b=-a;out=np.array([-7.])
        self.assertEqual(_kernel('sp_grid_difference')(1,a.ctypes.data,b.ctypes.data,out.ctypes.data),4)
        np.testing.assert_array_equal(out,[-7.])


class OrbitalMagneticReturnTests(unittest.TestCase):
    def test_yee_curl_and_electric_current_counted_once(self):
        from salmonpy.isolated_tddft import OrbitalMagneticResponse
        from salmonpy.errors import InputError
        origin=np.stack((np.eye(3)*.5,.5-np.eye(3)*.5))
        M=np.zeros((7,7,7,3));M[2:4,2:4,2:4,2]=.001
        result=OrbitalMagneticResponse(2,.006,origin,np.ones(3),M)
        Jm=result.current_density();self.assertGreater(np.max(np.abs(Jm)),0.)
        Je=np.full_like(M,.2);combined=result.total_current(Je)
        np.testing.assert_array_equal(combined,Je+Jm)
        np.testing.assert_array_equal(Je,.2);self.assertFalse(combined.flags.writeable)
        previous=OrbitalMagneticResponse(1,.003,origin,np.ones(3),M/3)
        np.testing.assert_allclose(result.current_density(previous=previous),Jm*(2./3.),rtol=1e-14,atol=1e-17)
        with self.assertRaises(InputError):result.current_density(previous=result)
        div=np.zeros(M.shape[:3])
        for d,axis in enumerate((2,1,0)):
            div+=Jm[...,d]-np.roll(Jm[...,d],1,axis=axis)
        np.testing.assert_allclose(div,0.,atol=1e-16)
        with self.assertRaises(ValueError):result.magnetization[0,0,0,0]=1.
        bad=OrbitalMagneticResponse(2,.006,np.zeros((2,3,3)),np.ones(3),M)
        with self.assertRaises(InputError):bad.current_density()
        M[-1,-1,-1]=1.
        edge=OrbitalMagneticResponse(2,.006,origin,np.ones(3),M)
        with self.assertRaises(InputError):edge.current_density()

    def test_public_readout_requires_runtime_evidence(self):
        from salmonpy.isolated_tddft import IsolatedTDDFTState
        from salmonpy.errors import InputError
        state=object.__new__(IsolatedTDDFTState)
        with self.assertRaises(InputError):state.orbital_magnetic_response()
