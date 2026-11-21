import hashlib
from pathlib import Path
import tempfile
import unittest
import re
import numpy as np
from unittest.mock import patch
from salmonpy.fdtd import FDTD, Boundaries, Dielectric, Impulse, SalmonExecutable, PlaneWave
from salmonpy.geometry import Grid, Geometry, Shape
from salmonpy import LorentzDrude
from salmonpy.fields import AcosPulse

class FDTDTests(unittest.TestCase):
    def calculation(self):
        return FDTD(Geometry(Grid((24,24,24),(2,2,2)),(Shape('ellipsoid',(8,8,8)),)),
                    materials=(LorentzDrude(plasma_frequency=1.,strength=[1.],resonance=[0.],damping=[.1]),),dt=.001,steps=10)

    def test_native_plane_tiles_colocation_and_incomplete_data(self):
        import struct
        from salmonpy.fdtd import PlaneMonitor,read_native_planes
        from dataclasses import replace
        for axis in range(3):
            with self.subTest(axis=axis),tempfile.TemporaryDirectory() as tmp:
                shape=[8,8,8];shape[axis]=24
                boundary=['periodic']*6;boundary[2*axis:2*axis+2]=['pml']*2
                pol=[0,0,0];pol[(axis+1)%3]=1
                calc=FDTD(Geometry(Grid(tuple(shape),(1.,1.,1.))),(),.001,2,
                    boundaries=Boundaries(tuple(boundary)),
                    source=PlaneWave(AcosPulse(.001,1.,.1,polarization=tuple(pol)),axis,-1.5),
                    planes=(PlaneMonitor(axis,.5),))
                self.assertIn('salmonpy-planes.dat',calc.input_assets())
                self.assertNotIn('salmonpy-planes.dat',replace(calc,planes=()).input_assets())
                tangent=[d for d in range(3) if d!=axis];dims=tuple(shape[d] for d in tangent)
                base=np.arange(3*np.prod(dims),dtype=float).reshape(*dims,3)
                e=np.stack([base*i for i in range(3)])
                h=np.stack([(base+10)*i for i in range(3)])
                ep=np.stack([(base+100)*i for i in range(3)])
                hp=np.stack([(base+200)*i for i in range(3)])
                lo=np.array([-n//2+1 for n in shape]);hi=lo+shape-1
                offset=int(np.argmin(abs(calc.geometry.grid.coordinates[axis]-.5)));idx=int(lo[axis]+offset)
                def write(root,ranks):
                    root.mkdir()
                    for rank in range(ranks):
                        local_lo=lo.copy();local_hi=hi.copy();a=rank*(dims[0]//ranks);b=(rank+1)*(dims[0]//ranks)
                        local_lo[tangent[0]]=lo[tangent[0]]+a;local_hi[tangent[0]]=lo[tangent[0]]+b-1
                        header=[16909060,1,axis+1,idx,idx-1,*lo,*hi,*local_lo,*local_hi,1,2,rank,0]
                        path=root/f'salmonpy_plane_0001_rank_{rank:08d}.bin'
                        with path.open('wb') as stream:
                            stream.write(b'SALMONPY_NATIVE_PLANES_V1'.ljust(32,b' '))
                            stream.write(struct.pack('<21i',*header));stream.write(struct.pack('<d',.001))
                            for step in range(3):
                                stream.write(struct.pack('<i',step))
                                for E,H in ((e,h),(ep,hp)):
                                    raw=np.concatenate((E[step,a:b],H[step,a:b]),axis=-1)
                                    stream.write(np.ascontiguousarray(raw.transpose(1,0,2),dtype='<f8').tobytes())
                root=Path(tmp);write(root/'serial',1);write(root/'split',2)
                full=read_native_planes(root/'serial',calc)[0]
                split=read_native_planes(root/'split',calc,ranks=2)[0]
                np.testing.assert_array_equal(full.E,split.E);np.testing.assert_array_equal(full.H,split.H)
                for k in range(3):
                    if k==axis:expected=.5*(e[:,:,:,k]+ep[:,:,:,k])
                    else:
                        d=tangent.index(k)+1
                        expected=.5*(e[:,:,:,k]+np.roll(e[:,:,:,k],1,axis=d))
                    np.testing.assert_array_equal(full.E[:,:,:,k],expected)
                    expected=np.zeros_like(h[:,:,:,k])
                    perpendicular=[d for d in range(3) if d!=k]
                    for i in (0,1):
                        for j in (0,1):
                            previous=(i if perpendicular[0]==axis else j if perpendicular[1]==axis else 0)
                            values=(hp if previous else h)[:,:,:,k]
                            for d,shift in zip(perpendicular,(i,j)):
                                if d!=axis:values=np.roll(values,shift,axis=tangent.index(d)+1)
                            expected+=values/4
                    np.testing.assert_array_equal(full.H[:,:,:,k],expected)
                np.testing.assert_array_equal(full.time,[0.,.001,.002])
                self.assertEqual(full.config['native_H_time'],'already centered at E time')
                self.assertFalse(full.E.flags.writeable)
                path=root/'split/salmonpy_plane_0001_rank_00000000.bin'
                path.write_bytes(path.read_bytes()[:-8])
                with self.assertRaisesRegex(ValueError,'truncated'):read_native_planes(root/'split',calc,ranks=2)
                serial=root/'serial/salmonpy_plane_0001_rank_00000000.bin'
                original=serial.read_bytes()
                for value in (float('nan'),float('inf')):
                    serial.write_bytes(original[:128]+struct.pack('<d',value)+original[136:])
                    with self.assertRaisesRegex(ValueError,'nonfinite'):
                        read_native_planes(root/'serial',calc)
                serial.write_bytes(original)
                invalid=list(boundary);invalid[2*tangent[0]:2*tangent[0]+2]=['pec','pec']
                with self.assertRaisesRegex(ValueError,'periodic'):
                    replace(calc,boundaries=Boundaries(tuple(invalid)))

    def test_plane_diagnostic_failure_preserves_completed_calculation(self):
        from dataclasses import replace
        from types import SimpleNamespace
        import json
        from salmonpy.fdtd import PlaneMonitor
        calculation=replace(self.calculation(),
            boundaries=Boundaries(('pml','pml','periodic','periodic','periodic','periodic')),
            planes=(PlaneMonitor(0,1.),))
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);exe=root/'native';exe.write_bytes(b'SALMONPY_NATIVE_PLANES_V1')
            engine=SalmonExecutable(exe,hashlib.sha256(exe.read_bytes()).hexdigest())
            def run(*args,**kw):
                kw['stdout'].write('FDTD end\nend SALMON\n');kw['stdout'].flush()
                (Path(kw['cwd'])/'obs1_at_point_rt.data').write_text(
                    '# 1:Time[a.u.] 2:E_z[a.u.]\n'+''.join(f'{n*.001} 0\n' for n in range(1,11)))
                return SimpleNamespace(returncode=0)
            with patch('salmonpy.fdtd.subprocess.run',side_effect=run),patch(
                    'salmonpy.fdtd.read_native_planes',side_effect=MemoryError('injected diagnostic resource failure')):
                with self.assertWarnsRegex(UserWarning,'FDTD completed; plane diagnostics unavailable'):
                    result=engine.run(calculation,root/'out')
            record=json.loads((result.directory/'salmonpy-run.json').read_text())
            self.assertEqual(record['status'],'complete')
            self.assertEqual(record['plane_observations_status'],'unavailable')
            self.assertEqual(len(result.observations[0].values),10)
            self.assertEqual(result.planes,())
            old=root/'old-runtime';old.write_bytes(b'old runtime without plane capability')
            old_engine=SalmonExecutable(old,hashlib.sha256(old.read_bytes()).hexdigest())
            with patch('salmonpy.fdtd.subprocess.run',side_effect=run),patch('salmonpy.fdtd.read_native_planes') as read:
                with self.assertWarnsRegex(UserWarning,'calculation continues'):
                    old_result=old_engine.run(calculation,root/'old-out')
                read.assert_not_called()
            self.assertEqual(old_result.planes,())
            self.assertEqual(len(old_result.observations[0].values),10)

    def test_input_and_fixed_settings(self):
        c=self.calculation();text=c.input_text()
        self.assertIn("theory = 'maxwell'",text)
        self.assertIn("unit_system = 'au'",text)
        self.assertIn('omega_p_ld(1) = 1',text)
        self.assertIn('obs_samp_em = 1',text)
        from dataclasses import replace
        periodic=Boundaries(('pml','pml','periodic','periodic','periodic','periodic'))
        mixed=replace(c,boundaries=periodic)
        self.assertIn("boundary_em(2,1) = 'periodic'",mixed.input_text())
        self.assertIn("yn_periodic = 'n'",mixed.input_text())
        with self.assertRaisesRegex(ValueError,'pairs'):
            Boundaries(('pml','pml','periodic','pec','pml','pml'))

        with self.assertRaises(AttributeError): c.dt=.5
        with self.assertRaises(ValueError): FDTD(c.geometry,c.materials,dt=1.,steps=10)

    def test_undefined_material_and_pml_contact(self):
        g=Geometry(Grid((24,24,24),(1,1,1)),(Shape('box',(23,23,23)),))
        with self.assertRaises(ValueError): FDTD(g,materials=(),dt=.001,steps=10)
        with self.assertRaises(ValueError): FDTD(g,materials=(LorentzDrude(plasma_frequency=1.,strength=[1.],resonance=[0.],damping=[.1]),),dt=.001,steps=10)
        FDTD(g,materials=(LorentzDrude(plasma_frequency=1.,strength=[1.],resonance=[0.],damping=[.1]),),dt=.001,steps=10,boundaries=Boundaries(('pec',)*6))

    def test_runner_failure_and_no_overwrite(self):
        with tempfile.TemporaryDirectory() as d:
            exe=Path(d)/'salmon';exe.write_bytes(b'not executed')
            runner=SalmonExecutable(exe,hashlib.sha256(exe.read_bytes()).hexdigest())
            target=Path(d)/'run'
            with patch('salmonpy.fdtd.subprocess.run') as run:
                run.return_value.returncode=0
                with self.assertRaises(RuntimeError): runner.run(self.calculation(),target)
                self.assertFalse(run.call_args.kwargs.get('shell',False))
            with self.assertRaisesRegex(FileExistsError,'Choose a new directory') as caught:
                runner.run(self.calculation(),target)
            self.assertEqual(caught.exception.filename,str(target.resolve()))
            self.assertEqual(caught.exception.errno,17)
            exe.write_bytes(b'changed')
            with self.assertRaises(ValueError): runner.run(self.calculation(),Path(d)/'other')
            self.assertFalse((Path(d)/'other').exists())

    def test_plane_wave_matches_common_A_convention(self):
        g=Geometry(Grid((32,32,32),(2,2,2)),(Shape('ellipsoid',(8,8,8)),))
        pulse=AcosPulse(.1,2.,4.,cep=.12)
        c=FDTD(g,(Dielectric(2.),),.001,10,source=PlaneWave(pulse,0,-12.))
        cep=float(re.search(r'phi_cep1 = (\S+)',c.input_text()).group(1))
        t=np.array([.4,1.3,2.7,3.5]);h=1e-6
        derivative=-(pulse.sample(t+h).values[:,2]-pulse.sample(t-h).values[:,2])/(2*h)
        theta1=np.pi/pulse.duration*(t-.5*pulse.duration)
        theta2=pulse.omega*(t-.5*pulse.duration)+2*np.pi*cep
        # Frozen formula from fdtd_eh::eh_calc_e_inc.
        legacy=-pulse.amplitude*(-np.pi/pulse.duration*np.sin(2*theta1)*np.cos(theta2)
                   -pulse.omega*np.cos(theta1)**2*np.sin(theta2))/pulse.omega
        np.testing.assert_allclose(legacy,derivative,atol=1e-10,rtol=1e-8)

    def test_direct_LD_config_is_copied(self):
        from salmonpy.ld import LDModelConfig
        f=[1.];m=LDModelConfig(1.,f,[0.],[.1]);base=self.calculation()
        c=FDTD(base.geometry,(m,),.001,10)
        f[0]=2.
        self.assertEqual(c.materials[0].strength,(1.,))

    def test_pml_interface_is_clear(self):
        g=Geometry(Grid((24,24,24),(1,1,1)),(Shape('box',(1,1,1),center=(-3.5,0,0)),))
        with self.assertRaisesRegex(ValueError,'intersects'):
            FDTD(g,(self.calculation().materials[0],),.001,10)

    def test_multimaterial_rejects_old_runtime(self):
        from dataclasses import replace
        import types
        with tempfile.TemporaryDirectory() as tmp:
            exe=Path(tmp)/'salmon';exe.write_bytes(b'fixture')
            engine=SalmonExecutable(exe,hashlib.sha256(exe.read_bytes()).hexdigest())
            c=self.calculation();c=replace(c,materials=c.materials*2)
            def old_run(*args,**kwargs):
                kwargs['stdout'].write('FDTD end\nend SALMON\n');kwargs['stdout'].flush()
                return types.SimpleNamespace(returncode=0)
            with patch('salmonpy.fdtd.subprocess.run',side_effect=old_run):
                with self.assertRaisesRegex(RuntimeError,'material-ID-independent'):
                    engine.run(c,Path(tmp)/'result')

    def test_ld_background_and_fastest_cell_cfl(self):
        from dataclasses import replace
        from salmonpy.ld import LDModelConfig
        base=self.calculation()
        cfg=LDModelConfig(1.,(1.,),(0.,),(.1,),epsilon_infinity=2.25)
        calc=replace(base,materials=(cfg,),source=PlaneWave(AcosPulse(.001,10.,2.),0,-12.),boundaries=Boundaries(('pec',)*6))
        self.assertEqual(calc.materials[0].epsilon_infinity,2.25)
        with self.assertWarnsRegex(RuntimeWarning,'raw data continue'):
            kicked=replace(base,materials=(cfg,))
        self.assertIsInstance(kicked.source,Impulse)
        with self.assertWarnsRegex(RuntimeWarning,'calibrated'):
            extrapolated=replace(calc,materials=(replace(cfg,valid_frequency_range=(.1,.3)),))
        self.assertEqual(extrapolated.materials[0].valid_frequency_range,(.1,.3))


        self.assertIn('epsilon_em(1) = 2.25',calc.input_text())
        # Slower matter cannot loosen the vacuum/PML timestep limit.
        with self.assertRaisesRegex(ValueError,'CFL'):
            replace(calc,dt=.009)
        fast=replace(calc,materials=(replace(cfg,epsilon_infinity=.25),))
        with self.assertRaisesRegex(ValueError,'CFL'):
            replace(fast,dt=.005)

    def test_mpi_grid_configuration(self):
        import salmonpy.fdtd as f
        self.assertTrue(hasattr(f,'FDTDMPI'),'FDTD MPI launch configuration missing')
        mpi=f.FDTDMPI((1,2,2))
        self.assertEqual(mpi.ranks,4)
        self.assertEqual(mpi.command('/site/salmon'),('mpiexec','-n','4','/site/salmon'))
        site=f.FDTDMPI((1,2,2),launcher=('site-launch',),rank_arguments=('--tasks={ranks}',))
        self.assertEqual(site.command('/site path/salmon'),('site-launch','--tasks=4','/site path/salmon'))
        implicit=f.FDTDMPI((1,2,2),launcher=('allocation-launch',),rank_arguments=())
        self.assertEqual(implicit.command('/site/salmon'),('allocation-launch','/site/salmon'))
        for args in ('-n',('bad\x00arg',)):
            with self.assertRaises(ValueError):f.FDTDMPI((1,1,1),rank_arguments=args)

        self.assertIn('nproc_rgrid = 1, 2, 2',mpi.input_text(self.calculation()))
        with self.assertRaises(ValueError): f.FDTDMPI((1,0,2))
        with self.assertRaises(ValueError): f.FDTDMPI((True,1,2))
        with self.assertRaises(ValueError): f.FDTDMPI((1,1))
        with self.assertRaises(ValueError): f.FDTDMPI((1,1,100)).input_text(self.calculation())

    def test_mpi_rejects_serial_and_nested_launch_before_output(self):
        import os,json
        from salmonpy.fdtd import FDTDMPI
        with tempfile.TemporaryDirectory() as tmp:
            exe=Path(tmp)/'salmon';exe.write_bytes(b'fixture')
            engine=SalmonExecutable(exe,hashlib.sha256(exe.read_bytes()).hexdigest())
            dest=Path(tmp)/'run'
            with self.assertRaisesRegex(ValueError,'MPI-enabled'):
                engine.run(self.calculation(),dest,mpi=FDTDMPI((1,2,2)))
            self.assertFalse(dest.exists())
            if os.name=='posix':
                engine=SalmonExecutable(exe,engine.sha256,mpi_enabled=True)
                with patch.dict(os.environ,{'OMPI_COMM_WORLD_RANK':'0'}):
                    with self.assertRaisesRegex(ValueError,'inside an MPI rank'):
                        engine.run(self.calculation(),dest,mpi=FDTDMPI((1,2,2)))
                self.assertFalse(dest.exists())

    def test_mpi_build_record_checks_capability_and_hash(self):
        import json
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);build=root/'build';build.mkdir();exe=build/'salmon';exe.write_bytes(b'fixture')
            record={'status':'built-and-audited','executable':str(exe),'sha256':hashlib.sha256(exe.read_bytes()).hexdigest(),'build_configuration':{'mpi_enabled':False}}
            (root/'build.json').write_text(json.dumps(record));cache=build/'CMakeCache.txt'
            cache.write_text('USE_MPI:BOOL=OFF\n')
            with self.assertRaisesRegex(ValueError,'not built with MPI'): SalmonExecutable.from_mpi_build(root/'build.json')
            cache.write_text('USE_MPI:BOOL=ON\n')
            with self.assertRaisesRegex(ValueError,'not built with MPI'): SalmonExecutable.from_mpi_build(root/'build.json')
            record['build_configuration']['mpi_enabled']=True
            (root/'build.json').write_text(json.dumps(record))
            self.assertTrue(SalmonExecutable.from_mpi_build(root/'build.json').mpi_enabled)
            exe.write_bytes(b'changed')
            with self.assertRaisesRegex(ValueError,'differs'): SalmonExecutable.from_mpi_build(root/'build.json')

    @unittest.skipUnless(__import__('os').name=='posix','POSIX MPI launcher')
    def test_mpi_interrupt_records_cancelled(self):
        import json
        from salmonpy.fdtd import FDTDMPI
        with tempfile.TemporaryDirectory() as tmp:
            exe=Path(tmp)/'salmon';exe.write_bytes(b'fixture')
            engine=SalmonExecutable(exe,hashlib.sha256(exe.read_bytes()).hexdigest(),mpi_enabled=True)
            dest=Path(tmp)/'run'
            with patch('salmonpy.fdtd._run_fdtd_mpi',side_effect=KeyboardInterrupt):
                with self.assertRaises(KeyboardInterrupt): engine.run(self.calculation(),dest,mpi=FDTDMPI((1,1,1)))
            self.assertEqual(json.loads((dest/'salmonpy-run.json').read_text())['status'],'cancelled')


    @unittest.skipUnless(__import__('os').name=='posix','POSIX MPI launcher')
    def test_initialized_mpi_rejected_even_without_launcher(self):
        import sys,types
        fake=types.SimpleNamespace(Is_initialized=lambda:True,Is_finalized=lambda:False)
        with tempfile.TemporaryDirectory() as tmp:
            exe=Path(tmp)/'salmon';exe.write_bytes(b'fixture')
            engine=SalmonExecutable(exe,hashlib.sha256(exe.read_bytes()).hexdigest(),mpi_enabled=True)
            with patch.dict(sys.modules,{'mpi4py.MPI':fake}), patch('salmonpy.fdtd.subprocess.run') as run:
                with self.assertRaisesRegex(ValueError,'MPI'):
                    engine.run(self.calculation(),Path(tmp)/'run')
                run.assert_not_called()


    def test_single_component_ld_tensor_uses_same_eh_input(self):
        from dataclasses import replace
        from salmonpy.ld import LDTensorModelConfig
        original=self.calculation()
        scalar=original.materials[0]
        tensor=scalar.tensor_config
        self.assertIsInstance(tensor,LDTensorModelConfig)
        self.assertEqual(tensor.independent_components,1)
        self.assertEqual(tensor.representation,'isotropic_diagonal')
        # Frozen slotted dataclasses use TypeError on some supported Python versions.
        with self.assertRaises((AttributeError,TypeError)):tensor.representation='full'
        self.assertEqual(tensor.representation,'isotropic_diagonal')
        converted=replace(original,materials=(tensor,))
        self.assertEqual(original.input_text(),converted.input_text())
        self.assertEqual(converted.tensor_materials,(tensor,))
        self.assertEqual(original.tensor_materials,(tensor,))
        from salmonpy.material_response import LDMaterialProvider
        standard=LorentzDrude.from_config(scalar)
        converted_model=LorentzDrude.from_config(tensor)
        provider=LDMaterialProvider(tensor)
        self.assertEqual(provider.identity,LDMaterialProvider(standard).identity)
        fields=np.array([[0.,0.,0.],[.01,.02,.03],[.02,0.,-.01]])
        first=standard.respond(fields,dt=.001)
        second=converted_model.respond(fields,dt=.001)
        np.testing.assert_array_equal(first.polarization,second.polarization)
        np.testing.assert_array_equal(first.current,second.current)


    def test_diagonal_tensor_asset_and_explicit_capability(self):
        from dataclasses import replace
        from salmonpy.ld import TensorLorentzDrude
        base=self.calculation()
        model=TensorLorentzDrude(plasma_frequency=.5,strength=[np.diag([.1,.2,.3])],resonance=[.4],damping=[.03],epsilon_infinity=np.diag([2.,3.,4.]))
        with self.assertWarns(RuntimeWarning):calculation=replace(base,materials=(model,))
        text=calculation.input_text()
        self.assertIn("media_type(1) = 'ld-diagonal'",text)
        assets=calculation.input_assets()
        self.assertEqual(set(assets),{'salmonpy-ld-diagonal.dat'})
        self.assertIn(hashlib.sha256(assets['salmonpy-ld-diagonal.dat'].encode()).hexdigest(),text)
        self.assertEqual(base.input_assets(),{})
        with tempfile.TemporaryDirectory() as directory:
            exe=Path(directory)/'old';exe.write_bytes(b'not an anisotropic runtime')
            runner=SalmonExecutable(exe,hashlib.sha256(exe.read_bytes()).hexdigest())
            target=Path(directory)/'not-created'
            with self.assertRaisesRegex(ValueError,'diagonal'):
                runner.run(calculation,target)
            self.assertFalse(target.exists())
        offdiag=TensorLorentzDrude(plasma_frequency=1.,strength=[[[.2,.05,0.],[.05,.3,0.],[0.,0.,.4]]],resonance=[.4],damping=[.03])
        full=replace(base,materials=(offdiag,))
        self.assertIn('ld-tensor',full.input_text())
        self.assertIn('salmonpy-ld-tensor.dat',full.input_assets())
        self.assertIn('SALMONPY_LD_TENSOR 1',full.input_assets()['salmonpy-ld-tensor.dat'])
        for settings in ({'tensor_rtol':1.},{'tensor_rtol':float('nan')},{'tensor_maxiter':0}):
            with self.assertRaises(ValueError):replace(full,**settings)
