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

    def test_input_and_fixed_settings(self):
        c=self.calculation();text=c.input_text()
        self.assertIn("theory = 'maxwell'",text)
        self.assertIn("unit_system = 'au'",text)
        self.assertIn('omega_p_ld(1) = 1',text)
        self.assertIn('obs_samp_em = 1',text)
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
            with self.assertRaises(FileExistsError): runner.run(self.calculation(),target)
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

    def test_mpi_grid_configuration(self):
        import salmonpy.fdtd as f
        self.assertTrue(hasattr(f,'FDTDMPI'),'FDTD MPI launch configuration missing')
        mpi=f.FDTDMPI((1,2,2))
        self.assertEqual(mpi.ranks,4)
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
