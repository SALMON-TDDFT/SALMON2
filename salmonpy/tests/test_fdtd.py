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
