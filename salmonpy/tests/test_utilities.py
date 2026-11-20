import tempfile
import unittest
from pathlib import Path
import numpy as np
from salmonpy.fields import AcosPulse, sample_pulses
from salmonpy.geometry import Grid, Shape, Geometry
from salmonpy.analysis import read_table, spectrum, optical_conductivity, dielectric_function, optical_linearity, optical_sum_rule

class UtilityTests(unittest.TestCase):
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
