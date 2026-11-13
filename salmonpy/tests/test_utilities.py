import tempfile
import unittest
from pathlib import Path
import numpy as np
from salmonpy.fields import AcosPulse, sample_pulses
from salmonpy.geometry import Grid, Shape, Geometry
from salmonpy.analysis import read_table, spectrum

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
