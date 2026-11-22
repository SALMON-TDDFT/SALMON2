import unittest
import numpy as np
from salmonpy.errors import InputError


class DipoleSourceTests(unittest.TestCase):
    def model(self):
        from salmonpy.dipole_source import DipoleMap
        return DipoleMap(('mol-a','mol-b'), (10,20), [[.25,1.],[.75,0.]], [2.,5.])

    def test_conservation_superposition_and_static_offset(self):
        m = self.model()
        from salmonpy._sbe_handles import default_api
        self.assertTrue(callable(default_api().library.sp_molecular_map))
        p = np.array([[1.,2.,3.],[-1.,4.,2.]])
        P = m.polarization(p)
        np.testing.assert_allclose(np.sum(P*np.array([2,5])[:,None],axis=0),p.sum(axis=0))
        j = m.current(p,p*2,dt=.25)
        np.testing.assert_allclose(np.sum(j*np.array([2,5])[:,None],axis=0),p.sum(axis=0)/.25)
        np.testing.assert_array_equal(m.current(p,p,dt=.25),np.zeros((2,3)))
        np.testing.assert_allclose(m.current(p+10,2*p+10,dt=.25),j)
        np.testing.assert_allclose(P,m.polarization([p[0],[0,0,0]])+m.polarization([[0,0,0],p[1]]))

    def test_ownership(self):
        from salmonpy.dipole_source import DipoleMap
        w=np.ones((1,1));v=np.array([2.])
        m=DipoleMap(('m',),(1,),w,v);w[:]=0;v[:]=10
        np.testing.assert_array_equal(m.polarization([[2,0,0]]),[[1,0,0]])
        with self.assertRaises(ValueError):m.weights[0,0]=0
        with self.assertRaises(ValueError):m.weights.setflags(write=True)

    def test_rejects_bad_mapping_and_values(self):
        from salmonpy.dipole_source import DipoleMap
        for w,v in (([[.9]],[1]),([[-1]],[1]),([[1]],[0]),([[np.nan]],[1]),([[True]],[1])):
            with self.subTest(w=w,v=v),self.assertRaises(InputError):DipoleMap(('m',),(1,),w,v)
        with self.assertRaises(InputError):DipoleMap(('m','m'),(1,),[[.5,.5]],[1])
        m=self.model()
        for p in ([[1,2,3]],[[1j,0,0],[0,0,0]],[[np.nan,0,0],[0,0,0]]):
            with self.assertRaises(InputError):m.polarization(p)
        for dt in (0,-1,np.inf,True):
            with self.assertRaises(InputError):m.current(np.zeros((2,3)),np.ones((2,3)),dt=dt)
        from salmonpy._sbe_handles import default_api
        # Bind via the public map, then check native finite rejection directly.
        m.polarization(np.zeros((2,3)))
        fn=default_api().library.sp_molecular_map
        weights=np.ones((1,1));vol=np.ones(1);p=np.zeros((1,3));out=np.full((1,3),-7.)
        for bad in (float('nan'),float('inf')):
            p[0,0]=bad
            self.assertEqual(fn(0,1,1,1.,weights.ctypes.data,vol.ctypes.data,p.ctypes.data,p.ctypes.data,out.ctypes.data),1)
            np.testing.assert_array_equal(out,np.full((1,3),-7.))

    def test_fortran_centered_oscillator_work_and_rejection(self):
        from salmonpy.dipole_source import _midpoint_oscillator
        from salmonpy.errors import NumericalError
        p=np.array([.2,.1,0.]);v=np.array([.3,0.,0.]);left=np.array([.5,.2,0.]);right=np.array([.7,-.1,0.])
        old_p=p.copy();old_v=v.copy()
        pn,vn,j,stats=_midpoint_oscillator(p,v,left,right,dt=.01,omega=2.,volume=8.)
        np.testing.assert_array_equal(p,old_p);np.testing.assert_array_equal(v,old_v)
        np.testing.assert_allclose(j,(pn-p)/.08,rtol=1e-13,atol=1e-15)
        self.assertLess(abs(stats[5]),1e-12);self.assertLess(abs(stats[6]),1e-12)
        self.assertGreater(abs(stats[7]),1e-8)
        own=np.array([.1,-.2,.3])
        q,w,k,other=_midpoint_oscillator(p,v,left+own,right+own,dt=.01,omega=2.,volume=8.,own_left=own,own_right=own)
        np.testing.assert_allclose(q,pn,atol=1e-15,rtol=1e-13)
        for invalid in (float('nan'),float('inf')):
            with self.assertRaises((InputError,NumericalError)):
                _midpoint_oscillator(p,v,[invalid,0,0],right,dt=.01,omega=2.,volume=8.)
        np.testing.assert_array_equal(p,old_p)
        # Exercise the Fortran IEEE rejection directly, bypassing Python validation.
        from salmonpy._sbe_handles import default_api
        fn=default_api().library.sp_molecular_midpoint
        vectors=[old_p.copy(),old_v.copy(),left.copy(),right.copy(),np.zeros(3),np.zeros(3)]
        outputs=[np.full(3,-7.),np.full(3,-7.),np.full(3,-7.),np.full(8,-7.)]
        for bad in (float('nan'),float('inf')):
            vectors[2][0]=bad
            self.assertEqual(fn(.01,2.,8.,*(a.ctypes.data for a in vectors+outputs)),1)
            for out in outputs:np.testing.assert_array_equal(out,np.full(out.shape,-7.))
        vectors[2][:]=left
        vectors[0][:]=1e300
        self.assertEqual(fn(.01,2.,8.,*(a.ctypes.data for a in vectors+outputs)),4)
        for out in outputs:np.testing.assert_array_equal(out,np.full(out.shape,-7.))


class MolecularEnergyWorkTests(unittest.TestCase):
    def test_held_work_fresh_energy_and_sign(self):
        from salmonpy.dipole_source import _energy_work
        fields=[[0.,0.,0.],[2.,0.,0.],[4.,0.,0.]]
        dipoles=[[0.,0.,0.],[1.,0.,0.],[2.,0.,0.]]
        result=_energy_work(fields,dipoles,[10.,999.,16.],[1,0,1])
        np.testing.assert_array_equal(result[-1],[6.,4.,6.,0.,2.,1.])
        np.testing.assert_array_equal(result[1,2:],[0.,0.,0.,0.])
        reversed_result=_energy_work(-np.array(fields),dipoles,[10.,999.,16.],[1,0,1])
        self.assertEqual(reversed_result[-1,3],12.)

        import ctypes as c
        from salmonpy.dipole_source import _molecular_function
        fn=_molecular_function('sp_molecular_energy_work',[c.c_int]+[c.c_void_p]*5)
        E=np.array(fields,dtype=float);p=np.array(dipoles,dtype=float)
        energy=np.array([10.,999.,16.]);fresh=np.array([1.,0.,1.])
        out=np.full((3,6),123.)
        for bad in (float('nan'),float('inf'),-float('inf')):
            energy[1]=bad
            self.assertNotEqual(fn(3,E.ctypes.data,p.ctypes.data,energy.ctypes.data,fresh.ctypes.data,out.ctypes.data),0)
            np.testing.assert_array_equal(out,123.)
        energy[1]=999.;E[-1,0]=1.e308;p[-1,0]=1.e308
        self.assertEqual(fn(3,E.ctypes.data,p.ctypes.data,energy.ctypes.data,fresh.ctypes.data,out.ctypes.data),4)
        np.testing.assert_array_equal(out,123.)

        mean=_molecular_function('sp_molecular_interval_field',[c.c_void_p]*3)
        left=np.full(3,1.e308);right=np.full(3,1.e308);selected=np.empty(3)
        self.assertEqual(mean(left.ctypes.data,right.ctypes.data,selected.ctypes.data),0)
        np.testing.assert_array_equal(selected,1.e308)
        for bad in (float('nan'),float('inf'),-float('inf')):
            left[-1]=bad;selected[:]=123.
            self.assertNotEqual(mean(left.ctypes.data,right.ctypes.data,selected.ctypes.data),0)
            np.testing.assert_array_equal(selected,123.)

        total=_energy_work([[1.,0.,0.],[3.,0.,0.],[5.,0.,0.]],dipoles,[10.,999.,16.],[1,0,1])
        own=_energy_work([[1.,0.,0.]]*3,dipoles,[10.,999.,16.],[1,0,1])
        self.assertEqual(total[-1,1],6.)
        self.assertEqual(own[-1,1],2.)
        self.assertEqual(result[-1,1],4.)
