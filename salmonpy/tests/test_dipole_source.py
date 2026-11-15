import unittest
import numpy as np
from salmonpy.errors import InputError


class DipoleSourceTests(unittest.TestCase):
    def model(self):
        from salmonpy.dipole_source import DipoleMap
        return DipoleMap(('mol-a','mol-b'), (10,20), [[.25,1.],[.75,0.]], [2.,5.])

    def test_conservation_superposition_and_static_offset(self):
        m = self.model()
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
