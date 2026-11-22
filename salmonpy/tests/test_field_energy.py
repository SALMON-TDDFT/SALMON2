"""Independent manufactured Yee step, not a Maxwell update implementation."""
import ctypes as c
import unittest
import numpy as np


class FieldBoxEnergyTests(unittest.TestCase):
    def test_manufactured_balance_source_sign_partition_and_finite(self):
        from salmonpy._sbe_handles import default_api
        fn=default_api().library.sp_fdtd_box_energy
        fn.argtypes=[c.c_void_p]*6+[c.c_double]*2+[c.c_int]+[c.c_void_p]*5
        fn.restype=c.c_int
        shape=(8,8,8);dt=.1;speed=1.;beta=.2
        old=[np.zeros(shape,order='F') for _ in range(12)]
        new=[np.zeros(shape,order='F') for _ in range(12)]
        z=np.arange(8)[None,None,:];eta=.03
        old[6][:]=beta*z*z;old[2][:]=eta*z*z
        new[6][:]=old[6]+dt*(eta*(2*z+1)+2*dt*beta)
        source=np.zeros((3,*shape),order='F');source[1]=.01
        poles=[np.zeros((*shape,1),order='F') for _ in range(6)]
        poles[1][:]=.002;poles[4][:]=.002
        new[2][:]=old[2]+dt*beta*(2*z-1)-4*np.pi*dt*.012
        lo=np.array([0,0,0],dtype=np.int32);hi=np.array([7,7,7],dtype=np.int32)
        box=np.array([[2,2,2],[5,5,5]],dtype=np.int32);dx=np.ones(3)
        ptr=lambda arrays:(c.c_void_p*len(arrays))(*[a.ctypes.data for a in arrays])
        oldp,newp,jp=ptr(old),ptr(new),ptr(poles)
        def sample(ownedlo=(1,1,1),ownedhi=(6,6,6),out=None):
            a=np.array(ownedlo,dtype=np.int32);b=np.array(ownedhi,dtype=np.int32)
            if out is None:out=np.full(17,123.)
            status=fn(lo.ctypes.data,hi.ctypes.data,a.ctypes.data,b.ctypes.data,
                box.ctypes.data,dx.ctypes.data,dt,speed,1,oldp,newp,jp,source.ctypes.data,out.ctypes.data)
            return status,out
        status,full=sample();self.assertEqual(status,0)
        self.assertLess(abs(full[15]),1.e-13)
        self.assertGreater(abs(full[16]),1.e-7)
        self.assertGreater(full[6],0.);self.assertGreater(full[7],0.)
        _,a=sample(ownedhi=(3,6,6));_,b=sample(ownedlo=(4,1,1))

        # Cancellation of ~26 Hartree storage leaves ~1e-14 roundoff residuals.
        np.testing.assert_allclose(a+b,full,rtol=1.e-12,atol=1.e-13)
        source[1]*=-1
        status,bad=sample();self.assertEqual(status,0);self.assertGreater(abs(bad[15]),1.e-5)
        source[1]*=-1
        for value in (np.nan,np.inf,-np.inf):
            new[2][2,2,2]=value
            status,out=sample()
            self.assertNotEqual(status,0);np.testing.assert_array_equal(out,123.)
