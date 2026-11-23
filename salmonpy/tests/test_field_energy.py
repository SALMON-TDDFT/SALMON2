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

    def test_diagonal_constitutive_storage_balance_partition_and_bad_coefficients(self):
        from salmonpy._sbe_handles import default_api
        fn=default_api().library.sp_fdtd_box_energy_diagonal
        fn.argtypes=[c.c_void_p]*6+[c.c_double]*2+[c.c_int]+[c.c_void_p]*6
        fn.restype=c.c_int
        shape=(8,8,8);dt=.1;speed=1.;beta=.2;eta=.03
        old=[np.zeros(shape,order='F') for _ in range(12)]
        new=[np.zeros(shape,order='F') for _ in range(12)]
        z=np.arange(8)[None,None,:]
        eps=np.array([2.,3.,4.]);mu=np.array([1.5,2.,2.5]);dx=np.array([.8,1.2,1.5])
        old[6][:]=beta*z*z;old[2][:]=eta*z*z
        new[2][:]=old[2]+dt*(beta*(2*z-1)/dx[2]-4*np.pi*.012)/eps[1]
        new[6][:]=old[6]+dt*(eta*(2*z+1)+2*dt*beta/(eps[1]*dx[2]))/(mu[0]*dx[2])
        coeff=[np.ones(shape,order='F') for _ in range(9)]
        for d in range(3):
            c4=eps[d]+2*np.pi*dt*.03
            coeff[d][:]=(eps[d]-2*np.pi*dt*.03)/c4
            coeff[d+3][:]=-4*np.pi*dt/c4
            coeff[d+6][:]=-speed*dt/(mu[d]*dx[(d+1)%3])
        ptr=lambda arrays:(c.c_void_p*len(arrays))(*[a.ctypes.data for a in arrays])
        oldp,newp,cp=ptr(old),ptr(new),ptr(coeff)
        poles=[np.zeros((*shape,1),order='F') for _ in range(6)]
        poles[1][:]=.002;poles[4][:]=.002
        current=ptr(poles)
        source=np.zeros((3,*shape),order='F');source[1]=.01
        lo=np.zeros(3,dtype=np.int32);hi=np.full(3,7,dtype=np.int32)
        box=np.array([[2,2,2],[5,5,5]],dtype=np.int32)
        def sample(a=(1,1,1),b=(6,6,6)):
            a=np.array(a,dtype=np.int32);b=np.array(b,dtype=np.int32);out=np.full(17,123.)
            code=fn(lo.ctypes.data,hi.ctypes.data,a.ctypes.data,b.ctypes.data,box.ctypes.data,
                    dx.ctypes.data,dt,speed,1,oldp,newp,current,source.ctypes.data,cp,out.ctypes.data)
            return code,out
        code,whole=sample();self.assertEqual(code,0);self.assertLess(abs(whole[15]),2e-13)
        _,left=sample(b=(3,6,6));_,right=sample(a=(4,1,1))
        np.testing.assert_allclose(left+right,whole,atol=2e-13,rtol=1e-12)
        self.assertGreater(whole[6],0.);self.assertGreater(whole[7],0.)
        e=np.array([.13,.27,.31]);h=np.array([.17,.21,.53])
        for arrays in (old,new):
            for d in range(12):arrays[d][:]=0.
            for d in range(3):arrays[2*d][:]=e[d];arrays[6+2*d][:]=h[d]
        for array in poles:array[:]=0.
        source[:]=0.
        code,static=sample();self.assertEqual(code,0)
        expected=64*np.prod(dx)*(np.dot(eps*e,e)+np.dot(mu*h,h))/(8*np.pi)
        for index in (0,1,3,4):self.assertAlmostEqual(static[index],expected,places=12)
        for value in (0.,np.nan,np.inf):
            coeff[3][3,3,3]=value
            code,out=sample();self.assertNotEqual(code,0);np.testing.assert_array_equal(out,123.)

    def test_constitutive_cross_energy_boundary_partition_and_rejection(self):
        from salmonpy._sbe_handles import default_api
        fn=default_api().library.sp_fdtd_cross_energy
        fn.argtypes=[c.c_void_p]*11;fn.restype=c.c_int
        shape=(3,4,2,2);size=np.prod(shape)
        rng=np.random.default_rng(419)
        # A globally symmetric nonlocal operator, independent of production interpolation.
        raw=rng.normal(size=(size,size));operator=(raw+raw.T)*.01
        old=np.asfortranarray(rng.normal(size=shape));new=np.asfortranarray(rng.normal(size=shape))
        ko=np.asfortranarray((operator@old.ravel(order='F')).reshape(shape,order='F'))
        kn=np.asfortranarray((operator@new.ravel(order='F')).reshape(shape,order='F'))
        lo=np.zeros(3,dtype=np.int32);hi=np.array(shape[1:],dtype=np.int32)-1
        dx=np.array([.8,1.2,1.5]);scale=np.prod(dx)/(8*np.pi)
        def sample(lower,upper,ownedlo=lo,ownedhi=hi):
            box=np.array([lower,upper],dtype=np.int32);out=np.full(3,123.)
            code=fn(*[x.ctypes.data for x in (lo,hi,ownedlo,ownedhi,box,dx,old,new,ko,kn,out)])
            return code,out
        code,whole=sample(lo,hi);self.assertEqual(code,0)
        expected=scale*np.array([np.sum(old*ko),np.sum(new*kn),np.sum(old*kn-new*ko)])
        np.testing.assert_allclose(whole,expected,atol=1e-14)
        self.assertLess(abs(whole[2]),1e-14)
        code,left=sample(lo,[1,1,1]);self.assertEqual(code,0)
        code,right=sample([2,0,0],hi);self.assertEqual(code,0)
        np.testing.assert_allclose(left+right,whole,atol=1e-14)
        self.assertGreater(abs(left[2]),1e-5)
        # Local energy change + outward constitutive exchange equals midpoint work.
        work=2*scale*np.sum((.5*(old+new)*(kn-ko))[:,:2])
        self.assertAlmostEqual(left[1]-left[0]+left[2],work,places=13)
        _,partition=sample(lo,hi,ownedhi=np.array([1,1,1],dtype=np.int32))
        np.testing.assert_allclose(partition,left,atol=1e-14)
        for value in (np.nan,np.inf,-np.inf):
            ko[0,0,0,0]=value
            code,out=sample(lo,hi);self.assertNotEqual(code,0)
            np.testing.assert_array_equal(out,123.)
        summary=default_api().library.sp_fdtd_box_summary
        summary.argtypes=[c.c_int,c.c_void_p,c.c_double,c.c_double,c.c_double,c.c_void_p]
        summary.restype=c.c_int
        steps=np.zeros((1,19));steps[0,[2,6,7,14,15,5,17,18]]=[1,2,3,4,15,9,5,7]
        out=np.full(12,123.)
        self.assertEqual(summary(1,steps.ctypes.data,6.,7.,8.,out.ctypes.data),0)
        np.testing.assert_array_equal(out,[1,2,3,4,15,15,32,-5,9,40,5,7])

    def test_ohmic_loss_excludes_ld_implicit_term_and_preserves_invalid_output(self):
        from salmonpy._sbe_handles import default_api
        fn=default_api().library.sp_fdtd_box_energy_conductive
        fn.argtypes=[c.c_void_p]*6+[c.c_double]*2+[c.c_int]+[c.c_void_p]*6
        fn.restype=c.c_int
        shape=(8,8,8);dt=.1;dx=np.array([.8,1.2,1.5]);eps=np.array([2.,3.,4.])
        sigma=np.array([.1,.2,.3]);implicit=np.array([.03,.05,.07]);e0=np.array([.2,.3,.4])
        old=[np.zeros(shape,order='F') for _ in range(12)]
        new=[np.zeros(shape,order='F') for _ in range(12)]
        coeff=[np.ones(shape,order='F') for _ in range(9)]+[np.zeros((*shape,1),order='F') for _ in range(3)]
        currents=[np.zeros((*shape,1),order='F') for _ in range(6)]
        ptr=lambda arrays:(c.c_void_p*len(arrays))(*[v.ctypes.data for v in arrays])
        lo=np.zeros(3,dtype=np.int32);hi=np.full(3,7,dtype=np.int32)
        box=np.array([[2,2,2],[5,5,5]],dtype=np.int32)
        def prepare():
            e1=np.empty(3)
            for d in range(3):
                c4=eps[d]+2*np.pi*dt*(sigma[d]+implicit[d])
                ratio=(eps[d]-2*np.pi*dt*(sigma[d]+implicit[d]))/c4
                e1[d]=ratio*e0[d];old[2*d][:]=e0[d];new[2*d][:]=e1[d]
                coeff[d][:]=ratio;coeff[d+3][:]=-4*np.pi*dt/c4
                coeff[d+6][:]=-dt/dx[(d+1)%3];coeff[d+9][:]=implicit[d]
                currents[d+3][:]=implicit[d]*(e0[d]+e1[d])
            return e1
        def sample(ownedlo=lo,ownedhi=hi):
            out=np.full(18,123.)
            code=fn(*[x.ctypes.data for x in (lo,hi,ownedlo,ownedhi,box,dx)],dt,1.,1,
                    ptr(old),ptr(new),ptr(currents),None,ptr(coeff),out.ctypes.data)
            return code,out
        e1=prepare();code,out=sample();self.assertEqual(code,0)
        expected=64*np.prod(dx)*dt*np.dot(sigma,(.5*(e0+e1))**2)
        self.assertAlmostEqual(out[17],expected,places=13)
        self.assertGreater(out[17],0);self.assertLess(abs(out[15]),1e-13)
        self.assertGreater(out[6],0)
        _,left=sample(ownedhi=np.array([3,7,7],dtype=np.int32))
        _,right=sample(ownedlo=np.array([4,0,0],dtype=np.int32))
        np.testing.assert_allclose(left+right,out,atol=1e-13,rtol=1e-12)
        sigma[:]=0;prepare();code,out=sample();self.assertEqual(code,0)
        self.assertLess(abs(out[17]),1e-15);self.assertLess(abs(out[15]),1e-13)
        for bad in (np.nan,np.inf,-np.inf):
            prepare();coeff[3][3,3,3]=bad;code,out=sample()
            self.assertNotEqual(code,0);np.testing.assert_array_equal(out,123.)
            prepare();coeff[9][3,3,3,0]=bad;code,out=sample()
            self.assertNotEqual(code,0);np.testing.assert_array_equal(out,123.)
        prepare();coeff[0][3,3,3]=1.1;code,out=sample()
        self.assertNotEqual(code,0);np.testing.assert_array_equal(out,123.)
