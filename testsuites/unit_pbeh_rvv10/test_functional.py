"""Compiled kernels checked against independent quadrature and variational derivatives."""
import ctypes as ct
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent

class FunctionalTest(unittest.TestCase):
    def compile(self, sources, probe):
        tmp=tempfile.TemporaryDirectory();self.addCleanup(tmp.cleanup)
        exe=Path(tmp.name)/'probe'
        cmd=['gfortran','-O0','-g','-fcheck=all','-ffree-line-length-none','-I/opt/homebrew/include']
        cmd += [str(ROOT/'src/xc'/s) for s in sources]+[str(HERE/probe)]
        cmd += ['-L/opt/homebrew/lib','-lxc','-lfftw3','-o',str(exe)]
        p=subprocess.run(cmd,cwd=tmp.name,capture_output=True,text=True)
        self.assertEqual(p.returncode,0,p.stderr)
        return exe

    def test_rvv10_variational_and_kernel(self):
        exe=self.compile(['rvv10.f90'],'rvv10_probe.f90')
        p=subprocess.run([exe],capture_output=True,text=True,check=True)
        lines=p.stdout.splitlines()
        derivative=np.fromstring(lines[0],sep=' ')
        self.assertLess(max(derivative),2e-8)
        convergence=np.fromstring(lines[1],sep=' ')
        differences=np.abs(convergence[:-1]-convergence[-1])
        self.assertLess(differences[2],differences[1])
        self.assertLess(differences[1],differences[0])
        self.assertLess(differences[1],1e-6)
        uniform=np.fromstring(lines[2],sep=' ')
        self.assertLess(max(uniform[:2]),2e-6)
        actual=np.array([float(x) for x in lines[3:]])
        expected=[]
        for i in range(1,4):
            a=.02;b=.02*i
            def phi(r):return -1.5/((1+a*r*r)*(1+b*r*r)*(2+(a+b)*r*r))
            for k in np.arange(4)*.1:
                if k==0:
                    t,weight=np.polynomial.legendre.leggauss(200)
                    u=(t+1)/2;r=u/(1-u)
                    val=np.sum(weight/2*4*np.pi*r*r*phi(r)/(1-u)**2)
                else:
                    r=np.linspace(0,4000,400001)
                    val=4*np.pi/k*np.trapezoid(r*phi(r)*np.sin(k*r),r)
                expected.append(val)
        np.testing.assert_allclose(actual,expected,rtol=1e-9,atol=1e-7)

    def test_periodic_density_derivative(self):
        exe=self.compile(["rvv10.f90"],"periodic_probe.f90")
        subprocess.run([exe],capture_output=True,text=True,check=True)

    def test_pbeh_semilocal(self):
        self.check_semilocal(.4)

    def test_pbe0_semilocal(self):
        self.check_semilocal(.25)

    def check_semilocal(self, fraction):
        exe=self.compile(['hse_semilocal.f90'],'semilocal_probe.f90')
        p=subprocess.run([exe,str(fraction)],capture_output=True,text=True,check=True)
        actual=np.loadtxt(p.stdout.splitlines())
        lib=ct.CDLL('/opt/homebrew/lib/libxc.dylib');ptr=ct.POINTER(ct.c_double)
        lib.xc_func_alloc.restype=ct.c_void_p
        lib.xc_func_init.argtypes=[ct.c_void_p,ct.c_int,ct.c_int]
        lib.xc_gga_exc_vxc.argtypes=[ct.c_void_p,ct.c_size_t]+[ptr]*5
        lib.xc_func_end.argtypes=lib.xc_func_free.argtypes=[ct.c_void_p]
        r=np.array([1e-7,.001,.1,1.]);s=np.array([1e-12,.0001,.03,.2])
        expected=np.zeros((4,3))
        for id,weight in [(101,1-fraction),(130,1.)]:
            f=lib.xc_func_alloc();self.assertEqual(lib.xc_func_init(f,id,1),0)
            arrays=[np.zeros(4) for _ in range(3)]
            lib.xc_gga_exc_vxc(f,4,*[a.ctypes.data_as(ptr) for a in [r,s]+arrays])
            expected += weight*np.array(arrays).T
            lib.xc_func_end(f);lib.xc_func_free(f)
        np.testing.assert_allclose(actual,expected,atol=1e-12,rtol=1e-12)

if __name__=='__main__':unittest.main()
