"""Compiled rectangular Wannier exchange versus an independent Bloch sum."""
import os
import sys
from pathlib import Path
import subprocess
import tempfile
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[2]


def multiplier(shape, spacing, omega, shift):
    axes = [2*np.pi*np.fft.fftfreq(n, d=h) for n, h in zip(shape, spacing)]
    q = np.stack(np.meshgrid(*axes, indexing='ij'), axis=-1)
    nyquist = np.pi/np.array(spacing)
    q = (q+shift+nyquist) % (2*nyquist)-nyquist
    q2 = np.sum(q*q, axis=-1)
    result = np.full(shape, np.pi/omega**2)
    np.divide(4*np.pi*(-np.expm1(-q2/(4*omega**2))), q2, out=result, where=q2>1e-24)
    return result


def reference(source, target, k, spacing, occ, omega):
    result = np.zeros_like(target)
    for ik, kv in enumerate(k):
        for iq, qv in enumerate(k):
            v = multiplier(source.shape[2:], spacing, omega, kv-qv)
            for j, phi in enumerate(source[iq]):
                pairs = phi.conj()*target[ik]
                potential = np.fft.ifftn(np.fft.fftn(pairs, axes=(-3,-2,-1))*v, axes=(-3,-2,-1))
                result[ik] -= occ[iq,j]/(2*len(k))*phi*potential
    return result


def wire(a):
    return a.transpose(2,3,4,1,0).tobytes(order='F')


class WannierTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.path = Path(cls.tmp.name)
        cls.exe = cls.path/'probe'
        module = ROOT/'src/xc/hse_wannier.f90'
        if not module.exists():
            raise AssertionError('Rectangular fractional-occupation Wannier backend is missing')
        fftw = Path(os.environ.get('FFTW_ROOT','/opt/homebrew/opt/fftw'))
        blas = Path(os.environ.get('OPENBLAS_ROOT','/opt/homebrew/opt/openblas'))
        cmd = [os.environ.get('FC','gfortran'),'-O0','-g','-fcheck=all','-fopenmp',
               '-I'+str(fftw/'include'),str(ROOT/'src/xc/hse_wannier_gauge.f90'),str(module),
               str(ROOT/'src/xc/hse_ace.f90'),str(Path(__file__).with_name('probe.f90')),
               '-L'+str(fftw/'lib'),'-lfftw3','-L'+str(blas/'lib'),'-lopenblas','-o',str(cls.exe)]
        p = subprocess.run(cmd,cwd=cls.path,capture_output=True,text=True)
        if p.returncode:
            raise AssertionError(p.stderr)

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def run_case(self, mesh, fractional=True, empty=False):
        shape = (4,3,2); h=np.array([.6,.8,.9]); omega=.11; no=2; nt=3
        rng=np.random.default_rng(573)
        k=np.array(list(np.ndindex(*mesh)))*2*np.pi/(np.array(mesh)*shape*h)+[.017,-.021,.003]
        k=k[rng.permutation(len(k))]
        source=np.empty((len(k),no)+shape,complex)
        for ik in range(len(k)):
            a=rng.normal(size=(np.prod(shape),no))+1j*rng.normal(size=(np.prod(shape),no))
            source[ik]=(np.linalg.qr(a)[0].T/np.sqrt(np.prod(h))).reshape((no,)+shape)
        target=rng.normal(size=(len(k),nt)+shape)+1j*rng.normal(size=(len(k),nt)+shape)
        occ=np.tile([2.,.37 if fractional else 2.],(len(k),1))
        if fractional and (len(k)>1 or empty):
            occ[-1,1]=0.
        inp=self.path/'in.bin';out=self.path/'out.bin'
        with inp.open('wb') as f:
            f.write(np.array([*shape,*mesh,no,nt],np.int32).tobytes())
            f.write(np.array([*h,omega],np.float64).tobytes())
            f.write(k.T.tobytes(order='F'));f.write(occ.T.tobytes(order='F'))
            f.write(wire(source));f.write(wire(target))
        p=subprocess.run([str(self.exe),str(inp),str(out)],capture_output=True,text=True,
                         env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS=str(getattr(self,'threads',1)),WANNIER_EXPECT_WORKERS=str(getattr(self,'threads',1)),WANNIER_TEST_BATCH=str(getattr(self,'batch',1))))
        self.assertEqual(p.returncode,0,p.stdout+p.stderr)
        a=np.fromfile(out,np.complex128).reshape(shape+(nt,len(k)),order='F').transpose(4,3,0,1,2)
        np.testing.assert_allclose(a,reference(source,target,k,h,occ,omega),rtol=2e-11,atol=2e-11)
        localized=np.fromfile(str(out)+'.localized',np.complex128).reshape(a.transpose(2,3,4,1,0).shape,order='F').transpose(4,3,0,1,2)
        np.testing.assert_allclose(localized,a,rtol=2e-11,atol=2e-11)
        sys.path.insert(0,str(ROOT/'samples/dc_hse'))
        from read_snapshot import read_snapshot
        snap=read_snapshot(str(out)+'.snapshot')
        np.testing.assert_allclose(snap['occupation'],occ.T)
        self.assertTrue(snap['converged'])
        self.assertEqual(snap['scf_iterations'],7)
        for ik in range(len(k)):
            recovered=snap['phi'][:,:,ik]@snap['u'][:,:,ik].conj().T
            expected=source[ik].transpose(1,2,3,0).reshape((np.prod(shape),no),order='F')
            np.testing.assert_allclose(recovered,expected,atol=1e-11)
        checks=np.loadtxt(str(out)+'.checks')
        self.assertLess(np.max(checks[:4]),1e-10,p.stdout)
        self.assertLessEqual(checks[4],1e-10)  # largest exchange-metric eigenvalue

    def test_post_scf_localization(self):
        self.run_case((1,1,1),empty=True)
        exe=self.path/'localize'
        fftw=Path(os.environ.get('FFTW_ROOT','/opt/homebrew/opt/fftw'))
        blas=Path(os.environ.get('OPENBLAS_ROOT','/opt/homebrew/opt/openblas'))
        cmd=[os.environ.get('FC','gfortran'),'-O2','-fcheck=all','-I'+str(fftw/'include'),
             str(ROOT/'src/xc/hse_wannier_gauge.f90'),str(ROOT/'src/xc/hse_wannier.f90'),
             str(ROOT/'samples/dc_hse/localize_snapshot.f90'),'-L'+str(fftw/'lib'),'-lfftw3',
             '-L'+str(blas/'lib'),'-lopenblas','-o',str(exe)]
        p=subprocess.run(cmd,cwd=self.path,capture_output=True,text=True)
        self.assertEqual(p.returncode,0,p.stderr)
        source=self.path/'out.bin.snapshot';dest=self.path/'reloc.bin'
        p=subprocess.run([str(exe),str(source),str(dest),'0','100','1e-6'],capture_output=True,text=True)
        self.assertEqual(p.returncode,0,p.stdout+p.stderr)
        from read_snapshot import read_snapshot
        a=read_snapshot(source);b=read_snapshot(dest)
        np.testing.assert_allclose(a['q']@a['q'].conj().T,b['q']@b['q'].conj().T,atol=1e-11)
        self.assertEqual(b['scf_iterations'],7)
        self.assertEqual(b['q'].shape[1],1)
        self.assertTrue(b['converged'])
        # A positive cutoff changes the density and cannot inherit certification.
        self.run_case((1,1,1))
        p=subprocess.run([str(exe),str(source),str(dest),'0.5','100','1e-6'],capture_output=True,text=True)
        self.assertEqual(p.returncode,0,p.stdout+p.stderr)
        self.assertFalse(read_snapshot(dest)['converged'])

    def test_rectangular_fractional_shifted_mesh(self):
        self.run_case((1,2,2))

    def test_gamma_fractional(self):
        self.run_case((1,1,1))

    def test_threaded_action(self):
        self.threads=1;self.run_case((1,2,2))
        serial=np.fromfile(self.path/'out.bin',np.complex128)
        self.threads=2;self.run_case((1,2,2))
        parallel=np.fromfile(self.path/'out.bin',np.complex128)
        np.testing.assert_allclose(parallel,serial,atol=1e-13,rtol=1e-13)

    def test_batched_shifted_mesh(self):
        self.batch=1;self.run_case((1,2,2))
        serial=np.fromfile(self.path/'out.bin',np.complex128)
        self.batch=4;self.run_case((1,2,2))
        batched=np.fromfile(self.path/'out.bin',np.complex128)
        np.testing.assert_allclose(batched,serial,atol=1e-13,rtol=1e-13)

    def test_equal_occupations(self):
        self.run_case((2,1,1),False)

if __name__=='__main__':
    unittest.main()
