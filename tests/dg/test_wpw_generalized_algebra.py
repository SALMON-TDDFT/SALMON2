from pathlib import Path
import shutil, subprocess, tempfile, unittest
import numpy as np

ROOT=Path(__file__).resolve().parents[2]

class GeneralizedAlgebraTest(unittest.TestCase):
 def test_driver(self):
  fc=shutil.which('gfortran'); self.assertIsNotNone(fc)
  with tempfile.TemporaryDirectory(prefix='wpw-generalized-') as td:
   td=Path(td)
   (td/'CMakeLists.txt').write_text(f'''cmake_minimum_required(VERSION 3.16)
project(wpw_generalized LANGUAGES Fortran)
find_package(LAPACK REQUIRED)
add_executable(check "{ROOT/'src/common/dg_generalized_algebra.f90'}" "{ROOT/'tests/dg/fixtures/wpw_generalized_driver.F90'}")
target_link_libraries(check PRIVATE ${{LAPACK_LIBRARIES}})
target_compile_options(check PRIVATE -fcheck=all -fbacktrace)
''')
   for cmd in (["cmake","-S",str(td),"-B",str(td/'b'),f"-DCMAKE_Fortran_COMPILER={fc}"],
               ["cmake","--build",str(td/'b'),"-j","2"]):
    r=subprocess.run(cmd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
    self.assertEqual(r.returncode,0,r.stdout)
   r=subprocess.run([str(td/'b/check')],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
   self.assertEqual(r.returncode,0,r.stdout); self.assertIn('PASS generalized algebra',r.stdout)
   lines={line.split()[0]:[float(v) for v in line.split()[1:]] for line in r.stdout.splitlines() if line.split() and line.split()[0] in {'EVAL','X','C','EXP','EXPNF'}}
   s=np.array([[1.4,.12+.04j,.03-.02j],[.12-.04j,1.1,.08+.01j],[.03+.02j,.08-.01j,.9]],complex)
   h0=np.array([[.8,.11-.03j,.02+.04j],[.11+.03j,1.3,.07-.02j],[.02-.04j,.07+.02j,1.8]],complex)
   hf=np.array([[0,-.19+.06j,0],[-.19-.06j,.13,0],[0,0,0]],complex);h=h0+hf
   se,sv=np.linalg.eigh(s);x=(sv*se**-.5)@sv.conj().T
   he,hv=np.linalg.eigh(x.conj().T@h@x)
   np.testing.assert_allclose(lines['EVAL'],he,atol=2e-13)
   xf=np.array(lines['X'][0::2])+1j*np.array(lines['X'][1::2]);xf=xf.reshape((3,3),order='F')
   np.testing.assert_allclose(xf,x,atol=2e-13);np.testing.assert_allclose(xf.conj().T@s@xf,np.eye(3),atol=2e-13)
   cf_eig=np.array(lines['C'][0::2])+1j*np.array(lines['C'][1::2]);cf_eig=cf_eig.reshape((3,3),order='F')
   for i in range(3):
    phase=np.vdot(cf_eig[:,i],x@hv[:,i]);cf_eig[:,i]*=phase/abs(phase)
   np.testing.assert_allclose(cf_eig,x@hv,atol=4e-13)
   np.testing.assert_allclose(h@cf_eig,s@cf_eig@np.diag(he),atol=5e-13)
   c0=np.array([.7+.1j,-.2+.3j,.1-.4j]);c0/=np.sqrt(np.vdot(c0,s@c0).real)
   d0=x.conj().T@s@c0;cref=x@hv@(np.exp(-1j*he*.17)*(hv.conj().T@d0))
   cf=np.array(lines['EXP'][0::2])+1j*np.array(lines['EXP'][1::2])
   np.testing.assert_allclose(cf,cref,atol=3e-13)
   np.testing.assert_allclose(np.vdot(cf,s@cf),1.0,atol=3e-13)
   hne,hne_v=np.linalg.eigh(x.conj().T@h0@x);cref_nf=x@hne_v@(np.exp(-1j*hne*.17)*(hne_v.conj().T@d0))
   cnf=np.array(lines['EXPNF'][0::2])+1j*np.array(lines['EXPNF'][1::2])
   np.testing.assert_allclose(cnf,cref_nf,atol=3e-13);self.assertGreater(np.linalg.norm(cnf-cref),1e-4)

if __name__=='__main__': unittest.main()
