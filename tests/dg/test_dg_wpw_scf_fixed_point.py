from pathlib import Path
import shutil, subprocess, tempfile, unittest
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
class WpwScfFixedPointTest(unittest.TestCase):
 def test_driver(self):
  fc=shutil.which('gfortran');self.assertIsNotNone(fc)
  source=ROOT/'src/gs/dc/dg_wannier_pw_scf.f90';self.assertTrue(source.exists(),'missing fixed-basis SCF solver')
  with tempfile.TemporaryDirectory(prefix='wpw-scf-') as td:
   td=Path(td);(td/'CMakeLists.txt').write_text(f'''cmake_minimum_required(VERSION 3.16)
project(wpw_scf LANGUAGES Fortran)
find_package(LAPACK REQUIRED)
add_executable(check "{ROOT/'src/common/dg_wpw_fixed_operator.f90'}" "{ROOT/'src/common/dg_generalized_algebra.f90'}" "{source}" "{ROOT/'tests/dg/fixtures/dg_wpw_scf_driver.F90'}")
target_link_libraries(check PRIVATE LAPACK::LAPACK)
''')
   for cmd in (["cmake","-S",str(td),"-B",str(td/'b'),f"-DCMAKE_Fortran_COMPILER={fc}"],["cmake","--build",str(td/'b'),"-j","2"]):
    r=subprocess.run(cmd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True);self.assertEqual(r.returncode,0,r.stdout)
   r=subprocess.run([str(td/'b/check')],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True);self.assertEqual(r.returncode,0,r.stdout);self.assertIn('PASS DG WPW fixed-point SCF',r.stdout)
   line=next(x for x in r.stdout.splitlines() if x.startswith('SCF_REF'))
   got=np.array([float(x) for x in line.split()[1:]])
   h0=np.array([[-.8,.04,.08],[.04,.15,-.03+.01j],[.08,-.03-.01j,.9]],complex)
   occ=np.array([2.,1.]);d=np.array([1.5,1.,.5]);v=np.zeros((3,3),complex)
   for _ in range(20000):
    ev,c=np.linalg.eigh(h0+v);p=(c[:,:2]*occ)@c[:,:2].conj().T;dn=np.real(np.diag(p));vo=np.diag(.15*dn)
    if np.linalg.norm(vo-v)/max(1e-30,np.linalg.norm(vo))<1e-13: break
    v=.65*v+.35*vo;d=dn
   ev,c=np.linalg.eigh(h0+vo);p=(c[:,:2]*occ)@c[:,:2].conj().T;d=np.real(np.diag(p));eh=.5*np.sum(d*(.15*d));energy=np.dot(occ,ev[:2])-eh+.2
   ref=np.r_[ev,d,energy]
   np.testing.assert_allclose(got,ref,rtol=2e-10,atol=2e-11)
if __name__=='__main__':unittest.main()
