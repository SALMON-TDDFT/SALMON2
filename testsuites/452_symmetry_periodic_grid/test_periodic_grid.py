"""Exercise both real production density-symmetry paths with checked arrays."""
from pathlib import Path
import os,subprocess,tempfile,unittest
ROOT=Path(__file__).resolve().parents[2]
STUBS='''module sym_sub
 implicit none
 real(8),allocatable :: SymMatA(:,:,:),SymMatB(:,:,:)
 logical :: use_symmetry=.true.
end module
module communication
 implicit none
contains
 subroutine comm_summation(a,b,n,comm)
  implicit none
  integer,intent(in)::n,comm
  real(8),intent(in)::a(:,:,:)
  real(8),intent(out)::b(:,:,:)
  b=a
 end subroutine
end module
module nvtx_wrapper
 implicit none
contains
 subroutine nvtxStartRange(name,line)
  implicit none
  character(*),intent(in)::name
  integer,intent(in)::line
 end subroutine
 subroutine nvtxEndRange
  implicit none
 end subroutine
end module
'''
DRIVER='''program periodic_grid_test
 use sym_sub
 use sym_rho_sub,only:init_sym_rho,sym_rho
 use,intrinsic::ieee_arithmetic,only:ieee_value,ieee_quiet_nan
 implicit none
 integer::n(3),i,j,k,ix,iy,iz,mode,shift(3),nsym
 real(8),allocatable::rho(:,:,:),ref(:,:,:)
 real(8)::v
 character(32)::arg
 call get_command_argument(1,arg);read(arg,*)mode
 n=[10,10,10];shift=0
 if(mode==2)n=[11,11,11]
 if(mode==3)n=[17,17,18]
 nsym=1;if(mode==9)nsym=2
 allocate(SymMatA(3,4,nsym),SymMatB(3,4,nsym));SymMatA=0d0;SymMatB=0d0
 do i=1,3
  SymMatA(i,i,1)=1d0;SymMatB(i,i,1)=1d0
 enddo
 select case(mode)
 case(1,2,3)
  ! Actual W/Se failure: a rotated zero acquires a tiny negative value.
  SymMatB(1,3,1)=-1.78975993575726293d-17
  SymMatB(1,2,1)=-5.55111512312578270d-17
 case(4)
  ! Values just inside the upper boundary still round to N.
  SymMatA(:,4,1)=1d-17
 case(5)
  ! Integer translations spanning multiple periods are equivalent.
  SymMatA(:,4,1)=[2d0,-3d0,1d0]
 case(6)
  ! A commensurate half-cell translation must keep the right mapping.
  SymMatA(1,4,1)=.5d0;shift(1)=5
 case(7)
  ! Retain the original1e-3 grid-commensurability rejection.
  SymMatA(1,4,1)=.125d0
 case(9)
  SymMatA(:,:,2)=SymMatA(:,:,1);SymMatB(:,:,2)=SymMatB(:,:,1)
  SymMatA(1,1,2)=-1d0;SymMatB(1,1,2)=-1d0
 case(10)
  SymMatA(1,4,1)=.000099d0
 case(11)
  SymMatA(1,4,1)=.000101d0
 case(8)
  SymMatB(1,3,1)=ieee_value(0d0,ieee_quiet_nan)
 end select
 allocate(rho(n(1),n(2),n(3)),ref(n(1),n(2),n(3)))
 do k=1,n(3)
 do j=1,n(2)
 do i=1,n(1)
  rho(i,j,k)=dble(i+100*j+10000*k)
 enddo
 enddo
 enddo
 ref=rho
 call init_sym_rho(n,[1,1,1],n,0)
 call sym_rho(rho)
 if(mode==7.or.mode==8.or.mode==11)error stop 'INVALID_INPUT_WAS_NOT_REJECTED'
 do k=1,n(3)
 do j=1,n(2)
 do i=1,n(1)
  ix=modulo(i-1-shift(1),n(1))+1
  iy=modulo(j-1-shift(2),n(2))+1
  iz=modulo(k-1-shift(3),n(3))+1
  v=ref(ix,iy,iz)
  if(mode==9)v=.5d0*(v+ref(modulo(-(i-1),n(1))+1,j,k))
  if(abs(rho(i,j,k)-v)>1d-10)error stop 'DENSITY_MAPPING_MISMATCH'
 enddo
 enddo
 enddo
 if(abs(sum(rho)-sum(ref))>1d-8)error stop 'CHARGE_NOT_PRESERVED'
 print *, 'PERIODIC_GRID_PASS'
end program
'''
class PeriodicGrid(unittest.TestCase):
 @classmethod
 def setUpClass(cls):
  cls.tmp=tempfile.TemporaryDirectory();p=Path(cls.tmp.name);(p/'stubs.f90').write_text(STUBS);(p/'driver.f90').write_text(DRIVER);cls.exe=p/'probe'
  q=subprocess.run([os.environ.get('FC','gfortran'),'-cpp','-fopenmp','-fcheck=all','-fbacktrace',str(p/'stubs.f90'),str(ROOT/'src/symmetry/symmetry_rho.f90'),str(p/'driver.f90'),'-o',str(cls.exe)],cwd=p,capture_output=True,text=True)
  if q.returncode:raise RuntimeError(q.stderr)
 @classmethod
 def tearDownClass(cls):cls.tmp.cleanup()
 def test_periodic_boundaries_and_density(self):
  for mode in list(range(1,7))+[9,10]:
   for threads in [1,2]:
    with self.subTest(mode=mode,threads=threads):
     q=subprocess.run([str(self.exe),str(mode)],capture_output=True,text=True,env=dict(os.environ,OMP_NUM_THREADS=str(threads)),timeout=15)
     self.assertEqual(q.returncode,0,q.stdout+q.stderr);self.assertIn('PERIODIC_GRID_PASS',q.stdout)
 def test_noncommensurate_and_nonfinite_rejected(self):
  for mode in [7,8,11]:
   q=subprocess.run([str(self.exe),str(mode)],capture_output=True,text=True,env=dict(os.environ,OMP_NUM_THREADS='1'),timeout=15)
   self.assertNotIn('INVALID_INPUT_WAS_NOT_REJECTED',q.stdout+q.stderr);self.assertNotIn('PERIODIC_GRID_PASS',q.stdout)
   self.assertIn('density symmetry',q.stdout+q.stderr)
if __name__=='__main__':unittest.main(verbosity=2)
