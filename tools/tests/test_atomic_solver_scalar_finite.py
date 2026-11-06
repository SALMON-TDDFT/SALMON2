from pathlib import Path
import re,subprocess,tempfile,os
root=Path(__file__).resolve().parents[2]
s=(root/'src/atom/pp/atomic_solver.f90').read_text()
routine=re.search(r'function interp_rho_in\(.*?end function interp_rho_in',s,re.S).group()
with tempfile.TemporaryDirectory() as d:
 p=Path(d)
 stub='''module structures
 type s_pp_info
 integer :: mr(1)=3
 real(8) :: rad(3,1),rho_pp_tbl(3,1)
 end type
end module
module math_constants
 real(8),parameter :: pi=acos(-1d0)
end module
module salmon_math
contains
 real(8) function interp_linear(n,x,y,r) result(v)
 integer :: n,i
 real(8) :: x(n),y(n),r
 do i=1,n-1
 if(r<=x(i+1))then
 v=y(i)+(y(i+1)-y(i))*(r-x(i))/(x(i+1)-x(i));return
 endif
 enddo
 v=y(n)
 end function
end module
module test_atomic
contains
'''
 driver='''program test
 use structures
 use test_atomic
 use ieee_arithmetic,only:ieee_value,ieee_quiet_nan,ieee_positive_inf
 implicit none
 type(s_pp_info)::pp
 real(8)::rho(2)
 integer::mode,pos
 character(12)::a
 call get_command_argument(1,a);read(a,*)mode
 call get_command_argument(2,a);read(a,*)pos
 pp%rad(:,1)=[1d0,2d0,3d0];pp%rho_pp_tbl(:,1)=4d0*acos(-1d0)
 if(mode==1)pp%rho_pp_tbl(pos,1)=ieee_value(1d0,ieee_quiet_nan)
 if(mode==2)pp%rho_pp_tbl(pos,1)=ieee_value(1d0,ieee_positive_inf)
 if(mode==3)pp%rho_pp_tbl(pos,1)=-1d0
 rho=interp_rho_in(pp,1,2,[1d0,2d0])
 if(maxval(abs(rho-[1d0,.25d0]))>1d-14)error stop 'wrong interpolation'
 end program
'''
 (p/'probe.f90').write_text(stub+routine+'\nend module\n'+driver)
 q=subprocess.run([os.environ.get('FC','gfortran'),'-fcheck=all',str(root/'src/common/ieee_finite_check.f90'),str(p/'probe.f90'),'-o',str(p/'probe')],cwd=p,capture_output=True,text=True);assert q.returncode==0,q.stderr
 count=0
 for mode in range(4):
  for pos in range(1,4):
   q=subprocess.run([str(p/'probe'),str(mode),str(pos)],capture_output=True,text=True)
   assert (q.returncode==0)==(mode==0),(mode,pos,q.stdout,q.stderr)
   if mode:assert 'invalid valence density' in q.stderr
   count+=1
 (p/'bad.f90').write_text('program bad\nuse ieee_finite_check\nlogical::ok(2)\nok=is_finite([1d0,2d0])\nend program\n')
 q=subprocess.run([os.environ.get('FC','gfortran'),'-c',str(p/'bad.f90')],cwd=p,capture_output=True,text=True);assert q.returncode!=0
 print('PASS:',count,'finite/NaN/Inf/negative density cases; array wrapper compile rejection')
