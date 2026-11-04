program gaussian
  use kli_self_hartree
  use ieee_arithmetic
  implicit none
  integer,parameter::n=24
  type(isolated_box)::box
  real(8)::step(3,3),rho(0:n-1,0:n-1,0:n-1),v(0:n-1,0:n-1,0:n-1),r(3),eh,b,pi
  integer::x,y,z,s,origin,i
  step=0d0
  do i=1,3
    step(i,i)=.4d0
  enddo
  pi=acos(-1d0)
  call prepare_box(box,[n,n,n],step,s)
  if(s/=0)error stop 1
  do origin=0,n/2,n/2
    do z=0,n-1;do y=0,n-1;do x=0,n-1
      r=.4d0*(modulo([x,y,z]-origin+n/2,n)-n/2)
      rho(x,y,z)=exp(-sum(r*r))/pi**1.5d0
    enddo;enddo;enddo
    call self_hartree(box,rho,v,eh,b,s)
    if(s/=0.or.abs(eh-1/sqrt(2*pi))>1d-7.or.abs(v(origin,origin,origin)-2/sqrt(pi))>1d-6)error stop 2
    print *, origin,eh,v(origin,origin,origin),b
  enddo
  rho(0,0,0)=ieee_value(0d0,ieee_quiet_nan)
  call self_hartree(box,rho,v,eh,b,s)
  if(s==0)error stop 3
  rho(0,0,0)=ieee_value(0d0,ieee_positive_inf)
  call self_hartree(box,rho,v,eh,b,s)
  if(s==0)error stop 4
  call release_box(box)
end program
