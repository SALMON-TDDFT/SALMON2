program stub_probe
 use exx_k_backend,only:s_exx_k_backend
 use exx_k_cufft,only:exx_k_cufft_create,s_exx_k_cufft
 use ieee_arithmetic,only:ieee_value,ieee_quiet_nan
 implicit none
 class(s_exx_k_backend),allocatable :: backend
 real(8) :: kernel(0:3,0:3,0:3)
 integer :: point(3,8),shift(3,8),slot(8),status,x,y,z,i
 complex(8) :: buffer(3,8,10),action(3,8,10)
 kernel=1d0;buffer=(1d0,2d0);slot=[10,1,3,4,5,6,7,8];i=0
 do z=0,1;do y=0,1;do x=0,1
  i=i+1;point(:,i)=[x,y,z];shift(:,i)=2*[x,y,z]
 enddo;enddo;enddo
 call exx_k_cufft_create(backend)
 call backend%release(status)
 if(status/=0)error stop 'initial release'
 call prepare(-1)
 call backend%apply(99,0,buffer,action,status)
 if(status/=0.or.any(action/=(0d0,0d0)))error stop 'zero rows beyond mesh'
 call backend%apply(8,1,buffer,action,status)
 if(status/=-1.or.any(action/=(0d0,0d0)))error stop 'disabled nonempty'
 call backend%apply(8,2,buffer,action,status)
 if(status/=-2)error stop 'rows exceed global mesh'
 call prepare(-1)
 call backend%apply(0,0,buffer,action,status)
 if(status/=-2)error stop 'invalid lower bound'
 slot(2)=10;call prepare(-2);slot(2)=1
 shift(1,1)=1;call prepare(-2);shift(1,1)=0
 point(1,1)=2;call prepare(-2);point(1,1)=0
 kernel(0,0,0)=ieee_value(0d0,ieee_quiet_nan);call prepare(-2);kernel=1d0
 call backend%prepare([huge(0),2,2],[2,2,2],3,kernel,point,shift,slot,10,status)
 if(status/=-2)error stop 'dimension overflow'
 call prepare(-1)
 buffer(1,1,1)=cmplx(ieee_value(0d0,ieee_quiet_nan),0d0,8)
 call backend%apply(1,1,buffer,action,status)
 if(status/=-2)error stop 'nonfinite tile'
 call backend%release(status)
 if(status/=0)error stop 'release'
 call backend%release(status)
 if(status/=0)error stop 'repeated release'
 select type(backend)
 type is(s_exx_k_cufft)
  if(backend%plan_builds/=0.or.backend%constant_uploads/=0.or.backend%tile_uploads/=0)error stop 'CPU counters'
 class default
  error stop 'Wrong factory dynamic type'
 end select
 buffer=(1d0,2d0)
 call prepare(-1)
 deallocate(backend) ! Finalize a validated, unreleased object.
 call exx_k_cufft_create(backend)
 call backend%apply(1,0,buffer,action,status)
 if(status/=-2)error stop 'fresh backend must not be ready'
 deallocate(backend)
 print *, 'k-cuFFT CPU stub validation passed'
contains
 subroutine prepare(expected)
  implicit none
  integer,intent(in) :: expected
  call backend%prepare([2,2,2],[2,2,2],3,kernel,point,shift,slot,10,status)
  if(status/=expected)error stop 'prepare status'
 end subroutine
end program
