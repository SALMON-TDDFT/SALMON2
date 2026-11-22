program test_frozen_kinetic
  use hsex_frozen_kinetic, only: frozen_kinetic_setup,frozen_kinetic_apply
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan,ieee_positive_inf
  implicit none
  integer,parameter :: n(3)=[4,4,4],ng=64
  real(8) :: cell(3,3),denominator(ng),k(3),energy
  complex(8) :: x(ng),y(ng)
  integer :: i,status
  cell=0d0
  do i=1,3
    cell(i,i)=8d0
  enddo
  k=[.2d0,-.1d0,.05d0]
  call frozen_kinetic_setup(n,cell,k,denominator,status)
  if(status/=0)error stop 'kinetic setup'
  do i=0,ng-1
    x(i+1)=exp(cmplx(0d0,2d0*acos(-1d0)*(mod(i,4)-mod(i/4,4))/4d0,8))
  enddo
  energy=1d0+.5d0*sum((k+2d0*acos(-1d0)/8d0*[1d0,-1d0,0d0])**2)
  call frozen_kinetic_apply(n,denominator,x,y,status)
  if(status/=0.or.maxval(abs(y-x/energy))>1d-13)error stop 'shifted Fourier kinetic inverse'
  denominator(1)=ieee_value(0d0,ieee_quiet_nan)
  call frozen_kinetic_apply(n,denominator,x,y,status)
  if(status==0)error stop 'NaN precondition accepted'
  denominator(1)=ieee_value(0d0,ieee_positive_inf)
  call frozen_kinetic_apply(n,denominator,x,y,status)
  if(status==0)error stop 'infinite precondition accepted'
  denominator(1)=-1d0
  call frozen_kinetic_apply(n,denominator,x,y,status)
  if(status==0)error stop 'negative precondition accepted'
  call frozen_kinetic_setup([3,4,4],cell,k,denominator(:48),status)
  if(status/=0)error stop 'non power of two grid rejected'
  call general_grid()
  print *,'FROZEN_KINETIC_PHASE_IEEE_PASS'
contains
  subroutine general_grid()
    use hsex_screen_math, only: screen_geometry
    implicit none
    integer,parameter :: dims(3)=[12,12,12],m=1728
    real(8) :: primitive(3,3),rec(3,3),vol,d(m),momentum(3),value
    complex(8) :: plane(m),result(m)
    integer :: j,code
    primitive(:,1)=[5.13d0,5.13d0,0d0]
    primitive(:,2)=[0d0,5.13d0,5.13d0]
    primitive(:,3)=[5.13d0,0d0,5.13d0]
    call screen_geometry(primitive,rec,vol,code)
    if(code/=0)error stop 'general grid geometry'
    momentum=matmul(rec,[.2d0,-.1d0,.15d0])
    call frozen_kinetic_setup(dims,primitive,momentum,d,code)
    if(code/=0)error stop 'Nr12 setup rejected'
    do j=0,m-1
      plane(j+1)=exp(cmplx(0d0,2d0*acos(-1d0)*(mod(j,12)-mod(j/12,12))/12d0,8))
    enddo
    value=1d0+.5d0*sum((momentum+matmul(rec,[1d0,-1d0,0d0]))**2)
    call frozen_kinetic_apply(dims,d,plane,result,code)
    if(code/=0.or.maxval(abs(result-plane/value))>1d-12)error stop 'Nr12 nonorthogonal kinetic inverse'
  end subroutine
end program
