program test_longwave
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan,ieee_positive_inf
  use hsex_longwave, only: longwave_band,longwave_cell
  use hsex_screen_math, only: screen_assemble,screen_inverse
  implicit none
  complex(8) :: occupied(4,1),velocity(4,3),body(4,2),response(5,5),expected(5,5),vertices(4,5)
  complex(8) :: chi(5,5),raw(5,5),point(1,1),average(1,1),chif(3,3),mat(1,1),isotropic(3,3)
  real(8) :: energy(4),macro(3,3),rec(3,3),q(3,3),qvec(3),head,qhead,vol,pi,factor,qc,wanted,error
  integer :: i,j,l,status,g(3,2),empty_g(3,0),reverse(5)
  block
    complex(8) :: a(2,2),inverse(2,2),identity(2,2)
    ! A two-by-two indefinite pivot also exercises the Hermitian solve.
    a=0d0;a(1,2)=cmplx(1d0,2d0,8);a(2,1)=conjg(a(1,2))
    identity=0d0;identity(1,1)=1d0;identity(2,2)=1d0
    call screen_inverse(a,inverse,status)
    if(status/=0)error stop 'Hermitian indefinite inverse rejected'
    if(maxval(abs(matmul(a,inverse)-identity))>1d-13)error stop 'Hermitian inverse residual'
    a=0d0
    call screen_inverse(a,inverse,status)
    if(status==0)error stop 'singular dielectric accepted'
    a=identity;a(1,2)=1d0
    call screen_inverse(a,inverse,status)
    if(status==0)error stop 'non-Hermitian dielectric accepted'
  end block
  occupied=0d0;occupied(1,1)=1d0;energy=[0d0,1d0,2d0,3d0]
  do j=1,3
    do i=1,4
      velocity(i,j)=cmplx(sin(dble(i+j)),cos(dble(2*i+j)),8)
    enddo
  enddo
  body(:,1)=[(1d0,0d0),(0d0,1d0),(.2d0,.3d0),(.1d0,-.4d0)]
  body(:,2)=conjg(body(:,1))
  call longwave_band(apply,occupied,0d0,velocity,body,response,status)
  if(status/=0)error stop 'longwave solve rejected finite input'
  vertices(:,4:)=body
  do i=2,4
    vertices(i,:3)=velocity(i,:)/energy(i)
  enddo
  vertices(1,:)=0d0;expected=0d0
  do j=1,5
    do i=1,5
      do l=2,4
        expected(i,j)=expected(i,j)+conjg(vertices(l,i))*vertices(l,j)/energy(l)
      enddo
    enddo
  enddo
  error=maxval(abs(response-expected))
  if(error>1d-11)error stop 'Sternheimer differs from explicit conduction-state reference'
  raw=response
  reverse=[1,2,3,5,4]
  do j=1,5
    do i=1,5
      chi(i,j)=raw(i,j)+merge(-1d0,1d0,(i<=3).neqv.(j<=3))*conjg(raw(reverse(i),reverse(j)))
    enddo
  enddo
  pi=acos(-1d0);vol=1000d0;rec=0d0
  do i=1,3
    rec(i,i)=2d0*pi/10d0
  enddo
  g(:,1)=[1,0,0];g(:,2)=-g(:,1)
  call longwave_cell([4,4,4],1,g,rec,vol,8,chi,macro,point,average,status)
  if(status/=0)error stop 'local-field optical closure failed'
  qvec=[1d-7,2d-7,3d-7]
  chif(1,1)=sum(qvec*matmul(chi(:3,:3),qvec))
  do i=1,2
    chif(i+1,1)=sum(chi(i+3,:3)*qvec)
    chif(1,i+1)=sum(chi(:3,i+3)*qvec)
  enddo
  chif(2:,2:)=chi(4:,4:)
  q(:,1)=qvec/rec(1,1)
  q(:,2)=g(:,1)+q(:,1);q(:,3)=g(:,2)+q(:,1)
  call screen_assemble([4,4,4],1,q,rec,vol,8,chif,.false.,mat,head,qhead,status)
  wanted=sum(qvec*qvec)/dot_product(qvec,matmul(macro,qvec))
  if(status/=0.or.abs(head-wanted)>1d-7)error stop 'optical dielectric limit differs from finite q'
  factor=2d0/(8*vol);isotropic=0d0
  do i=1,3
    isotropic(i,i)=11d0/(4d0*pi*factor)
  enddo
  ! One-point box form factor is exactly one, isolating the Coulomb integral.
  call longwave_cell([1,1,1],1,empty_g,rec,vol,8,isotropic,macro,point,average,status)
  qc=(6d0*pi*pi/(vol*8))**(1d0/3d0)
  wanted=12d0*pi/(12d0*qc*qc*vol)
  if(status/=0.or.abs(real(average(1,1),8)-wanted)>1d-12)error stop 'isotropic mini-sphere normalization'
  isotropic=0d0;isotropic(1,1)=ieee_value(0d0,ieee_quiet_nan)
  call screen_inverse(isotropic,expected(:3,:3),status)
  if(status==0)error stop 'NaN dielectric accepted'
  isotropic(1,1)=ieee_value(0d0,ieee_positive_inf)
  call screen_inverse(isotropic,expected(:3,:3),status)
  if(status==0)error stop 'infinite dielectric accepted'
  velocity(4,3)=ieee_value(0d0,ieee_quiet_nan)
  call longwave_band(apply,occupied,0d0,velocity,body,response,status)
  if(status==0)error stop 'NaN optical vertex accepted'
  velocity(4,3)=ieee_value(0d0,ieee_positive_inf)
  call longwave_band(apply,occupied,0d0,velocity,body,response,status)
  if(status==0)error stop 'infinite optical vertex accepted'
  print *, 'LONGWAVE_STERNHEIMER_LOCAL_FIELD_PASS',error
contains
  subroutine apply(x,y,status)
    implicit none
    complex(8),intent(in) :: x(:)
    complex(8),intent(out) :: y(:)
    integer,intent(out) :: status
    y=energy*x
    status=0
  end subroutine
end program
