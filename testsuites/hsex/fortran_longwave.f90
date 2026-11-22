program test_longwave
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan,ieee_positive_inf
  use hsex_longwave, only: longwave_band,longwave_cell
  use hsex_screen_math, only: screen_assemble,screen_inverse,screen_finite
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
  block
    real(8) :: u(4,4),position(4,3),direction(3),qmag,head_coefficient,finite_head,previous_error,ratio
    complex(8) :: optical(5,5),shifted(5,5),v(4,3),density(4,2),rho(4,2),wing(2),finite_wing(2)
    complex(8) :: uniform(4),phase(4),x(4,3)
    real(8) :: errors(4),head_error,wing_error,hsite(4,4),diamagnetic(3,3),paramagnetic(3,3)
    integer :: site,other,second_axis
    integer :: axis,target,step,sign_index,geometry
    ! Independent density vertices: diagonal position in a four-site basis,
    ! transformed into the Hamiltonian eigenbasis. No optical chi enters them.
    u(:,1)=[1d0,1d0,1d0,1d0]/2d0
    u(:,2)=[1d0,-1d0,1d0,-1d0]/2d0
    u(:,3)=[1d0,1d0,-1d0,-1d0]/2d0
    u(:,4)=[1d0,-1d0,-1d0,1d0]/2d0
    position(:,1)=[0d0,1d0,.3d0,-.7d0]
    position(:,2)=[.2d0,-.4d0,1.1d0,.6d0]
    position(:,3)=[-.3d0,.8d0,-.9d0,1.2d0]
    do axis=1,3
      x(:,axis)=matmul(transpose(u),position(:,axis)*u(:,1))
      ! dH_k/dk at k=0 for H_k=exp(-ik.r) H exp(ik.r).
      v(:,axis)=cmplx(0d0,1d0,8)*(energy-energy(1))*x(:,axis)
    enddo
    ! Independent current sum rule: <d2H/dk_a dk_b> + Pi_ab(0) = 0.
    ! The contact term uses real-space commutators, not velocity matrix elements.
    hsite=0d0
    do target=1,4
      do site=1,4
        hsite(:,site)=hsite(:,site)+u(:,target)*energy(target)*u(site,target)
      enddo
    enddo
    diamagnetic=0d0;paramagnetic=0d0
    do axis=1,3
      do second_axis=1,3
        do site=1,4
          do other=1,4
            diamagnetic(axis,second_axis)=diamagnetic(axis,second_axis)- &
              u(site,1)*hsite(site,other)*u(other,1)*(position(site,axis)-position(other,axis))* &
              (position(site,second_axis)-position(other,second_axis))
          enddo
        enddo
        do target=2,4
          paramagnetic(axis,second_axis)=paramagnetic(axis,second_axis)- &
            2d0*real(conjg(v(target,axis))*v(target,second_axis),8)/(energy(target)-energy(1))
        enddo
      enddo
    enddo
    if(.not.screen_finite(diamagnetic).or..not.screen_finite(paramagnetic))error stop 'nonfinite DC sum rule'
    if(maxval(abs(diamagnetic+paramagnetic))>1d-12)error stop 'insulating current DC cancellation failed'
    print *, 'INSULATING_DC_CONTACT_CANCELLATION',maxval(abs(diamagnetic+paramagnetic))
    do axis=1,2
      phase=exp(cmplx(0d0,1d0,8)*position(:,axis))
      density(:,axis)=matmul(transpose(u),phase*u(:,1))
    enddo
    call longwave_band(apply,occupied,energy(1),v,density,optical,status)
    if(status/=0.or..not.screen_finite(optical))error stop 'zero-limit optical solve failed'
    ! A uniform scalar potential cannot excite a gapped system.
    uniform=matmul(transpose(u),u(:,1))
    if(maxval(abs(uniform(2:)))>1d-14)error stop 'uniform density has interband response'
    ! Occupied components of velocity/density vertices must have no effect.
    v(1,:)=v(1,:)+cmplx(100d0,-77d0,8)
    density(1,:)=density(1,:)+cmplx(-91d0,123d0,8)
    call longwave_band(apply,occupied,energy(1),v,density,shifted,status)
    if(status/=0.or..not.screen_finite(shifted))error stop 'occupied-vertex invariance solve failed'
    if(maxval(abs(shifted-optical))>1d-12)error stop 'occupied vertex changes optical response'
    do geometry=1,2
      direction=[1d0,2d0,-1d0]
      if(geometry==2)direction=[-2d0,.5d0,3d0]
      direction=direction/sqrt(sum(direction**2))
      head_coefficient=real(dot_product(cmplx(direction,0d0,8),matmul(optical(:3,:3),direction)),8)
      do axis=1,2
        wing(axis)=sum(direction*optical(:3,axis+3))
      enddo
      previous_error=huge(1d0)
      do step=1,4
        qmag=.02d0/2**(step-1)
        do sign_index=1,2
          phase=exp(cmplx(0d0,(-1d0)**(sign_index-1)*qmag,8)*matmul(position,direction))
          rho(:,sign_index)=matmul(transpose(u),phase*u(:,1))
        enddo
        finite_head=0d0;finite_wing=0d0
        do target=2,4
          finite_head=finite_head+abs(rho(target,1))**2/(energy(target)-energy(1))/qmag**2
          do axis=1,2
            ! Rebuild the unshifted body vertex independently of the optical result.
            phase=exp(cmplx(0d0,1d0,8)*position(:,axis))
            uniform=matmul(transpose(u),phase*u(:,1))
            finite_wing(axis)=finite_wing(axis)+(conjg(rho(target,1))-conjg(rho(target,2)))* &
              uniform(target)/(2d0*qmag*(energy(target)-energy(1)))
          enddo
        enddo
        head_error=abs(finite_head-head_coefficient)
        wing_error=maxval(abs(finite_wing-wing))
        errors(step)=max(head_error,wing_error)
        if(.not.screen_finite(errors(:step)))error stop 'nonfinite zero-limit error'
        if(step>1)then
          ratio=previous_error/errors(step)
          if(ratio<3.9d0.or.ratio>4.1d0)error stop 'head/wing zero-limit convergence is not quadratic'
        endif
        previous_error=errors(step)
        print *, 'INSULATING_ZERO_LIMIT q/head_error/wing_error',qmag,head_error,wing_error
      enddo
      if(errors(4)>2d-6)error stop 'finite-q head/wings differ from Sternheimer optical limit'
    enddo
    print *, 'INSULATING_UNIFORM_DENSITY_AND_OCCUPIED_GAUGE_PASS'
  end block
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
