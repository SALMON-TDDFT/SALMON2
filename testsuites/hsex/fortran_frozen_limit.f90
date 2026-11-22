program test_frozen_limit
  use hsex_frozen_limit, only: frozen_optical_tensor,frozen_eta_limit,frozen_static_part
  use hsex_screen_math, only: screen_assemble
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan,ieee_positive_inf
  implicit none
  complex(8) :: psi(5,2),vertices(5,2),chi(5,5),reference(5,5),rotation(2,2),rotated(5,2)
  complex(8) :: vel(5,3),b(5,5),samples(2,2,3),limit(2,2),wanted(2,2)
  complex(8) :: raw_static(2,2),kept(2,2),expected_static(2,2),field(2)
  real(8) :: q_test(3,2),rec_test(3,3),head_test,qhead_test
  complex(8) :: mat_test(1,1),pole_h(2,2)
  real(8) :: e(5),rate,error,eta,q,chi_fixed,eta_first,q_first
  integer :: i,j,a,c,status,fail
  e=[-.8d0,-.3d0,.4d0,.7d0,1.1d0];psi=0d0;psi(1,1)=1d0;psi(2,2)=1d0;fail=0
  do j=1,2
    do i=1,5
      vertices(i,j)=exp(cmplx(0d0,dble(i*j),8))
    enddo
  enddo
  call frozen_optical_tensor(action,velocity,psi,vertices,1d-9,chi,rate,status)
  if(status/=0)error stop 'stationary occupied projector rejected'
  reference=0d0
  do a=1,2
    call velocity(psi(:,a),vel,status)
    b(:,:3)=vel
    b(:,4:)=vertices*spread(psi(:,a),2,2)
    do i=3,5
      b(i,:3)=b(i,:3)/(e(i)-e(a))
    enddo
    do j=1,5
      do c=1,5
        do i=3,5
          reference(c,j)=reference(c,j)+conjg(b(i,c))*b(i,j)/(e(i)-e(a))
        enddo
      enddo
    enddo
  enddo
  if(maxval(abs(chi-reference))>1d-10)error stop 'optical head/wing/body spectral reference mismatch'
  rotation(:,1)=[cmplx(cos(.4d0),0d0,8),cmplx(0d0,sin(.4d0),8)]
  rotation(:,2)=[cmplx(0d0,sin(.4d0),8),cmplx(cos(.4d0),0d0,8)]
  rotated=matmul(psi,rotation)
  call frozen_optical_tensor(action,velocity,rotated,vertices,1d-9,chi,rate,status)
  if(status/=0.or.maxval(abs(chi-reference))>1d-10)error stop 'occupied gauge dependence'
  rotated=psi;rotated(:,2)=0d0;rotated(2,2)=cos(.02d0);rotated(3,2)=sin(.02d0)
  call frozen_optical_tensor(action,velocity,rotated,vertices,1d-9,chi,rate,status)
  if(status/=3.or.rate<1d-3)error stop 'coherent nonstationary projector accepted as insulating equilibrium'
  do fail=1,3
    call frozen_optical_tensor(action,velocity,psi,vertices,1d-9,chi,rate,status)
    if(status==0)error stop 'failed/nonfinite callback accepted'
  enddo
  fail=0
  vertices(1,1)=cmplx(ieee_value(0d0,ieee_positive_inf),0d0,8)
  call frozen_optical_tensor(action,velocity,psi,vertices,1d-9,chi,rate,status)
  if(status==0)error stop 'infinite optical vertex accepted'
  wanted=reshape([cmplx(2d0,0d0,8),cmplx(.2d0,-.1d0,8), &
    cmplx(.2d0,.1d0,8),cmplx(1d0,0d0,8)],[2,2])
  do j=1,3
    eta=.01d0/2**(j-1)
    samples(:,:,j)=wanted*(1d0+.2d0*eta+3d0*eta**2)
  enddo
  call frozen_eta_limit(samples,.01d0,limit,error,status)
  if(status/=0.or.maxval(abs(limit-wanted))>1d-12)error stop 'eta linear/quadratic extrapolation'
  do j=1,3
    samples(:,:,j)=wanted*2d0**(j-1)
  enddo
  call frozen_eta_limit(samples,.01d0,limit,error,status)
  if(status==0)error stop 'eta pole accepted as static limit'
  ! Hydrodynamic finite-q response chi=-q^2/(q^2+eta^2). At fixed q the
  ! static limit is -1; q->0 first at finite eta instead gives zero.
  do i=1,3
    q=.2d0/2**(i-1)
    do j=1,3
      eta=q*.01d0/2**(j-1)
      samples(:,:,j)=-q*q/(q*q+eta*eta)
    enddo
    call frozen_eta_limit(samples,.001d0,limit,error,status)
    ! The retained eta^4 remainder is (7/32)*(.01)^4 to leading order.
    if(status/=0.or.maxval(abs(limit+1d0))>3d-9)error stop 'fixed q static extrapolation'
    eta_first=real(limit(1,1),8)
    eta=.02d0
    chi_fixed=-q*q/(q*q+eta*eta)
    print '(a,4es20.10)','FROZEN_Q_ETA_ORDER q/eta_fixed/chi_fixed/chi_eta0 ',q,eta,chi_fixed,eta_first
  enddo
  q=1d-8;eta=.02d0;q_first=-q*q/(q*q+eta*eta)
  if(abs(q_first)>1d-10.or.abs(eta_first-q_first)<.99d0)error stop 'noncommuting point limits'
  ! A fixed eta series that has not resolved v*q must be rejected, rather
  ! than silently treating suppressed screening as a converged static limit.
  q=.005d0
  do j=1,3
    eta=.02d0/2**(j-1);samples(:,:,j)=-q*q/(q*q+eta*eta)
  enddo
  call frozen_eta_limit(samples,.01d0,limit,error,status)
  if(status==0)error stop 'unresolved finite q eta accepted'
  samples=1d0;samples(2,1,3)=cmplx(ieee_value(0d0,ieee_quiet_nan),0d0,8)
  call frozen_eta_limit(samples,.01d0,limit,error,status)
  if(status==0)error stop 'NaN eta sample accepted'
  raw_static=reshape([cmplx(-2d0,.3d0,8),cmplx(.8d0,-.2d0,8), &
    cmplx(-.4d0,.5d0,8),cmplx(-1d0,-.1d0,8)],[2,2])
  expected_static=.5d0*(raw_static+conjg(transpose(raw_static)))
  call frozen_static_part(raw_static,'strict',.02d0,kept,error,status)
  if(status/=2)error stop 'strict closure accepted non-Hermitian response'
  call frozen_static_part(raw_static,'hermitian',.02d0,kept,error,status)
  if(status/=0.or.maxval(abs(kept-expected_static))>1d-14.or.error<.1d0) &
    error stop 'defined static projection mismatch'
  field=[cmplx(.7d0,-.2d0,8),cmplx(.3d0,.8d0,8)]
  if(abs(real(dot_product(field,matmul(raw_static,field)),8)- &
      real(dot_product(field,matmul(kept,field)),8))>1d-14)error stop 'real quadratic form changed'
  call frozen_static_part(kept,'strict',.02d0,expected_static,error,status)
  if(status/=0.or.error>1d-14)error stop 'projected closure not idempotent'
  raw_static(1,2)=cmplx(ieee_value(0d0,ieee_quiet_nan),0d0,8)
  call frozen_static_part(raw_static,'hermitian',.02d0,kept,error,status)
  if(status==0)error stop 'NaN response projected'
  raw_static(1,2)=cmplx(ieee_value(0d0,ieee_positive_inf),0d0,8)
  call frozen_static_part(raw_static,'hermitian',.02d0,kept,error,status)
  if(status==0)error stop 'infinite response projected'
  call frozen_static_part(expected_static,'invalid',.02d0,kept,error,status)
  if(status==0)error stop 'unknown closure accepted'
  rec_test=0d0
  do i=1,3
    rec_test(i,i)=2d0*acos(-1d0)/12d0
  enddo
  q_test=0d0;q_test(1,:)=[.25d0,1.25d0]
  raw_static=0d0;raw_static(1,1)=1d4;raw_static(2,2)=1d4
  raw_static(1,2)=cmplx(0d0,3d0,8);raw_static(2,1)=cmplx(0d0,3d0,8)
  call frozen_static_part(raw_static,'hermitian',.02d0,kept,error,status)
  if(status/=0)error stop 'projected positivity fixture rejected early'
  call screen_assemble([4,4,4],1,q_test,rec_test,12d0**3,2,-.5d0*kept,.true.,mat_test,head_test,qhead_test,status)
  if(status==0)error stop 'negative dielectric spectrum accepted after static projection'
  raw_static=cmplx(huge(1d0)/4d0,0d0,8)
  call frozen_static_part(raw_static,'hermitian',.02d0,kept,error,status)
  if(status==0)error stop 'overflowing norm accepted after static projection'
  call coherent_static_pole()
  print *, 'FROZEN_STATIC_HERMITIAN_QUADRATIC_IEEE_PASS'
  print *, 'FROZEN_OPTICAL_ETA_LIMIT_PASS'
contains
  subroutine coherent_static_pole()
    use hsex_frozen_response, only: frozen_source,frozen_solve,frozen_trace,s_frozen_solution,s_frozen_stats
    implicit none
    real(8) :: rotation(2,2),levels(2),eta,scaled(3),err
    complex(8) :: h(2,2),state(2,1),vertex(2),m(2,2),rho(2,2),source(2,2),exact,z(1),value(1)
    complex(8) :: samples(1,1,3),limit(1,1),kept(1,1)
    complex(8),allocatable :: u(:,:),v(:,:),c(:,:)
    type(s_frozen_solution) :: solution
    type(s_frozen_stats) :: stats
    integer :: j,a,b,code
    rotation(:,1)=[cos(.4d0),sin(.4d0)];rotation(:,2)=[-sin(.4d0),cos(.4d0)]
    levels=[-.4d0,.6d0]
    h=matmul(rotation*spread(levels,1,2),transpose(rotation));pole_h=h
    state(:,1)=cos(.3d0)*rotation(:,1)+cmplx(0d0,sin(.3d0),8)*rotation(:,2)
    vertex=[cmplx(1d0,0d0,8),cmplx(-1d0,0d0,8)]
    m=0d0;m(1,1)=1d0;m(2,2)=-1d0;m=matmul(transpose(rotation),matmul(m,rotation))
    rho=2d0*matmul(state,conjg(transpose(state)))
    rho=matmul(transpose(rotation),matmul(rho,rotation))
    source=matmul(m,rho)-matmul(rho,m)
    call frozen_source(state,state,[2d0],[2d0],vertex,u,v,c,code)
    if(code/=0)error stop 'coherent pole source'
    do j=1,3
      eta=.01d0/2**(j-1);z(1)=cmplx(0d0,eta,8)
      call frozen_solve(action_local,action_local,u,c,v,z,2,1d-8,solution,stats,code)
      if(code/=0)error stop 'coherent pole solve'
      call frozen_trace(solution,vertex,value,code)
      if(code/=0)error stop 'coherent pole trace'
      exact=0d0
      do b=1,2
        do a=1,2
          exact=exact+m(b,a)*source(a,b)/(z(1)-levels(a)+levels(b))
        enddo
      enddo
      if(abs(value(1)-exact)>1d-8)error stop 'coherent pole spectral mismatch'
      call frozen_static_part(reshape(value,[1,1]),'hermitian',.02d0,kept,err,code)
      if(code/=0)error stop 'coherent pole static projection'
      samples(:,:,j)=kept;scaled(j)=eta*real(kept(1,1),8)
      print '(a,3es20.10)','FROZEN_COHERENT_REAL_POLE eta/chiH/eta_chiH ',eta,real(kept(1,1),8),scaled(j)
    enddo
    if(abs(scaled(3))<.1d0.or.abs(scaled(1)-scaled(3))>.02d0*abs(scaled(3))) &
      error stop 'physical coherent pole fixture lost secular response'
    call frozen_eta_limit(samples,.02d0,limit,err,code)
    if(code==0)error stop 'real Hermitian coherent pole accepted as a static limit'
    print *,'FROZEN_COHERENT_HERMITIAN_POLE_REJECT_PASS'
  end subroutine
  subroutine action_local(x,y,status)
      implicit none
      complex(8),intent(in) :: x(:)
      complex(8),intent(out) :: y(:)
      integer,intent(out) :: status
      y=matmul(pole_h,x);status=0
  end subroutine
  subroutine action(x,y,code)
    implicit none
    complex(8),intent(in) :: x(:)
    complex(8),intent(out) :: y(:)
    integer,intent(out) :: code
    code=1
    if(fail==1)return
    y=e*x;code=0
    if(fail==2)y(1)=cmplx(ieee_value(0d0,ieee_quiet_nan),0d0,8)
  end subroutine
  subroutine velocity(x,y,code)
    implicit none
    complex(8),intent(in) :: x(:)
    complex(8),intent(out) :: y(:,:)
    integer,intent(out) :: code
    integer :: d,l,m
    complex(8) :: v
    code=0;y=0d0
    do d=1,3
      do l=1,5
        do m=1,5
          v=cmplx(cos(dble(d*(l+m))),sin(dble(d*(l-m))),8)
          y(l,d)=y(l,d)+v*x(m)
        enddo
      enddo
    enddo
    if(fail==3)y(1,1)=cmplx(ieee_value(0d0,ieee_positive_inf),0d0,8)
  end subroutine
end program
