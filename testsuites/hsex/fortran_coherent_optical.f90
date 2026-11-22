program test_coherent_optical
  use hsex_coherent_optical, only: coherent_optical,stationary_current_reference
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan,ieee_positive_inf
  implicit none
  integer,parameter :: n=4,nz=3
  complex(8) :: h(n,n),v(n,n),mass(n,n),position(n,n),psi(n,1),rho(n,n),z(nz),sigma(nz),expected(nz)
  complex(8) :: charge(1,nz),vertices(n,1),density(n,n),first(n,n),second(n,n),source(n,n),td(nz),fine(nz),seed(n,n)
  real(8) :: e(n),residual(2,nz),occ(1),volume,norm_error,dc,alpha,wanted_alpha
  integer :: i,j,s,coherent,status,actions,fail
  e=[-.6d0,-.1d0,.5d0,1d0];h=0d0;volume=2.3d0;occ=2d0;vertices=1d0;fail=0
  do i=1,n
    h(i,i)=e(i)
    do j=1,n
      position(i,j)=cmplx(cos(dble(i+j)),sin(dble(i-j)),8)/3d0
    enddo
  enddo
  v=cmplx(0d0,1d0,8)*(matmul(h,position)-matmul(position,h))
  mass=2d0*matmul(position,matmul(h,position))-matmul(matmul(position,position),h)- &
    matmul(h,matmul(position,position))
  z=[cmplx(.3d0,.12d0,8),cmplx(.3d0,.08d0,8),cmplx(0d0,.08d0,8)]
  do coherent=0,1
    psi=0d0;psi(1,1)=1d0
    if(coherent==1)then
      psi(1,1)=cos(.3d0);psi(3,1)=cmplx(0d0,sin(.3d0),8)
    endif
    rho=occ(1)*matmul(psi,conjg(transpose(psi)))
    call coherent_optical(apply,velocity,acceleration,psi,occ,vertices,z,n,1d-10,volume, &
      sigma,charge,residual,actions,status)
    if(status/=0)error stop 'coherent optical response rejected'
    if(maxval(abs(charge))>1d-10)error stop 'uniform optical charge not conserved'
    do s=1,nz
      do j=1,n
        do i=1,n
          first(i,j)=rho(i,j)/(z(s)-e(i)+e(j))
        enddo
      enddo
      source=matmul(v,first)-matmul(first,v)
      do j=1,n
        do i=1,n
          second(i,j)=source(i,j)/(z(s)-e(i)+e(j))
        enddo
      enddo
      density=matmul(v,second)+matmul(mass,first)
      expected(s)=cmplx(0d0,1d0,8)*sum([(density(i,i),i=1,n)])/volume
      ! Independent finite-dimensional Ward identity. The position impulse
      ! uses one resolvent and contains no velocity/contact matrix elements.
      source=matmul(position,rho)-matmul(rho,position)
      do j=1,n
        do i=1,n
          density(i,j)=source(i,j)/(z(s)-e(i)+e(j))
        enddo
      enddo
      density=matmul(position,density)
      if(abs(expected(s)-cmplx(0d0,1d0,8)*z(s)*sum([(density(i,i),i=1,n)])/volume)>1d-10) &
        error stop 'nonstationary density/current Ward identity'
    enddo
    if(maxval(abs(sigma-expected))>1d-9)error stop 'two-resolvent optical spectral mismatch'
    call stationary_current_reference(apply,velocity,acceleration,psi,occ(1),volume,1d-9,dc,alpha,status)
    if(coherent==0)then
      wanted_alpha=0d0
      do i=2,n
        wanted_alpha=wanted_alpha+2d0*occ(1)*abs(v(i,1))**2/(volume*(e(i)-e(1))**3)
      enddo
      if(status/=0.or.abs(dc)>1d-10.or.abs(alpha-wanted_alpha)>1d-10) &
        error stop 'stationary conductivity sum rule/reference mismatch'
    else
      if(status/=2)error stop 'coherent state accepted as stationary DC reference'
    endif
    seed=0d0
    do i=1,n
      seed(i,i)=1d0
    enddo
    call coherent_optical(apply,velocity,acceleration,psi,occ,vertices,z,n,1d-10,volume, &
      sigma,charge,residual,actions,status,seed)
    if(status/=0.or.maxval(abs(sigma-expected))>1d-9)error stop 'zero-source seeding changes optical response'
    call propagate(.01d0,td,norm_error)
    call propagate(.005d0,fine,norm_error)
    if(maxval(abs(sigma-fine))>2d-7.or.maxval(abs(td-fine))>2d-7)error stop 'independent optical time response'
    if(norm_error>1d-8)error stop 'optical reference norm drift'
    print '(a,i2,3es17.8)','COHERENT_OPTICAL_MODEL state/spectral/time/residual ',coherent, &
      maxval(abs(sigma-expected)),maxval(abs(sigma-fine)),maxval(residual)
  enddo
  do fail=1,4
    call coherent_optical(apply,velocity,acceleration,psi,occ,vertices,z,n,1d-10,volume, &
      sigma,charge,residual,actions,status)
    if(status==0)error stop 'invalid optical operator accepted'
  enddo
  fail=0
  call coherent_optical(apply,velocity,acceleration,psi,occ,vertices,z,1,1d-10,volume, &
    sigma,charge,residual,actions,status)
  if(status==0)error stop 'unconverged optical cap accepted'
  v=0d0;mass=0d0;occ=.37d0
  do i=1,n
    mass(i,i)=1d0
  enddo
  call coherent_optical(apply,velocity,acceleration,psi,occ,vertices,z,n,1d-10,volume, &
    sigma,charge,residual,actions,status)
  if(status/=0)error stop 'pure contact response rejected'
  if(maxval(abs(sigma-cmplx(0d0,1d0,8)*sum(occ)/(volume*z)))>1d-10)error stop 'Drude contact term lost'
  print *, 'COHERENT_OPTICAL_TIME_SPECTRAL_WARD_DRUDE_IEEE_PASS'
contains
  subroutine apply(x,y,code)
    implicit none
    complex(8),intent(in) :: x(:)
    complex(8),intent(out) :: y(:)
    integer,intent(out) :: code
    code=1
    if(fail==1)return
    y=matmul(h,x);code=0
    if(fail==2)y(1)=cmplx(ieee_value(0d0,ieee_quiet_nan),0d0,8)
  end subroutine
  subroutine velocity(x,y,code)
    implicit none
    complex(8),intent(in) :: x(:)
    complex(8),intent(out) :: y(:)
    integer,intent(out) :: code
    y=matmul(v,x);code=0
    if(fail==3)y(1)=cmplx(ieee_value(0d0,ieee_positive_inf),0d0,8)
  end subroutine
  subroutine acceleration(x,y,code)
    implicit none
    complex(8),intent(in) :: x(:)
    complex(8),intent(out) :: y(:)
    integer,intent(out) :: code
    code=1
    if(fail==4)return
    y=matmul(mass,x);code=0
  end subroutine
  subroutine derivative(state,rhs)
    implicit none
    complex(8),intent(in) :: state(n,2)
    complex(8),intent(out) :: rhs(n,2)
    rhs(:,1)=cmplx(0d0,-1d0,8)*matmul(h,state(:,1))
    rhs(:,2)=cmplx(0d0,-1d0,8)*(matmul(h,state(:,2))+matmul(v,state(:,1)))
  end subroutine
  subroutine propagate(dt,value,drift)
    implicit none
    real(8),intent(in) :: dt
    complex(8),intent(out) :: value(nz)
    real(8),intent(out) :: drift
    complex(8) :: state(n,2),k1(n,2),k2(n,2),k3(n,2),k4(n,2)
    real(8) :: time,weight,current
    integer :: step,steps
    steps=nint(400d0/dt);state(:,1)=psi(:,1);state(:,2)=0d0;value=0d0
    do step=0,steps
      time=step*dt;weight=2d0
      if(mod(step,2)==1)weight=4d0
      if(step==0.or.step==steps)weight=1d0
      current=occ(1)*(2d0*real(dot_product(state(:,1),matmul(v,state(:,2))),8)+ &
        real(dot_product(state(:,1),matmul(mass,state(:,1))),8))/volume
      value=value+weight*exp(cmplx(0d0,1d0,8)*z*time)*current
      if(step==steps)exit
      call derivative(state,k1);call derivative(state+dt*k1/2d0,k2)
      call derivative(state+dt*k2/2d0,k3);call derivative(state+dt*k3,k4)
      state=state+dt*(k1+2d0*k2+2d0*k3+k4)/6d0
    enddo
    value=value*dt/3d0;drift=abs(sum(abs(state(:,1))**2)-1d0)
  end subroutine
end program
