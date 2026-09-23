program test_update
  use tdcdft_lrc, only: advance_xc_field,proca_coefficients,instant_screening,advance_polarization_field,elf_value,elf_alpha
  implicit none
  real(8) :: old(3),now(3),next(3),j(3),dt,t,err(2),exact,alpha,gamma
  real(8) :: a(3),e(3),pol(3),response,coupling
  integer :: n,k,steps
  ! Current-corrected local ELF, HEG endpoint and unbounded normalized coupling.
  exact=3d0/5d0*(6d0*acos(-1d0)**2)**(2d0/3d0)
  if(abs(elf_value(1d0,exact,[0d0,0d0,0d0],[0d0,0d0,0d0])-.5d0)>1d-14) error stop 'ELF HEG'
  if(abs(elf_value(1d0,exact+4d0,[0d0,0d0,0d0],[2d0,0d0,0d0])-.5d0)>1d-14) error stop 'ELF boost'
  if(abs(elf_value(1d0,5d0,[1d0,0d0,0d0],[2d0,0d0,0d0])-1d0)>1d-14) error stop 'ELF single orbital'
  if(abs(elf_alpha(.2d0,.07d0,.07d0)-.2d0)>1d-14) error stop 'ELF initial'
  if(elf_alpha(.2d0,0d0,.07d0)/=0d0) error stop 'ELF zero'
  if(abs(elf_alpha(.2d0,.14d0,.07d0)-.4d0)>1d-14) error stop 'ELF no upper clip'
  ! Known integral of a'=-(0.2+0.1*t)*sin(t), j=-cos(t).
  do k=1,2
    steps=100*2**(k-1); dt=1d0/steps
    now=0.2d0*(cos(dt)-1d0)+0.1d0*(dt*cos(dt)-sin(dt))
    do n=1,steps-1
      t=n*dt; pol=sin(t); j=-cos(t)
      call advance_polarization_field(dt,0.2d0+0.1d0*(t-dt),0.2d0+0.1d0*t,pol,j,now,next)
      now=next
    end do
    exact=0.2d0*(cos(1d0)-1d0)+0.1d0*(cos(1d0)-sin(1d0))
    err(k)=maxval(abs(now-exact))
  end do
  if(err(1)/err(2)<3.8d0.or.err(1)/err(2)>4.2d0) error stop 'varying alpha second order'
  ! Once the new alpha is held, the centered electric field is alpha*P, without an offset.
  dt=0.01d0; pol=0.3d0; j=0d0; now=1d0
  call advance_polarization_field(dt,0.2d0,0.01d0,pol,j,now,next)
  old=next
  call advance_polarization_field(dt,0.01d0,0.01d0,pol,j,old,now)
  call advance_polarization_field(dt,0.01d0,0.01d0,pol,j,now,next)
  if(maxval(abs(-(next-old)/(2d0*dt)-0.01d0*pol))>1d-12) error stop 'no residual field offset'
  ! Constant alpha is identical to the second-order acceleration update.
  pol=0d0; j=[1d0,-2d0,0d0]; old=0d0; now=0.5d0*0.2d0*dt**2*j
  do n=1,100
    pol=pol-dt*j
    call advance_polarization_field(dt,0.2d0,0.2d0,pol,j,now,a)
    call advance_xc_field(dt,0.2d0,0d0,0d0,j,old,now,next)
    if(maxval(abs(a-next))>1d-13) error stop 'polarization fixed-alpha limit'
    old=now; now=next
  end do
  response=0d0
  coupling=0.2d0
  do n=0,100
    t=n*acos(-1d0)/50d0
    a=[sin(t),0d0,0d0]
    e=[-0.2d0*cos(t),0d0,0d0]
    j=0.03d0*a
    pol=-0.03d0*e/0.2d0**2
    call instant_screening(a,e,j,pol,0.2d0,0.01d0,1d0,1d-10,0.2d0,response,coupling)
    if(abs(response-0.03d0)>1d-14) error stop 'instant Drude quadratures'
    if(abs(coupling-0.2d0/(1d0+4d0*acos(-1d0)*0.02d0/0.04d0))>1d-14) &
      error stop 'screening closure'
  end do
  a=0d0; e=0d0
  call instant_screening(a,e,j,pol,0.2d0,0.01d0,1d0,1d-10,0.2d0,response,coupling)
  if(abs(response-0.03d0)>1d-14) error stop 'hold without field'
  a=[1d0,0d0,0d0]; j=-a
  call instant_screening(a,e,j,pol,0.2d0,0.01d0,1d0,1d-10,0.2d0,response,coupling)
  if(coupling/=0.2d0) error stop 'no antiscreening'
  j=0.01d0*a
  call instant_screening(a,e,j,pol,0.2d0,0.01d0,1d0,1d-10,0.2d0,response,coupling)
  if(coupling/=0.2d0) error stop 'calibrated weak reference'
  j=a
  call instant_screening(a,e,j,pol,0.2d0,0.01d0,0d0,1d-10,0.2d0,response,coupling)
  if(coupling/=0.2d0) error stop 'zero screening strength'
  ! A pulse tail with remanent P is not a stationary Drude response.
  ! Confirm bounded alpha and floor hold, not physical accuracy of the tail estimate.
  a=0d0; e=[-1d-6,0d0,0d0]; pol=[1d-3,0d0,0d0]; j=0d0
  call instant_screening(a,e,j,pol,0.2d0,0d0,1d0,1d-8,0.2d0,response,coupling)
  if(abs(response-40d0)>1d-10.or.coupling<=0d0.or.coupling>=0.2d0) error stop 'residual polarization tail'
  e=0d0
  call instant_screening(a,e,j,pol,0.2d0,0d0,1d0,1d-8,0.2d0,response,coupling)
  if(abs(response-40d0)>1d-10) error stop 'tail freeze'
  call proca_coefficients(-20d0*acos(-1d0),-0.2d0,alpha,gamma)
  if(abs(alpha-0.2d0)>1d-14) error stop 'Si sign and normalization'
  if(abs(gamma-0.01d0/acos(-1d0))>1d-14) error stop 'Si restoring term'
  call proca_coefficients(4d0*acos(-1d0),1d0,alpha,gamma)
  if(abs(alpha+1d0)>1d-14.or.gamma<=0d0) error stop 'positive a2 convention'
  call proca_coefficients(-20d0*acos(-1d0),0d0,alpha,gamma)
  if(abs(gamma)>tiny(1d0)) error stop 'massless reference'
  dt=0.01d0
  j=[1d0,-2d0,0d0]
  old=0d0
  now=0.5d0*0.2d0*j*dt**2
  do n=1,100
    call advance_xc_field(dt,0.2d0,0d0,0d0,j,old,now,next)
    old=now
    now=next
  end do
  if(maxval(abs(now-0.5d0*0.2d0*j*(101*dt)**2))>1d-12) error stop 'constant drive'
  old=0d0
  now=0d0
  call advance_xc_field(dt,0d0,0d0,0d0,j,old,now,next)
  if(maxval(abs(next))>tiny(1d0)) error stop 'zero coupling'
  ! Homogeneous damped oscillator: A=exp(-0.1*t)*cos(t), beta=0.2, gamma=1.01.
  do k=1,2
    steps=100*2**(k-1)
    dt=1d0/steps
    old=1d0
    now=exp(-0.1d0*dt)*cos(dt)
    j=0d0
    do n=1,steps-1
      call advance_xc_field(dt,0.4d0,0.2d0,1.01d0,j,old,now,next)
      old=now
      now=next
    end do
    exact=exp(-0.1d0)*cos(1d0)
    err(k)=maxval(abs(now-exact))
  end do
  if(err(1)/err(2)<3.8d0.or.err(1)/err(2)>4.2d0) error stop 'second order convergence'
  ! Smooth inhomogeneous LRC: A(t)=alpha*(t-sin(t)).
  dt=0.001d0
  old=0d0
  now=0.2d0*(dt-sin(dt))
  do n=1,999
    t=n*dt
    j=sin(t)
    call advance_xc_field(dt,0.2d0,0d0,0d0,j,old,now,next)
    old=now
    now=next
  end do
  if(maxval(abs(now-0.2d0*(1d0-sin(1d0))))>1d-8) error stop 'sinusoidal drive'
  print *, 'analytic TDCDFT update tests passed'
end program test_update
