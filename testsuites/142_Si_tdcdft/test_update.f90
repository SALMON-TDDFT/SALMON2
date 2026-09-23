program test_update
  use tdcdft_lrc, only: advance_xc_field
  implicit none
  real(8) :: old(3),now(3),next(3),j(3),dt,t,err(2),exact
  integer :: n,k,steps
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
