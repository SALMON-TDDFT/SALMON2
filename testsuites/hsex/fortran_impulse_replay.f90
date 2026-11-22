program replay_impulse_operator
  use hsex_impulse_resolvent, only: impulse_model,impulse_resolve,impulse_trajectory
  implicit none
  type(impulse_model) :: model
  character(1024) :: path,argument
  character(4) :: magic
  integer :: n,nt,unit,status,j,has_factor
  real(8) :: dt,contact
  real(8),allocatable :: current(:)
  complex(8) :: shifts(3),sigma(3)
  real(8) :: residual(3)
  call get_command_argument(1,path)
  call get_command_argument(2,argument);read(argument,*)nt
  call get_command_argument(3,argument);read(argument,*)dt
  open(newunit=unit,file=trim(path),access='stream',form='unformatted',status='old')
  read(unit)magic
  if(magic=='HIO2')then
    read(unit)n
  else
    rewind(unit);read(unit)n
  endif
  if(n<1.or.n>8192.or.nt<1.or.dt<=0d0)error stop 'invalid operator record or time grid'
  allocate(model%h(n,n),model%g(n),current(0:nt))
  read(unit)model%source_norm,model%next_norm,contact,model%h,model%g
  if(magic=='HIO2')then
    read(unit)model%validation_residual,has_factor
    if(has_factor==1)then
      allocate(model%residual_factor(n,n));read(unit)model%residual_factor
    else if(has_factor/=0)then
      error stop 'invalid residual factor in operator record'
    endif
  endif
  close(unit)
  shifts=[cmplx(.0007349864435d0,0d0,8),cmplx(.0007349864435d0,-.03d0,8),cmplx(.0007349864435d0,-.07d0,8)]
  call impulse_resolve(model,shifts,sigma,residual,status)
  print '(a,i6,3es20.10)','REPLAY_FREQUENCY_STATUS/residual ',status,residual
  if(status/=0)error stop 'frequency replay failed'
  call impulse_trajectory(model,contact,dt,current,status)
  print '(a,i6)','REPLAY_WINDOW_STATUS ',status
  if(status/=0)error stop 'window replay failed'
  open(newunit=unit,file='operator-time.dat',status='new')
  do j=0,nt
    write(unit,'(2es24.15)')j*dt,current(j)
  enddo
  close(unit)
end program
