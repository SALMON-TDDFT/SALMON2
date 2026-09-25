program native_semilocal_probe
  use hse_semilocal
  implicit none
  real(8) :: rho(4),sigma(4),eps(4),vrho(4),vsigma(4)
  integer :: ierr,i
  real(8) :: omega
  character(80) :: arg
  rho=[0d0,1d-4,.1d0,2d0];sigma=[0d0,.001d0,.2d0,3d0]
  omega=.11d0
  if(command_argument_count()>0)then
    call get_command_argument(1,arg)
    read(arg,*)omega
  endif
  call hse_semilocal_evaluate(rho,sigma,eps,vrho,vsigma,ierr,omega)
  if(ierr/=0)error stop 'HSE semilocal failed'
  do i=1,4
    write(*,'(3es25.16)')eps(i),vrho(i),vsigma(i)
  enddo
end program
