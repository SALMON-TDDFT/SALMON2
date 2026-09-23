program native_semilocal_probe
  use hse_semilocal
  implicit none
  real(8) :: rho(4),sigma(4),eps(4),vrho(4),vsigma(4)
  integer :: ierr,i
  rho=[0d0,1d-4,.1d0,2d0];sigma=[0d0,.001d0,.2d0,3d0]
  call hse_semilocal_evaluate(rho,sigma,eps,vrho,vsigma,ierr)
  if(ierr/=0)error stop 'HSE semilocal failed'
  do i=1,4
    write(*,'(3es25.16)')eps(i),vrho(i),vsigma(i)
  enddo
end program
