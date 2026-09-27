program gamma_memory_probe
  use iso_c_binding, only: c_int64_t
  use hse_wannier_gauge, only: gauge_minimize_gamma, gauge_minimize_gamma_inplace
  implicit none
  interface
    function peak_rss_bytes() bind(C) result(bytes)
      import c_int64_t
      integer(c_int64_t) :: bytes
    end function
  end interface
  complex(8),allocatable :: u(:,:,:),links(:,:,:,:)
  real(8) :: b(3,6),weights(6),spread,gradient
  integer :: n,i,l,iterations,status
  character(32) :: mode,arg
  call get_command_argument(1,mode)
  call get_command_argument(2,arg);read(arg,*) n
  allocate(u(n,n,1),links(n,n,6,1));u=0d0;links=0d0
  b=0d0;weights=0.5d0
  do l=1,3
    b(l,l)=1d0;b(l,l+3)=-1d0
  end do
  do i=1,n
    u(i,i,1)=1d0
    do l=1,3
      links(i,i,l,1)=exp(cmplx(0d0,0.2d0*i/n*l,8))
      links(i,i,l+3,1)=conjg(links(i,i,l,1))
    end do
  end do
  ! Zero sweeps isolates initial rotation and objective/gradient memory.
  if(trim(mode)=='reference')then
    call gauge_minimize_gamma(u,links,b,weights,0,1d-7,spread,gradient,iterations,status)
  else if(trim(mode)=='inplace')then
    call gauge_minimize_gamma_inplace(u,links,b,weights,0,1d-7,spread,gradient,iterations,status)
  else
    stop 2
  end if
  if(status/=1.or.iterations/=0.or.abs(spread)>1d-8.or.gradient>1d-8) stop 3
  write(*,'(a,1x,i0,1x,i0,2(1x,es24.16))') trim(mode),n,peak_rss_bytes(),spread,gradient
end program
