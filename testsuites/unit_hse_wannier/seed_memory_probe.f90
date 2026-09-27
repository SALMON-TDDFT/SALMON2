program seed_memory_probe
 use iso_c_binding, only: c_int64_t
 use hse_wannier_gauge, only: gauge_seed,gauge_seed_gamma
 implicit none
 interface
  function peak_rss_bytes() bind(C) result(bytes)
   import c_int64_t
   integer(c_int64_t) :: bytes
  end function
 end interface
 complex(8),allocatable :: coeff(:,:),u(:,:,:),gram(:,:)
 real(8),allocatable :: position(:,:)
 real(8) :: err
 integer :: n,ng,i,j,status
 character(32) :: mode,arg
 call get_command_argument(1,mode);call get_command_argument(2,arg);read(arg,*) n
 ng=16*n
 allocate(coeff(ng,n),u(n,n,1))
 do j=1,n;do i=1,ng
  coeff(i,j)=cmplx(sin(0.731d0*i*j),cos(0.219d0*i*(j+1)),8)/sqrt(real(ng,8))
 enddo;enddo
 select case(trim(mode))
 case('reference')
  allocate(position(3,ng));position=0d0
  call gauge_seed(reshape(coeff,[ng,n,1]),position,reshape([0d0,0d0,0d0],[3,1]),u,status)
 case('gamma2d')
  call gauge_seed_gamma(coeff,u(:,:,1),status)
 case default
  stop 2
 end select
 if(status/=0)stop 3
 allocate(gram(n,n));gram=matmul(conjg(transpose(u(:,:,1))),u(:,:,1))
 do j=1,n;gram(j,j)=gram(j,j)-1d0;enddo
 err=maxval(abs(gram));if(err>1d-10)stop 4
 write(*,'(a,1x,i0,1x,i0,1x,es24.16)')trim(mode),n,peak_rss_bytes(),err
 open(unit=17,file='seed-u.bin',access='stream',form='unformatted',status='replace')
 write(17)u
 close(17)
end program
