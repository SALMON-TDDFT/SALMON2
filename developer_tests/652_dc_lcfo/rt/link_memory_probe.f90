program link_memory_probe
 use mpi
 use iso_c_binding,only:c_int64_t
 use lcfo_mlwf_links
 implicit none
 interface
  function peak_rss_bytes() bind(C) result(bytes)
   import c_int64_t
   integer(c_int64_t)::bytes
  end function
 end interface
 complex(8),allocatable :: grid(:,:),shifted(:,:),local(:,:,:,:),raw(:,:,:,:)
 real(8),allocatable :: position(:,:)
 real(8) :: length(3),delta,checksum
 integer :: rank,np,ierr,n,ng,i,j,a
 integer(c_int64_t)::rss
 character(16)::mode,arg
 call MPI_Init(ierr);call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr);call MPI_Comm_size(MPI_COMM_WORLD,np,ierr)
 call get_command_argument(1,mode);call get_command_argument(2,arg);read(arg,*)n
 ng=4096;length=[20d0,23d0,26d0]
 allocate(grid(ng,n),position(3,ng))
 do i=1,ng
  position(:,i)=[mod(i,16),mod(i/16,16),mod(i/256,16)]*.42d0
  do j=1,n
   grid(i,j)=cmplx(sin(.13d0*(i+j+rank)),cos(.07d0*(i-j+rank)),8)/sqrt(dble(ng*np))
  enddo
 enddo
 if(trim(mode)=='dense')then
  allocate(local(n,n,6,1),raw(n,n,6,1),shifted(ng,n))
  do a=1,3
   delta=2*acos(-1d0)/length(a)
   do j=1,n;shifted(:,j)=grid(:,j)*exp(cmplx(0d0,-delta*position(a,:),8));enddo
   local(:,:,a,1)=matmul(conjg(transpose(grid)),shifted)*.037d0
   local(:,:,a+3,1)=conjg(transpose(local(:,:,a,1)))
  enddo
  call MPI_Allreduce(local,raw,size(raw),MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
 else if(trim(mode)=='tiled')then
  call lcfo_initial_links(grid,position,length,.037d0,MPI_COMM_WORLD,raw)
 else
  error stop 'unknown mode'
 endif
 checksum=0d0;if(rank==0)checksum=sum(abs(raw))
 rss=peak_rss_bytes()
 write(*,'(a,1x,a,2i8,i18,es24.14)')'LINK_MEMORY',trim(mode),n,rank,rss,checksum
 call MPI_Finalize(ierr)
end program
