program seed_stream_memory_probe
 use mpi
 use iso_c_binding, only: c_int64_t
 use lcfo_dist_rows, only: lcfo_gather_root
 use lcfo_seed, only: lcfo_seed_gamma
 use hse_wannier_gauge, only: gauge_seed_gamma
 implicit none
 interface
  function peak_rss_bytes() bind(C) result(bytes)
   import c_int64_t
   integer(c_int64_t) :: bytes
  end function
 end interface
 integer :: rank,np,ierr,n,nb,lo,i,j,iu,status,p
 integer,allocatable :: counts(:)
 integer(c_int64_t) :: rss
 integer(c_int64_t),allocatable :: peaks(:)
 complex(8),allocatable :: local(:,:),full(:,:),u(:,:),gram(:,:)
 real(8) :: error,started,seed_seconds
 character(32) :: mode,arg,backend
 call MPI_Init(ierr);call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr);call MPI_Comm_size(MPI_COMM_WORLD,np,ierr)
 call get_command_argument(1,mode);call get_command_argument(2,arg);read(arg,*) n
 call get_command_argument(3,backend)
 nb=16*n;allocate(counts(np),peaks(np));counts=nb/np;counts(np)=nb-sum(counts(:np-1))
 lo=sum(counts(:rank));allocate(local(counts(rank+1),n),u(n,n));u=0d0
 do j=1,n;do i=1,size(local,1)
  local(i,j)=cmplx(sin(0.731d0*(lo+i)*j),cos(0.219d0*(lo+i)*(j+1)),8)/sqrt(real(nb,8))
 enddo;enddo
 iu=0;status=0
 if(rank==0)open(newunit=iu,file='seed-coeff.bin',access='stream',form='unformatted',status='replace')
 started=MPI_Wtime()
 select case(trim(mode))
 case('reference')
  call lcfo_gather_root(local,counts,MPI_COMM_WORLD,full)
  if(rank==0)then
   call gauge_seed_gamma(full,u,status)
   write(iu)full
  endif
  deallocate(full)
 case('streamed')
  call lcfo_seed_gamma(local,counts,MPI_COMM_WORLD,u,status,snapshot_unit=iu,distributed=backend=='distributed')
 case default
  stop 2
 end select
 seed_seconds=MPI_Wtime()-started
 if(status/=0)stop 3
 if(rank==0)then
  close(iu);allocate(gram(n,n));gram=matmul(conjg(transpose(u)),u)
  do j=1,n;gram(j,j)=gram(j,j)-1d0;enddo
  error=maxval(abs(gram));if(error>1d-10)stop 4
  open(newunit=iu,file='seed-u.bin',access='stream',form='unformatted',status='replace');write(iu)u;close(iu)
 endif
 rss=peak_rss_bytes()
 call MPI_Gather(rss,1,MPI_INTEGER8,peaks,1,MPI_INTEGER8,0,MPI_COMM_WORLD,ierr)
 if(rank==0)then
  do p=1,np
   write(*,'(a,3(1x,i0),2(1x,es24.16))')trim(mode),n,p-1,peaks(p),error,seed_seconds
  enddo
 endif
 call MPI_Finalize(ierr)
end program
