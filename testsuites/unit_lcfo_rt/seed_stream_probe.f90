program seed_stream_probe
 use mpi
 use lcfo_dist_rows, only: lcfo_gather_root
 use lcfo_seed, only: lcfo_seed_gamma
 use hse_wannier_gauge, only: gauge_seed_gamma
 implicit none
 integer :: rank,np,ierr,n,lo,i,j,iu,status,reference,a,lowest,highest
 integer,allocatable :: counts(:)
 complex(8),allocatable :: local(:,:),full(:,:),u(:,:),v(:,:),stream(:,:),saved(:,:)
 call MPI_Init(ierr);call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr);call MPI_Comm_size(MPI_COMM_WORLD,np,ierr)
 allocate(counts(np))
 do a=1,5
  n=5;if(a==2)n=65
  if(a==5)n=0
  counts=n+3;counts(1)=0;counts(np)=n+5
  if(a==4)then
   counts=1;counts(1)=0
  endif
  lo=sum(counts(:rank))
  allocate(local(counts(rank+1),n),u(n,n),v(n,n),saved(counts(rank+1),n))
  do j=1,n;do i=1,size(local,1)
   local(i,j)=cmplx(sin(0.731d0*(lo+i)*j),cos(0.219d0*(lo+i)*(j+1)),8)
  enddo;enddo
  if(a==3)local(:,n)=local(:,1)
  saved=local;iu=0
  call lcfo_gather_root(local,counts,MPI_COMM_WORLD,full)
  if(rank==0)then
   call gauge_seed_gamma(full,v,reference)
   open(newunit=iu,file='seed-stream.bin',access='stream',form='unformatted',status='replace')
  endif
  call lcfo_seed_gamma(local,counts,MPI_COMM_WORLD,u,status,snapshot_unit=iu)
  call MPI_Allreduce(status,lowest,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,ierr)
  call MPI_Allreduce(status,highest,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
  if(lowest/=highest)error stop 'Inconsistent seed status across ranks'
  if(a>=3.and.status==0)error stop 'Invalid seed accepted'
  if(any(local/=saved))error stop 'Changed local coefficients'
  if(rank==0)then
   close(iu)
   if(status/=reference)error stop 'Seed status mismatch'
   if(a==3.and.status==0)error stop 'Singular seed accepted'
   if(status==0)then
    if(maxval(abs(u-v))>1d-11)error stop 'Distributed seed mismatch'
   endif
   allocate(stream(sum(counts),n))
   open(newunit=iu,file='seed-stream.bin',access='stream',form='unformatted',status='old')
   read(iu)stream;close(iu)
   if(any(stream/=full))error stop 'Coefficient snapshot order changed'
   deallocate(stream)
  endif
  if(a==1)then
   call lcfo_seed_gamma(local,counts,MPI_COMM_WORLD,u,status)
   if(status/=0)error stop 'Seed without snapshot failed'
   if(rank==0)then
    if(maxval(abs(u-v))>1d-11)error stop 'Seed without snapshot mismatch'
   endif
  endif
  deallocate(local,full,u,v,saved)
 enddo
 if(rank==0)print *,'Streamed seed, snapshot, empty root, tile tail and singular fallback passed'
 call MPI_Finalize(ierr)
end program
