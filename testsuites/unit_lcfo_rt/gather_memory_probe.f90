program gather_memory_probe
 use mpi
 use lcfo_dist_rows
 implicit none
 integer :: rank,np,ierr,j,k,lo,buffer_elements
 integer,allocatable :: counts(:)
 complex(8),allocatable :: local(:,:),global(:,:)
 call MPI_Init(ierr);call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr);call MPI_Comm_size(MPI_COMM_WORLD,np,ierr)
 allocate(counts(np));counts=3;counts(1)=0
 lo=sum(counts(:rank));allocate(local(counts(rank+1),130))
 do j=1,130;do k=1,size(local,1);local(k,j)=cmplx(lo+k,j,8);enddo;enddo
 call lcfo_gather_root(local,counts,MPI_COMM_WORLD,global,buffer_elements)
 if(rank==0)then
  if(buffer_elements>sum(counts)*64)error stop 'unbounded gather buffer'
  do j=1,130;do k=1,sum(counts)
   if(global(k,j)/=cmplx(k,j,8))error stop 'gather transpose/tail error'
  enddo;enddo
 else
  if(size(global)/=0.or.buffer_elements/=1)error stop 'nonroot global storage'
 endif
 if(rank==0)print *,'Tiled root gather passed'
 call MPI_Finalize(ierr)
end program
