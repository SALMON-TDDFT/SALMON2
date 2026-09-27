program gram_benchmark
 use mpi
 use lcfo_gram,only:lcfo_gram_error
 implicit none
 integer :: rank,np,ierr,n,j,k,r,pass,method,order(2),repeat
 complex(8),allocatable :: c(:,:)
 real(8) :: error,start,elapsed,maximum,reference
 call MPI_Init(ierr);call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr);call MPI_Comm_size(MPI_COMM_WORLD,np,ierr)
 n=16*np;repeat=200
 allocate(c(64,n))
 do j=1,n;do k=1,64
  c(k,j)=cmplx(sin(.13d0*(k+j+rank)),cos(.07d0*(k-j+rank)),8)/sqrt(dble(64*np))
 enddo;enddo
 call dense_error(reference);call lcfo_gram_error(c,MPI_COMM_WORLD,error)
 if(abs(error-reference)>1d-12)error stop 'benchmark mismatch'
 do pass=1,2
  order=[1,2];if(pass==2)order=[2,1]
  do j=1,2
   method=order(j)
   call MPI_Barrier(MPI_COMM_WORLD,ierr);start=MPI_Wtime()
   do r=1,repeat
    if(method==1)then
     call dense_error(error)
    else
     call lcfo_gram_error(c,MPI_COMM_WORLD,error)
    endif
   enddo
   elapsed=(MPI_Wtime()-start)/repeat
   call MPI_Reduce(elapsed,maximum,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
   if(rank==0)write(*,'(a,4i8,2es20.10)')'GRAM_TIMING ',np,n,pass,method,maximum,error
  enddo
 enddo
 call MPI_Finalize(ierr)
contains
 subroutine dense_error(value)
  implicit none
  real(8),intent(out) :: value
  integer :: jj
  complex(8),allocatable :: local(:,:),total(:,:)
  allocate(local(n,n),total(n,n))
  local=matmul(conjg(transpose(c)),c)
  call MPI_Allreduce(local,total,n*n,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  do jj=1,n;total(jj,jj)=total(jj,jj)-1d0;enddo
  value=maxval(abs(total))
 end subroutine
end program
