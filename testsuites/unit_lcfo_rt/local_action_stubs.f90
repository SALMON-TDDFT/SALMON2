module communication
 use mpi
 implicit none
 interface comm_summation
  module procedure sum_c2
 end interface
contains
 subroutine sum_c2(a,b,n,comm)
 complex(8),intent(in) :: a(:,:)
 complex(8),intent(out) :: b(:,:)
 integer,intent(in) :: n,comm
 integer :: ierr
 call MPI_Allreduce(a,b,n,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
 if(ierr/=MPI_SUCCESS)error stop 'test collective failed'
 end subroutine
end module
