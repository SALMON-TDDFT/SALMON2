module communication
 use mpi
 implicit none
 interface comm_summation
  module procedure sum_c2
 end interface
contains
 subroutine sum_c2(a,b,n,comm,dest)
 complex(8),intent(in) :: a(:,:)
 complex(8),intent(out) :: b(:,:)
 integer,intent(in) :: n,comm
 integer,optional,intent(in) :: dest
 integer :: ierr
 if(present(dest))then
 call MPI_Reduce(a,b,n,MPI_DOUBLE_COMPLEX,MPI_SUM,dest,comm,ierr)
 else
 call MPI_Allreduce(a,b,n,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
 endif
 if(ierr/=MPI_SUCCESS)error stop 'test collective failed'
 end subroutine
 subroutine comm_get_groupinfo(comm,rank,np)
 integer,intent(in)::comm
 integer,intent(out)::rank,np
 integer::ierr
 call MPI_Comm_rank(comm,rank,ierr)
 call MPI_Comm_size(comm,np,ierr)
 end subroutine
end module
