module misc_routines
 implicit none
contains
 function get_wtime() result(t)
 use mpi
 implicit none
 real(8)::t
 t=MPI_Wtime()
 end function
end module
module backup_routines
 implicit none
contains
 subroutine load_value(fd,a)
 implicit none
 integer,intent(in)::fd
 real(8),allocatable,intent(out)::a(:,:)
 error stop 'unused checkpoint reader'
 end subroutine
 subroutine save_value(fd,a)
 implicit none
 integer,intent(in)::fd
 real(8),intent(in)::a(:,:)
 error stop 'unused checkpoint writer'
 end subroutine
end module
module parallelization
 implicit none
 integer::nproc_group_global,nproc_id_global,nproc_size_global
end module
module communication
 use mpi
 implicit none
contains
 subroutine comm_get_min(a,b,n,c)
 implicit none
 integer,intent(in)::n,c
 real(8),intent(in)::a(n)
 real(8),intent(out)::b(n)
 integer::ierr
 call MPI_Allreduce(a,b,n,MPI_DOUBLE_PRECISION,MPI_MIN,c,ierr)
 end subroutine
 subroutine comm_get_max(a,b,n,c)
 implicit none
 integer,intent(in)::n,c
 real(8),intent(in)::a(n)
 real(8),intent(out)::b(n)
 integer::ierr
 call MPI_Allreduce(a,b,n,MPI_DOUBLE_PRECISION,MPI_MAX,c,ierr)
 end subroutine
 subroutine comm_summation(a,b,n,c)
 implicit none
 integer,intent(in)::n,c
 real(8),intent(in)::a(n)
 real(8),intent(out)::b(n)
 integer::ierr
 call MPI_Allreduce(a,b,n,MPI_DOUBLE_PRECISION,MPI_SUM,c,ierr)
 end subroutine
 logical function comm_is_root(rank)
 implicit none
 integer,intent(in)::rank
 comm_is_root=rank==0
 end function
end module
