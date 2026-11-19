! SALMONpy MPI attachment ABI v1. Explicitly built against the host MPI.
module salmonpy_parallel_bridge
  use iso_c_binding
  use communication
  implicit none
contains
  integer(c_int) function sp_attach(handle) bind(C)
    integer(c_int64_t),value :: handle
    integer :: ierr
    ! Current SALMON uses default INTEGER MPI handles. Reject incompatible kinds.
    if(kind(0)/=c_int.or.handle>huge(0).or.handle<(-int(huge(0),c_int64_t)-1)) then
      sp_attach=-1
      return
    end if
    call comm_attach(int(handle),ierr)
    sp_attach=ierr
  end function
  integer(c_int) function sp_detach() bind(C)
    integer :: ierr
    call comm_detach(ierr)
    sp_detach=ierr
  end function
  subroutine sp_info(rank,n) bind(C)
    integer(c_int),intent(out) :: rank,n
    integer :: comm,r,s
    call comm_get_globalinfo(comm,r,s)
    rank=r; n=s
  end subroutine
  subroutine sp_barrier() bind(C)
    call comm_sync_all
  end subroutine
  real(c_double) function sp_sum(value) bind(C)
    real(c_double),value :: value
    integer :: comm,rank,n
    real(8) :: output
    call comm_get_globalinfo(comm,rank,n)
    call comm_summation(value,output,comm)
    sp_sum=output
  end function
  subroutine sp_library(buffer,capacity,initialized) bind(C)
    use mpi, only: MPI_MAX_LIBRARY_VERSION_STRING
    integer(c_int),value :: capacity
    character(c_char),intent(out) :: buffer(capacity)
    integer(c_int),intent(out) :: initialized
    character(MPI_MAX_LIBRARY_VERSION_STRING) :: version
    integer :: n,i,ierr
    logical :: active,finished
    call MPI_Initialized(active,ierr)
    call MPI_Finalized(finished,ierr)
    initialized=0
    if(active.and..not.finished) initialized=1
    call MPI_Get_library_version(version,n,ierr)
    do i=1,min(n,capacity-1)
      buffer(i)=version(i:i)
    end do
    buffer(min(n,capacity-1)+1)=c_null_char
  end subroutine
end module
