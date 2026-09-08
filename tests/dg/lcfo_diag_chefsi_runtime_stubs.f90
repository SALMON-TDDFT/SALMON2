module structures
  implicit none
  type :: s_dcdft
    integer :: i_frag=1,icomm_frag=0,icomm_tot=0,id_frag=0,id_tot=0
    integer :: isize_frag=1,isize_tot=1,n_frag=1,nstate_tot=1
  end type s_dcdft
end module structures

module timer
  implicit none
  integer, parameter :: LOG_CHEFSI_SETUP=1,LOG_CHEFSI_TOTAL=2
  integer, parameter :: LOG_CHEFSI_H_APPLY=3,LOG_CHEFSI_H_COMM_POST=4
  integer, parameter :: LOG_CHEFSI_H_DIAG=5,LOG_CHEFSI_H_HALO=6
  integer, parameter :: LOG_CHEFSI_H_RECV_WAIT=7,LOG_CHEFSI_H_SEND_WAIT=8
  integer, parameter :: LOG_CHEFSI_FILTER=9,LOG_CHEFSI_ORTHO=10
  integer, parameter :: LOG_CHEFSI_PROJECT=11,LOG_CHEFSI_PROJECT_EIGEN=12
  integer, parameter :: LOG_CHEFSI_RAYLEIGH_RITZ=13,LOG_CHEFSI_REDISTRIBUTE=14
  integer, parameter :: LOG_CHEFSI_RESIDUAL=15,LOG_CHEFSI_ROTATE=16
  integer, parameter :: LOG_CHEFSI_LANCZOS=17,LOG_CHEFSI_EXPORT=18
contains
  subroutine timer_begin(id)
    integer, intent(in) :: id
  end subroutine timer_begin

  subroutine timer_end(id)
    integer, intent(in) :: id
  end subroutine timer_end
end module timer

module eigen_subdiag_sub
  implicit none
contains
  subroutine eigen_dsyev(matrix,eigenvalue,eigenvector)
    real(8), intent(inout) :: matrix(:,:)
    real(8), intent(out) :: eigenvalue(:),eigenvector(:,:)
    real(8), allocatable :: work(:)
    real(8) :: work_query(1)
    integer :: info,n,lwork

    n=size(matrix,1)
    eigenvector=matrix
    call dsyev('V','U',n,eigenvector,n,eigenvalue,work_query,-1,info)
    if(info/=0) error stop 'test eigen_dsyev workspace query failed'
    lwork=max(1,nint(work_query(1)))
    allocate(work(lwork))
    call dsyev('V','U',n,eigenvector,n,eigenvalue,work,lwork,info)
    if(info/=0) error stop 'test eigen_dsyev failed'
    deallocate(work)
  end subroutine eigen_dsyev
end module eigen_subdiag_sub

module communication
  use mpi
  implicit none
  interface comm_bcast
    module procedure comm_bcast_real4
  end interface
  interface comm_isend
    module procedure comm_isend_real2,comm_isend_real3
  end interface
  interface comm_irecv
    module procedure comm_irecv_real2,comm_irecv_real3
  end interface
  interface comm_get_max
    module procedure comm_get_max_integer,comm_get_max_real1
  end interface
  interface comm_summation
    module procedure comm_sum_real0,comm_sum_real1,comm_sum_real2
    module procedure comm_sum_integer1,comm_sum_integer2
  end interface
contains
  integer function comm_create_group(comm,color,key)
    integer, intent(in) :: comm,color,key
    integer :: ierr
    call MPI_Comm_split(comm,color,key,comm_create_group,ierr)
    if(ierr/=MPI_SUCCESS) error stop 'test MPI_Comm_split failed'
  end function comm_create_group

  subroutine comm_free_group(comm)
    integer, intent(inout) :: comm
    integer :: ierr
    call MPI_Comm_free(comm,ierr)
    if(ierr/=MPI_SUCCESS) error stop 'test MPI_Comm_free failed'
  end subroutine comm_free_group

  subroutine comm_bcast_real4(value,comm,root)
    real(8), intent(inout) :: value(:,:,:,:)
    integer, intent(in) :: comm,root
    integer :: ierr
    call MPI_Bcast(value,size(value),MPI_DOUBLE_PRECISION,root,comm,ierr)
    if(ierr/=MPI_SUCCESS) error stop 'test MPI_Bcast failed'
  end subroutine comm_bcast_real4

  integer function comm_isend_real2(value,dest,tag,comm)
    real(8), intent(in) :: value(:,:)
    integer, intent(in) :: dest,tag,comm
    integer :: ierr
    call MPI_Isend(value,size(value),MPI_DOUBLE_PRECISION,dest,tag,comm,comm_isend_real2,ierr)
    if(ierr/=MPI_SUCCESS) error stop 'test MPI_Isend failed'
  end function comm_isend_real2

  integer function comm_isend_real3(value,dest,tag,comm)
    real(8), intent(in) :: value(:,:,:)
    integer, intent(in) :: dest,tag,comm
    integer :: ierr
    call MPI_Isend(value,size(value),MPI_DOUBLE_PRECISION,dest,tag,comm,comm_isend_real3,ierr)
    if(ierr/=MPI_SUCCESS) error stop 'test MPI_Isend failed'
  end function comm_isend_real3

  integer function comm_irecv_real2(value,source,tag,comm)
    real(8), intent(out) :: value(:,:)
    integer, intent(in) :: source,tag,comm
    integer :: ierr
    call MPI_Irecv(value,size(value),MPI_DOUBLE_PRECISION,source,tag,comm,comm_irecv_real2,ierr)
    if(ierr/=MPI_SUCCESS) error stop 'test MPI_Irecv failed'
  end function comm_irecv_real2

  integer function comm_irecv_real3(value,source,tag,comm)
    real(8), intent(out) :: value(:,:,:)
    integer, intent(in) :: source,tag,comm
    integer :: ierr
    call MPI_Irecv(value,size(value),MPI_DOUBLE_PRECISION,source,tag,comm,comm_irecv_real3,ierr)
    if(ierr/=MPI_SUCCESS) error stop 'test MPI_Irecv failed'
  end function comm_irecv_real3

  subroutine comm_wait_all(request)
    integer, intent(inout) :: request(:)
    integer :: ierr
    integer, allocatable :: statuses(:,:)
    if(size(request)==0) return
    allocate(statuses(MPI_STATUS_SIZE,size(request)))
    call MPI_Waitall(size(request),request,statuses,ierr)
    if(ierr/=MPI_SUCCESS) error stop 'test MPI_Waitall failed'
    deallocate(statuses)
  end subroutine comm_wait_all

  subroutine comm_get_max_integer(value,comm)
    integer, intent(inout) :: value
    integer, intent(in) :: comm
    integer :: ierr,input
    input=value
    call MPI_Allreduce(input,value,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS) error stop 'test integer MPI_Allreduce failed'
  end subroutine comm_get_max_integer

  subroutine comm_get_max_real1(input,output,n,comm)
    real(8), intent(in) :: input(:)
    real(8), intent(out) :: output(:)
    integer, intent(in) :: n,comm
    integer :: ierr
    call MPI_Allreduce(input,output,n,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS) error stop 'test real max MPI_Allreduce failed'
  end subroutine comm_get_max_real1

  subroutine comm_sum_real0(input,output,comm)
    real(8), intent(in) :: input
    real(8), intent(out) :: output
    integer, intent(in) :: comm
    integer :: ierr
    call MPI_Allreduce(input,output,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS) error stop 'test scalar MPI_Allreduce failed'
  end subroutine comm_sum_real0

  subroutine comm_sum_real1(input,output,n,comm)
    real(8), intent(in) :: input(:)
    real(8), intent(out) :: output(:)
    integer, intent(in) :: n,comm
    integer :: ierr
    call MPI_Allreduce(input,output,n,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS) error stop 'test real vector MPI_Allreduce failed'
  end subroutine comm_sum_real1

  subroutine comm_sum_real2(input,output,n,comm)
    real(8), intent(in) :: input(:,:)
    real(8), intent(out) :: output(:,:)
    integer, intent(in) :: n,comm
    integer :: ierr
    call MPI_Allreduce(input,output,n,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS) error stop 'test real matrix MPI_Allreduce failed'
  end subroutine comm_sum_real2

  subroutine comm_sum_integer1(input,output,n,comm)
    integer, intent(in) :: input(:)
    integer, intent(out) :: output(:)
    integer, intent(in) :: n,comm
    integer :: ierr
    call MPI_Allreduce(input,output,n,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS) error stop 'test integer vector MPI_Allreduce failed'
  end subroutine comm_sum_integer1

  subroutine comm_sum_integer2(input,output,n,comm)
    integer, intent(in) :: input(:,:)
    integer, intent(out) :: output(:,:)
    integer, intent(in) :: n,comm
    integer :: ierr
    call MPI_Allreduce(input,output,n,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS) error stop 'test integer matrix MPI_Allreduce failed'
  end subroutine comm_sum_integer2
end module communication
