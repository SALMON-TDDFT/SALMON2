program test_dg_dc_local_basis_layout_mpi
  use mpi
  use dg_dc_local_basis_ground_state
  implicit none
  type(s_dg_dc_local_basis_layout) :: layout
  integer :: ierr,rank,nproc
  logical :: ok
  character(256) :: message

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nproc,ierr)
  if(nproc/=2) error stop 'test requires two ranks'

  call initialize_dg_dc_local_basis_layout(layout,rank+1,rank+2,3,101_8,202_8, &
    MPI_COMM_WORLD,ok,message)
  call require(ok,'unequal local basis layout')
  call require(layout%local_basis_count==rank+2,'local basis count preserved')
  call require(layout%global_basis_count==5,'global basis sum')
  call require(layout%global_band_count==3,'global band count independent of basis sum')
  call require(all(layout%basis_offsets==[0,2,5]),'exclusive basis offsets')
  call require(all(layout%fragment_ids==[1,2]),'fragment identity ordering')
  call require(layout%geometry_fingerprint==101_8 .and. layout%operator_fingerprint==202_8, &
    'layout fingerprints')
  call release_dg_dc_local_basis_layout(layout)

  call initialize_dg_dc_local_basis_layout(layout,rank+1,merge(0,3,rank==0),3,101_8,202_8, &
    MPI_COMM_WORLD,ok,message)
  call require(.not.ok,'zero local basis rejected collectively')

  call initialize_dg_dc_local_basis_layout(layout,rank+1,rank+2,3+rank,101_8,202_8, &
    MPI_COMM_WORLD,ok,message)
  call require(.not.ok,'rank-disagreeing global band count rejected')

  call initialize_dg_dc_local_basis_layout(layout,1,rank+2,3,101_8,202_8, &
    MPI_COMM_WORLD,ok,message)
  call require(.not.ok,'duplicate fragment identity rejected')

  call initialize_dg_dc_local_basis_layout(layout,rank+1,rank+2,3,101_8+rank,202_8, &
    MPI_COMM_WORLD,ok,message)
  call require(.not.ok,'rank-disagreeing geometry fingerprint rejected')

  if(rank==0) print '(a)','PASS DG-DC local basis layout keeps global bands independent'
  call MPI_Finalize(ierr)

contains

  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(*), intent(in) :: label
    logical :: global
    integer :: mpi_error
    call MPI_Allreduce(condition,global,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,mpi_error)
    if(.not.global) then
      if(rank==0) write(*,'(a,1x,a)') 'FAIL',trim(label)
      call MPI_Abort(MPI_COMM_WORLD,1,mpi_error)
    end if
  end subroutine require
end program test_dg_dc_local_basis_layout_mpi
