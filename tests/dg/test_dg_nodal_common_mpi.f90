program test_dg_nodal_common_mpi
  use mpi
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  use dg_nodal_state
  use dg_nodal_interfaces
  implicit none

  type(s_dg_nodal_common_state) :: state, snapshot
  integer :: ierr, rank, nproc
  integer(8) :: owner(2,1,1)
  logical :: ok
  character(len=256) :: message

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)
  if (nproc /= 2) error stop 'test requires two MPI ranks'

  owner(:,1,1) = [int(2*rank+1,8), int(2*rank+2,8)]
  call initialize_dg_nodal_common_state(state, rank+1, [2,1,1], [4,3,3], [1,1,1], 1, 2, 1, &
                                        owner, 101_8, 202_8, MPI_COMM_WORLD, ok, message)
  call require_collective(ok, 'valid common-state allocation')
  call require_collective(state%initialized .and. .not.state%ground_state_ready, 'readiness initialization')
  call require_collective(all(shape(state%psi_core) == [2,1,1,2,1]), 'core orbital allocation')
  call require_collective(all(shape(state%occupations) == [2,1]), 'occupation allocation')
  call require_collective(size(state%faces) == nodal_face_count, 'face allocation')
  call require_collective(all(state%core_owner == owner), 'explicit core ownership')
  call require_collective(state%geometry_fingerprint == 101_8 .and. &
                          state%operator_fingerprint == 202_8, 'provenance fingerprints')

  state%psi_core = cmplx(1.0d0 + rank, -0.25d0, 8)
  state%occupations = reshape([2.0d0, 0.0d0], [2,1])
  call validate_dg_nodal_common_state_mpi(state, MPI_COMM_WORLD, ok, message)
  call require_collective(ok, 'finite state validation')

  snapshot = state
  if (rank == 0) owner(:,1,1) = [1_8,2_8]
  if (rank == 1) owner(:,1,1) = [2_8,3_8]
  call initialize_dg_nodal_common_state(state, rank+1, [2,1,1], [4,3,3], [1,1,1], 1, 2, 1, &
                                        owner, 303_8, 404_8, MPI_COMM_WORLD, ok, message)
  call require_collective(.not.ok, 'duplicate ownership rejection')
  call require_collective(state%geometry_fingerprint == snapshot%geometry_fingerprint .and. &
                          all(state%core_owner == snapshot%core_owner), 'duplicate failure is transactional')

  call initialize_dg_nodal_common_state(state, rank+1, [2,1,1], [3,3,3], [1,1,1], 1, 2, 1, &
                                        snapshot%core_owner, 303_8, 404_8, MPI_COMM_WORLD, ok, message)
  call require_collective(.not.ok, 'malformed core-buffer geometry rejection')
  call require_collective(state%geometry_fingerprint == snapshot%geometry_fingerprint, &
                          'geometry failure is transactional')

  owner = snapshot%core_owner
  call initialize_dg_nodal_common_state(state, rank+1, [2,1,1], [4,3,3], [1,1,1], 1, 2, 1, &
                                        owner, 101_8, 202_8+int(rank,8), MPI_COMM_WORLD, ok, message)
  call require_collective(.not.ok, 'rank-disagreeing provenance rejection')
  call require_collective(state%operator_fingerprint == snapshot%operator_fingerprint, &
                          'provenance disagreement is transactional')

  state%psi_core = snapshot%psi_core
  if (rank == 1) state%psi_core(1,1,1,1,1) = cmplx(ieee_value(0.0d0,ieee_quiet_nan),0.0d0,8)
  call validate_dg_nodal_common_state_mpi(state, MPI_COMM_WORLD, ok, message)
  call require_collective(.not.ok, 'nonfinite state rejection agrees on all ranks')
  state = snapshot

  call release_dg_nodal_common_state(state)
  call require_collective(.not.state%initialized .and. .not.allocated(state%psi_core) .and. &
                          .not.allocated(state%occupations) .and. .not.allocated(state%faces) .and. &
                          .not.allocated(state%core_owner), 'complete release')

  if (rank == 0) print '(a)', 'PASS GS/RT-neutral nodal common state'
  call MPI_Finalize(ierr)

contains

  subroutine require_collective(condition, label)
    logical, intent(in) :: condition
    character(*), intent(in) :: label
    logical :: all_condition
    integer :: mpi_error
    call MPI_Allreduce(condition, all_condition, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, mpi_error)
    if (.not.all_condition) then
      if (rank == 0) write(*,'(a,1x,a)') 'FAIL', trim(label)
      call MPI_Abort(MPI_COMM_WORLD, 1, mpi_error)
    end if
  end subroutine require_collective

end program test_dg_nodal_common_mpi
