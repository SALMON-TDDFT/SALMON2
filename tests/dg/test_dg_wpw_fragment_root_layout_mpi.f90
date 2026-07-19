program test_dg_wpw_fragment_root_layout_mpi
  use mpi
  use rt_dg_wpw_column_layout, only: s_dg_wpw_column_layout, initialize_wpw_fragment_root_layout, &
    wpw_fragment_root_owner
  use dg_wpw_production_context, only: s_dg_wpw_production_context, &
    initialize_dg_wpw_fragment_root_context
  implicit none
  type(s_dg_wpw_column_layout) :: layout
  type(s_dg_wpw_production_context) :: context
  integer :: ierr, world_rank, fragment_id, fragment_rank
  integer :: production_comm, production_rank, production_size, info, owner, i
  logical :: is_root

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, world_rank, ierr)
  if (ierr /= MPI_SUCCESS) error stop 1
  if (world_rank < 0 .or. world_rank > 3) error stop 2

  fragment_id = world_rank / 2 + 1
  fragment_rank = modulo(world_rank, 2)
  is_root = fragment_rank == 0
  call MPI_Comm_split(MPI_COMM_WORLD, merge(1, MPI_UNDEFINED, is_root), fragment_id, production_comm, ierr)
  if (ierr /= MPI_SUCCESS) error stop 3

  if (is_root) then
    call MPI_Comm_rank(production_comm, production_rank, ierr)
    call MPI_Comm_size(production_comm, production_size, ierr)
    if (ierr /= MPI_SUCCESS) error stop 4
    call initialize_wpw_fragment_root_layout(layout, 2, 3, fragment_id, production_rank, production_size, info)
    if (info /= 0) error stop 5
    if (production_rank /= fragment_id - 1) error stop 6
    if (trim(layout%ownership_kind) /= 'fragment_root') error stop 7
    if (any(layout%owned_column_ids /= [(3*(fragment_id-1)+i, i=1,3)])) error stop 8
    owner = wpw_fragment_root_owner(layout%owned_column_ids(2), 3, 2, info)
    if (info /= 0 .or. owner /= fragment_id - 1) error stop 9
    owner = wpw_fragment_root_owner(7, 3, 2, info)
    if (info == 0 .or. owner /= -1) error stop 10
    call initialize_dg_wpw_fragment_root_context(context, production_comm, 2, 3, 4, fragment_id, &
      [1,2], [fragment_id], [1,2], info)
    if (info /= 0) error stop 12
    if (context%owned_fragment_id /= fragment_id) error stop 13
    if (trim(context%ownership_kind) /= 'fragment_root') error stop 14
    if (any(context%owned_column_ids /= layout%owned_column_ids)) error stop 15
    if (any(context%support_column_ids /= [1,2,3,4,5,6])) error stop 16
    call MPI_Comm_free(production_comm, ierr)
  else
    if (production_comm /= MPI_COMM_NULL) error stop 11
  end if

  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  if (world_rank == 0) print '(a)', 'PASS fragment-root WPW stable-column layout'
  call MPI_Finalize(ierr)
end program test_dg_wpw_fragment_root_layout_mpi
