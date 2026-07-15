program test_dg_wpw_owner_exchange_mpi
  use mpi
  use rt_dg_wpw_owner_exchange, only: reduce_w_partial_to_owners, fetch_rows_from_owners
  implicit none

  integer, parameter :: nrow = 7, nvec = 2
  integer :: ierr, rank_id, nrank, i, nowned, info, min_info, max_info
  integer :: row_ids(nrow)
  integer, allocatable :: owned_ids(:)
  complex(8) :: partial(nrow,nvec)
  complex(8), allocatable :: result(:,:)
  complex(8), allocatable :: owned_x(:,:), fetched_x(:,:)
  integer :: requested_ids(2)
  real(8) :: expected

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank_id, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nrank, ierr)
  if (nrank /= 3) then
    if (rank_id == 0) write(*,'(a)') 'FAIL owner exchange fixture requires 3 ranks'
    call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
  end if

  do i = 1, nrow
    row_ids(i) = i
    partial(i,1) = cmplx(dble((rank_id + 1) * i), dble(rank_id), kind=8)
    partial(i,2) = cmplx(-dble((rank_id + 1) * i), dble(2 * rank_id), kind=8)
  end do
  nowned = 0
  do i = 1, nrow
    if (cyclic_owner(i, nrank, nrow, info) == rank_id) nowned = nowned + 1
  end do
  allocate(owned_ids(nowned), result(nowned,nvec), owned_x(nowned,nvec))
  nowned = 0
  do i = 1, nrow
    if (cyclic_owner(i, nrank, nrow, info) /= rank_id) cycle
    nowned = nowned + 1
    owned_ids(nowned) = i
    owned_x(nowned,1) = cmplx(10.0d0 * dble(i) + 1.0d0, -dble(i), kind=8)
    owned_x(nowned,2) = cmplx(10.0d0 * dble(i) + 2.0d0, dble(i), kind=8)
  end do

  call reduce_w_partial_to_owners(MPI_COMM_WORLD, rank_id, nrank, cyclic_owner, nrow, row_ids, partial, &
                                  owned_ids, result, info)
  if (info /= 0) call MPI_Abort(MPI_COMM_WORLD, 10 + info, ierr)
  do i = 1, size(owned_ids)
    expected = 6.0d0 * dble(owned_ids(i))
    if (abs(result(i,1) - cmplx(expected, 3.0d0, kind=8)) > 1.0d-12) &
      call MPI_Abort(MPI_COMM_WORLD, 30, ierr)
    if (abs(result(i,2) - cmplx(-expected, 6.0d0, kind=8)) > 1.0d-12) &
      call MPI_Abort(MPI_COMM_WORLD, 31, ierr)
  end do

  ! A missing owner-local row must fail collectively with the same nonzero code.
  if (rank_id == 0) then
    call reduce_w_partial_to_owners(MPI_COMM_WORLD, rank_id, nrank, cyclic_owner, nrow, row_ids, partial, &
                                    owned_ids(2:), result(2:,:), info)
  else
    call reduce_w_partial_to_owners(MPI_COMM_WORLD, rank_id, nrank, cyclic_owner, nrow, row_ids, partial, &
                                    owned_ids, result, info)
  end if
  call MPI_Allreduce(info, min_info, 1, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierr)
  call MPI_Allreduce(info, max_info, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
  if (min_info == 0 .or. min_info /= max_info) call MPI_Abort(MPI_COMM_WORLD, 32, ierr)

  select case (rank_id)
  case (0)
    requested_ids = [2, 6]
  case (1)
    requested_ids = [1, 7]
  case default
    requested_ids = [3, 5]
  end select
  allocate(fetched_x(2,nvec))
  call fetch_rows_from_owners(MPI_COMM_WORLD, rank_id, nrank, cyclic_owner, nrow, owned_ids, owned_x, &
                              requested_ids, fetched_x, info)
  if (info /= 0) call MPI_Abort(MPI_COMM_WORLD, 40 + info, ierr)
  do i = 1, 2
    if (abs(fetched_x(i,1) - cmplx(10.0d0 * dble(requested_ids(i)) + 1.0d0, &
                                   -dble(requested_ids(i)), kind=8)) > 1.0d-12) &
      call MPI_Abort(MPI_COMM_WORLD, 60, ierr)
    if (abs(fetched_x(i,2) - cmplx(10.0d0 * dble(requested_ids(i)) + 2.0d0, &
                                    dble(requested_ids(i)), kind=8)) > 1.0d-12) &
      call MPI_Abort(MPI_COMM_WORLD, 61, ierr)
  end do

  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  if (rank_id == 0) write(*,'(a)') 'PASS owner-targeted WP reduce/fetch MPI fixture'
  call MPI_Finalize(ierr)

contains

  integer function cyclic_owner(row_id, nrank_arg, context, owner_info) result(owner)
    integer, intent(in) :: row_id, nrank_arg
    class(*), intent(in) :: context
    integer, intent(out) :: owner_info

    owner = -1
    owner_info = 0
    select type (nrows => context)
    type is (integer)
      if (row_id < 1 .or. row_id > nrows .or. nrank_arg <= 0) then
        owner_info = 1
      else
        owner = modulo(row_id - 1, nrank_arg)
      end if
    class default
      owner_info = 2
    end select
  end function cyclic_owner
end program test_dg_wpw_owner_exchange_mpi
