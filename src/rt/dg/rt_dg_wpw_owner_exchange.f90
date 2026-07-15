module rt_dg_wpw_owner_exchange
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use, intrinsic :: iso_fortran_env, only: error_unit, int64
  use mpi, only: MPI_Abort, MPI_Allreduce, MPI_Alltoall, MPI_Alltoallv, MPI_DOUBLE_COMPLEX, &
                 MPI_INTEGER, MPI_MAX, MPI_SUCCESS
  implicit none
  private

  abstract interface
    integer function row_owner_function(row_id, nrank, context, info) result(owner)
      integer, intent(in) :: row_id, nrank
      class(*), intent(in) :: context
      integer, intent(out) :: info
    end function row_owner_function
  end interface

  public :: reduce_w_partial_to_owners
  public :: fetch_rows_from_owners
  public :: row_owner_function

contains

  subroutine reduce_w_partial_to_owners(icomm, rank_id, nrank, row_owner, owner_context, row_ids, yw_partial, &
                                        owned_row_ids, owned_sum, info)
    integer, intent(in) :: icomm, rank_id, nrank
    procedure(row_owner_function) :: row_owner
    class(*), intent(in) :: owner_context
    integer, intent(in) :: row_ids(:), owned_row_ids(:)
    complex(8), intent(in) :: yw_partial(:,:)
    complex(8), intent(out) :: owned_sum(:,:)
    integer, intent(out) :: info

    integer :: local_info, global_info, ierr, irow, ivec, owner, owner_info, position
    integer :: nsend, nrecv, nvec, slot
    integer, allocatable :: send_counts(:), recv_counts(:), send_displs(:), recv_displs(:), fill(:)
    integer, allocatable :: send_ids(:), recv_ids(:), send_value_counts(:), recv_value_counts(:)
    integer, allocatable :: send_value_displs(:), recv_value_displs(:)
    complex(8), allocatable :: send_values(:,:), recv_values(:,:)

    owned_sum = (0.0d0, 0.0d0)
    local_info = 0
    if (nrank <= 0 .or. rank_id < 0 .or. rank_id >= nrank) local_info = 1
    if (size(yw_partial, 1) /= size(row_ids) .or. &
        size(owned_sum, 1) /= size(owned_row_ids) .or. &
        size(owned_sum, 2) /= size(yw_partial, 2)) local_info = max(local_info, 2)
    if (.not. strictly_increasing(row_ids) .or. .not. strictly_increasing(owned_row_ids)) &
      local_info = max(local_info, 3)
    do irow = 1, size(row_ids)
      owner = row_owner(row_ids(irow), nrank, owner_context, owner_info)
      if (owner_info /= 0) then
        local_info = max(local_info, 4)
      else if (owner < 0 .or. owner >= nrank) then
        local_info = max(local_info, 5)
      end if
    end do
    do irow = 1, size(owned_row_ids)
      owner = row_owner(owned_row_ids(irow), nrank, owner_context, owner_info)
      if (owner_info /= 0) then
        local_info = max(local_info, 6)
      else if (owner /= rank_id) then
        local_info = max(local_info, 7)
      end if
    end do
    do ivec = 1, size(yw_partial, 2)
      do irow = 1, size(yw_partial, 1)
        if (.not. ieee_is_finite(real(yw_partial(irow, ivec), 8)) .or. &
            .not. ieee_is_finite(aimag(yw_partial(irow, ivec)))) local_info = max(local_info, 8)
      end do
    end do

    if (local_info /= 0) then
      write(error_unit,'(a,i0,a,i0)') 'DG WPW owner exchange: rank-local validation failed on rank ', &
                                     rank_id, ', info=', local_info
    end if
    call MPI_Allreduce(local_info, global_info, 1, MPI_INTEGER, MPI_MAX, icomm, ierr)
    if (ierr /= MPI_SUCCESS) then
      call abort_mpi_collective(icomm, rank_id, 'partial reduction input', ierr)
      info = 9
      return
    end if
    if (global_info /= 0) then
      info = global_info
      return
    end if

    nvec = size(yw_partial, 2)
    allocate(send_counts(0:nrank-1), recv_counts(0:nrank-1), send_displs(0:nrank-1), &
             recv_displs(0:nrank-1), fill(0:nrank-1))
    send_counts = 0
    do irow = 1, size(row_ids)
      owner = row_owner(row_ids(irow), nrank, owner_context, owner_info)
      send_counts(owner) = send_counts(owner) + 1
    end do
    call MPI_Alltoall(send_counts, 1, MPI_INTEGER, recv_counts, 1, MPI_INTEGER, icomm, ierr)
    if (ierr /= MPI_SUCCESS) then
      call abort_mpi_collective(icomm, rank_id, 'partial reduction counts exchange', ierr)
      info = 10
      return
    end if
    call validate_exchange_sizes(send_counts, recv_counts, nvec, local_info)
    call collective_validation_handshake(icomm, rank_id, local_info, global_info, ierr, 'partial reduction counts')
    if (ierr /= MPI_SUCCESS) then
      call abort_mpi_collective(icomm, rank_id, 'partial reduction count validation', ierr)
      info = 14
      return
    end if
    if (global_info /= 0) then
      info = global_info
      return
    end if
    send_displs(0) = 0
    recv_displs(0) = 0
    do owner = 1, nrank-1
      send_displs(owner) = send_displs(owner-1) + send_counts(owner-1)
      recv_displs(owner) = recv_displs(owner-1) + recv_counts(owner-1)
    end do
    nsend = sum(send_counts)
    nrecv = sum(recv_counts)
    allocate(send_ids(max(1, nsend)), recv_ids(max(1, nrecv)))
    allocate(send_values(max(1, nvec), max(1, nsend)), recv_values(max(1, nvec), max(1, nrecv)))
    fill = 0
    do irow = 1, size(row_ids)
      owner = row_owner(row_ids(irow), nrank, owner_context, owner_info)
      slot = send_displs(owner) + fill(owner) + 1
      fill(owner) = fill(owner) + 1
      send_ids(slot) = row_ids(irow)
      if (nvec > 0) send_values(1:nvec, slot) = yw_partial(irow, 1:nvec)
    end do
    call MPI_Alltoallv(send_ids, send_counts, send_displs, MPI_INTEGER, recv_ids, recv_counts, &
                       recv_displs, MPI_INTEGER, icomm, ierr)
    if (ierr /= MPI_SUCCESS) then
      call abort_mpi_collective(icomm, rank_id, 'partial reduction id exchange', ierr)
      info = 11
      return
    end if

    allocate(send_value_counts(0:nrank-1), recv_value_counts(0:nrank-1), &
             send_value_displs(0:nrank-1), recv_value_displs(0:nrank-1))
    send_value_counts = send_counts * nvec
    recv_value_counts = recv_counts * nvec
    send_value_displs = send_displs * nvec
    recv_value_displs = recv_displs * nvec
    call MPI_Alltoallv(send_values, send_value_counts, send_value_displs, MPI_DOUBLE_COMPLEX, &
                       recv_values, recv_value_counts, recv_value_displs, MPI_DOUBLE_COMPLEX, icomm, ierr)
    if (ierr /= MPI_SUCCESS) then
      call abort_mpi_collective(icomm, rank_id, 'partial reduction value exchange', ierr)
      info = 12
      return
    end if

    local_info = 0
    do slot = 1, nrecv
      if (find_sorted_id(owned_row_ids, recv_ids(slot)) <= 0) then
        local_info = 13
        exit
      end if
    end do
    call collective_validation_handshake(icomm, rank_id, local_info, global_info, ierr, 'partial reduction ownership')
    if (ierr /= MPI_SUCCESS) then
      call abort_mpi_collective(icomm, rank_id, 'partial reduction ownership validation', ierr)
      info = 13
      return
    end if
    if (global_info /= 0) then
      info = global_info
      return
    end if

    do slot = 1, nrecv
      position = find_sorted_id(owned_row_ids, recv_ids(slot))
      do ivec = 1, nvec
        owned_sum(position, ivec) = owned_sum(position, ivec) + recv_values(ivec, slot)
      end do
    end do
    info = 0
  end subroutine reduce_w_partial_to_owners

  subroutine fetch_rows_from_owners(icomm, rank_id, nrank, row_owner, owner_context, owned_row_ids, owned_x, &
                                    requested_row_ids, fetched_x, info)
    integer, intent(in) :: icomm, rank_id, nrank
    procedure(row_owner_function) :: row_owner
    class(*), intent(in) :: owner_context
    integer, intent(in) :: owned_row_ids(:), requested_row_ids(:)
    complex(8), intent(in) :: owned_x(:,:)
    complex(8), intent(out) :: fetched_x(:,:)
    integer, intent(out) :: info

    integer :: local_info, global_info, ierr, irow, ivec, owner, owner_info, slot, position
    integer :: nsend, nrecv, nvec
    integer, allocatable :: send_counts(:), recv_counts(:), send_displs(:), recv_displs(:), fill(:)
    integer, allocatable :: send_ids(:), recv_ids(:), send_slots(:)
    integer, allocatable :: response_send_counts(:), response_recv_counts(:)
    integer, allocatable :: response_send_displs(:), response_recv_displs(:)
    complex(8), allocatable :: response_send(:,:), response_recv(:,:)

    fetched_x = (0.0d0, 0.0d0)
    local_info = 0
    if (nrank <= 0 .or. rank_id < 0 .or. rank_id >= nrank) local_info = 1
    if (size(owned_x, 1) /= size(owned_row_ids) .or. &
        size(fetched_x, 1) /= size(requested_row_ids) .or. &
        size(fetched_x, 2) /= size(owned_x, 2)) local_info = max(local_info, 2)
    if (.not. strictly_increasing(owned_row_ids) .or. .not. strictly_increasing(requested_row_ids)) &
      local_info = max(local_info, 3)
    do irow = 1, size(owned_row_ids)
      owner = row_owner(owned_row_ids(irow), nrank, owner_context, owner_info)
      if (owner_info /= 0) then
        local_info = max(local_info, 4)
      else if (owner /= rank_id) then
        local_info = max(local_info, 5)
      end if
    end do
    do irow = 1, size(requested_row_ids)
      owner = row_owner(requested_row_ids(irow), nrank, owner_context, owner_info)
      if (owner_info /= 0) then
        local_info = max(local_info, 6)
      else if (owner < 0 .or. owner >= nrank) then
        local_info = max(local_info, 7)
      end if
    end do
    do ivec = 1, size(owned_x, 2)
      do irow = 1, size(owned_x, 1)
        if (.not. ieee_is_finite(real(owned_x(irow, ivec), 8)) .or. &
            .not. ieee_is_finite(aimag(owned_x(irow, ivec)))) local_info = max(local_info, 8)
      end do
    end do
    call collective_validation_handshake(icomm, rank_id, local_info, global_info, ierr, 'row fetch input')
    if (ierr /= MPI_SUCCESS) then
      call abort_mpi_collective(icomm, rank_id, 'row fetch input validation', ierr)
      info = 9
      return
    end if
    if (global_info /= 0) then
      info = global_info
      return
    end if

    nvec = size(owned_x, 2)
    allocate(send_counts(0:nrank-1), recv_counts(0:nrank-1), send_displs(0:nrank-1), &
             recv_displs(0:nrank-1), fill(0:nrank-1))
    send_counts = 0
    do irow = 1, size(requested_row_ids)
      owner = row_owner(requested_row_ids(irow), nrank, owner_context, owner_info)
      send_counts(owner) = send_counts(owner) + 1
    end do
    call MPI_Alltoall(send_counts, 1, MPI_INTEGER, recv_counts, 1, MPI_INTEGER, icomm, ierr)
    if (ierr /= MPI_SUCCESS) then
      call abort_mpi_collective(icomm, rank_id, 'row fetch counts exchange', ierr)
      info = 10
      return
    end if
    call validate_exchange_sizes(send_counts, recv_counts, nvec, local_info)
    call collective_validation_handshake(icomm, rank_id, local_info, global_info, ierr, 'row fetch counts')
    if (ierr /= MPI_SUCCESS) then
      call abort_mpi_collective(icomm, rank_id, 'row fetch count validation', ierr)
      info = 15
      return
    end if
    if (global_info /= 0) then
      info = global_info
      return
    end if
    send_displs(0) = 0
    recv_displs(0) = 0
    do owner = 1, nrank-1
      send_displs(owner) = send_displs(owner-1) + send_counts(owner-1)
      recv_displs(owner) = recv_displs(owner-1) + recv_counts(owner-1)
    end do
    nsend = sum(send_counts)
    nrecv = sum(recv_counts)
    allocate(send_ids(max(1, nsend)), recv_ids(max(1, nrecv)), send_slots(max(1, nsend)))
    fill = 0
    do irow = 1, size(requested_row_ids)
      owner = row_owner(requested_row_ids(irow), nrank, owner_context, owner_info)
      slot = send_displs(owner) + fill(owner) + 1
      fill(owner) = fill(owner) + 1
      send_ids(slot) = requested_row_ids(irow)
      send_slots(slot) = irow
    end do
    call MPI_Alltoallv(send_ids, send_counts, send_displs, MPI_INTEGER, recv_ids, recv_counts, &
                       recv_displs, MPI_INTEGER, icomm, ierr)
    if (ierr /= MPI_SUCCESS) then
      call abort_mpi_collective(icomm, rank_id, 'row fetch id exchange', ierr)
      info = 11
      return
    end if

    local_info = 0
    do slot = 1, nrecv
      if (find_sorted_id(owned_row_ids, recv_ids(slot)) <= 0) then
        local_info = 12
        exit
      end if
    end do
    call collective_validation_handshake(icomm, rank_id, local_info, global_info, ierr, 'row fetch ownership')
    if (ierr /= MPI_SUCCESS) then
      call abort_mpi_collective(icomm, rank_id, 'row fetch ownership validation', ierr)
      info = 13
      return
    end if
    if (global_info /= 0) then
      info = global_info
      return
    end if

    allocate(response_send(max(1, nvec), max(1, nrecv)), &
             response_recv(max(1, nvec), max(1, nsend)))
    do slot = 1, nrecv
      position = find_sorted_id(owned_row_ids, recv_ids(slot))
      if (nvec > 0) response_send(1:nvec, slot) = owned_x(position, 1:nvec)
    end do
    allocate(response_send_counts(0:nrank-1), response_recv_counts(0:nrank-1), &
             response_send_displs(0:nrank-1), response_recv_displs(0:nrank-1))
    response_send_counts = recv_counts * nvec
    response_recv_counts = send_counts * nvec
    response_send_displs = recv_displs * nvec
    response_recv_displs = send_displs * nvec
    call MPI_Alltoallv(response_send, response_send_counts, response_send_displs, MPI_DOUBLE_COMPLEX, &
                       response_recv, response_recv_counts, response_recv_displs, MPI_DOUBLE_COMPLEX, icomm, ierr)
    if (ierr /= MPI_SUCCESS) then
      call abort_mpi_collective(icomm, rank_id, 'row fetch value exchange', ierr)
      info = 14
      return
    end if
    do slot = 1, nsend
      irow = send_slots(slot)
      if (nvec > 0) fetched_x(irow, 1:nvec) = response_recv(1:nvec, slot)
    end do
    info = 0
  end subroutine fetch_rows_from_owners

  subroutine abort_mpi_collective(icomm, rank_id, label, mpi_error)
    integer, intent(in) :: icomm, rank_id, mpi_error
    character(*), intent(in) :: label
    integer :: abort_ierr
    write(error_unit,'(a,a,a,i0,a,i0)') 'DG WPW MPI collective failed: ', trim(label), &
      ' on rank ', rank_id, ', ierr=', mpi_error
    call MPI_Abort(icomm, mpi_error, abort_ierr)
  end subroutine abort_mpi_collective

  subroutine collective_validation_handshake(icomm, rank_id, local_info, global_info, ierr, label)
    integer, intent(in) :: icomm, rank_id, local_info
    integer, intent(out) :: global_info, ierr
    character(*), intent(in) :: label
    if (local_info /= 0) then
      write(error_unit,'(a,a,a,i0,a,i0)') 'DG WPW ', trim(label), &
        ': rank-local validation failed on rank ', rank_id, ', info=', local_info
    end if
    call MPI_Allreduce(local_info, global_info, 1, MPI_INTEGER, MPI_MAX, icomm, ierr)
  end subroutine collective_validation_handshake

  subroutine validate_exchange_sizes(send_counts, recv_counts, nvec, info)
    integer, intent(in) :: send_counts(0:), recv_counts(0:), nvec
    integer, intent(out) :: info
    integer(int64) :: send_total, recv_total, int_limit

    info = 0
    if (nvec < 0 .or. any(send_counts < 0) .or. any(recv_counts < 0)) then
      info = 14
      return
    end if
    send_total = sum(int(send_counts, int64))
    recv_total = sum(int(recv_counts, int64))
    int_limit = int(huge(0), int64)
    if (send_total > int_limit .or. recv_total > int_limit) then
      info = 14
      return
    end if
    if (int(nvec, int64) > 0_int64) then
      if (send_total > int_limit / int(nvec, int64) .or. &
          recv_total > int_limit / int(nvec, int64)) info = 14
    end if
  end subroutine validate_exchange_sizes

  logical function strictly_increasing(ids) result(ok)
    integer, intent(in) :: ids(:)
    integer :: i
    ok = .true.
    do i = 2, size(ids)
      if (ids(i) <= ids(i-1)) then
        ok = .false.
        return
      end if
    end do
  end function strictly_increasing

  integer function find_sorted_id(ids, target) result(position)
    integer, intent(in) :: ids(:), target
    integer :: left, middle, right
    position = 0
    left = 1
    right = size(ids)
    do while (left <= right)
      middle = left + (right - left) / 2
      if (ids(middle) == target) then
        position = middle
        return
      else if (ids(middle) < target) then
        left = middle + 1
      else
        right = middle - 1
      end if
    end do
  end function find_sorted_id

end module rt_dg_wpw_owner_exchange
