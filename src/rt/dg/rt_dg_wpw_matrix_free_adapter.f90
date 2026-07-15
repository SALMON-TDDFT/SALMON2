module rt_dg_wpw_matrix_free_adapter
  use mpi, only: MPI_Allreduce, MPI_INTEGER, MPI_MAX, MPI_SUCCESS
  use rt_dg_fragment_types, only: s_dg_fragment_rt, matrix_block_info, complex_matrix_block_info
  use rt_dg_fragment_ops, only: apply_matrix_blocks_batch_compact, apply_complex_matrix_blocks_batch_compact
  use rt_dg_wpw_sparse_blocks, only: s_dg_wpw_sparse_blocks, apply_wp_owned_columns, apply_pp_owned_rows
  use rt_dg_wpw_owner_exchange, only: fetch_rows_from_owners, reduce_w_partial_to_owners, row_owner_function
  implicit none
  private

  integer, parameter, public :: ww_real_local_sipg = 1
  integer, parameter, public :: ww_complex_nonlocal = 2

  type, public :: s_rt_dg_wpw_operator_context
    type(s_dg_fragment_rt), pointer :: dg_frag => null()
    integer :: ispin = 0
    procedure(row_owner_function), pointer, nopass :: w_row_owner => null()
    procedure(row_owner_function), pointer, nopass :: pw_row_owner => null()
    class(*), pointer :: w_owner_context => null()
    class(*), pointer :: pw_owner_context => null()
    integer, pointer :: owned_w_row_ids(:) => null(), owned_pw_row_ids(:) => null()
    integer, pointer :: w_input_row_ids(:) => null(), pw_input_row_ids(:) => null()
    integer, pointer :: ww_block_ids(:) => null(), ww_nl_block_ids(:) => null()
    type(matrix_block_info), pointer :: ww_blocks(:) => null()
    type(complex_matrix_block_info), pointer :: ww_nl_blocks(:) => null()
    type(s_dg_wpw_sparse_blocks), pointer :: mixed_blocks => null()
    logical :: ww_identity = .false.
    integer :: ww_real_component = 0, ww_complex_component = 0
    logical :: context_bound = .false.
  end type s_rt_dg_wpw_operator_context

  public :: apply_h_wpw_callback
  public :: apply_s_wpw_callback
  public :: bind_rt_dg_wpw_operator_context

contains

  subroutine bind_rt_dg_wpw_operator_context(ctx, real_component, complex_component, info)
    type(s_rt_dg_wpw_operator_context), intent(inout) :: ctx
    integer, intent(in) :: real_component, complex_component
    integer, intent(out) :: info
    ! All pointer targets associated by the caller must have TARGET storage and
    ! outlive every callback using ctx.  Binding records that this lifetime and
    ! component-provenance contract has been checked by the construction site.
    ctx%context_bound = .false.
    ctx%ww_real_component = real_component
    ctx%ww_complex_component = complex_component
    info = 0
    if (real_component /= ww_real_local_sipg .or. &
        complex_component /= ww_complex_nonlocal) info = 1
    if (.not. context_ready(ctx, .true., .false.)) info = max(info, 2)
    if (info == 0) ctx%context_bound = .true.
  end subroutine bind_rt_dg_wpw_operator_context

  subroutine apply_h_wpw_callback(context, xw_owned, xp_owned, yw_owned, yp_owned, info)
    class(*), intent(inout) :: context
    complex(8), intent(in) :: xw_owned(:,:), xp_owned(:,:)
    complex(8), intent(out) :: yw_owned(:,:), yp_owned(:,:)
    integer, intent(out) :: info

    select type (ctx => context)
    type is (s_rt_dg_wpw_operator_context)
      if (.not. ctx%context_bound) then
        info = 79
        return
      end if
      if (.not. context_ready(ctx, .true., .true.)) then
        info = 80
        return
      end if
      call apply_h_wpw_distributed(ctx%dg_frag, ctx%ispin, ctx%w_row_owner, ctx%w_owner_context, &
        ctx%pw_row_owner, ctx%pw_owner_context, ctx%owned_w_row_ids, ctx%owned_pw_row_ids, &
        ctx%w_input_row_ids, ctx%pw_input_row_ids, xw_owned, xp_owned, ctx%ww_blocks, &
        ctx%ww_block_ids, ctx%ww_nl_blocks, ctx%ww_nl_block_ids, ctx%mixed_blocks, yw_owned, yp_owned, info)
    class default
      info = 81
    end select
  end subroutine apply_h_wpw_callback

  subroutine apply_s_wpw_callback(context, xw_owned, xp_owned, yw_owned, yp_owned, info)
    class(*), intent(inout) :: context
    complex(8), intent(in) :: xw_owned(:,:), xp_owned(:,:)
    complex(8), intent(out) :: yw_owned(:,:), yp_owned(:,:)
    integer, intent(out) :: info

    select type (ctx => context)
    type is (s_rt_dg_wpw_operator_context)
      if (.not. ctx%context_bound) then
        info = 78
        return
      end if
      if (.not. context_ready(ctx, .false., .true.)) then
        info = 82
        return
      end if
      call apply_s_wpw_distributed(ctx%dg_frag, ctx%ispin, ctx%w_row_owner, ctx%w_owner_context, &
        ctx%pw_row_owner, ctx%pw_owner_context, ctx%owned_w_row_ids, ctx%owned_pw_row_ids, &
        ctx%w_input_row_ids, ctx%pw_input_row_ids, xw_owned, xp_owned, ctx%ww_blocks, &
        ctx%ww_block_ids, ctx%mixed_blocks, ctx%ww_identity, yw_owned, yp_owned, info)
    class default
      info = 83
    end select
  end subroutine apply_s_wpw_callback

  logical function context_ready(ctx, require_nonlocal, require_bound) result(ready)
    type(s_rt_dg_wpw_operator_context), intent(in) :: ctx
    logical, intent(in) :: require_nonlocal
    logical, intent(in) :: require_bound
    ready = associated(ctx%dg_frag) .and. associated(ctx%w_row_owner) .and. &
      associated(ctx%pw_row_owner) .and. associated(ctx%w_owner_context) .and. &
      associated(ctx%pw_owner_context) .and. associated(ctx%owned_w_row_ids) .and. &
      associated(ctx%owned_pw_row_ids) .and. associated(ctx%w_input_row_ids) .and. &
      associated(ctx%pw_input_row_ids) .and. associated(ctx%ww_blocks) .and. &
      associated(ctx%ww_block_ids) .and. associated(ctx%mixed_blocks)
    if (require_nonlocal) ready = ready .and. associated(ctx%ww_nl_blocks) .and. associated(ctx%ww_nl_block_ids)
    if (require_bound) ready = ready .and. ctx%context_bound .and. &
      ctx%ww_real_component == ww_real_local_sipg .and. ctx%ww_complex_component == ww_complex_nonlocal
  end function context_ready

  subroutine apply_h_wpw_distributed(dg_frag, ispin, w_row_owner, w_owner_context, pw_row_owner, pw_owner_context, &
                                     owned_w_row_ids, owned_pw_row_ids, w_input_row_ids, pw_input_row_ids, &
                                     xw_owned, xp_owned, ww_blocks, ww_block_ids, ww_nl_blocks, ww_nl_block_ids, mixed_blocks, &
                                     yw_owned, yp_owned, info)
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ispin
    procedure(row_owner_function) :: w_row_owner, pw_row_owner
    class(*), intent(in) :: w_owner_context, pw_owner_context
    integer, intent(in) :: owned_w_row_ids(:), owned_pw_row_ids(:)
    integer, intent(in) :: w_input_row_ids(:), pw_input_row_ids(:), ww_block_ids(:)
    complex(8), intent(in) :: xw_owned(:,:), xp_owned(:,:)
    type(matrix_block_info), intent(in) :: ww_blocks(:)
    type(complex_matrix_block_info), intent(in) :: ww_nl_blocks(:)
    integer, intent(in) :: ww_nl_block_ids(:)
    type(s_dg_wpw_sparse_blocks), intent(in) :: mixed_blocks
    complex(8), intent(out) :: yw_owned(:,:), yp_owned(:,:)
    integer, intent(out) :: info

    call apply_wpw_distributed(dg_frag, ispin, w_row_owner, w_owner_context, pw_row_owner, pw_owner_context, owned_w_row_ids, &
      owned_pw_row_ids, w_input_row_ids, pw_input_row_ids, xw_owned, xp_owned, ww_blocks, &
      ww_block_ids, mixed_blocks, .false., yw_owned, yp_owned, info, ww_nl_blocks, ww_nl_block_ids)
  end subroutine apply_h_wpw_distributed

  subroutine apply_s_wpw_distributed(dg_frag, ispin, w_row_owner, w_owner_context, pw_row_owner, pw_owner_context, &
                                     owned_w_row_ids, owned_pw_row_ids, w_input_row_ids, pw_input_row_ids, &
                                     xw_owned, xp_owned, ww_blocks, ww_block_ids, mixed_blocks, ww_identity, &
                                     yw_owned, yp_owned, info)
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ispin
    procedure(row_owner_function) :: w_row_owner, pw_row_owner
    class(*), intent(in) :: w_owner_context, pw_owner_context
    integer, intent(in) :: owned_w_row_ids(:), owned_pw_row_ids(:)
    integer, intent(in) :: w_input_row_ids(:), pw_input_row_ids(:), ww_block_ids(:)
    complex(8), intent(in) :: xw_owned(:,:), xp_owned(:,:)
    type(matrix_block_info), intent(in) :: ww_blocks(:)
    type(s_dg_wpw_sparse_blocks), intent(in) :: mixed_blocks
    logical, intent(in) :: ww_identity
    complex(8), intent(out) :: yw_owned(:,:), yp_owned(:,:)
    integer, intent(out) :: info

    call apply_wpw_distributed(dg_frag, ispin, w_row_owner, w_owner_context, pw_row_owner, pw_owner_context, owned_w_row_ids, &
      owned_pw_row_ids, w_input_row_ids, pw_input_row_ids, xw_owned, xp_owned, ww_blocks, &
      ww_block_ids, mixed_blocks, ww_identity, yw_owned, yp_owned, info)
  end subroutine apply_s_wpw_distributed

  subroutine apply_wpw_distributed(dg_frag, ispin, w_row_owner, w_owner_context, pw_row_owner, pw_owner_context, &
                                   owned_w_row_ids, owned_pw_row_ids, w_input_row_ids, pw_input_row_ids, &
                                   xw_owned, xp_owned, ww_blocks, ww_block_ids, mixed_blocks, ww_identity, &
                                   yw_owned, yp_owned, info, ww_nl_blocks, ww_nl_block_ids)
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ispin
    procedure(row_owner_function) :: w_row_owner, pw_row_owner
    class(*), intent(in) :: w_owner_context, pw_owner_context
    integer, intent(in) :: owned_w_row_ids(:), owned_pw_row_ids(:)
    integer, intent(in) :: w_input_row_ids(:), pw_input_row_ids(:), ww_block_ids(:)
    complex(8), intent(in) :: xw_owned(:,:), xp_owned(:,:)
    type(matrix_block_info), intent(in) :: ww_blocks(:)
    type(s_dg_wpw_sparse_blocks), intent(in) :: mixed_blocks
    logical, intent(in) :: ww_identity
    complex(8), intent(out) :: yw_owned(:,:), yp_owned(:,:)
    integer, intent(out) :: info
    type(complex_matrix_block_info), intent(in), optional :: ww_nl_blocks(:)
    integer, intent(in), optional :: ww_nl_block_ids(:)

    complex(8), allocatable :: xw_input(:,:), xp_input(:,:), yw_input(:,:), yw_wp_partial(:,:)
    complex(8), allocatable :: yw_wp_owned(:,:)
    integer :: local_info, global_info, ierr, irow, position, nvec

    ! The adapter stores only rank-local owned rows and support-neighbor halo
    ! rows.  Rank-local validation is completed before each collective phase.
    yw_owned = (0.0d0, 0.0d0)
    yp_owned = (0.0d0, 0.0d0)
    local_info = 0
    if (size(xw_owned, 1) /= size(owned_w_row_ids) .or. &
        size(xp_owned, 1) /= size(owned_pw_row_ids) .or. &
        size(xw_owned, 2) /= size(xp_owned, 2) .or. &
        size(yw_owned, 1) /= size(owned_w_row_ids) .or. &
        size(yp_owned, 1) /= size(owned_pw_row_ids) .or. &
        size(yw_owned, 2) /= size(xw_owned, 2) .or. size(yp_owned, 2) /= size(xp_owned, 2)) local_info = 1
    call collective_info(dg_frag, local_info, global_info, ierr)
    if (ierr /= MPI_SUCCESS) then
      info = 90
      return
    end if
    if (global_info /= 0) then
      info = global_info
      return
    end if

    nvec = size(xw_owned, 2)
    allocate(xw_input(size(w_input_row_ids), nvec), xp_input(size(pw_input_row_ids), nvec))
    call fetch_rows_from_owners(dg_frag%icomm, dg_frag%id, dg_frag%isize, w_row_owner, w_owner_context, &
      owned_w_row_ids, xw_owned, w_input_row_ids, xw_input, info)
    if (info /= 0) return
    call fetch_rows_from_owners(dg_frag%icomm, dg_frag%id, dg_frag%isize, pw_row_owner, pw_owner_context, &
      owned_pw_row_ids, xp_owned, pw_input_row_ids, xp_input, info)
    if (info /= 0) return

    allocate(yw_input(size(w_input_row_ids), nvec), yw_wp_partial(size(w_input_row_ids), nvec), &
             yw_wp_owned(size(owned_w_row_ids), nvec))
    yw_input = (0.0d0, 0.0d0)
    yw_wp_partial = (0.0d0, 0.0d0)
    if (.not. ww_identity) then
      call apply_matrix_blocks_batch_compact(dg_frag, ww_blocks, ispin, w_input_row_ids, &
        xw_input, yw_input, ww_block_ids, local_info)
      call collective_info(dg_frag, local_info, global_info, ierr)
      if (ierr /= MPI_SUCCESS) then
        info = 91
        return
      end if
      if (global_info /= 0) then
        info = global_info
        return
      end if
      if (present(ww_nl_blocks) .neqv. present(ww_nl_block_ids)) then
        local_info = 21
      else if (present(ww_nl_blocks)) then
        call apply_complex_matrix_blocks_batch_compact(dg_frag, ww_nl_blocks, ispin, w_input_row_ids, &
          xw_input, yw_input, ww_nl_block_ids, local_info)
      end if
      call collective_info(dg_frag, local_info, global_info, ierr)
      if (ierr /= MPI_SUCCESS) then
        info = 95
        return
      end if
      if (global_info /= 0) then
        info = global_info
        return
      end if
      local_info = 0
      do irow = 1, size(owned_w_row_ids)
        position = find_sorted_id(w_input_row_ids, owned_w_row_ids(irow))
        if (position <= 0) then
          local_info = 20
          exit
        end if
        yw_owned(irow,:) = yw_input(position,:)
      end do
      call collective_info(dg_frag, local_info, global_info, ierr)
      if (ierr /= MPI_SUCCESS) then
        info = 94
        return
      end if
      if (global_info /= 0) then
        info = global_info
        return
      end if
    else
      yw_owned = xw_owned
    end if

    call apply_wp_owned_columns(mixed_blocks, ispin, w_input_row_ids, owned_pw_row_ids, &
      xw_input, xp_owned, yw_wp_partial, yp_owned, local_info)
    call collective_info(dg_frag, local_info, global_info, ierr)
    if (ierr /= MPI_SUCCESS) then
      info = 92
      return
    end if
    if (global_info /= 0) then
      info = global_info
      return
    end if
    call apply_pp_owned_rows(mixed_blocks, ispin, pw_input_row_ids, owned_pw_row_ids, &
      xp_input, yp_owned, local_info)
    call collective_info(dg_frag, local_info, global_info, ierr)
    if (ierr /= MPI_SUCCESS) then
      info = 93
      return
    end if
    if (global_info /= 0) then
      info = global_info
      return
    end if
    call reduce_w_partial_to_owners(dg_frag%icomm, dg_frag%id, dg_frag%isize, w_row_owner, w_owner_context, &
      w_input_row_ids, yw_wp_partial, owned_w_row_ids, yw_wp_owned, info)
    if (info /= 0) return
    yw_owned = yw_owned + yw_wp_owned
    info = 0
  end subroutine apply_wpw_distributed

  subroutine collective_info(dg_frag, local_info, global_info, ierr)
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: local_info
    integer, intent(out) :: global_info, ierr
    call MPI_Allreduce(local_info, global_info, 1, MPI_INTEGER, MPI_MAX, dg_frag%icomm, ierr)
  end subroutine collective_info

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

end module rt_dg_wpw_matrix_free_adapter
