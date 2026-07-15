program test_dg_wpw_adapter_mpi
  use mpi
  use rt_dg_fragment_types, only: s_dg_fragment_rt, matrix_block_info, complex_matrix_block_info
  use rt_dg_wpw_sparse_blocks, only: s_dg_wpw_sparse_blocks
  use rt_dg_wpw_matrix_free_adapter, only: s_rt_dg_wpw_operator_context, &
    apply_h_wpw_callback, apply_s_wpw_callback, bind_rt_dg_wpw_operator_context, &
    ww_real_local_sipg, ww_complex_nonlocal
  implicit none

  type(s_dg_fragment_rt), target :: dg
  type(matrix_block_info), target :: ww_blocks(4)
  type(complex_matrix_block_info), target :: ww_nl_blocks(4)
  type(s_dg_wpw_sparse_blocks), target :: mixed_h, mixed_s
  type(s_rt_dg_wpw_operator_context) :: ctx_h, ctx_s
  integer, target :: owner_extent = 2
  integer, target :: owned_w(1), owned_p(1), input_w(2) = [1,2], input_p(2) = [1,2]
  integer, target :: block_ids(4) = [1,2,3,4], nl_block_ids(4) = [1,2,3,4]
  integer :: ierr, rank_id, nrank, info, iblk, irow, icol
  complex(8) :: xw(1,1), xp(1,1), yh_w(1,1), yh_p(1,1), ys_w(1,1), ys_p(1,1)
  complex(8) :: xall(4), expected_h(4), expected_s(4)
  complex(8) :: dense_h(4,4), dense_s(4,4), ww(2,2), nl(2,2), wp(2,2), pp(2,2)
  complex(8) :: swp(2,2), spp(2,2)
  real(8) :: local_error, global_error

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank_id, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nrank, ierr)
  if (nrank /= 2) call MPI_Abort(MPI_COMM_WORLD, 1, ierr)

  dg%icomm = MPI_COMM_WORLD
  dg%id = rank_id
  dg%isize = nrank
  dg%nspin = 1
  dg%n_frag = 2
  allocate(dg%n_basis(2,1), dg%index_basis(1,2,1))
  dg%n_basis(:,1) = 1
  dg%index_basis(1,:,1) = [1,2]

  ww = reshape([cmplx(2d0,0d0,8), cmplx(0.5d0,0d0,8), &
                cmplx(0.5d0,0d0,8), cmplx(3d0,0d0,8)], [2,2])
  nl = reshape([cmplx(0d0,0d0,8), cmplx(0d0,-0.2d0,8), &
                cmplx(0d0,0.2d0,8), cmplx(0.1d0,0d0,8)], [2,2])
  wp = reshape([cmplx(1d0,0.1d0,8), cmplx(0d0,0.3d0,8), &
                cmplx(0.2d0,0d0,8), cmplx(0.4d0,-0.1d0,8)], [2,2])
  pp = reshape([cmplx(4d0,0d0,8), cmplx(0d0,-0.1d0,8), &
                cmplx(0d0,0.1d0,8), cmplx(5d0,0d0,8)], [2,2])
  swp = 0.05d0 * wp
  spp = reshape([cmplx(1.1d0,0d0,8), cmplx(0d0,-0.01d0,8), &
                 cmplx(0d0,0.01d0,8), cmplx(1.2d0,0d0,8)], [2,2])

  iblk = 0
  do icol = 1, 2
    do irow = 1, 2
      iblk = iblk + 1
      ww_blocks(iblk)%ifrag_row = irow
      ww_blocks(iblk)%ifrag_col = icol
      ww_nl_blocks(iblk)%ifrag_row = irow
      ww_nl_blocks(iblk)%ifrag_col = icol
      allocate(ww_blocks(iblk)%val(1,1,1), ww_nl_blocks(iblk)%val(1,1,1))
      ww_blocks(iblk)%val(1,1,1) = real(ww(irow,icol),8)
      ww_nl_blocks(iblk)%val(1,1,1) = nl(irow,icol)
    end do
  end do

  owned_w = rank_id + 1
  owned_p = rank_id + 1
  call init_mixed(mixed_h, owned_p(1), wp, pp)
  call init_mixed(mixed_s, owned_p(1), swp, spp)
  call bind_context(ctx_h, mixed_h, .false.)
  call bind_context(ctx_s, mixed_s, .true.)
  call bind_rt_dg_wpw_operator_context(ctx_h, ww_real_local_sipg, ww_complex_nonlocal, info)
  if (info /= 0) call MPI_Abort(MPI_COMM_WORLD, 5 + info, ierr)
  call bind_rt_dg_wpw_operator_context(ctx_s, ww_real_local_sipg, ww_complex_nonlocal, info)
  if (info /= 0) call MPI_Abort(MPI_COMM_WORLD, 7 + info, ierr)

  xall = [cmplx(1d0,0.2d0,8), cmplx(-0.3d0,0.4d0,8), &
          cmplx(0.7d0,-0.1d0,8), cmplx(-0.2d0,0.5d0,8)]
  xw(1,1) = xall(owned_w(1))
  xp(1,1) = xall(2 + owned_p(1))
  call apply_h_wpw_callback(ctx_h, xw, xp, yh_w, yh_p, info)
  if (info /= 0) call MPI_Abort(MPI_COMM_WORLD, 10 + info, ierr)
  call apply_s_wpw_callback(ctx_s, xw, xp, ys_w, ys_p, info)
  if (info /= 0) call MPI_Abort(MPI_COMM_WORLD, 20 + info, ierr)

  dense_h = (0d0,0d0)
  dense_h(1:2,1:2) = ww + nl
  dense_h(1:2,3:4) = wp
  dense_h(3:4,1:2) = transpose(conjg(wp))
  dense_h(3:4,3:4) = pp
  dense_s = (0d0,0d0)
  dense_s(1,1) = (1d0,0d0)
  dense_s(2,2) = (1d0,0d0)
  dense_s(1:2,3:4) = swp
  dense_s(3:4,1:2) = transpose(conjg(swp))
  dense_s(3:4,3:4) = spp
  expected_h = matmul(dense_h, xall)
  expected_s = matmul(dense_s, xall)
  local_error = max(abs(yh_w(1,1)-expected_h(owned_w(1))), abs(yh_p(1,1)-expected_h(2+owned_p(1))), &
                    abs(ys_w(1,1)-expected_s(owned_w(1))), abs(ys_p(1,1)-expected_s(2+owned_p(1))))
  call MPI_Allreduce(local_error, global_error, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
  if (global_error > 1d-12) call MPI_Abort(MPI_COMM_WORLD, 90, ierr)
  if (rank_id == 0) write(*,'(a,es12.4)') 'PASS end-to-end WPW H/S callback MPI equivalence, error=', global_error
  call MPI_Finalize(ierr)

contains

  subroutine init_mixed(mixed, owned_col, wp_dense, pp_dense)
    type(s_dg_wpw_sparse_blocks), intent(out) :: mixed
    integer, intent(in) :: owned_col
    complex(8), intent(in) :: wp_dense(2,2), pp_dense(2,2)
    allocate(mixed%wp_w_row_ids(2), mixed%wp_pw_col_ids(2), mixed%wp_values(2,1))
    allocate(mixed%pp_pw_row_ids(2), mixed%pp_pw_col_ids(2), mixed%pp_values(2,1))
    mixed%wp_w_row_ids = [1,2]
    mixed%wp_pw_col_ids = owned_col
    mixed%wp_values(:,1) = wp_dense(:,owned_col)
    mixed%pp_pw_row_ids = owned_col
    mixed%pp_pw_col_ids = [1,2]
    mixed%pp_values(:,1) = pp_dense(owned_col,:)
  end subroutine init_mixed

  subroutine bind_context(ctx, mixed, identity_ww)
    type(s_rt_dg_wpw_operator_context), intent(out) :: ctx
    type(s_dg_wpw_sparse_blocks), target, intent(inout) :: mixed
    logical, intent(in) :: identity_ww
    ctx%dg_frag => dg
    ctx%ispin = 1
    ctx%w_row_owner => cyclic_owner
    ctx%pw_row_owner => cyclic_owner
    ctx%w_owner_context => owner_extent
    ctx%pw_owner_context => owner_extent
    ctx%owned_w_row_ids => owned_w
    ctx%owned_pw_row_ids => owned_p
    ctx%w_input_row_ids => input_w
    ctx%pw_input_row_ids => input_p
    ctx%ww_blocks => ww_blocks
    ctx%ww_block_ids => block_ids
    ctx%ww_nl_blocks => ww_nl_blocks
    ctx%ww_nl_block_ids => nl_block_ids
    ctx%mixed_blocks => mixed
    ctx%ww_identity = identity_ww
  end subroutine bind_context

  integer function cyclic_owner(row_id, nrank_arg, context, owner_info) result(owner)
    integer, intent(in) :: row_id, nrank_arg
    class(*), intent(in) :: context
    integer, intent(out) :: owner_info
    owner = -1
    owner_info = 0
    select type (nrow => context)
    type is (integer)
      if (row_id < 1 .or. row_id > nrow) then
        owner_info = 1
      else
        owner = modulo(row_id - 1, nrank_arg)
      end if
    class default
      owner_info = 2
    end select
  end function cyclic_owner
end program test_dg_wpw_adapter_mpi
