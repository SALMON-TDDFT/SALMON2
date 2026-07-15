program test_dg_wpw_sparse_builder
  use rt_dg_wpw_column_layout, only: s_dg_wpw_column_layout, initialize_wpw_column_layout
  use rt_dg_wpw_sparse_blocks, only: s_dg_wpw_sparse_blocks
  use rt_dg_wpw_sparse_builder, only: s_dg_wpw_sparse_candidates, &
    build_windowed_sparse_wpw_operators, wpw_normalized_window_at_grid, wpw_grad_chi, &
    wpw_candidate_volume, wpw_candidate_face
  implicit none
  type(s_dg_wpw_column_layout) :: layout
  type(s_dg_wpw_sparse_candidates) :: candidates
  type(s_dg_wpw_sparse_blocks) :: h_blocks, s_blocks
  integer :: info, rank_id
  complex(8) :: value, grad_p(3), phase
  real(8), parameter :: tol = 1.0d-13

  allocate(candidates%wp_w_row_ids(2), candidates%wp_pw_col_ids(2), candidates%wp_origin(2))
  allocate(candidates%wp_h_values(2,1), candidates%wp_s_values(2,1))
  allocate(candidates%pp_pw_row_ids(2), candidates%pp_pw_col_ids(2), candidates%pp_origin(2))
  allocate(candidates%pp_h_values(2,1), candidates%pp_s_values(2,1))

  do rank_id = 0, 1
    call initialize_wpw_column_layout(layout, 1, 2, rank_id, 2, info)
    if (info /= 0) error stop 1
    candidates%wp_w_row_ids = [1,2]
    candidates%wp_pw_col_ids = rank_id + 1
    candidates%wp_origin = [wpw_candidate_volume, wpw_candidate_face]
    candidates%wp_h_values(:,1) = [1d0,2d0] + 2d0 * rank_id
    candidates%wp_s_values(:,1) = [0.1d0,0.2d0] + 0.2d0 * rank_id
    candidates%pp_pw_row_ids = rank_id + 1
    candidates%pp_pw_col_ids = [1,2]
    candidates%pp_origin = wpw_candidate_volume
    candidates%pp_h_values(:,1) = [5d0,6d0] + 2d0 * rank_id
    candidates%pp_s_values(:,1) = [1d0,0.1d0]
    call build_windowed_sparse_wpw_operators(layout, rank_id, 2, candidates, h_blocks, s_blocks, info)
    if (info /= 0) error stop 2
    if (size(h_blocks%wp_values,1) /= 2 .or. size(h_blocks%pp_values,1) /= 2) error stop 3
    if (any(h_blocks%wp_pw_col_ids /= rank_id + 1)) error stop 4
    if (any(h_blocks%pp_pw_row_ids /= rank_id + 1)) error stop 5
    if (maxval(abs(h_blocks%wp_values - candidates%wp_h_values)) > tol) error stop 6
    candidates%pp_origin(1) = wpw_candidate_face
    call build_windowed_sparse_wpw_operators(layout, rank_id, 2, candidates, h_blocks, s_blocks, info)
    if (info /= 13) error stop 10
    candidates%pp_origin = wpw_candidate_volume
    candidates%wp_pw_col_ids = 2 - rank_id
    call build_windowed_sparse_wpw_operators(layout, rank_id, 2, candidates, h_blocks, s_blocks, info)
    if (info /= 12) error stop 11
  end do

  phase = cmplx(0.6d0,0.8d0,8)
  value = wpw_normalized_window_at_grid(2d0, phase, 4d0, info)
  if (info /= 0 .or. abs(value-phase) > tol) error stop 7
  call wpw_grad_chi(2d0, [1d0,2d0,3d0], [0.5d0,0d0,-0.5d0], phase, 4d0, grad_p, info)
  if (info /= 0) error stop 8
  if (abs(grad_p(1) - cmplx(0.5d0,0.5d0,8)*phase) > tol) error stop 9
  write(*,'(a)') 'PASS direct windowed sparse builder fixture'
end program test_dg_wpw_sparse_builder
