program test_dg_wpw_column_layout
  use rt_dg_wpw_column_layout, only: s_dg_wpw_column_layout, initialize_wpw_column_layout, &
                                     wpw_column_id, wpw_column_pair, wpw_column_owner
  implicit none

  integer, parameter :: nfrag = 5, ng = 7, nrank = 4
  type(s_dg_wpw_column_layout) :: layout
  integer :: rank_id, info, i, column_id, fragment_id, g_id, owner
  integer :: counts(0:nrank-1), seen(nfrag*ng)

  counts = 0
  seen = 0
  do rank_id = 0, nrank-1
    call initialize_wpw_column_layout(layout, nfrag, ng, rank_id, nrank, info)
    if (info /= 0) error stop 10
    if (trim(layout%basis_kind) /= 'windowed_kg') error stop 11
    if (layout%n_global_columns /= nfrag * ng .or. layout%n_g_modes /= ng) error stop 12
    counts(rank_id) = size(layout%owned_column_ids)
    if (size(layout%pw_fragment_ids) /= counts(rank_id) .or. &
        size(layout%pw_g_ids) /= counts(rank_id) .or. size(layout%pw_owner) /= counts(rank_id)) error stop 13
    do i = 1, counts(rank_id)
      column_id = layout%owned_column_ids(i)
      if (column_id < 1 .or. column_id > nfrag * ng) error stop 14
      seen(column_id) = seen(column_id) + 1
      call wpw_column_pair(column_id, ng, fragment_id, g_id, info)
      if (info /= 0) error stop 15
      if (fragment_id /= layout%pw_fragment_ids(i) .or. g_id /= layout%pw_g_ids(i)) error stop 16
      if (wpw_column_id(fragment_id, g_id, ng, info) /= column_id .or. info /= 0) error stop 17
      owner = wpw_column_owner(column_id, nfrag * ng, nrank, info)
      if (info /= 0 .or. owner /= rank_id .or. layout%pw_owner(i) /= rank_id) error stop 18
    end do
  end do
  if (any(seen /= 1)) error stop 19
  if (maxval(counts) - minval(counts) > 1) error stop 20

  call initialize_wpw_column_layout(layout, 0, ng, 0, nrank, info)
  if (info == 0) error stop 21
  call initialize_wpw_column_layout(layout, huge(1), 2, 0, nrank, info)
  if (info == 0) error stop 22
  write(*,'(a)') 'PASS production (K,G) column-layout fixture'
end program test_dg_wpw_column_layout
