!=======================================================================
  ! Distribute fragments across MPI ranks
  !=======================================================================
  subroutine distribute_fragments(n_frag, myrank, nprocs, ifrag_start, ifrag_end)
    implicit none
    integer, intent(in)  :: n_frag, myrank, nprocs
    integer, intent(out) :: ifrag_start, ifrag_end
    call get_fragment_range_for_rank(n_frag, myrank, nprocs, ifrag_start, ifrag_end)
    
  end subroutine distribute_fragments

  !=======================================================================
  ! Return [ifrag_start, ifrag_end] assigned to one rank by block distribution
  !=======================================================================
  subroutine get_fragment_range_for_rank(n_frag, myrank, nprocs, ifrag_start, ifrag_end)
    implicit none
    integer, intent(in)  :: n_frag, myrank, nprocs
    integer, intent(out) :: ifrag_start, ifrag_end
    integer :: n_per_proc, n_remainder

    ifrag_start = 1
    ifrag_end = 0
    if (n_frag <= 0 .or. nprocs <= 0) return
    if (myrank < 0 .or. myrank >= nprocs) return

    n_per_proc = n_frag / nprocs
    n_remainder = mod(n_frag, nprocs)

    if (myrank < n_remainder) then
      ifrag_start = myrank * (n_per_proc + 1) + 1
      ifrag_end = ifrag_start + n_per_proc
    else
      ifrag_start = n_remainder * (n_per_proc + 1) + (myrank - n_remainder) * n_per_proc + 1
      ifrag_end = ifrag_start + n_per_proc - 1
    end if

    ifrag_start = max(1, min(n_frag, ifrag_start))
    ifrag_end = max(0, min(n_frag, ifrag_end))
  end subroutine get_fragment_range_for_rank

  !=======================================================================
  ! Return owner rank of fragment index ifrag for block distribution
  !=======================================================================
  integer function get_fragment_owner_rank(ifrag, n_frag, nprocs) result(owner_rank)
    implicit none
    integer, intent(in) :: ifrag, n_frag, nprocs
    integer :: n_per_proc, n_remainder, cutoff

    owner_rank = 0
    if (nprocs <= 0) return
    if (ifrag < 1 .or. ifrag > n_frag) return

    n_per_proc = n_frag / nprocs
    n_remainder = mod(n_frag, nprocs)
    cutoff = n_remainder * (n_per_proc + 1)

    if (ifrag <= cutoff) then
      owner_rank = (ifrag - 1) / (n_per_proc + 1)
    else
      owner_rank = n_remainder + (ifrag - cutoff - 1) / n_per_proc
    end if
  end function get_fragment_owner_rank


  logical function halo_axis_matches_direction(lo_ref, hi_ref, lo_nei, hi_nei, ngrid, dir) result(matches)
    implicit none
    integer, intent(in) :: lo_ref, hi_ref, lo_nei, hi_nei, ngrid, dir
    integer :: shift, nei_lo, nei_hi

    matches = .false.
    do shift = -ngrid, ngrid, ngrid
      nei_lo = lo_nei + shift
      nei_hi = hi_nei + shift
      select case (dir)
      case (-1)
        if (nei_hi == lo_ref - 1) then
          matches = .true.
          return
        end if
      case (0)
        if (max(lo_ref, nei_lo) <= min(hi_ref, nei_hi)) then
          matches = .true.
          return
        end if
      case (1)
        if (nei_lo == hi_ref + 1) then
          matches = .true.
          return
        end if
      end select
    end do
  end function halo_axis_matches_direction

  logical function fragment_matches_direction(dg_frag, ifrag, jfrag, dvec) result(matches)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ifrag, jfrag, dvec(3)

    integer :: axis, lo_ref, hi_ref, lo_nei, hi_nei

    matches = .true.
    do axis = 1, 3
      lo_ref = dg_frag%ixyz_frag(axis, ifrag)
      hi_ref = lo_ref + dg_frag%nxyz_domain(axis, ifrag) - 1
      lo_nei = dg_frag%ixyz_frag(axis, jfrag)
      hi_nei = lo_nei + dg_frag%nxyz_domain(axis, jfrag) - 1
      if (.not. halo_axis_matches_direction(lo_ref, hi_ref, lo_nei, hi_nei, dg_frag%lgnum_total(axis), dvec(axis))) then
        matches = .false.
        return
      end if
    end do
  end function fragment_matches_direction

  integer function map_global_to_phi_box_index(gidx, lb, ub, ngrid) result(idx)
    implicit none
    integer, intent(in) :: gidx, lb, ub, ngrid

    idx = modulo(gidx - 1, ngrid) + 1
    do while (idx < lb)
      idx = idx + ngrid
    end do
    do while (idx > ub)
      idx = idx - ngrid
    end do
  end function map_global_to_phi_box_index

  subroutine map_fragment_local_to_phi_box(dg_frag, ifrag, lx, ly, lz, px, py, pz)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ifrag, lx, ly, lz
    integer, intent(out) :: px, py, pz

    px = map_global_to_phi_box_index(dg_frag%ixyz_frag(1, ifrag) + lx - 1, &
      lbound(dg_frag%phi_frag, 1), ubound(dg_frag%phi_frag, 1), dg_frag%lgnum_total(1))
    py = map_global_to_phi_box_index(dg_frag%ixyz_frag(2, ifrag) + ly - 1, &
      lbound(dg_frag%phi_frag, 2), ubound(dg_frag%phi_frag, 2), dg_frag%lgnum_total(2))
    pz = map_global_to_phi_box_index(dg_frag%ixyz_frag(3, ifrag) + lz - 1, &
      lbound(dg_frag%phi_frag, 3), ubound(dg_frag%phi_frag, 3), dg_frag%lgnum_total(3))
  end subroutine map_fragment_local_to_phi_box

  subroutine compute_halo_axis_segment(dg_frag, ifrag, axis, dir, start_send_g, start_recv_g, len_raw)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ifrag, axis, dir
    integer, intent(out) :: start_send_g, start_recv_g, len_raw
    integer :: lo_frag, nloc, nbuf

    lo_frag = dg_frag%ixyz_frag(axis, ifrag)
    nloc = dg_frag%nxyz_domain(axis, ifrag)
    nbuf = dg_frag%nxyz_buffer(axis)

    select case (dir)
    case (0)
      start_send_g = lo_frag
      start_recv_g = lo_frag
      len_raw = nloc
    case (1)
      start_send_g = lo_frag + nloc - nbuf
      start_recv_g = lo_frag + nloc
      len_raw = nbuf
    case (-1)
      start_send_g = lo_frag
      start_recv_g = lo_frag - nbuf
      len_raw = nbuf
    case default
      start_send_g = lo_frag
      start_recv_g = lo_frag
      len_raw = 0
    end select
  end subroutine compute_halo_axis_segment

  integer function halo_axis_local_index(dg_frag, ifrag, axis, dir, ioff, recv_side) result(idx)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ifrag, axis, dir, ioff
    logical, intent(in) :: recv_side
    integer :: start_send_g, start_recv_g, len_raw
    integer :: gidx, lb, ub

    call compute_halo_axis_segment(dg_frag, ifrag, axis, dir, start_send_g, start_recv_g, len_raw)
    if (ioff < 1 .or. ioff > len_raw) then
      idx = 0
      return
    end if
    if (recv_side) then
      gidx = start_recv_g + ioff - 1
    else
      gidx = start_send_g + ioff - 1
    end if
    lb = lbound(dg_frag%phi_frag, axis)
    ub = ubound(dg_frag%phi_frag, axis)
    idx = map_global_to_phi_box_index(gidx, lb, ub, dg_frag%lgnum_total(axis))
  end function halo_axis_local_index

  subroutine get_halo_point_indices(dg_frag, halo, ix_buf, iy_buf, iz_buf, send_idx, recv_idx)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    type(halo_info), intent(in) :: halo
    integer, intent(in) :: ix_buf, iy_buf, iz_buf
    integer, intent(out) :: send_idx(3), recv_idx(3)

    send_idx(1) = halo_axis_local_index(dg_frag, halo%ifrag_dst, 1, halo%dvec(1), ix_buf, .false.)
    send_idx(2) = halo_axis_local_index(dg_frag, halo%ifrag_dst, 2, halo%dvec(2), iy_buf, .false.)
    send_idx(3) = halo_axis_local_index(dg_frag, halo%ifrag_dst, 3, halo%dvec(3), iz_buf, .false.)
    recv_idx(1) = halo_axis_local_index(dg_frag, halo%ifrag_dst, 1, halo%dvec(1), ix_buf, .true.)
    recv_idx(2) = halo_axis_local_index(dg_frag, halo%ifrag_dst, 2, halo%dvec(2), iy_buf, .true.)
    recv_idx(3) = halo_axis_local_index(dg_frag, halo%ifrag_dst, 3, halo%dvec(3), iz_buf, .true.)
  end subroutine get_halo_point_indices

  subroutine compute_halo_axis_window(dg_frag, ifrag, axis, dir, dsp_send, dsp_recv, length)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ifrag, axis, dir
    integer, intent(out) :: dsp_send, dsp_recv, length
    integer :: lb, ub
    integer :: start_send_g, start_recv_g, len_raw
    integer :: ioff, idx_send, idx_recv

    lb = dg_frag%mg%is(axis) - dg_frag%nxyz_buffer(axis)
    ub = dg_frag%mg%ie(axis) + dg_frag%nxyz_buffer(axis)
    call compute_halo_axis_segment(dg_frag, ifrag, axis, dir, start_send_g, start_recv_g, len_raw)

    length = 0
    dsp_send = lb - 1
    dsp_recv = lb - 1
    do ioff = 0, max(0, len_raw - 1)
      idx_send = map_global_to_phi_box_index(start_send_g + ioff, lb, ub, dg_frag%lgnum_total(axis))
      idx_recv = map_global_to_phi_box_index(start_recv_g + ioff, lb, ub, dg_frag%lgnum_total(axis))
      if (idx_send < lb .or. idx_send > ub) cycle
      if (idx_recv < lb .or. idx_recv > ub) cycle
      if (length == 0) then
        dsp_send = idx_send - 1
        dsp_recv = idx_recv - 1
      end if
      length = length + 1
    end do
  end subroutine compute_halo_axis_window

  !=======================================================================
  ! Initialize halo communication structures for fragment boundaries
  ! Following lcfo.f90 halo exchange pattern with periodic boundaries
  !=======================================================================
  subroutine init_halo_communication(dg_frag, info)
    use structures
    use communication, only: comm_summation, comm_is_root
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_parallel_info),  intent(in)    :: info

    integer :: nh(3), lx, ly, lz, i, n, ifrag, jfrag, neighbor_root_rank
    integer :: d(3)
    integer, allocatable :: id_tmp(:)
    integer :: ifrag_count

    ! Build MPI rank array for all fragments (comm_summation across all ranks)
    allocate(id_tmp(dg_frag%n_frag))
    id_tmp = 0
    ifrag_count = dg_frag%ifrag_end - dg_frag%ifrag_start + 1
    if (dg_frag%is_frag_root) then
      do i = 1, ifrag_count
        ifrag = dg_frag%ifrag_start + i - 1
        id_tmp(ifrag) = dg_frag%id + 1  ! +1 to distinguish from unset (0)
      end do
    end if
    call comm_summation(id_tmp, dg_frag%id_array, dg_frag%n_frag, dg_frag%icomm)
    dg_frag%id_array = dg_frag%id_array - 1  ! Convert back to 0-based rank

    ! Debug-only diagnostics for halo setup
#ifdef DEBUG
    if (comm_is_root(dg_frag%id)) then
      write(*,'(1x,a)') "  Fragment-to-rank mapping:"
      do ifrag = 1, dg_frag%n_frag
        write(*,'(4x,a,i2,a,i2)') "Fragment ", ifrag, " -> Rank ", dg_frag%id_array(ifrag)
      end do
      write(*,'(1x,a)') "  Fragment geometry (ixyz_frag):"
      do ifrag = 1, dg_frag%n_frag
        write(*,'(4x,a,i2,a,3i6,a,3i5,a)') &
          "Fragment ", ifrag, ": origin=(", dg_frag%ixyz_frag(1:3, ifrag), &
          "), domain=(", dg_frag%nxyz_domain(1:3, ifrag), ")"
      end do
      write(*,'(4x,a,3i5)') "Total grid size (lgnum_total): ", dg_frag%lgnum_total(1:3)
    end if
#endif

    if (ifrag_count <= 0) then
      dg_frag%n_halo = 0
      dg_frag%has_halo_exchange = .false.
      return
    end if

    ! Determine which directions require halo communication
    nh = 0
    do n = 1, 3
      if (dg_frag%num_fragment(n) > 1) nh(n) = 1
    end do

    if (allocated(dg_frag%halo)) deallocate(dg_frag%halo)
    allocate(dg_frag%halo(max(1, 26 * ifrag_count)))

    i = 0
    do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
      do lx = -nh(1), nh(1)
      do ly = -nh(2), nh(2)
      do lz = -nh(3), nh(3)
        if (lx == 0 .and. ly == 0 .and. lz == 0) cycle

        i = i + 1
        dg_frag%halo(i)%dvec(1:3) = [lx, ly, lz]
        dg_frag%halo(i)%id_dst = -1
        dg_frag%halo(i)%id_src = -1
        dg_frag%halo(i)%ifrag_src = -1
        dg_frag%halo(i)%ifrag_dst = ifrag

        do jfrag = 1, dg_frag%n_frag
          if (.not. fragment_matches_direction(dg_frag, ifrag, jfrag, dg_frag%halo(i)%dvec)) cycle

          neighbor_root_rank = get_fragment_group_root_rank(jfrag, dg_frag%nproc_frag)
          if (dg_frag%id_array(jfrag) /= neighbor_root_rank) then
            write(*,'(a,i0,a,i0,a,i0,a,i0)') &
              "[ERROR] DG-Fragment RT: inconsistent fragment root rank for ifrag=", ifrag, &
              " jfrag=", jfrag, " stored=", dg_frag%id_array(jfrag), " expected=", neighbor_root_rank
            stop "DG-Fragment RT: inconsistent fragment-group root rank"
          end if

          if (dg_frag%halo(i)%id_dst < 0) then
            dg_frag%halo(i)%id_dst = neighbor_root_rank + dg_frag%id_frag
          end if
          if (dg_frag%halo(i)%id_src < 0) then
            dg_frag%halo(i)%id_src = neighbor_root_rank + dg_frag%id_frag
            dg_frag%halo(i)%ifrag_src = jfrag
          else if (dg_frag%halo(i)%ifrag_src /= jfrag) then
            write(*,'(a,i0,a,3(i0,a),a,i0,a,i0)') &
              "[ERROR] DG-Fragment RT: ambiguous halo source for ifrag=", ifrag, " dir=(", &
              dg_frag%halo(i)%dvec(1), ",", dg_frag%halo(i)%dvec(2), ",", dg_frag%halo(i)%dvec(3), &
              ") first=", dg_frag%halo(i)%ifrag_src, " second=", jfrag
            stop "DG-Fragment RT: ambiguous halo source fragment"
          end if
        end do

        if (dg_frag%halo(i)%id_dst < 0 .or. dg_frag%halo(i)%id_src < 0) then
          write(*,'(a,i2,a,i2,a,i2,a,i2,a,i3,a,i3)') &
            "[ERROR] DG-Fragment RT: invalid halo neighbors for dir=(", &
            dg_frag%halo(i)%dvec(1), ",", dg_frag%halo(i)%dvec(2), ",", &
            dg_frag%halo(i)%dvec(3), ") (id_dst=", dg_frag%halo(i)%id_dst, &
            ", id_src=", dg_frag%halo(i)%id_src, ")"
          stop "DG-Fragment RT: dst, src"
        end if

        do n = 1, 3
          call compute_halo_axis_window(dg_frag, ifrag, n, dg_frag%halo(i)%dvec(n), &
                                        dg_frag%halo(i)%dsp_send(n), dg_frag%halo(i)%dsp_recv(n), &
                                        dg_frag%halo(i)%length(n))
        end do

#ifdef DEBUG
        if (comm_is_root(dg_frag%id)) then
          write(*,'(a,i2,a,i2,a,i2,a,i2,a,i2,a,i2,a,i2,a,i2,a,i2,a)') &
            "  [Rank ", dg_frag%id, "] Halo ", i, " dir=(", dg_frag%halo(i)%dvec(1), ",", &
            dg_frag%halo(i)%dvec(2), ",", dg_frag%halo(i)%dvec(3), " ): frag ", ifrag, &
            " -> dst frag ", dg_frag%halo(i)%ifrag_dst, " (rank ", dg_frag%halo(i)%id_dst, ")"
        end if
#endif

        allocate(dg_frag%halo(i)%buf_send(dg_frag%halo(i)%length(1), dg_frag%halo(i)%length(2), &
                                          dg_frag%halo(i)%length(3), dg_frag%nstate_frag, 1))
        allocate(dg_frag%halo(i)%buf_recv(dg_frag%halo(i)%length(1), dg_frag%halo(i)%length(2), &
                                          dg_frag%halo(i)%length(3), dg_frag%nstate_frag, 1))
        allocate(dg_frag%halo(i)%buf_send_c(dg_frag%halo(i)%length(1), dg_frag%halo(i)%length(2), &
                                            dg_frag%halo(i)%length(3), dg_frag%nstate_frag, 1))
        allocate(dg_frag%halo(i)%buf_recv_c(dg_frag%halo(i)%length(1), dg_frag%halo(i)%length(2), &
                                            dg_frag%halo(i)%length(3), dg_frag%nstate_frag, 1))

      end do
      end do
      end do
    end do

    dg_frag%n_halo = i
    dg_frag%has_halo_exchange = (dg_frag%n_halo > 0)

    if (comm_is_root(dg_frag%id)) then
      write(*,'(1x,a,i0,a)') "Halo communication initialized: ", dg_frag%n_halo, " neighbor regions"
    end if

    deallocate(id_tmp)

  end subroutine init_halo_communication

  !=======================================================================
  ! Exchange halo regions for phi_frag between neighboring fragments
  ! Must be called after any modification to phi_frag interior
  ! System boundaries use PERIODIC boundary conditions
  !=======================================================================
  subroutine exchange_phi_frag_halo(dg_frag)
    use communication, only: comm_isend, comm_irecv, comm_wait_all
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag

    integer :: i_halo, ix, iy, iz, istate
    integer :: l(3), i_local
    integer :: ifrag_send, ifrag_recv, dir_code
    integer :: itag_send, itag_recv
    integer, allocatable :: ireq_send(:), ireq_recv(:)
    logical :: use_complex
    integer :: send_idx(3), recv_idx(3)

    if (.not. dg_frag%has_halo_exchange) return
    if (dg_frag%n_halo <= 0) return

    use_complex = allocated(dg_frag%phi_frag_c)
    allocate(ireq_send(dg_frag%n_halo), ireq_recv(dg_frag%n_halo))

    do i_halo = 1, dg_frag%n_halo
      i_local = dg_frag%halo(i_halo)%ifrag_dst - dg_frag%ifrag_start + 1
      if (i_local < 1 .or. i_local > (dg_frag%ifrag_end - dg_frag%ifrag_start + 1)) then
        stop "DG-Fragment RT: invalid i_local in exchange_phi_frag_halo"
      end if

      l = dg_frag%halo(i_halo)%length
      do istate = 1, dg_frag%nstate_frag
      do iz = 1, l(3)
      do iy = 1, l(2)
      do ix = 1, l(1)
        call get_halo_point_indices(dg_frag, dg_frag%halo(i_halo), ix, iy, iz, send_idx, recv_idx)
        if (use_complex) then
          dg_frag%halo(i_halo)%buf_send_c(ix, iy, iz, istate, 1) = &
            dg_frag%phi_frag_c(send_idx(1), send_idx(2), send_idx(3), istate, i_local)
          dg_frag%halo(i_halo)%buf_send(ix, iy, iz, istate, 1) = &
            real(dg_frag%halo(i_halo)%buf_send_c(ix, iy, iz, istate, 1), kind=8)
        else
          dg_frag%halo(i_halo)%buf_send(ix, iy, iz, istate, 1) = &
            dg_frag%phi_frag(send_idx(1), send_idx(2), send_idx(3), istate, i_local)
        end if
      end do
      end do
      end do
      end do

      ifrag_send = dg_frag%halo(i_halo)%ifrag_dst
      dir_code = (dg_frag%halo(i_halo)%dvec(1) + 1) * 9 + &
                 (dg_frag%halo(i_halo)%dvec(2) + 1) * 3 + &
                 (dg_frag%halo(i_halo)%dvec(3) + 1)
      itag_send = (ifrag_send - 1) * 27 + dir_code
      if (use_complex) then
        ireq_send(i_halo) = comm_isend(dg_frag%halo(i_halo)%buf_send_c, &
                                       dg_frag%halo(i_halo)%id_dst, &
                                       itag_send, dg_frag%icomm)
      else
        ireq_send(i_halo) = comm_isend(dg_frag%halo(i_halo)%buf_send, &
                                       dg_frag%halo(i_halo)%id_dst, &
                                       itag_send, dg_frag%icomm)
      end if

      ifrag_recv = dg_frag%halo(i_halo)%ifrag_src
      itag_recv = (ifrag_recv - 1) * 27 + (26 - dir_code)
      if (use_complex) then
        ireq_recv(i_halo) = comm_irecv(dg_frag%halo(i_halo)%buf_recv_c, &
                                       dg_frag%halo(i_halo)%id_src, &
                                       itag_recv, dg_frag%icomm)
      else
        ireq_recv(i_halo) = comm_irecv(dg_frag%halo(i_halo)%buf_recv, &
                                       dg_frag%halo(i_halo)%id_src, &
                                       itag_recv, dg_frag%icomm)
      end if
    end do

    call comm_wait_all(ireq_recv)
    call comm_wait_all(ireq_send)

    do i_halo = 1, dg_frag%n_halo
      i_local = dg_frag%halo(i_halo)%ifrag_dst - dg_frag%ifrag_start + 1
      l = dg_frag%halo(i_halo)%length
      do istate = 1, dg_frag%nstate_frag
      do iz = 1, l(3)
      do iy = 1, l(2)
      do ix = 1, l(1)
        call get_halo_point_indices(dg_frag, dg_frag%halo(i_halo), ix, iy, iz, send_idx, recv_idx)
        if (use_complex) then
          dg_frag%phi_frag_c(recv_idx(1), recv_idx(2), recv_idx(3), istate, i_local) = &
            dg_frag%halo(i_halo)%buf_recv_c(ix, iy, iz, istate, 1)
          dg_frag%phi_frag(recv_idx(1), recv_idx(2), recv_idx(3), istate, i_local) = &
            real(dg_frag%halo(i_halo)%buf_recv_c(ix, iy, iz, istate, 1), kind=8)
          dg_frag%halo(i_halo)%buf_recv(ix, iy, iz, istate, 1) = &
            real(dg_frag%halo(i_halo)%buf_recv_c(ix, iy, iz, istate, 1), kind=8)
        else
          dg_frag%phi_frag(recv_idx(1), recv_idx(2), recv_idx(3), istate, i_local) = &
            dg_frag%halo(i_halo)%buf_recv(ix, iy, iz, istate, 1)
        end if
      end do
      end do
      end do
      end do
    end do

    deallocate(ireq_send, ireq_recv)

  end subroutine exchange_phi_frag_halo
