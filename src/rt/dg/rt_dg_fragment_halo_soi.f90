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
    if (idx < lb) then
      idx = idx + ((lb - idx + ngrid - 1) / ngrid) * ngrid
    end if
    if (idx > ub) then
      idx = idx - ((idx - ub + ngrid - 1) / ngrid) * ngrid
    end if
    if (idx < lb .or. idx > ub) idx = 0
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

  subroutine compute_halo_axis_block(dg_frag, axis, dir, send_lo, send_hi, recv_lo, recv_hi, length)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: axis, dir
    integer, intent(out) :: send_lo, send_hi, recv_lo, recv_hi, length
    integer :: core_lo, core_hi, buf, phi_lo, phi_hi

    core_lo = dg_frag%mg%is(axis)
    core_hi = dg_frag%mg%ie(axis)
    buf = dg_frag%nxyz_buffer(axis)
    phi_lo = lbound(dg_frag%phi_frag, axis)
    phi_hi = ubound(dg_frag%phi_frag, axis)

    select case (dir)
    case (0)
      send_lo = phi_lo
      send_hi = phi_hi
      recv_lo = phi_lo
      recv_hi = phi_hi
    case (1)
      send_lo = core_hi - buf + 1
      send_hi = core_hi
      recv_lo = core_hi + 1
      recv_hi = core_hi + buf
    case (-1)
      send_lo = core_lo
      send_hi = core_lo + buf - 1
      recv_lo = core_lo - buf
      recv_hi = core_lo - 1
    case default
      send_lo = 1
      send_hi = 0
      recv_lo = 1
      recv_hi = 0
    end select

    if (send_lo < phi_lo .or. send_hi > phi_hi .or. recv_lo < phi_lo .or. recv_hi > phi_hi) then
      stop "DG-Fragment RT SOI: halo block outside local phi box"
    end if
    if ((send_hi - send_lo) /= (recv_hi - recv_lo)) then
      stop "DG-Fragment RT SOI: halo send/recv block mismatch"
    end if
    length = max(0, send_hi - send_lo + 1)
  end subroutine compute_halo_axis_block

  subroutine get_halo_block_point_indices(halo, ix_buf, iy_buf, iz_buf, send_idx, recv_idx)
    implicit none
    type(halo_info), intent(in) :: halo
    integer, intent(in) :: ix_buf, iy_buf, iz_buf
    integer, intent(out) :: send_idx(3), recv_idx(3)

    send_idx(1) = halo%send_lo(1) + ix_buf - 1
    send_idx(2) = halo%send_lo(2) + iy_buf - 1
    send_idx(3) = halo%send_lo(3) + iz_buf - 1
    recv_idx(1) = halo%recv_lo(1) + ix_buf - 1
    recv_idx(2) = halo%recv_lo(2) + iy_buf - 1
    recv_idx(3) = halo%recv_lo(3) + iz_buf - 1
  end subroutine get_halo_block_point_indices

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
          call compute_halo_axis_block(dg_frag, n, dg_frag%halo(i)%dvec(n), &
                                       dg_frag%halo(i)%send_lo(n), dg_frag%halo(i)%send_hi(n), &
                                       dg_frag%halo(i)%recv_lo(n), dg_frag%halo(i)%recv_hi(n), &
                                       dg_frag%halo(i)%length(n))
          dg_frag%halo(i)%dsp_send(n) = dg_frag%halo(i)%send_lo(n) - 1
          dg_frag%halo(i)%dsp_recv(n) = dg_frag%halo(i)%recv_lo(n) - 1
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
    integer :: lb1, ub1, lb2, ub2, lb3, ub3

    if (.not. dg_frag%has_halo_exchange) return
    if (dg_frag%n_halo <= 0) return

    use_complex = allocated(dg_frag%phi_frag_c)
    allocate(ireq_send(dg_frag%n_halo), ireq_recv(dg_frag%n_halo))
    lb1 = lbound(dg_frag%phi_frag, 1)
    ub1 = ubound(dg_frag%phi_frag, 1)
    lb2 = lbound(dg_frag%phi_frag, 2)
    ub2 = ubound(dg_frag%phi_frag, 2)
    lb3 = lbound(dg_frag%phi_frag, 3)
    ub3 = ubound(dg_frag%phi_frag, 3)

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
        send_idx(1) = dg_frag%halo(i_halo)%send_lo(1) + ix - 1
        send_idx(2) = dg_frag%halo(i_halo)%send_lo(2) + iy - 1
        send_idx(3) = dg_frag%halo(i_halo)%send_lo(3) + iz - 1
        recv_idx(1) = dg_frag%halo(i_halo)%recv_lo(1) + ix - 1
        recv_idx(2) = dg_frag%halo(i_halo)%recv_lo(2) + iy - 1
        recv_idx(3) = dg_frag%halo(i_halo)%recv_lo(3) + iz - 1
        if (send_idx(1) < lb1 .or. send_idx(1) > ub1 .or. &
            send_idx(2) < lb2 .or. send_idx(2) > ub2 .or. &
            send_idx(3) < lb3 .or. send_idx(3) > ub3) then
          write(*,'(1x,a,1x,3(i0,1x),a,1x,3(i0,1x),a,1x,3(i0,1x),a,i0,a,i0,a,3(i0,1x))') &
            "DG-Fragment RT SOI halo send index out of range:", send_idx(1), send_idx(2), send_idx(3), &
            "dsp=", dg_frag%halo(i_halo)%dsp_send(1), dg_frag%halo(i_halo)%dsp_send(2), dg_frag%halo(i_halo)%dsp_send(3), &
            "len=", l(1), l(2), l(3), " i_local=", i_local, " i_halo=", i_halo, &
            " dvec=", dg_frag%halo(i_halo)%dvec(1), dg_frag%halo(i_halo)%dvec(2), dg_frag%halo(i_halo)%dvec(3)
          write(*,'(1x,a,6(i0,1x),a,6(i0,1x))') &
            "        phi_bounds=", lb1, ub1, lb2, ub2, lb3, ub3, &
            "mg_box=", dg_frag%mg%is(1), dg_frag%mg%ie(1), dg_frag%mg%is(2), dg_frag%mg%ie(2), dg_frag%mg%is(3), dg_frag%mg%ie(3)
          write(*,'(1x,a,i0,a,i0,a,i0,a,i0)') "        halo send context: rank=", dg_frag%id, &
            " id_frag=", dg_frag%id_frag, " ifrag_src=", dg_frag%halo(i_halo)%ifrag_src, " ifrag_dst=", dg_frag%halo(i_halo)%ifrag_dst
          stop "DG-Fragment RT SOI: halo send index out of range"
        end if
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
        send_idx(1) = dg_frag%halo(i_halo)%send_lo(1) + ix - 1
        send_idx(2) = dg_frag%halo(i_halo)%send_lo(2) + iy - 1
        send_idx(3) = dg_frag%halo(i_halo)%send_lo(3) + iz - 1
        recv_idx(1) = dg_frag%halo(i_halo)%recv_lo(1) + ix - 1
        recv_idx(2) = dg_frag%halo(i_halo)%recv_lo(2) + iy - 1
        recv_idx(3) = dg_frag%halo(i_halo)%recv_lo(3) + iz - 1
        if (recv_idx(1) < lb1 .or. recv_idx(1) > ub1 .or. &
            recv_idx(2) < lb2 .or. recv_idx(2) > ub2 .or. &
            recv_idx(3) < lb3 .or. recv_idx(3) > ub3) then
          write(*,'(1x,a,1x,3(i0,1x),a,1x,3(i0,1x),a,1x,3(i0,1x),a,i0,a,i0,a,3(i0,1x))') &
            "DG-Fragment RT SOI halo recv index out of range:", recv_idx(1), recv_idx(2), recv_idx(3), &
            "dsp=", dg_frag%halo(i_halo)%dsp_recv(1), dg_frag%halo(i_halo)%dsp_recv(2), dg_frag%halo(i_halo)%dsp_recv(3), &
            "len=", l(1), l(2), l(3), " i_local=", i_local, " i_halo=", i_halo, &
            " dvec=", dg_frag%halo(i_halo)%dvec(1), dg_frag%halo(i_halo)%dvec(2), dg_frag%halo(i_halo)%dvec(3)
          write(*,'(1x,a,6(i0,1x),a,6(i0,1x))') &
            "        phi_bounds=", lb1, ub1, lb2, ub2, lb3, ub3, &
            "mg_box=", dg_frag%mg%is(1), dg_frag%mg%ie(1), dg_frag%mg%is(2), dg_frag%mg%ie(2), dg_frag%mg%is(3), dg_frag%mg%ie(3)
          write(*,'(1x,a,i0,a,i0,a,i0,a,i0)') "        halo recv context: rank=", dg_frag%id, &
            " id_frag=", dg_frag%id_frag, " ifrag_src=", dg_frag%halo(i_halo)%ifrag_src, " ifrag_dst=", dg_frag%halo(i_halo)%ifrag_dst
          stop "DG-Fragment RT SOI: halo recv index out of range"
        end if
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
