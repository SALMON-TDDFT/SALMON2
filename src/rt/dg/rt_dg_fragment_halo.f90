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
      stop "DG-Fragment RT: halo block outside local phi box"
    end if
    if ((send_hi - send_lo) /= (recv_hi - recv_lo)) then
      stop "DG-Fragment RT: halo send/recv block mismatch"
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
    use unusedvar_mod, only: salmon_unusedvar
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_parallel_info),  intent(in)    :: info

    ! Halo communication is disabled by design in this branch.
    call salmon_unusedvar(info)
    dg_frag%n_halo = 0
    dg_frag%has_halo_exchange = .false.

  end subroutine init_halo_communication

  !=======================================================================
  ! Exchange halo regions for phi_frag between neighboring fragments
  ! Must be called after any modification to phi_frag interior
  ! System boundaries use PERIODIC boundary conditions
  !=======================================================================
  subroutine exchange_phi_frag_halo(dg_frag)
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag

    ! Halo communication is disabled by design in this branch.
    dg_frag%n_halo = 0
    dg_frag%has_halo_exchange = .false.
  end subroutine exchange_phi_frag_halo
