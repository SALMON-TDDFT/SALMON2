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
    integer :: rank, ifrag_start, ifrag_end

    owner_rank = 0
    if (nprocs <= 0) return
    if (ifrag < 1 .or. ifrag > n_frag) return

    do rank = 0, nprocs - 1
      call get_fragment_range_for_rank(n_frag, rank, nprocs, ifrag_start, ifrag_end)
      if (ifrag >= ifrag_start .and. ifrag <= ifrag_end) then
        owner_rank = rank
        return
      end if
    end do

    owner_rank = nprocs - 1
  end function get_fragment_owner_rank

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

    integer :: nh(3), lx, ly, lz, i, n, ifrag, jfrag, i_local
    integer :: ir1(3), ir2(3), d(3)
    integer :: id_tmp(dg_frag%n_frag)
    integer :: ifrag_count

    ! Build MPI rank array for all fragments (comm_summation across all ranks)
    id_tmp = 0
    ifrag_count = dg_frag%ifrag_end - dg_frag%ifrag_start + 1
    do i = 1, ifrag_count
      ifrag = dg_frag%ifrag_start + i - 1
      id_tmp(ifrag) = dg_frag%id + 1  ! +1 to distinguish from unset (0)
    end do
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
      i_local = ifrag - dg_frag%ifrag_start + 1
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
          ir1(1:3) = dg_frag%ixyz_frag(1:3, jfrag)

          ! destination neighbor (+direction)
          ir2(1:3) = dg_frag%ixyz_frag(1:3, ifrag) + &
                     dg_frag%halo(i)%dvec(1:3) * dg_frag%nxyz_domain(1:3, ifrag)
          d(1:3) = mod(ir1(1:3) - ir2(1:3), dg_frag%lgnum_total(1:3))
          if (d(1) == 0 .and. d(2) == 0 .and. d(3) == 0 .and. dg_frag%halo(i)%id_dst < 0) then
            dg_frag%halo(i)%id_dst = dg_frag%id_array(jfrag)
          end if

          ! source neighbor (-direction)
          ir2(1:3) = dg_frag%ixyz_frag(1:3, ifrag) - &
                     dg_frag%halo(i)%dvec(1:3) * dg_frag%nxyz_domain(1:3, ifrag)
          d(1:3) = mod(ir1(1:3) - ir2(1:3), dg_frag%lgnum_total(1:3))
          if (d(1) == 0 .and. d(2) == 0 .and. d(3) == 0 .and. dg_frag%halo(i)%id_src < 0) then
            dg_frag%halo(i)%id_src = dg_frag%id_array(jfrag)
            dg_frag%halo(i)%ifrag_src = jfrag
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
          select case (dg_frag%halo(i)%dvec(n))
          case (0)
            dg_frag%halo(i)%length(n) = dg_frag%nxyz_domain(n, ifrag)
            dg_frag%halo(i)%dsp_send(n) = 0
            dg_frag%halo(i)%dsp_recv(n) = 0
          case (1)
            dg_frag%halo(i)%length(n) = dg_frag%nxyz_buffer(n)
            dg_frag%halo(i)%dsp_send(n) = dg_frag%nxyz_domain(n, ifrag) - dg_frag%nxyz_buffer(n)
            dg_frag%halo(i)%dsp_recv(n) = dg_frag%nxyz_domain(n, ifrag)
          case (-1)
            dg_frag%halo(i)%length(n) = dg_frag%nxyz_buffer(n)
            dg_frag%halo(i)%dsp_send(n) = 0
            dg_frag%halo(i)%dsp_recv(n) = -dg_frag%nxyz_buffer(n)
          end select
        end do

        if (comm_is_root(dg_frag%id)) then
          write(*,'(a,i2,a,i2,a,i2,a,i2,a,i2,a,i2,a,i2,a,i2,a,i2,a)') &
            "  [Rank ", dg_frag%id, "] Halo ", i, " dir=(", dg_frag%halo(i)%dvec(1), ",", &
            dg_frag%halo(i)%dvec(2), ",", dg_frag%halo(i)%dvec(3), " ): frag ", ifrag, &
            " -> dst frag ", dg_frag%halo(i)%ifrag_dst, " (rank ", dg_frag%halo(i)%id_dst, ")"
        end if

      end do
      end do
      end do
    end do

    dg_frag%n_halo = i
    dg_frag%has_halo_exchange = (dg_frag%n_halo > 0)

    if (comm_is_root(dg_frag%id)) then
      write(*,'(1x,a,i0,a)') "Halo communication initialized: ", dg_frag%n_halo, " neighbor regions"
    end if

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
    integer :: l(3), d(3), i_local
    integer :: ifrag_send, ifrag_recv, dir_code
    integer :: itag_send, itag_recv
    integer, allocatable :: ireq_send(:), ireq_recv(:)

    if (.not. dg_frag%has_halo_exchange) return
    if (dg_frag%n_halo <= 0) return

    allocate(ireq_send(dg_frag%n_halo), ireq_recv(dg_frag%n_halo))

    do i_halo = 1, dg_frag%n_halo
      i_local = dg_frag%halo(i_halo)%ifrag_dst - dg_frag%ifrag_start + 1
      if (i_local < 1 .or. i_local > (dg_frag%ifrag_end - dg_frag%ifrag_start + 1)) then
        stop "DG-Fragment RT: invalid i_local in exchange_phi_frag_halo"
      end if

      l = dg_frag%halo(i_halo)%length
      d = dg_frag%halo(i_halo)%dsp_send

      if (allocated(dg_frag%halo(i_halo)%buf_send)) deallocate(dg_frag%halo(i_halo)%buf_send)
      if (allocated(dg_frag%halo(i_halo)%buf_recv)) deallocate(dg_frag%halo(i_halo)%buf_recv)
      allocate(dg_frag%halo(i_halo)%buf_send(l(1), l(2), l(3), dg_frag%nstate_frag, 1))
      allocate(dg_frag%halo(i_halo)%buf_recv(l(1), l(2), l(3), dg_frag%nstate_frag, 1))

      do istate = 1, dg_frag%nstate_frag
      do iz = 1, l(3)
      do iy = 1, l(2)
      do ix = 1, l(1)
        dg_frag%halo(i_halo)%buf_send(ix, iy, iz, istate, 1) = &
          dg_frag%phi_frag(d(1) + ix, d(2) + iy, d(3) + iz, istate, i_local)
      end do
      end do
      end do
      end do

      ifrag_send = dg_frag%halo(i_halo)%ifrag_dst
      dir_code = (dg_frag%halo(i_halo)%dvec(1) + 1) * 9 + &
                 (dg_frag%halo(i_halo)%dvec(2) + 1) * 3 + &
                 (dg_frag%halo(i_halo)%dvec(3) + 1)
      itag_send = (ifrag_send - 1) * 27 + dir_code
      ireq_send(i_halo) = comm_isend(dg_frag%halo(i_halo)%buf_send, &
                                     dg_frag%halo(i_halo)%id_dst, &
                                     itag_send, dg_frag%icomm)

      ifrag_recv = dg_frag%halo(i_halo)%ifrag_src
      itag_recv = (ifrag_recv - 1) * 27 + dir_code
      ireq_recv(i_halo) = comm_irecv(dg_frag%halo(i_halo)%buf_recv, &
                                     dg_frag%halo(i_halo)%id_src, &
                                     itag_recv, dg_frag%icomm)
    end do

    call comm_wait_all(ireq_recv)
    call comm_wait_all(ireq_send)

    do i_halo = 1, dg_frag%n_halo
      i_local = dg_frag%halo(i_halo)%ifrag_dst - dg_frag%ifrag_start + 1
      l = dg_frag%halo(i_halo)%length
      d = dg_frag%halo(i_halo)%dsp_recv

      do istate = 1, dg_frag%nstate_frag
      do iz = 1, l(3)
      do iy = 1, l(2)
      do ix = 1, l(1)
        dg_frag%phi_frag(d(1) + ix, d(2) + iy, d(3) + iz, istate, i_local) = &
          dg_frag%halo(i_halo)%buf_recv(ix, iy, iz, istate, 1)
      end do
      end do
      end do
      end do
    end do

    deallocate(ireq_send, ireq_recv)

  end subroutine exchange_phi_frag_halo
