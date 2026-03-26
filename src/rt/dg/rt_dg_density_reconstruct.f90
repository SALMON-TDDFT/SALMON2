  subroutine calculate_density_from_fragments(dg_frag, system, mg, rho, rho_s)
    use structures
    use salmon_global, only: nelec, nelec_spin
    use communication, only: comm_summation, comm_isend, comm_irecv, comm_wait_all
    use rt_dg_fragment_ops, only: refresh_pw_coef_cache
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system),     intent(in)    :: system
    type(s_rgrid),          intent(in)    :: mg
    type(s_scalar),         intent(inout) :: rho
    type(s_scalar),         intent(inout) :: rho_s(system%nspin)

    integer, parameter :: rho_tag_base = 4101
    integer, parameter :: w_tag_base = 4201
    integer :: ifrag, io, i_local, ispin
    integer :: istate_frag
    integer :: ix, iy, iz, ixg, iyg, izg, owner_rank
    integer :: ig_i, nbf, ipw, n_pw
    integer :: nxyz(3), ixyz0(3)
    integer :: nocc_per_spin, nocc_spin, nocc_cache
    integer :: irank, nreq_send, nreq_recv, ireq, slot, npts
    integer :: igrid0, igrid, ngrid, npt_blk, io0, nbatch, tmp_idx, ipw0, npw_blk, ipw_loc
    integer, parameter :: grid_block_size = 512, state_block_size = 64, pw_block_size = 128
    logical, parameter :: enable_density_trace = .false.
    real(8) :: occ_factor
    real(8) :: phi_i, rho_contrib
    real(8) :: total_charge, total_charge_local, scale_rho, elec_num_scaled_local
    real(8) :: rx, ry, rz, boxL(3), inv_sqrt_vol, theta
    real(8) :: t_total0, t_total1, t_cache0, t_cache1
    real(8) :: t_project0, t_project1, t_comm0, t_comm1, t_norm0, t_norm1
    real(8) :: time_cache, time_project, time_comm, time_norm
    real(8), allocatable :: w_local(:,:,:)
    integer, allocatable :: req_send(:), req_recv(:)
    integer, allocatable :: ix_buf(:), iy_buf(:), iz_buf(:), owner_buf(:), ixg_buf(:), iyg_buf(:), izg_buf(:)
    type(s_scalar), allocatable :: rho_send(:), w_send(:), rho_recv(:), w_recv(:)
    real(8), allocatable :: phi_blk(:,:), rho_blk(:)
    complex(8), allocatable :: coef_blk(:,:), psi_blk(:,:), phase_blk(:,:), coef_pw_blk(:,:)

    rho%f = 0.0d0
    do ispin = 1, system%nspin
      rho_s(ispin)%f = 0.0d0
    end do
    call cpu_time(t_total0)
    time_cache = 0.0d0
    time_project = 0.0d0
    time_comm = 0.0d0
    time_norm = 0.0d0
    if (enable_density_trace .and. dg_frag%id == 0) then
      write(*,'(1x,a)') "        density trace: stage=entry"
      flush(6)
    end if

    if (.not. allocated(dg_frag%phi_frag)) return
    if ((.not. allocated(dg_frag%density_phi_cache)) .or. (.not. dg_frag%density_phi_cache_valid)) then
      call rebuild_density_phi_cache()
    end if

    allocate(w_local(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))
    allocate(ix_buf(grid_block_size), iy_buf(grid_block_size), iz_buf(grid_block_size))
    allocate(owner_buf(grid_block_size), ixg_buf(grid_block_size), iyg_buf(grid_block_size), izg_buf(grid_block_size))
    allocate(phi_blk(grid_block_size, dg_frag%nstate_frag))
    allocate(rho_blk(grid_block_size))
    allocate(coef_blk(dg_frag%nstate_frag, state_block_size))
    allocate(psi_blk(grid_block_size, state_block_size))
    allocate(rho_send(0:dg_frag%isize-1), w_send(0:dg_frag%isize-1))
    allocate(rho_recv(0:dg_frag%isize-1), w_recv(0:dg_frag%isize-1))

    rho%f = 0.0d0
    w_local = 0.0d0
    ! Closed-shell fallback: nelec = 2 * nocc_per_spin.
    nocc_per_spin = min(dg_frag%nstate_tot, int(nelec / 2.0d0 + 1.0d-12))
    nocc_cache = nocc_per_spin
    if (system%nspin == 2 .and. sum(nelec_spin(:)) > 0) then
      nocc_cache = min(dg_frag%nstate_tot, maxval(nelec_spin(1:system%nspin)))
    end if
    occ_factor = 2.0d0 / real(system%nspin, 8)
    n_pw = 0
    if (dg_frag%use_plane_wave_basis .and. allocated(dg_frag%coef_pw) .and. allocated(dg_frag%k_pw)) then
      n_pw = dg_frag%n_plane_waves
    end if
    if (n_pw > 0) then
      allocate(phase_blk(grid_block_size, pw_block_size))
      allocate(coef_pw_blk(pw_block_size, state_block_size))
    end if
    boxL(1) = dg_frag%hgs(1) * real(mg%num(1), 8)
    boxL(2) = dg_frag%hgs(2) * real(mg%num(2), 8)
    boxL(3) = dg_frag%hgs(3) * real(mg%num(3), 8)
    inv_sqrt_vol = 1.0d0 / sqrt(max(1.0d-16, boxL(1) * boxL(2) * boxL(3)))

    do irank = 0, dg_frag%isize - 1
      if (irank == dg_frag%id) cycle
      if (allocated(dg_frag%density_send_count)) then
        npts = dg_frag%density_send_count(irank)
      else
        npts = 0
      end if
      if (npts > 0) then
        allocate(rho_send(irank)%f(1:npts, 1:1, 1:1), w_send(irank)%f(1:npts, 1:1, 1:1))
        rho_send(irank)%f = 0.0d0
        w_send(irank)%f = 0.0d0
      end if
      if (allocated(dg_frag%density_recv_map)) then
        npts = dg_frag%density_recv_map(irank)%npts
      else
        npts = 0
      end if
      if (npts > 0) then
        allocate(rho_recv(irank)%f(1:npts, 1:1, 1:1), w_recv(irank)%f(1:npts, 1:1, 1:1))
        rho_recv(irank)%f = 0.0d0
        w_recv(irank)%f = 0.0d0
      end if
    end do

    if (n_pw > 0) then
      call cpu_time(t_cache0)
      if ((.not. allocated(dg_frag%coef_pw_full_cache)) .or. dg_frag%coef_pw_full_cache_nstate < nocc_cache) then
        call refresh_pw_coef_cache(dg_frag, nocc_cache)
      end if
      call cpu_time(t_cache1)
      time_cache = time_cache + (t_cache1 - t_cache0)
    end if
    if (enable_density_trace .and. dg_frag%id == 0) then
      write(*,'(1x,a,a,1pe12.4)') "        density trace: stage=after-pw-cache dt=", "", time_cache
      flush(6)
    end if

    if (dg_frag%is_frag_root) then
      if (enable_density_trace .and. dg_frag%id == 0) then
        write(*,'(1x,a)') "        density trace: stage=before-frag-loop"
        flush(6)
      end if
      i_local = 0
      do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
        i_local = i_local + 1

        nxyz(1:3) = dg_frag%nxyz_domain(1:3, ifrag)
        ixyz0(1:3) = dg_frag%ixyz_frag(1:3, ifrag)

        do iz = 1, nxyz(3)
          do iy = 1, nxyz(2)
            do ix = 1, nxyz(1)
              ixg = ixyz0(1) + ix - 1
              iyg = ixyz0(2) + iy - 1
              izg = ixyz0(3) + iz - 1

              ixg = mod(ixg - 1, mg%num(1)) + 1
              iyg = mod(iyg - 1, mg%num(2)) + 1
              izg = mod(izg - 1, mg%num(3)) + 1

              if (allocated(dg_frag%density_owner_map)) then
                owner_rank = dg_frag%density_owner_map(ix, iy, iz, i_local)
                ixg = dg_frag%density_ixg_map(ix, iy, iz, i_local)
                iyg = dg_frag%density_iyg_map(ix, iy, iz, i_local)
                izg = dg_frag%density_izg_map(ix, iy, iz, i_local)
              else
                owner_rank = find_grid_owner(ixg, iyg, izg)
              end if
                  if (owner_rank == dg_frag%id) then
                    w_local(ixg, iyg, izg) = w_local(ixg, iyg, izg) + 1.0d0
                  else if (allocated(dg_frag%density_send_slot_map)) then
                    slot = dg_frag%density_send_slot_map(ix, iy, iz, i_local)
                    if (slot > 0) then
                      w_send(owner_rank)%f(slot, 1, 1) = w_send(owner_rank)%f(slot, 1, 1) + 1.0d0
                    end if
                  end if
          end do
        end do
      end do

      call cpu_time(t_project0)
      ngrid = nxyz(1) * nxyz(2) * nxyz(3)
      do igrid0 = 1, ngrid, grid_block_size
        npt_blk = min(grid_block_size, ngrid - igrid0 + 1)
        do igrid = 1, npt_blk
          tmp_idx = igrid0 + igrid - 2
          ix = mod(tmp_idx, nxyz(1)) + 1
          tmp_idx = tmp_idx / nxyz(1)
          iy = mod(tmp_idx, nxyz(2)) + 1
          iz = tmp_idx / nxyz(2) + 1
          ix_buf(igrid) = ix
          iy_buf(igrid) = iy
          iz_buf(igrid) = iz
          if (allocated(dg_frag%density_owner_map)) then
            ixg = dg_frag%density_ixg_map(ix, iy, iz, i_local)
            iyg = dg_frag%density_iyg_map(ix, iy, iz, i_local)
            izg = dg_frag%density_izg_map(ix, iy, iz, i_local)
            owner_rank = dg_frag%density_owner_map(ix, iy, iz, i_local)
          else
            ixg = mod(ixyz0(1) + ix - 2, mg%num(1)) + 1
            iyg = mod(ixyz0(2) + iy - 2, mg%num(2)) + 1
            izg = mod(ixyz0(3) + iz - 2, mg%num(3)) + 1
            owner_rank = find_grid_owner(ixg, iyg, izg)
          end if
          ixg_buf(igrid) = ixg
          iyg_buf(igrid) = iyg
          izg_buf(igrid) = izg
          owner_buf(igrid) = owner_rank
        end do

        do ispin = 1, system%nspin
          nocc_spin = nocc_per_spin
          if (system%nspin == 2 .and. sum(nelec_spin(:)) > 0) then
            nocc_spin = min(dg_frag%nstate_tot, nelec_spin(ispin))
          end if
          nbf = dg_frag%n_basis(ifrag, ispin)
          if (nbf <= 0 .or. nocc_spin <= 0) cycle

          phi_blk(1:npt_blk, 1:nbf) = dg_frag%density_phi_cache(igrid0:igrid0+npt_blk-1, 1:nbf, i_local)

          do io0 = 1, nocc_spin, state_block_size
            nbatch = min(state_block_size, nocc_spin - io0 + 1)
            coef_blk(1:nbf, 1:nbatch) = (0.0d0, 0.0d0)
            do io = 1, nbatch
              do istate_frag = 1, nbf
                ig_i = dg_frag%index_basis(istate_frag, ifrag, ispin)
                if (ig_i < 1 .or. ig_i > dg_frag%n_mat_max) cycle
                coef_blk(istate_frag, io) = dg_frag%coef(ig_i, io0 + io - 1, ispin)
              end do
            end do

            psi_blk(1:npt_blk, 1:nbatch) = matmul(phi_blk(1:npt_blk, 1:nbf), coef_blk(1:nbf, 1:nbatch))

            if (n_pw > 0) then
              do ipw0 = 1, n_pw, pw_block_size
                npw_blk = min(pw_block_size, n_pw - ipw0 + 1)
                do igrid = 1, npt_blk
                  ixg = ixg_buf(igrid)
                  iyg = iyg_buf(igrid)
                  izg = izg_buf(igrid)
                  rx = real(ixg - 1, 8) * dg_frag%hgs(1)
                  ry = real(iyg - 1, 8) * dg_frag%hgs(2)
                  rz = real(izg - 1, 8) * dg_frag%hgs(3)
                  do ipw_loc = 1, npw_blk
                    ipw = ipw0 + ipw_loc - 1
                    theta = dg_frag%k_pw(1, ipw) * rx + dg_frag%k_pw(2, ipw) * ry + dg_frag%k_pw(3, ipw) * rz
                    phase_blk(igrid, ipw_loc) = cmplx(cos(theta), sin(theta), kind=8) * inv_sqrt_vol
                  end do
                end do
                coef_pw_blk(1:npw_blk, 1:nbatch) = dg_frag%coef_pw_full_cache(ipw0:ipw0+npw_blk-1, io0:io0+nbatch-1, ispin)
                psi_blk(1:npt_blk, 1:nbatch) = psi_blk(1:npt_blk, 1:nbatch) + &
                  matmul(phase_blk(1:npt_blk, 1:npw_blk), coef_pw_blk(1:npw_blk, 1:nbatch))
              end do
            end if

            rho_blk(1:npt_blk) = 0.0d0
            do io = 1, nbatch
              rho_blk(1:npt_blk) = rho_blk(1:npt_blk) + occ_factor * real(conjg(psi_blk(1:npt_blk, io)) * psi_blk(1:npt_blk, io), kind=8)
            end do

            do igrid = 1, npt_blk
              owner_rank = owner_buf(igrid)
              ixg = ixg_buf(igrid)
              iyg = iyg_buf(igrid)
              izg = izg_buf(igrid)
              rho_contrib = rho_blk(igrid)
              if (owner_rank == dg_frag%id) then
                rho%f(ixg, iyg, izg) = rho%f(ixg, iyg, izg) + rho_contrib
              else if (allocated(dg_frag%density_send_slot_map)) then
                slot = dg_frag%density_send_slot_map(ix_buf(igrid), iy_buf(igrid), iz_buf(igrid), i_local)
                if (slot > 0) then
                  rho_send(owner_rank)%f(slot, 1, 1) = rho_send(owner_rank)%f(slot, 1, 1) + rho_contrib
                end if
              end if
            end do
          end do
        end do
      end do
      call cpu_time(t_project1)
      time_project = time_project + (t_project1 - t_project0)
      if (enable_density_trace .and. dg_frag%id == 0) then
        write(*,'(1x,a,a,1pe12.4)') "        density trace: stage=after-project dt=", "", time_project
        flush(6)
      end if
      end do
    end if

    if (enable_density_trace .and. dg_frag%id == 0) then
      write(*,'(1x,a)') "        density trace: stage=before-comm"
      flush(6)
    end if
    call cpu_time(t_comm0)
    nreq_recv = 0
    do irank = 0, dg_frag%isize - 1
      if (allocated(rho_recv(irank)%f)) nreq_recv = nreq_recv + 2
    end do
    if (nreq_recv > 0) then
      allocate(req_recv(nreq_recv))
      ireq = 0
      do irank = 0, dg_frag%isize - 1
        if (.not. allocated(rho_recv(irank)%f)) cycle
        ireq = ireq + 1
        req_recv(ireq) = comm_irecv(rho_recv(irank)%f, irank, rho_tag_base + dg_frag%id, dg_frag%icomm)
        ireq = ireq + 1
        req_recv(ireq) = comm_irecv(w_recv(irank)%f, irank, w_tag_base + dg_frag%id, dg_frag%icomm)
      end do
    end if

    nreq_send = 0
    do irank = 0, dg_frag%isize - 1
      if (allocated(rho_send(irank)%f)) nreq_send = nreq_send + 2
    end do
    if (nreq_send > 0) then
      allocate(req_send(nreq_send))
      ireq = 0
      do irank = 0, dg_frag%isize - 1
        if (.not. allocated(rho_send(irank)%f)) cycle
        ireq = ireq + 1
        req_send(ireq) = comm_isend(rho_send(irank)%f, irank, rho_tag_base + irank, dg_frag%icomm)
        ireq = ireq + 1
        req_send(ireq) = comm_isend(w_send(irank)%f, irank, w_tag_base + irank, dg_frag%icomm)
      end do
    end if

    if (nreq_recv > 0) then
      call comm_wait_all(req_recv)
      do irank = 0, dg_frag%isize - 1
        if (.not. allocated(rho_recv(irank)%f)) cycle
        do slot = 1, dg_frag%density_recv_map(irank)%npts
          ixg = dg_frag%density_recv_map(irank)%ixg(slot)
          iyg = dg_frag%density_recv_map(irank)%iyg(slot)
          izg = dg_frag%density_recv_map(irank)%izg(slot)
          rho%f(ixg, iyg, izg) = rho%f(ixg, iyg, izg) + rho_recv(irank)%f(slot, 1, 1)
          w_local(ixg, iyg, izg) = w_local(ixg, iyg, izg) + w_recv(irank)%f(slot, 1, 1)
        end do
        deallocate(rho_recv(irank)%f, w_recv(irank)%f)
      end do
      deallocate(req_recv)
    end if
    if (nreq_send > 0) then
      call comm_wait_all(req_send)
      do irank = 0, dg_frag%isize - 1
        if (.not. allocated(rho_send(irank)%f)) cycle
        deallocate(rho_send(irank)%f, w_send(irank)%f)
      end do
      deallocate(req_send)
    end if
    call cpu_time(t_comm1)
    time_comm = time_comm + (t_comm1 - t_comm0)
    if (enable_density_trace .and. dg_frag%id == 0) then
      write(*,'(1x,a,a,1pe12.4)') "        density trace: stage=after-comm dt=", "", time_comm
      flush(6)
    end if

    if (enable_density_trace .and. dg_frag%id == 0) then
      write(*,'(1x,a)') "        density trace: stage=before-normalize"
      flush(6)
    end if
    call cpu_time(t_norm0)
    where (w_local > 0.5d0)
      rho%f = rho%f / w_local
    end where

    total_charge_local = sum(rho%f) * system%hvol
    call comm_summation(total_charge_local, total_charge, dg_frag%icomm)
    dg_frag%elec_num_raw = total_charge
    dg_frag%rho_scale_factor = 1.0d0
    if (total_charge > 1.0d-14 .and. total_charge == total_charge) then
      scale_rho = nelec / total_charge
      dg_frag%rho_scale_factor = scale_rho
      rho%f = rho%f * scale_rho
    end if
    elec_num_scaled_local = sum(rho%f) * system%hvol
    call comm_summation(elec_num_scaled_local, dg_frag%elec_num_scaled, dg_frag%icomm)

    do ispin = 1, system%nspin
      rho_s(ispin)%f = rho%f / real(system%nspin, 8)
    end do
    call cpu_time(t_norm1)
    time_norm = time_norm + (t_norm1 - t_norm0)
    if (enable_density_trace .and. dg_frag%id == 0) then
      write(*,'(1x,a,a,1pe12.4)') "        density trace: stage=after-normalize dt=", "", time_norm
      flush(6)
    end if

    deallocate(w_local, ix_buf, iy_buf, iz_buf, owner_buf, ixg_buf, iyg_buf, izg_buf)
    deallocate(phi_blk, rho_blk, coef_blk, psi_blk)
    if (allocated(phase_blk)) deallocate(phase_blk)
    if (allocated(coef_pw_blk)) deallocate(coef_pw_blk)
    deallocate(rho_send, w_send, rho_recv, w_recv)
    call cpu_time(t_total1)
    if (enable_density_trace .and. dg_frag%id == 0) then
      write(*,'(1x,a,1pe12.4,a,1pe12.4,a,1pe12.4,a,1pe12.4,a,1pe12.4)') &
        "        density timing: total=", t_total1 - t_total0, " cache=", time_cache, &
        " project=", time_project, " comm=", time_comm, " norm=", time_norm
      flush(6)
    end if

  contains

  subroutine rebuild_density_phi_cache()
    implicit none
    integer :: ifrag_cache, i_local_cache
    integer :: nxyz_cache(3), ngrid_cache, ngrid_max
    integer :: ix_cache, iy_cache, iz_cache, igrid_cache
    integer :: ib_cache

    ngrid_max = 0
    i_local_cache = 0
    do ifrag_cache = dg_frag%ifrag_start, dg_frag%ifrag_end
      i_local_cache = i_local_cache + 1
      nxyz_cache(1:3) = dg_frag%nxyz_domain(1:3, ifrag_cache)
      ngrid_cache = nxyz_cache(1) * nxyz_cache(2) * nxyz_cache(3)
      ngrid_max = max(ngrid_max, ngrid_cache)
    end do

    if (allocated(dg_frag%density_phi_cache)) deallocate(dg_frag%density_phi_cache)
    if (ngrid_max <= 0) then
      dg_frag%density_phi_cache_valid = .false.
      return
    end if

    allocate(dg_frag%density_phi_cache(ngrid_max, dg_frag%nstate_frag, max(1, dg_frag%ifrag_end - dg_frag%ifrag_start + 1)))
    dg_frag%density_phi_cache = 0.0d0

    i_local_cache = 0
    do ifrag_cache = dg_frag%ifrag_start, dg_frag%ifrag_end
      i_local_cache = i_local_cache + 1
      nxyz_cache(1:3) = dg_frag%nxyz_domain(1:3, ifrag_cache)
      igrid_cache = 0
      do iz_cache = 1, nxyz_cache(3)
        do iy_cache = 1, nxyz_cache(2)
          do ix_cache = 1, nxyz_cache(1)
            igrid_cache = igrid_cache + 1
            do ib_cache = 1, dg_frag%nstate_frag
              dg_frag%density_phi_cache(igrid_cache, ib_cache, i_local_cache) = &
                dg_frag%phi_frag(ix_cache, iy_cache, iz_cache, ib_cache, i_local_cache)
            end do
          end do
        end do
      end do
    end do

    dg_frag%density_phi_cache_valid = .true.
  end subroutine rebuild_density_phi_cache

  integer function find_grid_owner(ixg, iyg, izg) result(owner)
    implicit none
    integer, intent(in) :: ixg, iyg, izg
    integer :: jrank

    owner = dg_frag%id
    do jrank = 0, dg_frag%isize - 1
      if (ixg < mg%is_all(1, jrank) .or. ixg > mg%ie_all(1, jrank)) cycle
      if (iyg < mg%is_all(2, jrank) .or. iyg > mg%ie_all(2, jrank)) cycle
      if (izg < mg%is_all(3, jrank) .or. izg > mg%ie_all(3, jrank)) cycle
      owner = jrank
      return
    end do
  end function find_grid_owner

  logical function local_fragments_overlap_rank_box(dg_frag_local, mg_local, rank_box) result(overlap)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag_local
    type(s_rgrid), intent(in) :: mg_local
    integer, intent(in) :: rank_box
    integer :: jfrag

    overlap = .false.
    do jfrag = dg_frag_local%ifrag_start, dg_frag_local%ifrag_end
      if (fragment_overlaps_box(dg_frag_local%ixyz_frag(:, jfrag), dg_frag_local%nxyz_domain(:, jfrag), &
          mg_local%is_all(:, rank_box), mg_local%ie_all(:, rank_box), mg_local%num)) then
        overlap = .true.
        return
      end if
    end do
  end function local_fragments_overlap_rank_box

  logical function rank_fragments_overlap_local_box(dg_frag_local, mg_local, source_rank) result(overlap)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag_local
    type(s_rgrid), intent(in) :: mg_local
    integer, intent(in) :: source_rank
    integer :: jfrag

    overlap = .false.
    do jfrag = 1, dg_frag_local%n_frag
      if (dg_frag_local%id_array(jfrag) /= source_rank) cycle
      if (fragment_overlaps_box(dg_frag_local%ixyz_frag(:, jfrag), dg_frag_local%nxyz_domain(:, jfrag), &
          mg_local%is, mg_local%ie, mg_local%num)) then
        overlap = .true.
        return
      end if
    end do
  end function rank_fragments_overlap_local_box

  logical function fragment_overlaps_box(ixyz_frag, nxyz_domain, box_is, box_ie, num_grid) result(overlap)
    implicit none
    integer, intent(in) :: ixyz_frag(3), nxyz_domain(3), box_is(3), box_ie(3), num_grid(3)
    integer :: idir

    overlap = .true.
    do idir = 1, 3
      if (.not. periodic_range_overlaps(ixyz_frag(idir), nxyz_domain(idir), box_is(idir), box_ie(idir), num_grid(idir))) then
        overlap = .false.
        return
      end if
    end do
  end function fragment_overlaps_box

  logical function periodic_range_overlaps(start_idx, length_idx, box_is, box_ie, num_grid) result(overlap)
    implicit none
    integer, intent(in) :: start_idx, length_idx, box_is, box_ie, num_grid
    integer :: current_start, remaining_len, segment_len, seg_is, seg_ie

    overlap = .false.
    if (length_idx <= 0) return
    current_start = mod(start_idx - 1, num_grid) + 1
    remaining_len = length_idx
    do while (remaining_len > 0)
      segment_len = min(remaining_len, num_grid - current_start + 1)
      seg_is = current_start
      seg_ie = current_start + segment_len - 1
      if (max(seg_is, box_is) <= min(seg_ie, box_ie)) then
        overlap = .true.
        return
      end if
      remaining_len = remaining_len - segment_len
      current_start = 1
    end do
  end function periodic_range_overlaps

  end subroutine calculate_density_from_fragments
