  subroutine calculate_density_from_fragments(dg_frag, system, mg, rho, rho_s)
    use structures
    use salmon_global, only: nelec, nelec_spin
    use communication, only: comm_summation, comm_bcast, comm_isend, comm_irecv, comm_wait_all, COMM_GROUP_NULL
    use rt_dg_fragment_ops, only: refresh_pw_coef_cache
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system),     intent(in)    :: system
    type(s_rgrid),          intent(in)    :: mg
    type(s_scalar),         intent(inout) :: rho
    type(s_scalar),         intent(inout) :: rho_s(system%nspin)

    integer, parameter :: rho_tag_base = 4101
    integer :: ifrag, ifrag_pack, io, i_local, i_local_pack, ispin
    integer :: istate_frag
    integer :: ix, iy, iz, ixg, iyg, izg, owner_rank
    integer :: ig_i, nbf, nbf_max, ipw, n_pw, n_frag, n_tot, n_basis_mix, max_mixed_basis
    integer :: nxyz(3), nxyz_pack(3), ixyz0(3), ifrag_count, ngrid_max, igrid_cache
    integer :: nocc_per_spin, nocc_spin, nocc_cache
    integer :: irank, nreq_send, nreq_recv, ireq, slot, npts, idx_local, idx_remote, idx_subgroup
    integer :: local_grid_count, remote_grid_count, valid_remote_grid_count, valid_subgroup_grid_count
    integer :: igrid0, igrid, ngrid, npt_blk, io0, nbatch, tmp_idx, ipw0, npw_blk, ipw_loc
    integer :: total_send_pts, pack_offset, subgroup_root_rank, block_idx_global, self_slot_count, total_local_pts
    integer :: nblocks_ifrag, first_block_offset, block_step_blocks, block_offset
    integer :: valid_basis_count
    integer, parameter :: grid_block_size = 512, state_block_size = 64, pw_block_size = 128
    real(8) :: occ_factor, occ_scale
    real(8) :: phi_i, rho_contrib, rho_accum, rho_mix_accum
    real(8) :: total_charge, total_charge_local, scale_rho, elec_num_scaled_local
    real(8) :: rx, ry, rz, boxL(3), inv_sqrt_vol, theta
    real(8) :: t_total0, t_total1, t_cache0, t_cache1
    real(8) :: t_project0, t_project1, t_comm0, t_comm1, t_norm0, t_norm1
    real(8) :: t_setup0, t_setup1, t_psi0, t_psi1, t_rho0, t_rho1
    real(8) :: time_cache, time_project, time_comm, time_norm
    real(8) :: time_project_setup, time_project_psi, time_project_rho
    real(8) :: time_project_grid_prep, time_project_phi_pack, time_project_overhead
    real(8) :: time_project_dmat_build
    real(8) :: t_dmat0, t_dmat1
    logical :: use_mixed_density, distribute_project
    integer, allocatable :: req_send(:), req_recv(:)
    integer, allocatable :: ix_buf(:), iy_buf(:), iz_buf(:), owner_buf(:), ixg_buf(:), iyg_buf(:), izg_buf(:)
    integer, allocatable :: slot_buf(:), local_grid_ids(:), remote_grid_ids(:), valid_remote_grid_ids(:), valid_subgroup_grid_ids(:)
    integer, allocatable :: basis_gid(:), valid_basis_ids(:)
    type(s_scalar), allocatable :: rho_send(:), rho_recv(:)
    type(s_scalar), allocatable :: rho_s_send(:,:), rho_s_recv(:,:)
    type(s_scalar), allocatable :: rho_reduce(:), rho_s_reduce(:,:)
    real(8), allocatable :: phi_blk(:,:), rho_blk(:), rho_blk_accum(:), coef_blk_re(:,:), coef_blk_im(:,:), psi_blk_re(:,:), psi_blk_im(:,:)
    real(8), allocatable :: density_mat_re(:,:), density_tmp(:,:)
    real(8), allocatable :: D_frag_re(:,:,:)   ! (nbf_max, nbf_max, nspin) pre-computed D per fragment
    real(8), allocatable :: coef_re_frag(:,:)   ! (nbf_max, nocc_spin) real coef for current fragment
    real(8), allocatable :: D_partial_re(:,:)    ! (nbf_max, nbf_max) partial D per rank
    integer :: io_s_frag, io_e_frag, nocc_loc, nocc_per_rank_loc
    complex(8), allocatable :: psi_blk(:,:), phase_cache(:,:), coef_pw_blk(:,:), coef_occ_weighted(:,:)
    complex(8), allocatable :: density_mix(:,:,:), basis_mix_blk(:,:), density_mix_tmp(:,:)
    complex(8), allocatable :: transform_frag(:,:), transform_pw(:,:)
    real(8), allocatable :: send_pack(:), send_sum(:)
    integer, allocatable :: subgroup_send_offset(:,:)
    integer, allocatable :: subgroup_self_ixg_tmp(:), subgroup_self_iyg_tmp(:), subgroup_self_izg_tmp(:)

    rho%f = 0.0d0
    do ispin = 1, system%nspin
      rho_s(ispin)%f = 0.0d0
    end do
    call cpu_time(t_total0)
    time_cache = 0.0d0
    time_project = 0.0d0
    time_comm = 0.0d0
    time_norm = 0.0d0
    time_project_setup = 0.0d0
    time_project_psi = 0.0d0
    time_project_rho = 0.0d0
    time_project_grid_prep = 0.0d0
    time_project_phi_pack = 0.0d0
    time_project_overhead = 0.0d0
    time_project_dmat_build = 0.0d0
    if (dg_frag%id == 0) then
      write(*,'(1x,a)') "        density trace: stage=entry"
      flush(6)
    end if

    if (.not. allocated(dg_frag%phi_frag)) return

    ifrag_count = dg_frag%ifrag_end - dg_frag%ifrag_start + 1
    ngrid_max = 0
    if (ifrag_count > 0) then
      do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
        ngrid_max = max(ngrid_max, product(dg_frag%nxyz_domain(:, ifrag)))
      end do
    end if
    nbf_max = max(1, maxval(dg_frag%n_basis(:, 1:system%nspin)))

    allocate(ix_buf(grid_block_size), iy_buf(grid_block_size), iz_buf(grid_block_size))
    allocate(owner_buf(grid_block_size), ixg_buf(grid_block_size), iyg_buf(grid_block_size), izg_buf(grid_block_size))
    allocate(slot_buf(grid_block_size), local_grid_ids(grid_block_size), remote_grid_ids(grid_block_size))
    allocate(valid_remote_grid_ids(grid_block_size), valid_subgroup_grid_ids(grid_block_size))
    allocate(basis_gid(dg_frag%nstate_frag), valid_basis_ids(dg_frag%nstate_frag))
    allocate(phi_blk(grid_block_size, dg_frag%nstate_frag))
    allocate(rho_blk(grid_block_size))
    allocate(rho_blk_accum(grid_block_size))
    ! Closed-shell fallback: nelec = 2 * nocc_per_spin.
    nocc_per_spin = min(dg_frag%nstate_tot, int(nelec / 2.0d0 + 1.0d-12))
    nocc_cache = nocc_per_spin
    if (system%nspin == 2 .and. sum(nelec_spin(:)) > 0) then
      nocc_cache = min(dg_frag%nstate_tot, maxval(nelec_spin(1:system%nspin)))
    end if
    allocate(coef_blk_re(nbf_max, state_block_size), coef_blk_im(nbf_max, state_block_size))
    allocate(psi_blk_re(grid_block_size, state_block_size), psi_blk_im(grid_block_size, state_block_size))
    allocate(density_mat_re(nbf_max, nbf_max), density_tmp(grid_block_size, nbf_max))
    allocate(psi_blk(grid_block_size, state_block_size))
    allocate(coef_occ_weighted(nbf_max, max(1, nocc_cache)))
    allocate(coef_pw_blk(pw_block_size, state_block_size))
    allocate(rho_send(0:dg_frag%isize-1))
    allocate(rho_recv(0:dg_frag%isize-1))
    allocate(rho_s_send(0:dg_frag%isize-1, system%nspin), rho_s_recv(0:dg_frag%isize-1, system%nspin))
    allocate(rho_reduce(0:dg_frag%isize-1))
    allocate(rho_s_reduce(0:dg_frag%isize-1, system%nspin))

    rho%f = 0.0d0
    occ_factor = 2.0d0 / real(system%nspin, 8)
    occ_scale = sqrt(occ_factor)
    n_pw = 0
    if (dg_frag%use_plane_wave_basis .and. allocated(dg_frag%coef_pw) .and. allocated(dg_frag%k_pw)) then
      n_pw = dg_frag%n_plane_waves
    end if
    n_frag = dg_frag%n_mat_max
    n_tot = n_frag + n_pw
    if (n_pw > 0) allocate(phase_cache(grid_block_size, n_pw))
    use_mixed_density = (n_pw > 0 .and. dg_frag%mixed_basis_ready .and. allocated(dg_frag%mixed_transform) .and. &
      allocated(dg_frag%coef_mix) .and. allocated(dg_frag%mixed_basis_dim))
    distribute_project = (dg_frag%icomm_frag /= COMM_GROUP_NULL .and. dg_frag%isize_frag > 1)
    subgroup_root_rank = dg_frag%id - dg_frag%id_frag
    if (.not. allocated(dg_frag%density_matrix_frag) .or. &
        size(dg_frag%density_matrix_frag, 1) /= nbf_max .or. size(dg_frag%density_matrix_frag, 2) /= nbf_max .or. &
        size(dg_frag%density_matrix_frag, 3) /= system%nspin .or. size(dg_frag%density_matrix_frag, 4) /= ifrag_count) then
      if (allocated(dg_frag%density_matrix_frag)) deallocate(dg_frag%density_matrix_frag)
      allocate(dg_frag%density_matrix_frag(nbf_max, nbf_max, system%nspin, max(1, ifrag_count)))
    end if
    if (.not. allocated(dg_frag%density_matrix_frag_valid) .or. &
        size(dg_frag%density_matrix_frag_valid, 1) /= system%nspin .or. &
        size(dg_frag%density_matrix_frag_valid, 2) /= ifrag_count) then
      if (allocated(dg_frag%density_matrix_frag_valid)) deallocate(dg_frag%density_matrix_frag_valid)
      allocate(dg_frag%density_matrix_frag_valid(system%nspin, max(1, ifrag_count)))
    end if
    dg_frag%density_matrix_frag(:, :, :, :) = (0.0d0, 0.0d0)
    dg_frag%density_matrix_frag_valid(:, :) = .false.
    max_mixed_basis = 0
    if (use_mixed_density) then
      max_mixed_basis = maxval(dg_frag%mixed_basis_dim(1:system%nspin))
      use_mixed_density = (max_mixed_basis > 0)
    end if
    if (use_mixed_density) then
      allocate(density_mix(max_mixed_basis, max_mixed_basis, system%nspin))
      density_mix(:, :, :) = (0.0d0, 0.0d0)
      allocate(basis_mix_blk(grid_block_size, max_mixed_basis))
      allocate(density_mix_tmp(grid_block_size, max_mixed_basis))
      allocate(transform_frag(dg_frag%nstate_frag, max_mixed_basis))
      if (n_pw > 0) allocate(transform_pw(n_pw, max_mixed_basis))
      do ispin = 1, system%nspin
        nocc_spin = nocc_per_spin
        if (system%nspin == 2 .and. sum(nelec_spin(:)) > 0) then
          nocc_spin = min(dg_frag%nstate_tot, nelec_spin(ispin))
        end if
        n_basis_mix = min(dg_frag%mixed_basis_dim(ispin), max_mixed_basis, size(dg_frag%coef_mix, 1))
        if (n_basis_mix <= 0 .or. nocc_spin <= 0) cycle
        if (.not. distribute_project .or. dg_frag%is_frag_root) then
          density_mix(1:n_basis_mix, 1:n_basis_mix, ispin) = occ_factor * &
            matmul(dg_frag%coef_mix(1:n_basis_mix, 1:nocc_spin, ispin), &
                   conjg(transpose(dg_frag%coef_mix(1:n_basis_mix, 1:nocc_spin, ispin))))
        else
          density_mix(1:n_basis_mix, 1:n_basis_mix, ispin) = (0.0d0, 0.0d0)
        end if
        if (distribute_project) then
          call comm_bcast(density_mix(1:n_basis_mix, 1:n_basis_mix, ispin), dg_frag%icomm_frag, 0)
        end if
      end do
    end if
    boxL(1) = dg_frag%hgs(1) * real(mg%num(1), 8)
    boxL(2) = dg_frag%hgs(2) * real(mg%num(2), 8)
    boxL(3) = dg_frag%hgs(3) * real(mg%num(3), 8)
    inv_sqrt_vol = 1.0d0 / sqrt(max(1.0d-16, boxL(1) * boxL(2) * boxL(3)))

    call cpu_time(t_setup0)
    if (distribute_project .and. allocated(dg_frag%density_owner_map) .and. .not. allocated(dg_frag%density_subgroup_send_count)) then
      allocate(dg_frag%density_subgroup_send_count(0:dg_frag%isize-1))
      allocate(dg_frag%density_subgroup_send_slot_map(size(dg_frag%density_owner_map, 1), size(dg_frag%density_owner_map, 2), &
                                                     size(dg_frag%density_owner_map, 3), size(dg_frag%density_owner_map, 4)))
      dg_frag%density_subgroup_send_count = 0
      dg_frag%density_subgroup_send_slot_map = 0
      total_local_pts = 0
      do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
        total_local_pts = total_local_pts + product(dg_frag%nxyz_domain(1:3, ifrag))
      end do
      if (total_local_pts > 0) then
        allocate(subgroup_self_ixg_tmp(total_local_pts), subgroup_self_iyg_tmp(total_local_pts), subgroup_self_izg_tmp(total_local_pts))
      end if
      self_slot_count = 0
      i_local = 0
      do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
        i_local = i_local + 1
        nxyz(1:3) = dg_frag%nxyz_domain(1:3, ifrag)
        do iz = 1, nxyz(3)
          do iy = 1, nxyz(2)
            do ix = 1, nxyz(1)
              owner_rank = dg_frag%density_owner_map(ix, iy, iz, i_local)
              dg_frag%density_subgroup_send_count(owner_rank) = dg_frag%density_subgroup_send_count(owner_rank) + 1
              dg_frag%density_subgroup_send_slot_map(ix, iy, iz, i_local) = dg_frag%density_subgroup_send_count(owner_rank)
              if (owner_rank == subgroup_root_rank) then
                self_slot_count = self_slot_count + 1
                subgroup_self_ixg_tmp(self_slot_count) = dg_frag%density_ixg_map(ix, iy, iz, i_local)
                subgroup_self_iyg_tmp(self_slot_count) = dg_frag%density_iyg_map(ix, iy, iz, i_local)
                subgroup_self_izg_tmp(self_slot_count) = dg_frag%density_izg_map(ix, iy, iz, i_local)
              end if
            end do
          end do
        end do
      end do
      if (self_slot_count > 0) then
        allocate(dg_frag%density_subgroup_self_ixg(self_slot_count), dg_frag%density_subgroup_self_iyg(self_slot_count), &
                 dg_frag%density_subgroup_self_izg(self_slot_count))
        dg_frag%density_subgroup_self_ixg(:) = subgroup_self_ixg_tmp(1:self_slot_count)
        dg_frag%density_subgroup_self_iyg(:) = subgroup_self_iyg_tmp(1:self_slot_count)
        dg_frag%density_subgroup_self_izg(:) = subgroup_self_izg_tmp(1:self_slot_count)
      end if
      if (allocated(subgroup_self_ixg_tmp)) deallocate(subgroup_self_ixg_tmp)
      if (allocated(subgroup_self_iyg_tmp)) deallocate(subgroup_self_iyg_tmp)
      if (allocated(subgroup_self_izg_tmp)) deallocate(subgroup_self_izg_tmp)
    end if
    if (.not. allocated(dg_frag%density_block_nblocks)) then
      allocate(dg_frag%density_block_nblocks(ifrag_count), dg_frag%density_block_first_offset(ifrag_count), &
               dg_frag%density_block_step(ifrag_count))
      block_idx_global = 0
      i_local = 0
      do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
        i_local = i_local + 1
        nxyz(1:3) = dg_frag%nxyz_domain(1:3, ifrag)
        ngrid = nxyz(1) * nxyz(2) * nxyz(3)
        nblocks_ifrag = (ngrid + grid_block_size - 1) / grid_block_size
        dg_frag%density_block_nblocks(i_local) = nblocks_ifrag
        dg_frag%density_block_first_offset(i_local) = 0
        dg_frag%density_block_step(i_local) = 1
        if (distribute_project) then
          dg_frag%density_block_first_offset(i_local) = modulo(dg_frag%id_frag - modulo(block_idx_global, dg_frag%isize_frag), &
                                                               dg_frag%isize_frag)
          dg_frag%density_block_step(i_local) = dg_frag%isize_frag
        end if
        block_idx_global = block_idx_global + nblocks_ifrag
      end do
    end if
    call cpu_time(t_setup1)
    time_project_overhead = time_project_overhead + (t_setup1 - t_setup0)

    do irank = 0, dg_frag%isize - 1
      if (dg_frag%is_frag_root .and. allocated(dg_frag%density_send_count)) then
        if (irank == dg_frag%id) then
          npts = 0
        else
          npts = dg_frag%density_send_count(irank)
        end if
      else
        npts = 0
      end if
      if (npts > 0) then
        allocate(rho_send(irank)%f(1:npts, 1:1, 1:1))
        rho_send(irank)%f = 0.0d0
        do ispin = 1, system%nspin
          allocate(rho_s_send(irank, ispin)%f(1:npts, 1:1, 1:1))
          rho_s_send(irank, ispin)%f = 0.0d0
        end do
      end if
      if (distribute_project .and. allocated(dg_frag%density_subgroup_send_count)) then
        npts = dg_frag%density_subgroup_send_count(irank)
      else
        npts = 0
      end if
      if (npts > 0) then
        allocate(rho_reduce(irank)%f(1:npts, 1:1, 1:1))
        rho_reduce(irank)%f = 0.0d0
        do ispin = 1, system%nspin
          allocate(rho_s_reduce(irank, ispin)%f(1:npts, 1:1, 1:1))
          rho_s_reduce(irank, ispin)%f = 0.0d0
        end do
      end if
      if (irank == dg_frag%id) cycle
      if (allocated(dg_frag%density_recv_map)) then
        npts = dg_frag%density_recv_map(irank)%npts
      else
        npts = 0
      end if
      if (npts > 0) then
        allocate(rho_recv(irank)%f(1:npts, 1:1, 1:1))
        rho_recv(irank)%f = 0.0d0
        do ispin = 1, system%nspin
          allocate(rho_s_recv(irank, ispin)%f(1:npts, 1:1, 1:1))
          rho_s_recv(irank, ispin)%f = 0.0d0
        end do
      end if
    end do
    total_send_pts = 0
    if (distribute_project .and. allocated(dg_frag%density_subgroup_send_count)) then
      do irank = 0, dg_frag%isize - 1
        total_send_pts = total_send_pts + (1 + system%nspin) * dg_frag%density_subgroup_send_count(irank)
      end do
    end if
    if (total_send_pts > 0) then
      allocate(send_pack(total_send_pts), send_sum(total_send_pts))
      send_pack(:) = 0.0d0
      allocate(subgroup_send_offset(0:dg_frag%isize-1, 0:system%nspin))
      subgroup_send_offset(:, :) = 0
      pack_offset = 0
      do irank = 0, dg_frag%isize - 1
        npts = dg_frag%density_subgroup_send_count(irank)
        subgroup_send_offset(irank, 0) = pack_offset
        if (npts > 0) then
          pack_offset = pack_offset + npts
        end if
        do ispin = 1, system%nspin
          subgroup_send_offset(irank, ispin) = pack_offset
          if (npts > 0) then
            pack_offset = pack_offset + npts
          end if
        end do
      end do
    end if

    if (n_pw > 0) then
      call cpu_time(t_cache0)
      if ((.not. allocated(dg_frag%coef_pw_full_cache)) .or. dg_frag%coef_pw_full_cache_nstate < nocc_cache) then
        call refresh_pw_coef_cache(dg_frag, nocc_cache)
      end if
      call cpu_time(t_cache1)
      time_cache = time_cache + (t_cache1 - t_cache0)
    end if
    if (dg_frag%id == 0) then
      write(*,'(1x,a,a,1pe12.4)') "        density trace: stage=after-pw-cache dt=", "", time_cache
      flush(6)
    end if

    if (dg_frag%is_frag_root .or. distribute_project) then
      if (dg_frag%id == 0) then
        write(*,'(1x,a)') "        density trace: stage=before-frag-loop"
        flush(6)
        write(*,'(1x,a,i0,a,i0,a,i0,a,i0)') "        density trace: comm-state icomm_frag=", dg_frag%icomm_frag, &
          " id_frag=", dg_frag%id_frag, " isize_frag=", dg_frag%isize_frag, " ifrag_group=", dg_frag%ifrag_group
        flush(6)
      end if
      allocate(coef_re_frag(nbf_max, max(1, nocc_cache)))
      allocate(D_partial_re(nbf_max, nbf_max))
      allocate(D_frag_re(nbf_max, nbf_max, system%nspin))
      call cpu_time(t_project0)
      i_local = 0
      block_idx_global = 0
      do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
        i_local = i_local + 1

        nxyz(1:3) = dg_frag%nxyz_domain(1:3, ifrag)
        ixyz0(1:3) = dg_frag%ixyz_frag(1:3, ifrag)
        ngrid = nxyz(1) * nxyz(2) * nxyz(3)
        nblocks_ifrag = dg_frag%density_block_nblocks(i_local)
        first_block_offset = dg_frag%density_block_first_offset(i_local)
        block_step_blocks = dg_frag%density_block_step(i_local)
        ! --- D pre-pass: compute density matrix for all spins before block loop ---
        D_frag_re(:,:,:) = 0.0d0
        if (n_pw == 0) then
          do ispin = 1, system%nspin
            nbf = dg_frag%n_basis(ifrag, ispin)
            nocc_spin = nocc_per_spin
            if (system%nspin == 2 .and. sum(nelec_spin(:)) > 0) then
              nocc_spin = min(dg_frag%nstate_tot, nelec_spin(ispin))
            end if
            if (nbf <= 0 .or. nocc_spin <= 0) cycle
            valid_basis_count = 0
            do istate_frag = 1, nbf
              basis_gid(istate_frag) = dg_frag%index_basis(istate_frag, ifrag, ispin)
              if (basis_gid(istate_frag) < 1 .or. basis_gid(istate_frag) > dg_frag%n_mat_max) cycle
              valid_basis_count = valid_basis_count + 1
              valid_basis_ids(valid_basis_count) = istate_frag
            end do
            call cpu_time(t_dmat0)
            ! Step 3a: root fills coef_re_frag, bcasts to all ranks in icomm_frag
            coef_re_frag(1:nbf_max, 1:nocc_spin) = 0.0d0
            if (.not. distribute_project .or. dg_frag%is_frag_root) then
!$omp parallel do collapse(2) private(io, idx_local, istate_frag) schedule(static)
              do io = 1, nocc_spin
                do idx_local = 1, valid_basis_count
                  istate_frag = valid_basis_ids(idx_local)
                  coef_re_frag(istate_frag, io) = real(dg_frag%coef(basis_gid(istate_frag), io, ispin), kind=8)
                end do
              end do
!$omp end parallel do
            end if
            if (distribute_project) then
              call comm_bcast(coef_re_frag(1:nbf, 1:nocc_spin), dg_frag%icomm_frag, 0)
            end if

            ! Step 3b: each rank computes dsyrk on its state slice
            nocc_per_rank_loc = (nocc_spin + dg_frag%isize_frag - 1) / dg_frag%isize_frag
            io_s_frag = dg_frag%id_frag * nocc_per_rank_loc + 1
            io_e_frag = min((dg_frag%id_frag + 1) * nocc_per_rank_loc, nocc_spin)
            nocc_loc = max(0, io_e_frag - io_s_frag + 1)

            D_partial_re(1:nbf_max, 1:nbf_max) = 0.0d0
            if (nocc_loc > 0 .and. nbf > 0) then
              ! D_partial = occ_factor * coef_re_frag[:,io_s:io_e] * coef_re_frag[:,io_s:io_e]^T
              ! upper triangle only
              ! NOTE: use 'N' (not 'T') — A is (nbf_max x nocc_loc), LDA=nbf_max
              !   'N': C = alpha * A * A^T where A is (n x k) = (nbf x nocc_loc), LDA>=n
              call dsyrk('U', 'N', nbf, nocc_loc, occ_factor, &
                         coef_re_frag(1, io_s_frag), nbf_max, &
                         0.0d0, D_partial_re(1,1), nbf_max)
            end if

            ! Step 3c: AllReduce partial D across icomm_frag
            if (distribute_project) then
              call comm_summation(D_partial_re(1:nbf_max, 1:nbf_max), &
                                  D_frag_re(1:nbf_max, 1:nbf_max, ispin), &
                                  nbf_max * nbf_max, dg_frag%icomm_frag)
            else
              D_frag_re(1:nbf_max, 1:nbf_max, ispin) = D_partial_re(1:nbf_max, 1:nbf_max)
            end if

            ! Step 3d: symmetrize (copy upper triangle to lower)
!$omp parallel do private(io, istate_frag) schedule(static)
            do io = 1, nbf
              do istate_frag = io + 1, nbf
                D_frag_re(istate_frag, io, ispin) = D_frag_re(io, istate_frag, ispin)
              end do
            end do
!$omp end parallel do
            ! n_pw=0: D_frag_re is authoritative; density_matrix_frag (complex) not updated in this path
            dg_frag%density_matrix_frag_valid(ispin, i_local) = .true.
            call cpu_time(t_dmat1)
            time_project_dmat_build = time_project_dmat_build + (t_dmat1 - t_dmat0)
          end do
        end if
        ! --- end D pre-pass ---
        do block_offset = first_block_offset, nblocks_ifrag - 1, block_step_blocks
          igrid0 = 1 + block_offset * grid_block_size
          npt_blk = min(grid_block_size, ngrid - igrid0 + 1)
          local_grid_count = 0
          remote_grid_count = 0
          valid_remote_grid_count = 0
          valid_subgroup_grid_count = 0
          call cpu_time(t_setup0)
          if (allocated(dg_frag%density_owner_map)) then
            if (distribute_project .and. allocated(dg_frag%density_subgroup_send_slot_map)) then
              call prepare_grid_buffers_owner_map(i_local, igrid0, npt_blk, nxyz, .true.)
            else if (allocated(dg_frag%density_send_slot_map)) then
              call prepare_grid_buffers_owner_map(i_local, igrid0, npt_blk, nxyz, .false.)
            else
              call prepare_grid_buffers_owner_map_no_slot(i_local, igrid0, npt_blk, nxyz)
            end if
          else
            if (distribute_project .and. allocated(dg_frag%density_subgroup_send_slot_map)) then
              call prepare_grid_buffers_runtime_map(i_local, igrid0, npt_blk, nxyz, ixyz0, .true.)
            else if (allocated(dg_frag%density_send_slot_map)) then
              call prepare_grid_buffers_runtime_map(i_local, igrid0, npt_blk, nxyz, ixyz0, .false.)
            else
              call prepare_grid_buffers_runtime_map_no_slot(igrid0, npt_blk, nxyz, ixyz0)
            end if
          end if
          if (.not. distribute_project) then
            do igrid = 1, npt_blk
              if (owner_buf(igrid) == dg_frag%id) then
                local_grid_count = local_grid_count + 1
                local_grid_ids(local_grid_count) = igrid
              else
                remote_grid_count = remote_grid_count + 1
                remote_grid_ids(remote_grid_count) = igrid
                if (slot_buf(igrid) > 0) then
                  valid_remote_grid_count = valid_remote_grid_count + 1
                  valid_remote_grid_ids(valid_remote_grid_count) = igrid
                end if
              end if
            end do
          else
            do igrid = 1, npt_blk
              if (slot_buf(igrid) > 0) then
                valid_subgroup_grid_count = valid_subgroup_grid_count + 1
                valid_subgroup_grid_ids(valid_subgroup_grid_count) = igrid
              end if
            end do
          end if
          call cpu_time(t_setup1)
          time_project_grid_prep = time_project_grid_prep + (t_setup1 - t_setup0)

          if (n_pw > 0) then
!$omp parallel do private(ixg, iyg, izg, rx, ry, rz, ipw, theta) schedule(static)
            do igrid = 1, npt_blk
              ixg = ixg_buf(igrid)
              iyg = iyg_buf(igrid)
              izg = izg_buf(igrid)
              rx = real(ixg - 1, 8) * dg_frag%hgs(1)
              ry = real(iyg - 1, 8) * dg_frag%hgs(2)
              rz = real(izg - 1, 8) * dg_frag%hgs(3)
              do ipw = 1, n_pw
                theta = dg_frag%k_pw(1, ipw) * rx + dg_frag%k_pw(2, ipw) * ry + dg_frag%k_pw(3, ipw) * rz
                phase_cache(igrid, ipw) = cmplx(cos(theta), sin(theta), kind=8) * inv_sqrt_vol
              end do
            end do
!$omp end parallel do
          end if

          do ispin = 1, system%nspin
            nocc_spin = nocc_per_spin
            if (system%nspin == 2 .and. sum(nelec_spin(:)) > 0) then
            nocc_spin = min(dg_frag%nstate_tot, nelec_spin(ispin))
            end if
            nbf = dg_frag%n_basis(ifrag, ispin)
            if (nbf <= 0 .or. nocc_spin <= 0) cycle
            valid_basis_count = 0
            do istate_frag = 1, nbf
              basis_gid(istate_frag) = dg_frag%index_basis(istate_frag, ifrag, ispin)
              if (basis_gid(istate_frag) < 1 .or. basis_gid(istate_frag) > dg_frag%n_mat_max) cycle
              valid_basis_count = valid_basis_count + 1
              valid_basis_ids(valid_basis_count) = istate_frag
            end do

            if ((.not. allocated(dg_frag%density_phi_cache)) .or. (.not. dg_frag%density_phi_cache_valid)) then
              if (.not. allocated(dg_frag%density_phi_cache)) then
                allocate(dg_frag%density_phi_cache(ngrid_max, dg_frag%nstate_frag, ifrag_count))
              end if
              dg_frag%density_phi_cache(:, :, :) = 0.0d0
              do ifrag_pack = dg_frag%ifrag_start, dg_frag%ifrag_end
                i_local_pack = ifrag_pack - dg_frag%ifrag_start + 1
                nxyz_pack(:) = dg_frag%nxyz_domain(:, ifrag_pack)
!$omp parallel do collapse(4) private(ix, iy, iz, istate_frag, igrid_cache) schedule(static)
                do iz = 1, nxyz_pack(3)
                  do iy = 1, nxyz_pack(2)
                    do ix = 1, nxyz_pack(1)
                      do istate_frag = 1, dg_frag%nstate_frag
                        igrid_cache = ix + nxyz_pack(1) * ((iy - 1) + nxyz_pack(2) * (iz - 1))
                        dg_frag%density_phi_cache(igrid_cache, istate_frag, i_local_pack) = &
                          dg_frag%phi_frag(ix, iy, iz, istate_frag, i_local_pack)
                      end do
                    end do
                  end do
                end do
!$omp end parallel do
              end do
              dg_frag%density_phi_cache_valid = .true.
            end if

          call cpu_time(t_setup0)
!$omp parallel do private(igrid) schedule(static)
            do istate_frag = 1, nbf
!$omp simd
              do igrid = 1, npt_blk
                phi_blk(igrid, istate_frag) = dg_frag%density_phi_cache(igrid0 + igrid - 1, istate_frag, i_local)
              end do
            end do
!$omp end parallel do
          call cpu_time(t_setup1)
          time_project_phi_pack = time_project_phi_pack + (t_setup1 - t_setup0)

            if (use_mixed_density) then
              n_basis_mix = min(dg_frag%mixed_basis_dim(ispin), max_mixed_basis)
              if (n_basis_mix <= 0) cycle
              call cpu_time(t_setup0)
              transform_frag(1:nbf, 1:n_basis_mix) = (0.0d0, 0.0d0)
              do istate_frag = 1, nbf
                ig_i = dg_frag%index_basis(istate_frag, ifrag, ispin)
                if (ig_i < 1 .or. ig_i > n_frag) cycle
                transform_frag(istate_frag, 1:n_basis_mix) = dg_frag%mixed_transform(ig_i, 1:n_basis_mix, ispin)
              end do
              call cpu_time(t_setup1)
              time_project_setup = time_project_setup + (t_setup1 - t_setup0)
              call cpu_time(t_psi0)
              basis_mix_blk(1:npt_blk, 1:n_basis_mix) = matmul(phi_blk(1:npt_blk, 1:nbf), transform_frag(1:nbf, 1:n_basis_mix))
              if (n_pw > 0) then
                transform_pw(1:n_pw, 1:n_basis_mix) = dg_frag%mixed_transform(n_frag+1:n_tot, 1:n_basis_mix, ispin)
                basis_mix_blk(1:npt_blk, 1:n_basis_mix) = basis_mix_blk(1:npt_blk, 1:n_basis_mix) + &
                  matmul(phase_cache(1:npt_blk, 1:n_pw), transform_pw(1:n_pw, 1:n_basis_mix))
              end if
              density_mix_tmp(1:npt_blk, 1:n_basis_mix) = matmul(basis_mix_blk(1:npt_blk, 1:n_basis_mix), &
                density_mix(1:n_basis_mix, 1:n_basis_mix, ispin))
              call cpu_time(t_psi1)
              time_project_psi = time_project_psi + (t_psi1 - t_psi0)
              call cpu_time(t_rho0)

!$omp parallel private(io, igrid, owner_rank, ixg, iyg, izg, rho_contrib, slot, pack_offset, rho_mix_accum)
!$omp do schedule(static)
              do igrid = 1, npt_blk
                rho_mix_accum = 0.0d0
!$omp simd reduction(+:rho_mix_accum)
                do io = 1, n_basis_mix
                  rho_mix_accum = rho_mix_accum + real(conjg(basis_mix_blk(igrid, io)) * density_mix_tmp(igrid, io), kind=8)
                end do
                rho_blk(igrid) = rho_mix_accum
              end do
!$omp end do

              if (distribute_project .and. allocated(dg_frag%density_subgroup_send_slot_map)) then
!$omp do schedule(static)
                do idx_subgroup = 1, valid_subgroup_grid_count
                  igrid = valid_subgroup_grid_ids(idx_subgroup)
                  owner_rank = owner_buf(igrid)
                  slot = slot_buf(igrid)
                  rho_contrib = rho_blk(igrid)
                  pack_offset = subgroup_send_offset(owner_rank, 0)
                  send_pack(pack_offset + slot) = send_pack(pack_offset + slot) + rho_contrib
                  pack_offset = subgroup_send_offset(owner_rank, ispin)
                  send_pack(pack_offset + slot) = send_pack(pack_offset + slot) + rho_contrib
                end do
!$omp end do
              else
!$omp do schedule(static)
                do idx_local = 1, local_grid_count
                  igrid = local_grid_ids(idx_local)
                  ixg = ixg_buf(igrid)
                  iyg = iyg_buf(igrid)
                  izg = izg_buf(igrid)
                  rho_contrib = rho_blk(igrid)
                  rho%f(ixg, iyg, izg) = rho%f(ixg, iyg, izg) + rho_contrib
                  rho_s(ispin)%f(ixg, iyg, izg) = rho_s(ispin)%f(ixg, iyg, izg) + rho_contrib
                end do
!$omp end do nowait
!$omp do schedule(static)
                do idx_remote = 1, valid_remote_grid_count
                  igrid = valid_remote_grid_ids(idx_remote)
                  owner_rank = owner_buf(igrid)
                  slot = slot_buf(igrid)
                  rho_contrib = rho_blk(igrid)
                  rho_send(owner_rank)%f(slot, 1, 1) = rho_send(owner_rank)%f(slot, 1, 1) + rho_contrib
                  rho_s_send(owner_rank, ispin)%f(slot, 1, 1) = rho_s_send(owner_rank, ispin)%f(slot, 1, 1) + rho_contrib
                end do
!$omp end do
              end if
!$omp end parallel
              call cpu_time(t_rho1)
              time_project_rho = time_project_rho + (t_rho1 - t_rho0)
            else
              ! D already computed in pre-pass for n_pw == 0
              rho_blk_accum(1:npt_blk) = 0.0d0
              if (n_pw == 0) then
                call cpu_time(t_psi0)
!$omp parallel do private(io, istate_frag) schedule(static)
                do io = 1, nbf
!$omp simd
                  do istate_frag = 1, nbf
                    density_mat_re(istate_frag, io) = D_frag_re(istate_frag, io, ispin)
                  end do
                end do
!$omp end parallel do
                call dgemm('N', 'N', npt_blk, nbf, nbf, 1.0d0, phi_blk, grid_block_size, &
                           density_mat_re, nbf_max, 0.0d0, density_tmp, grid_block_size)
                call cpu_time(t_psi1)
                time_project_psi = time_project_psi + (t_psi1 - t_psi0)
                call cpu_time(t_rho0)
!$omp parallel do private(io, igrid, rho_accum) schedule(static)
                do igrid = 1, npt_blk
                  rho_accum = 0.0d0
!$omp simd reduction(+:rho_accum)
                  do io = 1, nbf
                    rho_accum = rho_accum + phi_blk(igrid, io) * density_tmp(igrid, io)
                  end do
                  rho_blk_accum(igrid) = rho_blk_accum(igrid) + rho_accum
                end do
!$omp end parallel do
                call cpu_time(t_rho1)
                time_project_rho = time_project_rho + (t_rho1 - t_rho0)
              else
                do io0 = 1, nocc_spin, state_block_size
                  nbatch = min(state_block_size, nocc_spin - io0 + 1)
                  call cpu_time(t_setup0)
                  if (.not. distribute_project .or. dg_frag%is_frag_root) then
                    if (valid_basis_count < nbf) then
                      coef_blk_re(1:nbf, 1:nbatch) = 0.0d0
                      coef_blk_im(1:nbf, 1:nbatch) = 0.0d0
                    end if
!$omp parallel do collapse(2) private(io, idx_local, istate_frag) schedule(static)
                    do io = 1, nbatch
                      do idx_local = 1, valid_basis_count
                        istate_frag = valid_basis_ids(idx_local)
                        coef_blk_re(istate_frag, io) = real(dg_frag%coef(basis_gid(istate_frag), io0 + io - 1, ispin), kind=8)
                        coef_blk_im(istate_frag, io) = aimag(dg_frag%coef(basis_gid(istate_frag), io0 + io - 1, ispin))
                      end do
                    end do
!$omp end parallel do
                  end if
                  if (distribute_project) then
                    call comm_bcast(coef_blk_re(1:nbf, 1:nbatch), dg_frag%icomm_frag, 0)
                    call comm_bcast(coef_blk_im(1:nbf, 1:nbatch), dg_frag%icomm_frag, 0)
                  end if
                  call cpu_time(t_setup1)
                  time_project_setup = time_project_setup + (t_setup1 - t_setup0)

                  call cpu_time(t_psi0)
                  call dgemm('N', 'N', npt_blk, nbatch, nbf, 1.0d0, phi_blk, grid_block_size, &
                             coef_blk_re, nbf_max, 0.0d0, psi_blk_re, grid_block_size)
                  call dgemm('N', 'N', npt_blk, nbatch, nbf, 1.0d0, phi_blk, grid_block_size, &
                             coef_blk_im, nbf_max, 0.0d0, psi_blk_im, grid_block_size)

!$omp parallel do collapse(2) private(io, igrid) schedule(static)
                  do io = 1, nbatch
                    do igrid = 1, npt_blk
                      psi_blk(igrid, io) = cmplx(psi_blk_re(igrid, io), psi_blk_im(igrid, io), kind=8)
                    end do
                  end do
!$omp end parallel do
                  do ipw0 = 1, n_pw, pw_block_size
                    npw_blk = min(pw_block_size, n_pw - ipw0 + 1)
                    coef_pw_blk(1:npw_blk, 1:nbatch) = (0.0d0, 0.0d0)
                    if (.not. distribute_project .or. dg_frag%is_frag_root) then
                      coef_pw_blk(1:npw_blk, 1:nbatch) = &
                        dg_frag%coef_pw_full_cache(ipw0:ipw0+npw_blk-1, io0:io0+nbatch-1, ispin)
                    end if
                    if (distribute_project) then
                      call comm_bcast(coef_pw_blk(1:npw_blk, 1:nbatch), dg_frag%icomm_frag, 0)
                    end if
                    psi_blk(1:npt_blk, 1:nbatch) = psi_blk(1:npt_blk, 1:nbatch) + &
                      matmul(phase_cache(1:npt_blk, ipw0:ipw0+npw_blk-1), coef_pw_blk(1:npw_blk, 1:nbatch))
                  end do
                  call cpu_time(t_psi1)
                  time_project_psi = time_project_psi + (t_psi1 - t_psi0)
                  call cpu_time(t_rho0)
!$omp parallel do private(io, igrid, rho_accum) schedule(static)
                  do igrid = 1, npt_blk
                    rho_accum = 0.0d0
!$omp simd reduction(+:rho_accum)
                    do io = 1, nbatch
                      rho_accum = rho_accum + occ_factor * real(conjg(psi_blk(igrid, io)) * psi_blk(igrid, io), kind=8)
                    end do
                    rho_blk_accum(igrid) = rho_blk_accum(igrid) + rho_accum
                  end do
!$omp end parallel do
                  call cpu_time(t_rho1)
                  time_project_rho = time_project_rho + (t_rho1 - t_rho0)
                end do
              end if
              call cpu_time(t_rho0)
!$omp parallel private(igrid, owner_rank, ixg, iyg, izg, rho_contrib, slot, pack_offset)
              if (distribute_project .and. allocated(dg_frag%density_subgroup_send_slot_map)) then
!$omp do schedule(static)
                  do idx_subgroup = 1, valid_subgroup_grid_count
                    igrid = valid_subgroup_grid_ids(idx_subgroup)
                    owner_rank = owner_buf(igrid)
                    slot = slot_buf(igrid)
                    rho_contrib = rho_blk_accum(igrid)
                    pack_offset = subgroup_send_offset(owner_rank, 0)
                    send_pack(pack_offset + slot) = send_pack(pack_offset + slot) + rho_contrib
                    pack_offset = subgroup_send_offset(owner_rank, ispin)
                    send_pack(pack_offset + slot) = send_pack(pack_offset + slot) + rho_contrib
                  end do
!$omp end do
              else
!$omp do schedule(static)
                  do idx_local = 1, local_grid_count
                    igrid = local_grid_ids(idx_local)
                    ixg = ixg_buf(igrid)
                    iyg = iyg_buf(igrid)
                    izg = izg_buf(igrid)
                    rho_contrib = rho_blk_accum(igrid)
                    rho%f(ixg, iyg, izg) = rho%f(ixg, iyg, izg) + rho_contrib
                    rho_s(ispin)%f(ixg, iyg, izg) = rho_s(ispin)%f(ixg, iyg, izg) + rho_contrib
                  end do
!$omp end do nowait
!$omp do schedule(static)
                  do idx_remote = 1, valid_remote_grid_count
                    igrid = valid_remote_grid_ids(idx_remote)
                    owner_rank = owner_buf(igrid)
                    slot = slot_buf(igrid)
                    rho_contrib = rho_blk_accum(igrid)
                    rho_send(owner_rank)%f(slot, 1, 1) = rho_send(owner_rank)%f(slot, 1, 1) + rho_contrib
                    rho_s_send(owner_rank, ispin)%f(slot, 1, 1) = rho_s_send(owner_rank, ispin)%f(slot, 1, 1) + rho_contrib
                  end do
!$omp end do
              end if
!$omp end parallel
                call cpu_time(t_rho1)
                time_project_rho = time_project_rho + (t_rho1 - t_rho0)
            end if
          end do
        end do
        block_idx_global = block_idx_global + nblocks_ifrag
      end do
      call cpu_time(t_project1)
      time_project = time_project + (t_project1 - t_project0)
      if (allocated(D_frag_re)) deallocate(D_frag_re)
      if (allocated(coef_re_frag)) deallocate(coef_re_frag)
      if (allocated(D_partial_re)) deallocate(D_partial_re)
      if (dg_frag%id == 0) then
        write(*,'(1x,a,a,1pe12.4)') "        density trace: stage=after-project dt=", "", time_project
        write(*,'(1x,a,7(a,1pe12.4))') "        density trace: project breakdown", &
          " setup=", time_project_setup, " psi=", time_project_psi, " rho=", time_project_rho, &
          " grid=", time_project_grid_prep, " phi=", time_project_phi_pack, &
          " over=", time_project_overhead, " dmat=", time_project_dmat_build
        flush(6)
      end if
    end if

    if (distribute_project .and. allocated(send_sum)) then
      call comm_summation(send_pack, send_sum, total_send_pts, dg_frag%icomm_frag)
      do irank = 0, dg_frag%isize - 1
        npts = dg_frag%density_subgroup_send_count(irank)
        if (npts <= 0) cycle
        pack_offset = subgroup_send_offset(irank, 0)
        if (dg_frag%is_frag_root) then
          if (irank == subgroup_root_rank) then
!$omp parallel do private(slot, ixg, iyg, izg) schedule(static)
            do slot = 1, npts
              ixg = dg_frag%density_subgroup_self_ixg(slot)
              iyg = dg_frag%density_subgroup_self_iyg(slot)
              izg = dg_frag%density_subgroup_self_izg(slot)
              rho%f(ixg, iyg, izg) = rho%f(ixg, iyg, izg) + send_sum(pack_offset+slot)
            end do
!$omp end parallel do
          else
            rho_send(irank)%f(:, 1, 1) = send_sum(pack_offset+1:pack_offset+npts)
          end if
        end if
        if (dg_frag%is_frag_root) then
          if (irank == subgroup_root_rank) then
!$omp parallel do collapse(2) private(ispin, slot, ixg, iyg, izg, pack_offset) schedule(static)
            do ispin = 1, system%nspin
              do slot = 1, npts
                pack_offset = subgroup_send_offset(irank, ispin)
                ixg = dg_frag%density_subgroup_self_ixg(slot)
                iyg = dg_frag%density_subgroup_self_iyg(slot)
                izg = dg_frag%density_subgroup_self_izg(slot)
                rho_s(ispin)%f(ixg, iyg, izg) = rho_s(ispin)%f(ixg, iyg, izg) + send_sum(pack_offset+slot)
              end do
            end do
!$omp end parallel do
          else
            do ispin = 1, system%nspin
              pack_offset = subgroup_send_offset(irank, ispin)
              rho_s_send(irank, ispin)%f(:, 1, 1) = send_sum(pack_offset+1:pack_offset+npts)
            end do
          end if
        end if
      end do
    end if

    if (dg_frag%id == 0) then
      write(*,'(1x,a)') "        density trace: stage=before-comm"
      flush(6)
    end if
    call cpu_time(t_comm0)
    nreq_recv = 0
    do irank = 0, dg_frag%isize - 1
      if (allocated(rho_recv(irank)%f)) nreq_recv = nreq_recv + 1 + system%nspin
    end do
    if (nreq_recv > 0) then
      allocate(req_recv(nreq_recv))
      ireq = 0
      do irank = 0, dg_frag%isize - 1
        if (.not. allocated(rho_recv(irank)%f)) cycle
        ireq = ireq + 1
        req_recv(ireq) = comm_irecv(rho_recv(irank)%f, irank, rho_tag_base + dg_frag%id, dg_frag%icomm)
        do ispin = 1, system%nspin
          ireq = ireq + 1
          req_recv(ireq) = comm_irecv(rho_s_recv(irank, ispin)%f, irank, rho_tag_base + 100 * ispin + dg_frag%id, dg_frag%icomm)
        end do
      end do
    end if

    nreq_send = 0
    if (dg_frag%is_frag_root) then
      do irank = 0, dg_frag%isize - 1
        if (allocated(rho_send(irank)%f)) nreq_send = nreq_send + 1 + system%nspin
      end do
    end if
    if (dg_frag%is_frag_root .and. nreq_send > 0) then
      allocate(req_send(nreq_send))
      ireq = 0
      do irank = 0, dg_frag%isize - 1
        if (.not. allocated(rho_send(irank)%f)) cycle
        ireq = ireq + 1
        req_send(ireq) = comm_isend(rho_send(irank)%f, irank, rho_tag_base + irank, dg_frag%icomm)
        do ispin = 1, system%nspin
          ireq = ireq + 1
          req_send(ireq) = comm_isend(rho_s_send(irank, ispin)%f, irank, rho_tag_base + 100 * ispin + irank, dg_frag%icomm)
        end do
      end do
    end if

    if (nreq_recv > 0) then
      write(*,'(1x,a,i0,a)') "        density collective: rank=", dg_frag%id, " stage=before-recv-wait"
      flush(6)
      call comm_wait_all(req_recv)
      write(*,'(1x,a,i0,a)') "        density collective: rank=", dg_frag%id, " stage=after-recv-wait"
      flush(6)
      do irank = 0, dg_frag%isize - 1
        if (.not. allocated(rho_recv(irank)%f)) cycle
!$omp parallel do private(slot, ixg, iyg, izg) schedule(static)
        do slot = 1, dg_frag%density_recv_map(irank)%npts
          ixg = dg_frag%density_recv_map(irank)%ixg(slot)
          iyg = dg_frag%density_recv_map(irank)%iyg(slot)
          izg = dg_frag%density_recv_map(irank)%izg(slot)
          rho%f(ixg, iyg, izg) = rho%f(ixg, iyg, izg) + rho_recv(irank)%f(slot, 1, 1)
        end do
!$omp end parallel do
!$omp parallel do collapse(2) private(ispin, slot, ixg, iyg, izg) schedule(static)
        do ispin = 1, system%nspin
          do slot = 1, dg_frag%density_recv_map(irank)%npts
            ixg = dg_frag%density_recv_map(irank)%ixg(slot)
            iyg = dg_frag%density_recv_map(irank)%iyg(slot)
            izg = dg_frag%density_recv_map(irank)%izg(slot)
            rho_s(ispin)%f(ixg, iyg, izg) = rho_s(ispin)%f(ixg, iyg, izg) + rho_s_recv(irank, ispin)%f(slot, 1, 1)
          end do
        end do
!$omp end parallel do
        deallocate(rho_recv(irank)%f)
        do ispin = 1, system%nspin
          if (allocated(rho_s_recv(irank, ispin)%f)) deallocate(rho_s_recv(irank, ispin)%f)
        end do
      end do
      deallocate(req_recv)
    end if
    if (dg_frag%is_frag_root .and. nreq_send > 0) then
      write(*,'(1x,a,i0,a)') "        density collective: rank=", dg_frag%id, " stage=before-send-wait"
      flush(6)
      call comm_wait_all(req_send)
      write(*,'(1x,a,i0,a)') "        density collective: rank=", dg_frag%id, " stage=after-send-wait"
      flush(6)
      do irank = 0, dg_frag%isize - 1
        if (.not. allocated(rho_send(irank)%f)) cycle
        deallocate(rho_send(irank)%f)
        do ispin = 1, system%nspin
          if (allocated(rho_s_send(irank, ispin)%f)) deallocate(rho_s_send(irank, ispin)%f)
        end do
      end do
      deallocate(req_send)
    end if
    do irank = 0, dg_frag%isize - 1
      if (allocated(rho_reduce(irank)%f)) deallocate(rho_reduce(irank)%f)
      do ispin = 1, system%nspin
        if (allocated(rho_s_reduce(irank, ispin)%f)) deallocate(rho_s_reduce(irank, ispin)%f)
      end do
      if (allocated(rho_send(irank)%f)) deallocate(rho_send(irank)%f)
      do ispin = 1, system%nspin
        if (allocated(rho_s_send(irank, ispin)%f)) deallocate(rho_s_send(irank, ispin)%f)
      end do
    end do
    if (allocated(subgroup_send_offset)) deallocate(subgroup_send_offset)
    call cpu_time(t_comm1)
    time_comm = time_comm + (t_comm1 - t_comm0)
    if (dg_frag%id == 0) then
      write(*,'(1x,a,a,1pe12.4)') "        density trace: stage=after-comm dt=", "", time_comm
      flush(6)
    end if

    if (dg_frag%id == 0) then
      write(*,'(1x,a)') "        density trace: stage=before-normalize"
      flush(6)
    end if
    call cpu_time(t_norm0)
    if (allocated(dg_frag%density_inv_weight_local)) then
      rho%f = rho%f * dg_frag%density_inv_weight_local
      do ispin = 1, system%nspin
        rho_s(ispin)%f = rho_s(ispin)%f * dg_frag%density_inv_weight_local
      end do
    end if

    total_charge_local = sum(rho%f) * system%hvol
    write(*,'(1x,a,i0,a)') "        density collective: rank=", dg_frag%id, " stage=before-total-charge-sum"
    flush(6)
    call comm_summation(total_charge_local, total_charge, dg_frag%icomm)
    write(*,'(1x,a,i0,a,1pe12.4)') "        density collective: rank=", dg_frag%id, &
      " stage=after-total-charge-sum total_charge=", total_charge
    flush(6)
    dg_frag%elec_num_raw = total_charge
    dg_frag%rho_scale_factor = 1.0d0
    if (total_charge > 1.0d-14 .and. total_charge == total_charge) then
      scale_rho = nelec / total_charge
      dg_frag%rho_scale_factor = scale_rho
      rho%f = rho%f * scale_rho
      do ispin = 1, system%nspin
        rho_s(ispin)%f = rho_s(ispin)%f * scale_rho
      end do
    end if
    elec_num_scaled_local = sum(rho%f) * system%hvol
    write(*,'(1x,a,i0,a)') "        density collective: rank=", dg_frag%id, " stage=before-scaled-charge-sum"
    flush(6)
    call comm_summation(elec_num_scaled_local, dg_frag%elec_num_scaled, dg_frag%icomm)
    write(*,'(1x,a,i0,a,1pe12.4)') "        density collective: rank=", dg_frag%id, &
      " stage=after-scaled-charge-sum elec_num_scaled=", dg_frag%elec_num_scaled
    flush(6)

    call cpu_time(t_norm1)
    time_norm = time_norm + (t_norm1 - t_norm0)
    if (dg_frag%id == 0) then
      write(*,'(1x,a,a,1pe12.4)') "        density trace: stage=after-normalize dt=", "", time_norm
      flush(6)
    end if

    deallocate(ix_buf, iy_buf, iz_buf, owner_buf, ixg_buf, iyg_buf, izg_buf)
    deallocate(slot_buf, local_grid_ids, remote_grid_ids, valid_remote_grid_ids, valid_subgroup_grid_ids)
    deallocate(basis_gid, valid_basis_ids)
    deallocate(phi_blk, rho_blk, rho_blk_accum, coef_blk_re, coef_blk_im, psi_blk_re, psi_blk_im, density_mat_re, density_tmp, &
      psi_blk, coef_occ_weighted)
    if (allocated(phase_cache)) deallocate(phase_cache)
    if (allocated(density_mix)) deallocate(density_mix)
    if (allocated(basis_mix_blk)) deallocate(basis_mix_blk)
    if (allocated(density_mix_tmp)) deallocate(density_mix_tmp)
    if (allocated(transform_frag)) deallocate(transform_frag)
    if (allocated(transform_pw)) deallocate(transform_pw)
    if (allocated(send_pack)) deallocate(send_pack)
    if (allocated(send_sum)) deallocate(send_sum)
    deallocate(rho_send, rho_recv, rho_s_send, rho_s_recv)
    call cpu_time(t_total1)
    if (dg_frag%id == 0) then
      write(*,'(1x,a,1pe12.4,a,1pe12.4,a,1pe12.4,a,1pe12.4,a,1pe12.4)') &
        "        density timing: total=", t_total1 - t_total0, " cache=", time_cache, &
        " project=", time_project, " comm=", time_comm, " norm=", time_norm
      write(*,'(1x,a,6(a,1pe12.4))') "        density timing detail: project", &
        " setup=", time_project_setup, " psi=", time_project_psi, " rho=", time_project_rho, &
        " grid=", time_project_grid_prep, " phi=", time_project_phi_pack, " over=", time_project_overhead
      flush(6)
    end if

  contains

    subroutine prepare_grid_buffers_owner_map(i_local_grid, igrid0_grid, npt_blk_grid, nxyz_grid, use_subgroup_slot)
      implicit none
      integer, intent(in) :: i_local_grid, igrid0_grid, npt_blk_grid
      integer, intent(in) :: nxyz_grid(3)
      logical, intent(in) :: use_subgroup_slot

!$omp parallel do private(igrid, tmp_idx, ix, iy, iz, ixg, iyg, izg, owner_rank) schedule(static)
      do igrid = 1, npt_blk_grid
        tmp_idx = igrid0_grid + igrid - 2
        ix = mod(tmp_idx, nxyz_grid(1)) + 1
        tmp_idx = tmp_idx / nxyz_grid(1)
        iy = mod(tmp_idx, nxyz_grid(2)) + 1
        iz = tmp_idx / nxyz_grid(2) + 1
        ix_buf(igrid) = ix
        iy_buf(igrid) = iy
        iz_buf(igrid) = iz
        ixg = dg_frag%density_ixg_map(ix, iy, iz, i_local_grid)
        iyg = dg_frag%density_iyg_map(ix, iy, iz, i_local_grid)
        izg = dg_frag%density_izg_map(ix, iy, iz, i_local_grid)
        owner_rank = dg_frag%density_owner_map(ix, iy, iz, i_local_grid)
        ixg_buf(igrid) = ixg
        iyg_buf(igrid) = iyg
        izg_buf(igrid) = izg
        owner_buf(igrid) = owner_rank
        if (use_subgroup_slot) then
          slot_buf(igrid) = dg_frag%density_subgroup_send_slot_map(ix, iy, iz, i_local_grid)
        else
          slot_buf(igrid) = dg_frag%density_send_slot_map(ix, iy, iz, i_local_grid)
        end if
      end do
!$omp end parallel do
    end subroutine prepare_grid_buffers_owner_map

    subroutine prepare_grid_buffers_owner_map_no_slot(i_local_grid, igrid0_grid, npt_blk_grid, nxyz_grid)
      implicit none
      integer, intent(in) :: i_local_grid, igrid0_grid, npt_blk_grid
      integer, intent(in) :: nxyz_grid(3)

!$omp parallel do private(igrid, tmp_idx, ix, iy, iz, ixg, iyg, izg, owner_rank) schedule(static)
      do igrid = 1, npt_blk_grid
        tmp_idx = igrid0_grid + igrid - 2
        ix = mod(tmp_idx, nxyz_grid(1)) + 1
        tmp_idx = tmp_idx / nxyz_grid(1)
        iy = mod(tmp_idx, nxyz_grid(2)) + 1
        iz = tmp_idx / nxyz_grid(2) + 1
        ix_buf(igrid) = ix
        iy_buf(igrid) = iy
        iz_buf(igrid) = iz
        ixg = dg_frag%density_ixg_map(ix, iy, iz, i_local_grid)
        iyg = dg_frag%density_iyg_map(ix, iy, iz, i_local_grid)
        izg = dg_frag%density_izg_map(ix, iy, iz, i_local_grid)
        owner_rank = dg_frag%density_owner_map(ix, iy, iz, i_local_grid)
        ixg_buf(igrid) = ixg
        iyg_buf(igrid) = iyg
        izg_buf(igrid) = izg
        owner_buf(igrid) = owner_rank
        slot_buf(igrid) = 0
      end do
!$omp end parallel do
    end subroutine prepare_grid_buffers_owner_map_no_slot

    subroutine prepare_grid_buffers_runtime_map(i_local_grid, igrid0_grid, npt_blk_grid, nxyz_grid, ixyz0_grid, use_subgroup_slot)
      implicit none
      integer, intent(in) :: i_local_grid, igrid0_grid, npt_blk_grid
      integer, intent(in) :: nxyz_grid(3), ixyz0_grid(3)
      logical, intent(in) :: use_subgroup_slot

!$omp parallel do private(igrid, tmp_idx, ix, iy, iz, ixg, iyg, izg, owner_rank) schedule(static)
      do igrid = 1, npt_blk_grid
        tmp_idx = igrid0_grid + igrid - 2
        ix = mod(tmp_idx, nxyz_grid(1)) + 1
        tmp_idx = tmp_idx / nxyz_grid(1)
        iy = mod(tmp_idx, nxyz_grid(2)) + 1
        iz = tmp_idx / nxyz_grid(2) + 1
        ix_buf(igrid) = ix
        iy_buf(igrid) = iy
        iz_buf(igrid) = iz
        ixg = mod(ixyz0_grid(1) + ix - 2, mg%num(1)) + 1
        iyg = mod(ixyz0_grid(2) + iy - 2, mg%num(2)) + 1
        izg = mod(ixyz0_grid(3) + iz - 2, mg%num(3)) + 1
        owner_rank = find_grid_owner(ixg, iyg, izg)
        ixg_buf(igrid) = ixg
        iyg_buf(igrid) = iyg
        izg_buf(igrid) = izg
        owner_buf(igrid) = owner_rank
        if (use_subgroup_slot) then
          slot_buf(igrid) = dg_frag%density_subgroup_send_slot_map(ix, iy, iz, i_local_grid)
        else
          slot_buf(igrid) = dg_frag%density_send_slot_map(ix, iy, iz, i_local_grid)
        end if
      end do
!$omp end parallel do
    end subroutine prepare_grid_buffers_runtime_map

    subroutine prepare_grid_buffers_runtime_map_no_slot(igrid0_grid, npt_blk_grid, nxyz_grid, ixyz0_grid)
      implicit none
      integer, intent(in) :: igrid0_grid, npt_blk_grid
      integer, intent(in) :: nxyz_grid(3), ixyz0_grid(3)

!$omp parallel do private(igrid, tmp_idx, ix, iy, iz, ixg, iyg, izg, owner_rank) schedule(static)
      do igrid = 1, npt_blk_grid
        tmp_idx = igrid0_grid + igrid - 2
        ix = mod(tmp_idx, nxyz_grid(1)) + 1
        tmp_idx = tmp_idx / nxyz_grid(1)
        iy = mod(tmp_idx, nxyz_grid(2)) + 1
        iz = tmp_idx / nxyz_grid(2) + 1
        ix_buf(igrid) = ix
        iy_buf(igrid) = iy
        iz_buf(igrid) = iz
        ixg = mod(ixyz0_grid(1) + ix - 2, mg%num(1)) + 1
        iyg = mod(ixyz0_grid(2) + iy - 2, mg%num(2)) + 1
        izg = mod(ixyz0_grid(3) + iz - 2, mg%num(3)) + 1
        owner_rank = find_grid_owner(ixg, iyg, izg)
        ixg_buf(igrid) = ixg
        iyg_buf(igrid) = iyg
        izg_buf(igrid) = izg
        owner_buf(igrid) = owner_rank
        slot_buf(igrid) = 0
      end do
!$omp end parallel do
    end subroutine prepare_grid_buffers_runtime_map_no_slot

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
