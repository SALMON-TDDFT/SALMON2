  subroutine calculate_density_from_fragments(dg_frag, system, mg, rho, rho_s)
    use structures
    use salmon_global, only: nelec, nelec_spin
    use communication, only: comm_summation, comm_bcast, comm_alltoallv, COMM_GROUP_NULL
    use rt_dg_fragment_ops, only: refresh_pw_coef_cache
    use rt_dg_fragment_types, only: density_grid_point_info
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system),     intent(in)    :: system
    type(s_rgrid),          intent(in)    :: mg
    type(s_scalar),         intent(inout) :: rho
    type(s_scalar),         intent(inout) :: rho_s(system%nspin)

    integer :: ifrag, ifrag_pack, io, i_local, i_local_pack, ispin
    integer :: istate_frag
    integer :: ix, iy, iz, ixg, iyg, izg, bx, by, bz, owner_rank
    integer :: ig_i, nbf, nbf_max, ipw, n_pw, n_frag, n_tot, n_basis_mix, max_mixed_basis
    integer :: nxyz(3), nxyz_pack(3), ifrag_count, ngrid_max, igrid_cache
    integer :: nocc_per_spin, nocc_spin, nocc_cache
    integer :: irank, slot, npts, idx_local, idx_remote
    integer :: local_grid_count, remote_grid_count, valid_remote_grid_count
    integer :: igrid0, igrid, ngrid, npt_blk, io0, nbatch, tmp_idx, ipw0, npw_blk, ipw_loc
    integer :: total_send_pts, subgroup_root_rank, block_idx_global, self_slot_count, total_local_pts
    integer :: nblocks_ifrag, first_block_offset, block_step_blocks, block_offset
    integer :: valid_basis_count
    integer :: handler_id_frag
    integer :: packed_npts, spin_offset
    integer, parameter :: grid_block_size = 512, state_block_size = 64, pw_block_size = 128
    real(8) :: occ_factor, occ_scale
    real(8) :: phi_i, rho_contrib, rho_raw_contrib, rho_accum, rho_mix_accum
    real(8) :: total_charge, total_charge_local, scale_rho, rho_sum_local
    real(8) :: charge_budget_local(6), charge_budget_total(6)
    real(8) :: rx, ry, rz, boxL(3), inv_sqrt_vol, theta
    real(8) :: t_total0, t_total1, t_cache0, t_cache1
    real(8) :: t_project0, t_project1, t_comm0, t_comm1, t_norm0, t_norm1
    real(8) :: t_setup0, t_setup1, t_psi0, t_psi1, t_rho0, t_rho1
    real(8) :: time_cache, time_project, time_comm, time_norm
    real(8) :: time_comm_subgroup_reduce, time_comm_pack, time_comm_exchange, time_comm_unpack
    real(8) :: time_project_setup, time_project_psi, time_project_rho
    real(8) :: time_project_grid_prep, time_project_phi_pack, time_project_overhead
    real(8) :: time_project_dmat_build
    real(8) :: t_dmat0, t_dmat1
    real(8) :: time_cache_pw_refresh, time_cache_phi_block_refresh
    logical :: use_mixed_density
    logical :: rebuilt_pw_cache, rebuilt_phi_block_cache
    logical :: need_pw_cache_alloc, need_pw_cache_expand
    logical :: need_phi_cache_alloc, need_phi_count_alloc, need_phi_cache_invalid, need_phi_cache_resize
    integer, allocatable :: ix_buf(:), iy_buf(:), iz_buf(:), owner_buf(:), ixg_buf(:), iyg_buf(:), izg_buf(:)
    integer, allocatable :: slot_buf(:), local_grid_ids(:), remote_grid_ids(:), valid_remote_grid_ids(:)
    integer, allocatable :: basis_gid(:), valid_basis_ids(:)
    type(s_scalar), allocatable :: rho_send(:), rho_recv(:)
    integer, allocatable :: send_counts(:), recv_counts(:), send_displs(:), recv_displs(:)
    real(8), allocatable :: send_flat(:), recv_flat(:)
    real(8), allocatable :: rho_bf(:,:,:), rho_s_bf(:,:,:,:)
    real(8), allocatable :: phi_blk(:,:), rho_blk(:), rho_blk_accum(:), coef_blk_re(:,:), coef_blk_im(:,:), psi_blk_re(:,:), psi_blk_im(:,:)
    real(8), allocatable :: density_mat_re(:,:), density_tmp(:,:)
    real(8), allocatable :: D_frag_re(:,:,:)   ! (nbf_max, nbf_max, nspin) pre-computed D per fragment
    real(8), allocatable :: coef_re_frag(:,:)   ! (nbf_max, nocc_spin) real coef for current fragment
    real(8), allocatable :: D_partial_re(:,:)    ! (nbf_max, nbf_max) partial D per rank
    real(8), allocatable :: coef_re_full(:,:,:)  ! (nbf_max, nocc_cache, nspin) upfront bcast coef (n_pw>0)
    real(8), allocatable :: coef_im_full(:,:,:)  ! (nbf_max, nocc_cache, nspin)
    real(8), allocatable :: rho_blk_partial(:)   ! (grid_block_size) partial rho for state slice
    real(8) :: time_project_rho_reduce, time_project_phi_block_build
    integer :: io_s_frag, io_e_frag, nocc_loc, nocc_per_rank_loc
    integer :: nblocks_max, block_cache_idx, npt_cache, rem_xy
    integer :: ibuf_x, ibuf_y, ibuf_z
    integer :: rho_x_lo, rho_x_hi, rho_y_lo, rho_y_hi, rho_z_lo, rho_z_hi
    integer :: rho_s_x_lo, rho_s_x_hi, rho_s_y_lo, rho_s_y_hi, rho_s_z_lo, rho_s_z_hi
    complex(8), allocatable :: psi_blk(:,:), phase_cache(:,:), coef_pw_blk(:,:), coef_occ_weighted(:,:)
    complex(8), allocatable :: density_mix(:,:,:), basis_mix_blk(:,:), density_mix_tmp(:,:)
    complex(8), allocatable :: transform_frag(:,:), transform_pw(:,:)
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
    time_comm_subgroup_reduce = 0.0d0
    time_comm_pack = 0.0d0
    time_comm_exchange = 0.0d0
    time_comm_unpack = 0.0d0
    time_project_setup = 0.0d0
    time_project_psi = 0.0d0
    time_project_rho = 0.0d0
    time_project_grid_prep = 0.0d0
    time_project_phi_pack = 0.0d0
    time_project_overhead = 0.0d0
    time_project_dmat_build = 0.0d0
    time_project_rho_reduce = 0.0d0
    time_project_phi_block_build = 0.0d0
    time_cache_pw_refresh = 0.0d0
    time_cache_phi_block_refresh = 0.0d0
    charge_budget_local(:) = 0.0d0
    charge_budget_total(:) = 0.0d0
    rebuilt_pw_cache = .false.
    rebuilt_phi_block_cache = .false.
    need_pw_cache_alloc = .false.
    need_pw_cache_expand = .false.
    need_phi_cache_alloc = .false.
    need_phi_count_alloc = .false.
    need_phi_cache_invalid = .false.
    need_phi_cache_resize = .false.
    if (dg_frag%is_frag_root) then
      write(*,'(1x,a,i0,a,i0,a,i0)') "        density trace: rank=", dg_frag%id, &
        " id_frag=", dg_frag%id_frag, " stage=entry_root root_world=", dg_frag%id - dg_frag%id_frag
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
    allocate(valid_remote_grid_ids(grid_block_size))
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
    ibuf_x = dg_frag%nxyz_buffer(1)
    ibuf_y = dg_frag%nxyz_buffer(2)
    ibuf_z = dg_frag%nxyz_buffer(3)
    rho_x_lo = mg%is(1)
    rho_x_hi = mg%ie(1)
    rho_y_lo = mg%is(2)
    rho_y_hi = mg%ie(2)
    rho_z_lo = mg%is(3)
    rho_z_hi = mg%ie(3)
    rho_s_x_lo = rho_x_lo
    rho_s_x_hi = rho_x_hi
    rho_s_y_lo = rho_y_lo
    rho_s_y_hi = rho_y_hi
    rho_s_z_lo = rho_z_lo
    rho_s_z_hi = rho_z_hi
    allocate(rho_bf(rho_x_lo-ibuf_x:rho_x_hi+ibuf_x, &
                    rho_y_lo-ibuf_y:rho_y_hi+ibuf_y, &
                    rho_z_lo-ibuf_z:rho_z_hi+ibuf_z))
    allocate(rho_s_bf(rho_s_x_lo:rho_s_x_hi, &
                      rho_s_y_lo:rho_s_y_hi, &
                      rho_s_z_lo:rho_s_z_hi, &
      system%nspin))
    allocate(rho_send(0:dg_frag%isize-1))
    allocate(rho_recv(0:dg_frag%isize-1))

    rho%f = 0.0d0
    rho_bf(:, :, :) = 0.0d0
    rho_s_bf(:, :, :, :) = 0.0d0
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
    ! Density reconstruction uses subgroup-distributed projection and collective reductions on icomm_frag.
    subgroup_root_rank = dg_frag%id - dg_frag%id_frag
    total_send_pts = 0
    do irank = 0, dg_frag%isize - 1
      npts = dg_frag%density_send_count(irank)
      if (npts <= 0) cycle
      allocate(rho_send(irank)%f((system%nspin + 1) * npts, 1, 1))
      rho_send(irank)%f(:, :, :) = 0.0d0
    end do
    do irank = 0, dg_frag%isize - 1
      npts = dg_frag%density_recv_map(irank)%npts
      if (npts <= 0) cycle
      allocate(rho_recv(irank)%f((system%nspin + 1) * npts, 1, 1))
      rho_recv(irank)%f(:, :, :) = 0.0d0
    end do
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
        density_mix(1:n_basis_mix, 1:n_basis_mix, ispin) = occ_factor * &
          matmul(dg_frag%coef_mix(1:n_basis_mix, 1:nocc_spin, ispin), &
                 conjg(transpose(dg_frag%coef_mix(1:n_basis_mix, 1:nocc_spin, ispin))))
      end do
    end if
    boxL(1) = dg_frag%hgs(1) * real(mg%num(1), 8)
    boxL(2) = dg_frag%hgs(2) * real(mg%num(2), 8)
    boxL(3) = dg_frag%hgs(3) * real(mg%num(3), 8)
    inv_sqrt_vol = 1.0d0 / sqrt(max(1.0d-16, boxL(1) * boxL(2) * boxL(3)))

    call cpu_time(t_setup0)
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
        block_idx_global = block_idx_global + nblocks_ifrag
      end do
    end if
    call cpu_time(t_setup1)
    time_project_overhead = time_project_overhead + (t_setup1 - t_setup0)

    if (n_pw > 0) then
      call cpu_time(t_cache0)
      need_pw_cache_alloc = (.not. allocated(dg_frag%coef_pw_full_cache))
      need_pw_cache_expand = (.not. need_pw_cache_alloc) .and. dg_frag%coef_pw_full_cache_nstate < nocc_cache
      if (need_pw_cache_alloc .or. need_pw_cache_expand) then
        call refresh_pw_coef_cache(dg_frag, nocc_cache)
        rebuilt_pw_cache = .true.
      end if
      call cpu_time(t_cache1)
      time_cache = time_cache + (t_cache1 - t_cache0)
      if (rebuilt_pw_cache) time_cache_pw_refresh = time_cache_pw_refresh + (t_cache1 - t_cache0)
    end if
    if (dg_frag%id == 0) then
      write(*,'(1x,a,a,1pe12.4)') "        density trace: stage=after-pw-cache dt=", "", time_cache
      write(*,'(1x,a,l1,a,l1,a,1pe12.4)') "        density cache: pw_refresh rebuilt=", rebuilt_pw_cache, &
        " expand_only=", need_pw_cache_expand, " dt=", time_cache_pw_refresh
      flush(6)
    end if

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
    call cpu_time(t_setup0)
    need_phi_cache_alloc = (.not. allocated(dg_frag%density_phi_block_cache))
    need_phi_count_alloc = (.not. allocated(dg_frag%density_phi_block_count))
    need_phi_cache_invalid = (.not. dg_frag%density_phi_block_cache_valid)
    need_phi_cache_resize = dg_frag%density_phi_block_size /= grid_block_size
    if (need_phi_cache_alloc .or. need_phi_count_alloc .or. need_phi_cache_invalid .or. need_phi_cache_resize) then
      if (allocated(dg_frag%density_phi_block_cache)) deallocate(dg_frag%density_phi_block_cache)
      if (allocated(dg_frag%density_phi_block_count)) deallocate(dg_frag%density_phi_block_count)
      nblocks_max = max(1, maxval(dg_frag%density_block_nblocks))
      allocate(dg_frag%density_phi_block_cache(grid_block_size, dg_frag%nstate_frag, nblocks_max, ifrag_count))
      allocate(dg_frag%density_phi_block_count(ifrag_count))
      dg_frag%density_phi_block_cache(:, :, :, :) = 0.0d0
      dg_frag%density_phi_block_count(:) = 0
      i_local = 0
      do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
        i_local = i_local + 1
        local_grid_count = dg_frag%density_grid_point_count(i_local)
        dg_frag%density_phi_block_count(i_local) = dg_frag%density_block_nblocks(i_local)
        do block_cache_idx = 1, dg_frag%density_phi_block_count(i_local)
          igrid0 = 1 + (block_cache_idx - 1) * grid_block_size
          npt_cache = min(grid_block_size, local_grid_count - igrid0 + 1)
!$omp parallel do private(istate_frag, igrid, ixg, iyg, izg, bx, by, bz) schedule(static)
          do istate_frag = 1, dg_frag%nstate_frag
            do igrid = 1, npt_cache
              ixg = dg_frag%density_grid_points(igrid0 + igrid - 1, i_local)%ixg
              iyg = dg_frag%density_grid_points(igrid0 + igrid - 1, i_local)%iyg
              izg = dg_frag%density_grid_points(igrid0 + igrid - 1, i_local)%izg
              bx = map_global_to_phi_box_coord_ham(ixg, lbound(dg_frag%phi_frag, 1), ubound(dg_frag%phi_frag, 1), dg_frag%lgnum_total(1))
              by = map_global_to_phi_box_coord_ham(iyg, lbound(dg_frag%phi_frag, 2), ubound(dg_frag%phi_frag, 2), dg_frag%lgnum_total(2))
              bz = map_global_to_phi_box_coord_ham(izg, lbound(dg_frag%phi_frag, 3), ubound(dg_frag%phi_frag, 3), dg_frag%lgnum_total(3))
              if (bx == 0 .or. by == 0 .or. bz == 0) cycle
              dg_frag%density_phi_block_cache(igrid, istate_frag, block_cache_idx, i_local) = &
                dg_frag%phi_frag(bx, by, bz, istate_frag, i_local)
            end do
          end do
!$omp end parallel do
        end do
      end do
      dg_frag%density_phi_block_size = grid_block_size
      dg_frag%density_phi_block_cache_valid = .true.
      rebuilt_phi_block_cache = .true.
    end if
    call cpu_time(t_setup1)
    time_project_phi_block_build = time_project_phi_block_build + (t_setup1 - t_setup0)
    if (rebuilt_phi_block_cache) time_cache_phi_block_refresh = time_cache_phi_block_refresh + (t_setup1 - t_setup0)
    if (dg_frag%id == 0) then
      write(*,'(1x,a,l1,a,l1,a,l1,a,l1,a,l1,a,1pe12.4)') "        density cache: phi_block rebuilt=", &
        rebuilt_phi_block_cache, " alloc=", need_phi_cache_alloc, " count_alloc=", need_phi_count_alloc, &
        " invalid=", need_phi_cache_invalid, " resize=", need_phi_cache_resize, " dt=", time_cache_phi_block_refresh
      flush(6)
    end if
    call cpu_time(t_project0)
      i_local = 0
      block_idx_global = 0
      do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
        i_local = i_local + 1

        nxyz(1:3) = dg_frag%nxyz_domain(1:3, ifrag)
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
            if (dg_frag%is_frag_root) then
!$omp parallel do collapse(2) private(io, idx_local, istate_frag) schedule(static)
              do io = 1, nocc_spin
                do idx_local = 1, valid_basis_count
                  istate_frag = valid_basis_ids(idx_local)
                  coef_re_frag(istate_frag, io) = real(dg_frag%coef(basis_gid(istate_frag), io, ispin), kind=8)
                end do
              end do
!$omp end parallel do
            end if
            ! bcast coef_re_frag from root to all ranks in icomm_frag
            if (dg_frag%isize_frag > 1 .and. dg_frag%icomm_frag /= COMM_GROUP_NULL) then
              call comm_bcast(coef_re_frag(1:nbf_max, 1:nocc_spin), dg_frag%icomm_frag, 0)
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
            if (dg_frag%isize_frag > 1 .and. dg_frag%icomm_frag /= COMM_GROUP_NULL) then
              call comm_summation(D_partial_re(1:nbf_max, 1:nbf_max), D_frag_re(1:nbf_max, 1:nbf_max, ispin), &
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
        else
          ! n_pw > 0: upfront bcast of coef_re/im to all icomm_frag ranks
          if (.not. allocated(coef_re_full)) then
            allocate(coef_re_full(nbf_max, max(1, nocc_cache), system%nspin))
            allocate(coef_im_full(nbf_max, max(1, nocc_cache), system%nspin))
          end if
          coef_re_full(1:nbf_max, 1:nocc_cache, 1:system%nspin) = 0.0d0
          coef_im_full(1:nbf_max, 1:nocc_cache, 1:system%nspin) = 0.0d0
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
            if (dg_frag%is_frag_root) then
!$omp parallel do collapse(2) private(io, idx_local, istate_frag) schedule(static)
              do io = 1, nocc_spin
                do idx_local = 1, valid_basis_count
                  istate_frag = valid_basis_ids(idx_local)
                  coef_re_full(istate_frag, io, ispin) = real(dg_frag%coef(basis_gid(istate_frag), io, ispin), kind=8)
                  coef_im_full(istate_frag, io, ispin) = aimag(dg_frag%coef(basis_gid(istate_frag), io, ispin))
                end do
              end do
!$omp end parallel do
            end if
          end do
          ! bcast coef_re/im_full from root to all ranks in icomm_frag
          if (dg_frag%isize_frag > 1 .and. dg_frag%icomm_frag /= COMM_GROUP_NULL) then
            call comm_bcast(coef_re_full(1:nbf_max, 1:nocc_cache, 1:system%nspin), dg_frag%icomm_frag, 0)
            call comm_bcast(coef_im_full(1:nbf_max, 1:nocc_cache, 1:system%nspin), dg_frag%icomm_frag, 0)
          end if
          ! NOTE: coef_pw_full_cache is already populated on all ranks by refresh_pw_coef_cache
          ! (called before the fragment loop), so no bcast is needed here.

          ! Compute state range for this rank
          nocc_per_rank_loc = (nocc_per_spin + dg_frag%isize_frag - 1) / dg_frag%isize_frag
          io_s_frag = dg_frag%id_frag * nocc_per_rank_loc + 1
          io_e_frag = min((dg_frag%id_frag + 1) * nocc_per_rank_loc, nocc_per_spin)
        end if
        ! --- end D pre-pass ---
        do block_offset = first_block_offset, nblocks_ifrag - 1, block_step_blocks
          igrid0 = 1 + block_offset * grid_block_size
          npt_blk = min(grid_block_size, ngrid - igrid0 + 1)
          local_grid_count = 0
          remote_grid_count = 0
          valid_remote_grid_count = 0
          call cpu_time(t_setup0)
          if (allocated(dg_frag%density_send_slot_map)) then
            call prepare_grid_buffers_owner_map(i_local, igrid0, npt_blk, nxyz, .false.)
          else
            call prepare_grid_buffers_owner_map_no_slot(i_local, igrid0, npt_blk, nxyz)
          end if
          do igrid = 1, npt_blk
            if (owner_buf(igrid) < 0) cycle
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

          call cpu_time(t_setup0)
          phi_blk(1:npt_blk, 1:nbf) = dg_frag%density_phi_block_cache(1:npt_blk, 1:nbf, block_offset + 1, i_local)
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

!$omp parallel private(io, igrid, owner_rank, ixg, iyg, izg, rho_contrib, rho_raw_contrib, slot, rho_mix_accum)
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

!$omp do schedule(static)
                do idx_local = 1, local_grid_count
                  igrid = local_grid_ids(idx_local)
                  ixg = ixg_buf(igrid)
                  iyg = iyg_buf(igrid)
                  izg = izg_buf(igrid)
                  rho_raw_contrib = rho_blk(igrid)
                  rho_contrib = rho_raw_contrib
                  if (allocated(dg_frag%density_inv_weight_local)) then
                    rho_contrib = rho_contrib * dg_frag%density_inv_weight_local(ixg, iyg, izg)
                  end if
!$omp atomic
                  charge_budget_local(1) = charge_budget_local(1) + rho_raw_contrib
!$omp atomic
                  charge_budget_local(2) = charge_budget_local(2) + rho_contrib
                  rho_bf(ixg, iyg, izg) = rho_bf(ixg, iyg, izg) + rho_contrib
                  rho_s_bf(ixg, iyg, izg, ispin) = rho_s_bf(ixg, iyg, izg, ispin) + rho_contrib
                end do
!$omp end do nowait
!$omp do schedule(static)
                do idx_remote = 1, valid_remote_grid_count
                  igrid = valid_remote_grid_ids(idx_remote)
                  owner_rank = owner_buf(igrid)
                  slot = slot_buf(igrid)
                  if (owner_rank < 0 .or. owner_rank > dg_frag%isize - 1) then
                    write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0)') &
                      "DG density remote owner out of range: rank=", dg_frag%id, &
                      " id_frag=", dg_frag%id_frag, " i_local=", i_local, &
                      " igrid=", igrid, " owner=", owner_rank, " valid_remote=", valid_remote_grid_count
                    flush(6)
                    stop "DG-Fragment RT: density remote owner out of range"
                  end if
                  if (.not. allocated(rho_send(owner_rank)%f)) then
                    write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0)') &
                      "DG density remote send buffer missing: rank=", dg_frag%id, &
                      " id_frag=", dg_frag%id_frag, " i_local=", i_local, &
                      " igrid=", igrid, " owner=", owner_rank, " nsend=", dg_frag%density_send_count(owner_rank)
                    flush(6)
                    stop "DG-Fragment RT: density remote send buffer missing"
                  end if
                  if (slot < 1 .or. slot > dg_frag%density_send_count(owner_rank)) then
                    write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0)') &
                      "DG density remote slot out of range: rank=", dg_frag%id, &
                      " id_frag=", dg_frag%id_frag, " i_local=", i_local, &
                      " igrid=", igrid, " owner=", owner_rank, " slot=", slot, &
                      " nsend=", dg_frag%density_send_count(owner_rank)
                    flush(6)
                    stop "DG-Fragment RT: density remote slot out of range"
                  end if
                  rho_contrib = rho_blk(igrid)
!$omp atomic update
                  rho_send(owner_rank)%f(slot, 1, 1) = rho_send(owner_rank)%f(slot, 1, 1) + rho_contrib
                  spin_offset = ispin * dg_frag%density_send_count(owner_rank)
                  if (spin_offset + slot < 1 .or. spin_offset + slot > size(rho_send(owner_rank)%f, 1)) then
                    write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0)') &
                      "DG density remote spin slot out of range: rank=", dg_frag%id, &
                      " id_frag=", dg_frag%id_frag, " i_local=", i_local, &
                      " igrid=", igrid, " owner=", owner_rank, " slot=", slot, &
                      " spin_slot=", spin_offset + slot, " send_size=", size(rho_send(owner_rank)%f, 1)
                    flush(6)
                    stop "DG-Fragment RT: density remote spin slot out of range"
                  end if
!$omp atomic update
                  rho_send(owner_rank)%f(spin_offset + slot, 1, 1) = rho_send(owner_rank)%f(spin_offset + slot, 1, 1) + rho_contrib
                end do
!$omp end do
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
                ! n_pw > 0: state-distributed loop, no per-batch bcast
                if (.not. allocated(rho_blk_partial)) allocate(rho_blk_partial(grid_block_size))
                rho_blk_partial(1:npt_blk) = 0.0d0

                do io0 = io_s_frag, io_e_frag, state_block_size
                  nbatch = min(state_block_size, io_e_frag - io0 + 1)

                  ! copy coef from upfront buffer (no bcast)
                  coef_blk_re(1:nbf, 1:nbatch) = coef_re_full(1:nbf, io0:io0+nbatch-1, ispin)
                  coef_blk_im(1:nbf, 1:nbatch) = coef_im_full(1:nbf, io0:io0+nbatch-1, ispin)

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
                    ! direct access from full_cache on all ranks (no bcast)
                    coef_pw_blk(1:npw_blk, 1:nbatch) = &
                      dg_frag%coef_pw_full_cache(ipw0:ipw0+npw_blk-1, io0:io0+nbatch-1, ispin)
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
                    rho_blk_partial(igrid) = rho_blk_partial(igrid) + rho_accum
                  end do
!$omp end parallel do
                  call cpu_time(t_rho1)
                  time_project_rho = time_project_rho + (t_rho1 - t_rho0)
                end do  ! io0

                ! AllReduce rho_blk_partial across icomm_frag → rho_blk_accum
                ! Note: comm_summation overwrites rho_blk_accum (does not accumulate into it)
                call cpu_time(t_rho0)
                rho_blk_accum(1:npt_blk) = rho_blk_partial(1:npt_blk)
                call cpu_time(t_rho1)
                time_project_rho_reduce = time_project_rho_reduce + (t_rho1 - t_rho0)
              end if
              ! rho_blk_accum: filled by dgemm-path (n_pw==0) or AllReduce (n_pw>0)
              call cpu_time(t_rho0)
!$omp parallel private(igrid, owner_rank, ixg, iyg, izg, rho_contrib, rho_raw_contrib, slot)
!$omp do schedule(static)
                  do idx_local = 1, local_grid_count
                    igrid = local_grid_ids(idx_local)
                    ixg = ixg_buf(igrid)
                    iyg = iyg_buf(igrid)
                    izg = izg_buf(igrid)
                    rho_raw_contrib = rho_blk_accum(igrid)
                    rho_contrib = rho_raw_contrib
                    if (allocated(dg_frag%density_inv_weight_local)) then
                      rho_contrib = rho_contrib * dg_frag%density_inv_weight_local(ixg, iyg, izg)
                    end if
!$omp atomic
                    charge_budget_local(1) = charge_budget_local(1) + rho_raw_contrib
!$omp atomic
                    charge_budget_local(2) = charge_budget_local(2) + rho_contrib
                    rho_bf(ixg, iyg, izg) = rho_bf(ixg, iyg, izg) + rho_contrib
                    rho_s_bf(ixg, iyg, izg, ispin) = rho_s_bf(ixg, iyg, izg, ispin) + rho_contrib
                  end do
!$omp end do nowait
!$omp do schedule(static)
                  do idx_remote = 1, valid_remote_grid_count
                    igrid = valid_remote_grid_ids(idx_remote)
                    owner_rank = owner_buf(igrid)
                    slot = slot_buf(igrid)
                    if (owner_rank < 0 .or. owner_rank > dg_frag%isize - 1) then
                      write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0)') &
                        "DG density accum owner out of range: rank=", dg_frag%id, &
                        " id_frag=", dg_frag%id_frag, " i_local=", i_local, &
                        " igrid=", igrid, " owner=", owner_rank, " valid_remote=", valid_remote_grid_count
                      flush(6)
                      stop "DG-Fragment RT: density accum owner out of range"
                    end if
                    if (.not. allocated(rho_send(owner_rank)%f)) then
                      write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0)') &
                        "DG density accum send buffer missing: rank=", dg_frag%id, &
                        " id_frag=", dg_frag%id_frag, " i_local=", i_local, &
                        " igrid=", igrid, " owner=", owner_rank, " nsend=", dg_frag%density_send_count(owner_rank)
                      flush(6)
                      stop "DG-Fragment RT: density accum send buffer missing"
                    end if
                    if (slot < 1 .or. slot > dg_frag%density_send_count(owner_rank)) then
                      write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0)') &
                        "DG density accum slot out of range: rank=", dg_frag%id, &
                        " id_frag=", dg_frag%id_frag, " i_local=", i_local, &
                        " igrid=", igrid, " owner=", owner_rank, " slot=", slot, &
                        " nsend=", dg_frag%density_send_count(owner_rank)
                      flush(6)
                      stop "DG-Fragment RT: density accum slot out of range"
                    end if
                    rho_contrib = rho_blk_accum(igrid)
!$omp atomic update
                    rho_send(owner_rank)%f(slot, 1, 1) = rho_send(owner_rank)%f(slot, 1, 1) + rho_contrib
                    spin_offset = ispin * dg_frag%density_send_count(owner_rank)
                    if (spin_offset + slot < 1 .or. spin_offset + slot > size(rho_send(owner_rank)%f, 1)) then
                      write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0)') &
                        "DG density accum spin slot out of range: rank=", dg_frag%id, &
                        " id_frag=", dg_frag%id_frag, " i_local=", i_local, &
                        " igrid=", igrid, " owner=", owner_rank, " slot=", slot, &
                        " spin_slot=", spin_offset + slot, " send_size=", size(rho_send(owner_rank)%f, 1)
                      flush(6)
                      stop "DG-Fragment RT: density accum spin slot out of range"
                    end if
!$omp atomic update
                    rho_send(owner_rank)%f(spin_offset + slot, 1, 1) = rho_send(owner_rank)%f(spin_offset + slot, 1, 1) + rho_contrib
                  end do
!$omp end do
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
    if (allocated(coef_re_full)) deallocate(coef_re_full)
    if (allocated(coef_im_full)) deallocate(coef_im_full)
    if (allocated(rho_blk_partial)) deallocate(rho_blk_partial)
    if (dg_frag%is_frag_root) then
      write(*,'(1x,a,i0,a,i0,a,1pe12.4)') "        density trace: rank=", dg_frag%id, &
        " id_frag=", dg_frag%id_frag, " stage=after-project dt=", time_project
      write(*,'(1x,a,i0,a,i0,8(a,1pe12.4))') "        density trace: rank=", dg_frag%id, &
        " id_frag=", dg_frag%id_frag, &
        " setup=", time_project_setup, " psi=", time_project_psi, " rho=", time_project_rho, &
        " grid=", time_project_grid_prep, " phi=", time_project_phi_pack, &
        " over=", time_project_overhead, " dmat=", time_project_dmat_build, &
        " rho_red=", time_project_rho_reduce
      flush(6)
    end if

    if (dg_frag%is_frag_root) then
      write(*,'(1x,a,i0,a,i0,a)') "        density trace: rank=", dg_frag%id, &
        " id_frag=", dg_frag%id_frag, " stage=before-comm"
      flush(6)
    end if
    call cpu_time(t_comm0)
    allocate(send_counts(0:dg_frag%isize-1), recv_counts(0:dg_frag%isize-1))
    allocate(send_displs(0:dg_frag%isize-1), recv_displs(0:dg_frag%isize-1))
    send_counts = 0
    recv_counts = 0
    do irank = 0, dg_frag%isize - 1
      if (allocated(rho_send(irank)%f)) send_counts(irank) = size(rho_send(irank)%f)
      if (allocated(rho_recv(irank)%f)) recv_counts(irank) = size(rho_recv(irank)%f)
    end do
    send_displs(0) = 0
    recv_displs(0) = 0
    do irank = 1, dg_frag%isize - 1
      send_displs(irank) = send_displs(irank-1) + send_counts(irank-1)
      recv_displs(irank) = recv_displs(irank-1) + recv_counts(irank-1)
    end do
    write(*,'(1x,a,i0,a,i0,a,i0,a,i0)') "        density alltoallv summary: rank=", dg_frag%id, &
      " id_frag=", dg_frag%id_frag, " send_total=", sum(send_counts), " recv_total=", sum(recv_counts)
    flush(6)
    do irank = 0, dg_frag%isize - 1
      if (send_counts(irank) > 0 .or. recv_counts(irank) > 0) then
        write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0)') &
          "        density alltoallv peer: rank=", dg_frag%id, " id_frag=", dg_frag%id_frag, &
          " peer=", irank, " send_count=", send_counts(irank), " recv_count=", recv_counts(irank), &
          " send_pts=", dg_frag%density_send_count(irank), " recv_pts=", dg_frag%density_recv_map(irank)%npts
        flush(6)
      end if
    end do
    allocate(send_flat(max(1, sum(send_counts))), recv_flat(max(1, sum(recv_counts))))
    send_flat = 0.0d0
    recv_flat = 0.0d0
    call cpu_time(t_setup0)
    do irank = 0, dg_frag%isize - 1
      if (.not. allocated(rho_send(irank)%f)) cycle
      if (send_counts(irank) <= 0) cycle
      send_flat(send_displs(irank)+1:send_displs(irank)+send_counts(irank)) = rho_send(irank)%f(:, 1, 1)
      deallocate(rho_send(irank)%f)
    end do
    call cpu_time(t_setup1)
    time_comm_pack = time_comm_pack + (t_setup1 - t_setup0)
    write(*,'(1x,a,i0,a)') "        density collective: rank=", dg_frag%id, " stage=before-alltoallv"
    flush(6)
    call cpu_time(t_setup0)
    call comm_alltoallv(send_flat, send_counts, send_displs, recv_flat, recv_counts, recv_displs, dg_frag%icomm)
    call cpu_time(t_setup1)
    time_comm_exchange = time_comm_exchange + (t_setup1 - t_setup0)
    write(*,'(1x,a,i0,a)') "        density collective: rank=", dg_frag%id, " stage=after-alltoallv"
    flush(6)
    call cpu_time(t_setup0)
    do irank = 0, dg_frag%isize - 1
      if (.not. allocated(rho_recv(irank)%f)) cycle
      if (recv_counts(irank) > 0) then
        rho_recv(irank)%f(:, 1, 1) = recv_flat(recv_displs(irank)+1:recv_displs(irank)+recv_counts(irank))
      end if
      npts = dg_frag%density_recv_map(irank)%npts
!$omp parallel do private(slot, ixg, iyg, izg, ispin, spin_offset, rho_contrib, rho_raw_contrib) schedule(static)
      do slot = 1, dg_frag%density_recv_map(irank)%npts
        ixg = dg_frag%density_recv_map(irank)%ixg(slot)
        iyg = dg_frag%density_recv_map(irank)%iyg(slot)
        izg = dg_frag%density_recv_map(irank)%izg(slot)
        rho_raw_contrib = rho_recv(irank)%f(slot, 1, 1)
        rho_contrib = rho_raw_contrib
        if (allocated(dg_frag%density_inv_weight_local)) then
          rho_contrib = rho_contrib * dg_frag%density_inv_weight_local(ixg, iyg, izg)
        end if
!$omp atomic
        charge_budget_local(5) = charge_budget_local(5) + rho_raw_contrib
!$omp atomic
        charge_budget_local(6) = charge_budget_local(6) + rho_contrib
        rho_bf(ixg, iyg, izg) = rho_bf(ixg, iyg, izg) + rho_contrib
        do ispin = 1, system%nspin
          spin_offset = ispin * npts
          rho_contrib = rho_recv(irank)%f(spin_offset + slot, 1, 1)
          if (allocated(dg_frag%density_inv_weight_local)) then
            rho_contrib = rho_contrib * dg_frag%density_inv_weight_local(ixg, iyg, izg)
          end if
          rho_s_bf(ixg, iyg, izg, ispin) = rho_s_bf(ixg, iyg, izg, ispin) + rho_contrib
        end do
      end do
!$omp end parallel do
      deallocate(rho_recv(irank)%f)
    end do
    call cpu_time(t_setup1)
    time_comm_unpack = time_comm_unpack + (t_setup1 - t_setup0)
    deallocate(send_flat, recv_flat, send_counts, recv_counts, send_displs, recv_displs)
    call cpu_time(t_comm1)
    time_comm = time_comm + (t_comm1 - t_comm0)
    if (dg_frag%is_frag_root) then
      write(*,'(1x,a,i0,a,i0,a,1pe12.4)') "        density trace: rank=", dg_frag%id, &
        " id_frag=", dg_frag%id_frag, " stage=after-comm dt=", time_comm
      write(*,'(1x,a,i0,a,i0,4(a,1pe12.4))') "        density trace: rank=", dg_frag%id, &
        " id_frag=", dg_frag%id_frag, &
        " subgroup_reduce=", time_comm_subgroup_reduce, " pack=", time_comm_pack, &
        " exchange=", time_comm_exchange, " unpack=", time_comm_unpack
      flush(6)
    end if

    if (dg_frag%id == 0) then
      write(*,'(1x,a)') "        density trace: stage=before-normalize"
      flush(6)
    end if
    ! rho_bf -> rho_s boundary: only the authoritative mg-local density is
    ! materialized below for downstream Hartree/XC/reconstruct consumers.
    rho%f(:, :, :) = rho_bf(:, :, :)
    do ispin = 1, system%nspin
      rho_s(ispin)%f(:, :, :) = rho_s_bf(:, :, :, ispin)
    end do
    call comm_summation(charge_budget_local, charge_budget_total, 6, dg_frag%icomm)
    call cpu_time(t_norm0)
    total_charge_local = 0.0d0
!$omp parallel do collapse(3) reduction(+:total_charge_local) private(ix,iy,iz) schedule(static)
    do iz = rho_z_lo, rho_z_hi
      do iy = rho_y_lo, rho_y_hi
        do ix = rho_x_lo, rho_x_hi
          total_charge_local = total_charge_local + rho%f(ix, iy, iz)
        end do
      end do
    end do
!$omp end parallel do
    total_charge_local = total_charge_local * system%hvol
    if (dg_frag%id == 0) then
      write(*,'(1x,a,8(a,1pe12.4))') "        density charge budget:", &
        " raw_local=", charge_budget_total(1) * system%hvol, &
        " weighted_local=", charge_budget_total(2) * system%hvol, &
        " raw_subgroup=", charge_budget_total(3) * system%hvol, &
        " weighted_subgroup=", charge_budget_total(4) * system%hvol, &
        " raw_recv=", charge_budget_total(5) * system%hvol, &
        " weighted_recv=", charge_budget_total(6) * system%hvol, &
        " raw_total=", (charge_budget_total(1) + charge_budget_total(3) + charge_budget_total(5)) * system%hvol, &
        " weighted_total=", (charge_budget_total(2) + charge_budget_total(4) + charge_budget_total(6)) * system%hvol
      flush(6)
      write(*,'(1x,a,i0,a)') "        density collective: rank=", dg_frag%id, " stage=before-total-charge-sum"
      flush(6)
    end if
    call comm_summation(total_charge_local, total_charge, dg_frag%icomm)
    if (dg_frag%id == 0) then
      write(*,'(1x,a,i0,a,1pe12.4)') "        density collective: rank=", dg_frag%id, &
        " stage=after-total-charge-sum total_charge=", total_charge
      flush(6)
    end if
    dg_frag%elec_num_raw = total_charge
    dg_frag%rho_scale_factor = 1.0d0
    if (total_charge > 1.0d-14 .and. total_charge == total_charge) then
      scale_rho = nelec / total_charge
      dg_frag%rho_scale_factor = scale_rho
!$omp parallel do collapse(3) private(ix,iy,iz,ispin,rho_sum_local) schedule(static)
      do iz = rho_z_lo, rho_z_hi
        do iy = rho_y_lo, rho_y_hi
          do ix = rho_x_lo, rho_x_hi
            rho_sum_local = 0.0d0
            do ispin = 1, system%nspin
              rho_s(ispin)%f(ix, iy, iz) = rho_s(ispin)%f(ix, iy, iz) * scale_rho
              rho_sum_local = rho_sum_local + rho_s(ispin)%f(ix, iy, iz)
            end do
            rho%f(ix, iy, iz) = rho_sum_local
          end do
        end do
      end do
!$omp end parallel do
      dg_frag%elec_num_scaled = total_charge * scale_rho
    else
      dg_frag%elec_num_scaled = total_charge
    end if
    if (dg_frag%id == 0) then
      write(*,'(1x,a,i0,a)') "        density collective: rank=", dg_frag%id, " stage=before-scaled-charge-sum"
      flush(6)
      write(*,'(1x,a,i0,a,1pe12.4)') "        density collective: rank=", dg_frag%id, &
        " stage=after-scaled-charge-sum elec_num_scaled=", dg_frag%elec_num_scaled
      flush(6)
    end if

    call cpu_time(t_norm1)
    time_norm = time_norm + (t_norm1 - t_norm0)
    if (dg_frag%id == 0) then
      write(*,'(1x,a,a,1pe12.4)') "        density trace: stage=after-normalize dt=", "", time_norm
      flush(6)
    end if

    deallocate(ix_buf, iy_buf, iz_buf, owner_buf, ixg_buf, iyg_buf, izg_buf)
    deallocate(slot_buf, local_grid_ids, remote_grid_ids, valid_remote_grid_ids)
    deallocate(basis_gid, valid_basis_ids)
    deallocate(phi_blk, rho_blk, rho_blk_accum, coef_blk_re, coef_blk_im, psi_blk_re, psi_blk_im, density_mat_re, density_tmp, &
      psi_blk, coef_occ_weighted)
    if (allocated(phase_cache)) deallocate(phase_cache)
    if (allocated(density_mix)) deallocate(density_mix)
    if (allocated(basis_mix_blk)) deallocate(basis_mix_blk)
    if (allocated(density_mix_tmp)) deallocate(density_mix_tmp)
    if (allocated(transform_frag)) deallocate(transform_frag)
    if (allocated(transform_pw)) deallocate(transform_pw)
    deallocate(rho_bf, rho_s_bf)
    deallocate(rho_send, rho_recv)
    call cpu_time(t_total1)
    if (dg_frag%id == 0) then
      write(*,'(1x,a,1pe12.4,a,1pe12.4,a,1pe12.4,a,1pe12.4,a,1pe12.4)') &
        "        density timing: total=", t_total1 - t_total0, " cache=", time_cache, &
        " project=", time_project, " comm=", time_comm, " norm=", time_norm
      write(*,'(1x,a,4(a,1pe12.4),2(a,l1))') "        density cache timing:", &
        " pw_refresh=", time_cache_pw_refresh, " phi_blk_build=", time_cache_phi_block_refresh, &
        " cache_total=", time_cache + time_cache_phi_block_refresh, " step_total=", time_cache_pw_refresh + time_cache_phi_block_refresh, &
        " pw_rebuilt=", rebuilt_pw_cache, " phi_rebuilt=", rebuilt_phi_block_cache
      write(*,'(1x,a,7(a,1pe12.4))') "        density timing detail: project", &
        " setup=", time_project_setup, " psi=", time_project_psi, " rho=", time_project_rho, &
        " grid=", time_project_grid_prep, " phi=", time_project_phi_pack, &
        " phi_blk_build=", time_project_phi_block_build, " over=", time_project_overhead
      flush(6)
    end if

  contains

    subroutine prepare_grid_buffers_owner_map(i_local_grid, igrid0_grid, npt_blk_grid, nxyz_grid, use_subgroup_slot)
      implicit none
      integer, intent(in) :: i_local_grid, igrid0_grid, npt_blk_grid
      integer, intent(in) :: nxyz_grid(3)
      logical, intent(in) :: use_subgroup_slot
      type(density_grid_point_info) :: point

!$omp parallel do private(igrid, point) schedule(static)
      do igrid = 1, npt_blk_grid
        point = dg_frag%density_grid_points(igrid0_grid + igrid - 1, i_local_grid)
        ix_buf(igrid) = point%ix
        iy_buf(igrid) = point%iy
        iz_buf(igrid) = point%iz
        if (.not. point%is_primary) then
          ixg_buf(igrid) = 1
          iyg_buf(igrid) = 1
          izg_buf(igrid) = 1
          owner_buf(igrid) = -1
          slot_buf(igrid) = 0
          cycle
        end if
        ixg_buf(igrid) = point%ixg
        iyg_buf(igrid) = point%iyg
        izg_buf(igrid) = point%izg
        owner_buf(igrid) = point%owner_rank
        if (use_subgroup_slot) then
          slot_buf(igrid) = point%subgroup_send_slot
        else
          slot_buf(igrid) = point%send_slot
        end if
      end do
!$omp end parallel do
    end subroutine prepare_grid_buffers_owner_map

    subroutine prepare_grid_buffers_owner_map_no_slot(i_local_grid, igrid0_grid, npt_blk_grid, nxyz_grid)
      implicit none
      integer, intent(in) :: i_local_grid, igrid0_grid, npt_blk_grid
      integer, intent(in) :: nxyz_grid(3)
      type(density_grid_point_info) :: point

!$omp parallel do private(igrid, point) schedule(static)
      do igrid = 1, npt_blk_grid
        point = dg_frag%density_grid_points(igrid0_grid + igrid - 1, i_local_grid)
        ix_buf(igrid) = point%ix
        iy_buf(igrid) = point%iy
        iz_buf(igrid) = point%iz
        if (.not. point%is_primary) then
          ixg_buf(igrid) = 1
          iyg_buf(igrid) = 1
          izg_buf(igrid) = 1
          owner_buf(igrid) = -1
          slot_buf(igrid) = 0
          cycle
        end if
        ixg_buf(igrid) = point%ixg
        iyg_buf(igrid) = point%iyg
        izg_buf(igrid) = point%izg
        owner_buf(igrid) = point%owner_rank
        slot_buf(igrid) = 0
      end do
!$omp end parallel do
    end subroutine prepare_grid_buffers_owner_map_no_slot

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

  logical function target_rank_owned_by_handler(target_rank) result(is_handler)
    implicit none
    integer, intent(in) :: target_rank
    integer :: handler_id_frag

    handler_id_frag = modulo(target_rank, dg_frag%isize_frag)
    is_handler = (dg_frag%id_frag == handler_id_frag)
  end function target_rank_owned_by_handler

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

  integer function map_global_to_phi_box_coord_ham(ig, lb, ub, lgtot) result(iloc)
    implicit none
    integer, intent(in) :: ig, lb, ub, lgtot

    iloc = modulo(ig - 1, lgtot) + 1
    if (iloc < lb) then
      iloc = iloc + ((lb - iloc + lgtot - 1) / lgtot) * lgtot
    end if
    if (iloc > ub) then
      iloc = iloc - ((iloc - ub + lgtot - 1) / lgtot) * lgtot
    end if
    if (iloc < lb .or. iloc > ub) iloc = 0
  end function map_global_to_phi_box_coord_ham

  end subroutine calculate_density_from_fragments
