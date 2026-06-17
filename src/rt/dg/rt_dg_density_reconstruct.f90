  subroutine calculate_density_from_fragments(dg_frag, system, mg, rho, rho_s, itt_debug)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use structures
    use communication, only: comm_summation, comm_bcast, COMM_GROUP_NULL
    use rt_dg_fragment_ops, only: refresh_pw_coef_cache, apply_overlap_operator_batch, fetch_remote_coef_rows
    use rt_dg_fragment_types, only: s_dg_fragment_rt, density_grid_point_info
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system),     intent(in)    :: system
    type(s_rgrid),          intent(in)    :: mg
    type(s_scalar),         intent(inout) :: rho
    type(s_scalar),         intent(inout) :: rho_s(system%nspin)
    integer, intent(in), optional :: itt_debug

    integer :: ifrag, io, i_local, ispin, iw
    integer :: istate_frag
    integer :: ix, iy, iz, ixg, iyg, izg, bx, by, bz, owner_rank
    integer :: im_pf, ipw_pf, ipw_max
    integer :: ig_i, nbf, nbf_max, ipw, n_pw, n_frag, n_tot, n_basis_mix, max_mixed_basis
    integer :: nxyz(3), ifrag_count, ngrid_max
    integer :: nkeep, nkeep_max_density, nw_owned
    integer :: nocc_spin, nocc_cache
    integer :: irank, slot, npts, idx_local, idx_remote, point_idx, local_coef_idx, local_state_col
    integer :: local_grid_count, remote_grid_count, valid_remote_grid_count
    integer :: igrid0, igrid, ngrid, npt_blk, io0, io1, nbatch, nstate, ipw0, npw_blk
    integer :: total_send_pts, subgroup_root_rank, block_idx_global
    integer :: send_total_count, recv_total_count
    integer :: nblocks_ifrag, first_block_offset, block_step_blocks, block_offset
    integer :: valid_basis_count
    integer :: spin_offset, density_payload_count
    integer, parameter :: grid_block_size = 8192
    integer, parameter :: state_block_size = 64
    integer, parameter :: rho_state_block_size = 16
    integer, parameter :: pw_block_size = 128
    integer, parameter :: mixed_io_block_size = 64
    real(8), parameter :: density_abs_limit = 1.0d12
    complex(8), parameter :: zzero = (0.0d0, 0.0d0), zone = (1.0d0, 0.0d0)
    real(8) :: rho_contrib, rho_raw_contrib, rho_accum
    real(8) :: total_charge, total_charge_local
    real(8) :: total_charge_reduce_in(1), total_charge_reduce_out(1)
    real(8) :: occ_factor, abs_max_pf
    real(8) :: rho_trace_min, rho_trace_max
    real(8) :: boxL(3), inv_sqrt_vol, theta, inv_lgnum1
    logical :: use_mixed_density
    logical :: enable_wannier_density
    logical :: use_buffer_wannier_density
    logical :: enable_density_phi_block_cache, enable_density_phase_block_cache
    logical :: rebuilt_pw_cache, rebuilt_phi_block_cache, rebuilt_phase_block_cache
    logical :: need_pw_cache_alloc, need_pw_cache_expand
    logical :: need_phi_cache_alloc, need_phi_count_alloc, need_phi_cache_invalid, need_phi_cache_resize
    logical :: need_phase_cache_alloc, need_phase_cache_invalid, need_phase_cache_resize, need_phase_cache_npw
    logical :: need_rhobf_box_cache
    logical :: env_phi_block_cache_seen, env_phase_block_cache_seen
    integer, allocatable :: ix_buf(:), iy_buf(:), iz_buf(:), owner_buf(:), ixg_buf(:), iyg_buf(:), izg_buf(:)
    integer, allocatable :: rhobf_bx_buf(:), rhobf_by_buf(:), rhobf_bz_buf(:)
    integer, allocatable :: slot_buf(:), local_grid_ids(:), remote_grid_ids(:), valid_remote_grid_ids(:)
    integer, allocatable :: basis_gid(:), valid_basis_ids(:)
    integer, allocatable :: basis_gid_spin(:,:), valid_basis_ids_spin(:,:), valid_basis_count_spin(:)
    integer, allocatable :: owned_basis_ids_spin(:,:), owned_coef_local_ids_spin(:,:), owned_basis_count_spin(:)
    type(s_scalar), allocatable :: rho_send(:), rho_recv(:)
    integer, allocatable :: send_counts(:), recv_counts(:), send_displs(:), recv_displs(:)
    real(8), allocatable :: send_flat(:), recv_flat(:)
    real(8), allocatable :: rho_bf(:,:,:), rho_s_bf(:,:,:,:)
    real(8), allocatable :: phi_blk(:,:), rho_blk(:)
    real(8), allocatable :: rho_blk_accum(:), rho_blk_reduced(:)
    real(8), allocatable :: coef_blk_re(:,:), coef_blk_im(:,:)
    real(8), allocatable :: psi_blk_re(:,:), psi_blk_im(:,:)
    real(8), allocatable :: coef_blk_ri(:,:), psi_blk_ri(:,:)
    real(8), allocatable :: D_frag_re(:,:,:)   ! (nbf_max, nbf_max, nspin) pre-computed D per fragment
    real(8), allocatable :: D_partial_re(:,:)    ! (nbf_max, nbf_max) partial D per rank
    real(8), allocatable :: rho_dmat_q_local(:,:), rho_dmat_q_full(:,:)
    real(8), allocatable :: D_wannier_re(:,:,:), D_wannier_tmp(:,:)
    real(8), allocatable :: wannier_coef_owned(:,:,:)
    real(8), allocatable :: wannier_blk_local(:,:), wannier_blk_full(:,:), wannier_dmat_q(:,:)
    integer, allocatable :: wannier_owned_count_spin(:)
    real(8), allocatable :: coef_re_full(:,:,:)  ! (nbf_max, nocc_cache, nspin) upfront bcast coef (n_pw>0)
    real(8), allocatable :: coef_im_full(:,:,:)  ! (nbf_max, nocc_cache, nspin)
    real(8), allocatable :: rho_blk_partial(:)   ! (grid_block_size) partial rho for state slice
    real(8), allocatable :: occ_cache(:), occ_sqrt_cache(:), occ_blk(:)
    complex(8), allocatable :: coef_c_full(:,:), coef_c_frag(:,:)
    integer :: io_s_frag, io_e_frag, io_loc, nocc_loc, nocc_mix_cols, nocc_per_rank_loc
    integer :: ib_s_frag, ib_e_frag, ib_loc, ib_global, coef_local
    integer :: nbf_frag_is, nbf_frag_ie, nbf_frag_count, nbf_frag_cap, nbf_per_rank_loc
    integer :: nblocks_max, block_cache_idx, npt_cache
    integer :: phi_lb1, phi_lb2, phi_lb3, phi_ub1, phi_ub2, phi_ub3
    integer :: phi_lg1, phi_lg2, phi_lg3
    integer :: ibuf_x, ibuf_y, ibuf_z
    integer :: rho_x_lo, rho_x_hi, rho_y_lo, rho_y_hi, rho_z_lo, rho_z_hi
    integer :: rho_s_x_lo, rho_s_x_hi, rho_s_y_lo, rho_s_y_hi, rho_s_z_lo, rho_s_z_hi
    complex(8), allocatable :: phase_cache(:,:), coef_pw_blk(:,:), pw_tmp_z(:,:)
    complex(8), allocatable :: density_mix(:,:,:), density_mix_partial(:,:), basis_mix_blk(:,:), density_mix_tmp(:,:)
    complex(8), allocatable :: basis_mix_blk_t(:,:), density_mix_tmp_t(:,:)
    complex(8), allocatable :: transform_frag_spin(:,:,:), transform_pw_spin(:,:,:)
    complex(8), allocatable :: mix_transform_spin(:,:), mix_overlap_spin(:,:), s_mix(:,:), s_mix_work(:,:)
    complex(8), allocatable :: coef_mix_eff(:,:), coef_mix_metric(:,:), coef_mix_spin(:,:,:)
    complex(8), allocatable :: d_raw_ff(:,:), d_raw_fp(:,:), d_raw_pp(:,:)
    complex(8) :: phase_fix
    integer, allocatable :: ipiv_mix(:)
    integer, allocatable :: n_basis_mix_spin(:)
    real(8), allocatable :: kpw_hx(:), kpw_hy(:), kpw_hz(:)
    character(32) :: env_phi_block_cache, env_phase_block_cache
    character(32) :: env_wannier_density
    character(32) :: env_rho_mix_mode
    character(32) :: env_trace_density_charge
    integer :: env_status
    integer :: rho_mix_mode_kind
    integer :: info_lapack
    integer :: itt_tag
    logical :: need_full_coef_mix_spin
    logical :: density_on_frag_root
    logical :: density_exchange_active
    logical :: trace_density_charge
    logical :: enable_fp_phase_fix
    character(32) :: env_fp_phase_fix
    logical, save :: density_env_initialized = .false.
    logical, save :: cfg_env_phi_block_cache_seen = .false.
    logical, save :: cfg_env_phase_block_cache_seen = .false.
    ! Keep the block-packed fragment basis cache enabled by default for RT
    ! throughput on large runs.  Set SALMON_DG_DENSITY_PHI_BLOCK_CACHE=0 to
    ! trade speed back for lower memory use.
    logical, save :: cfg_enable_density_phi_block_cache = .true.
    logical, save :: cfg_enable_density_phase_block_cache = .false.
    logical, save :: cfg_enable_wannier_density = .false.
    logical, save :: cfg_wannier_density_logged = .false.
    logical, save :: cfg_trace_density_charge = .false.
    character(32), save :: cfg_env_rho_mix_mode = 'legacy'
    integer, save :: cfg_rho_mix_mode_kind = 0
    logical, save :: cfg_enable_fp_phase_fix = .false.

    itt_tag = -1
    if (present(itt_debug)) itt_tag = itt_debug

    rho%f = 0.0d0
    do ispin = 1, system%nspin
      rho_s(ispin)%f = 0.0d0
    end do
    rebuilt_pw_cache = .false.
    rebuilt_phi_block_cache = .false.
    rebuilt_phase_block_cache = .false.
    enable_density_phi_block_cache = .false.
    enable_density_phase_block_cache = .false.
    env_phi_block_cache_seen = .false.
    env_phase_block_cache_seen = .false.
    env_phi_block_cache = ''
    env_phase_block_cache = ''
    env_wannier_density = ''
    env_rho_mix_mode = 'legacy'
    env_trace_density_charge = ''
    env_fp_phase_fix = ''
    rho_mix_mode_kind = 0
    need_full_coef_mix_spin = .false.
    enable_fp_phase_fix = .false.
    trace_density_charge = .false.
    ! These switches are process-level controls.  Cache them once because this
    ! routine is called many times per RT step in self-consistent propagation.
    if (.not. density_env_initialized) then
      call get_environment_variable('SALMON_DG_DENSITY_PHI_BLOCK_CACHE', env_phi_block_cache, status=env_status)
      if (env_status == 0) then
        cfg_env_phi_block_cache_seen = .true.
        select case (adjustl(trim(env_phi_block_cache)))
        case('1','y','Y','yes','YES','true','TRUE','on','ON')
          cfg_enable_density_phi_block_cache = .true.
        case('0','n','N','no','NO','false','FALSE','off','OFF')
          cfg_enable_density_phi_block_cache = .false.
        end select
      end if
      call get_environment_variable('SALMON_DG_DENSITY_PHASE_BLOCK_CACHE', env_phase_block_cache, status=env_status)
      if (env_status == 0) then
        cfg_env_phase_block_cache_seen = .true.
        select case (adjustl(trim(env_phase_block_cache)))
        case('1','y','Y','yes','YES','true','TRUE','on','ON')
          cfg_enable_density_phase_block_cache = .true.
        case('0','n','N','no','NO','false','FALSE','off','OFF')
          cfg_enable_density_phase_block_cache = .false.
        end select
      end if
      call get_environment_variable('SALMON_DG_DENSITY_WANNIER', env_wannier_density, status=env_status)
      if (env_status == 0) then
        select case (adjustl(trim(env_wannier_density)))
        case('1','y','Y','yes','YES','true','TRUE','on','ON')
          cfg_enable_wannier_density = .true.
        case('0','n','N','no','NO','false','FALSE','off','OFF')
          cfg_enable_wannier_density = .false.
        end select
      end if
      call get_environment_variable('SALMON_DG_DENSITY_TRACE_CHARGE', env_trace_density_charge, status=env_status)
      if (env_status == 0) then
        select case (adjustl(trim(env_trace_density_charge)))
        case('1','y','Y','yes','YES','true','TRUE','on','ON')
          cfg_trace_density_charge = .true.
        case('0','n','N','no','NO','false','FALSE','off','OFF')
          cfg_trace_density_charge = .false.
        end select
      end if
      call get_environment_variable('SALMON_DG_RHO_MIX_MODE', env_rho_mix_mode, status=env_status)
      if (env_status == 0) then
        select case (adjustl(trim(env_rho_mix_mode)))
        case('orthonormal_cc','ORTHONORMAL_CC')
          cfg_rho_mix_mode_kind = 1
          cfg_env_rho_mix_mode = env_rho_mix_mode
        case('metric_consistent','METRIC_CONSISTENT','overlap_metric','OVERLAP_METRIC')
          cfg_rho_mix_mode_kind = 2
          cfg_env_rho_mix_mode = env_rho_mix_mode
        case default
          cfg_rho_mix_mode_kind = 0
          cfg_env_rho_mix_mode = 'legacy'
        end select
      else
        cfg_env_rho_mix_mode = 'legacy'
      end if
      call get_environment_variable('SALMON_DG_FP_PHASE_FIX', env_fp_phase_fix, status=env_status)
      if (env_status == 0) then
        select case (adjustl(trim(env_fp_phase_fix)))
        case('1','y','Y','yes','YES','true','TRUE','on','ON')
          cfg_enable_fp_phase_fix = .true.
        case('0','n','N','no','NO','false','FALSE','off','OFF')
          cfg_enable_fp_phase_fix = .false.
        end select
      end if
      density_env_initialized = .true.
    end if
    env_phi_block_cache_seen = cfg_env_phi_block_cache_seen
    env_phase_block_cache_seen = cfg_env_phase_block_cache_seen
    enable_density_phi_block_cache = cfg_enable_density_phi_block_cache
    enable_density_phase_block_cache = cfg_enable_density_phase_block_cache
    enable_wannier_density = cfg_enable_wannier_density
    trace_density_charge = cfg_trace_density_charge
    use_buffer_wannier_density = enable_wannier_density .and. &
      dg_frag%buffer_wannier_flux_seed_applied .and. dg_frag%has_buffer_periodic_wannier_basis .and. &
      allocated(dg_frag%buffer_wannier_coef)
    if (enable_wannier_density) then
      if (use_buffer_wannier_density) then
        if (.not. allocated(dg_frag%buffer_wannier_nkeep)) &
          stop "DG density Wannier path requires buffer_periodic_wannier_basis.bin"
      else if (.not. dg_frag%has_local_wannier_basis) then
        stop "DG density Wannier path requires local_wannier_basis.bin"
      end if
    end if
    env_rho_mix_mode = cfg_env_rho_mix_mode
    rho_mix_mode_kind = cfg_rho_mix_mode_kind
    enable_fp_phase_fix = cfg_enable_fp_phase_fix
    need_full_coef_mix_spin = .false.
    density_payload_count = merge(system%nspin + 1, 1, system%nspin > 1)
    need_pw_cache_alloc = .false.
    need_pw_cache_expand = .false.
    need_phi_cache_alloc = .false.
    need_phi_count_alloc = .false.
    need_phi_cache_invalid = .false.
    need_phi_cache_resize = .false.
    need_phase_cache_alloc = .false.
    need_phase_cache_invalid = .false.
    need_phase_cache_resize = .false.
    need_phase_cache_npw = .false.
    need_rhobf_box_cache = .false.

    if (.not. allocated(dg_frag%phi_frag)) return
    inv_lgnum1 = 1.0d0 / dble(max(1, dg_frag%lgnum_total(1)))

    ifrag_count = dg_frag%ifrag_end - dg_frag%ifrag_start + 1
    ngrid_max = 0
    if (ifrag_count > 0) then
      if (allocated(dg_frag%density_grid_point_count)) then
        ngrid_max = max(1, maxval(dg_frag%density_grid_point_count(1:ifrag_count)))
      else
        do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
          ngrid_max = max(ngrid_max, product(dg_frag%nxyz_domain(:, ifrag)))
        end do
      end if
    end if
    nbf_max = max(1, maxval(dg_frag%n_basis(:, 1:system%nspin)))
    n_pw = max(0, dg_frag%n_plane_waves)
    if (enable_wannier_density) then
      if (n_pw > 0) &
        stop "DG density Wannier path is only implemented for pure fragment basis"
      if (.not. cfg_wannier_density_logged .and. dg_frag%id == 0) then
        if (use_buffer_wannier_density) then
          write(*,'(1x,a)') &
            "[DG-DENSITY-WANNIER] enabled: buffer-periodic Wannier projection on core density grid"
        else
          write(*,'(1x,a)') &
            "[DG-DENSITY-WANNIER] enabled: center-owned local Wannier projection on core density grid"
        end if
        cfg_wannier_density_logged = .true.
      end if
    end if
    if (n_pw > 0) then
      stop "DG density PW path requires row-local PW reconstruction; full PW coefficient cache is disabled"
    end if
    nbf_frag_cap = nbf_max
    if (dg_frag%parallel_mode_orbital .and. dg_frag%isize_frag > 1) then
      nbf_frag_cap = (nbf_max + dg_frag%isize_frag - 1) / dg_frag%isize_frag
    end if
    nbf_frag_cap = max(1, nbf_frag_cap)
    allocate(ix_buf(grid_block_size), iy_buf(grid_block_size), iz_buf(grid_block_size))
    allocate(owner_buf(grid_block_size), ixg_buf(grid_block_size), iyg_buf(grid_block_size), izg_buf(grid_block_size))
    allocate(rhobf_bx_buf(grid_block_size), rhobf_by_buf(grid_block_size), rhobf_bz_buf(grid_block_size))
    allocate(slot_buf(grid_block_size), local_grid_ids(grid_block_size), remote_grid_ids(grid_block_size))
    allocate(valid_remote_grid_ids(grid_block_size))
    allocate(basis_gid(dg_frag%nstate_frag), valid_basis_ids(dg_frag%nstate_frag))
    allocate(basis_gid_spin(nbf_max, system%nspin), &
      valid_basis_ids_spin(nbf_max, system%nspin), &
      valid_basis_count_spin(system%nspin))
    allocate(owned_basis_ids_spin(nbf_max, system%nspin), &
      owned_coef_local_ids_spin(nbf_max, system%nspin), &
      owned_basis_count_spin(system%nspin))
    allocate(phi_blk(grid_block_size, nbf_frag_cap))
    allocate(rho_blk(grid_block_size))
    allocate(rho_blk_accum(grid_block_size))
    allocate(rho_blk_reduced(grid_block_size))
    nocc_cache = 0
    if (allocated(dg_frag%nocc_spin)) then
      nocc_cache = maxval(dg_frag%nocc_spin(1:system%nspin))
    end if
    nocc_cache = max(1, nocc_cache)
    allocate(coef_blk_re(nbf_max, state_block_size), coef_blk_im(nbf_max, state_block_size))
    allocate(psi_blk_re(grid_block_size, state_block_size), psi_blk_im(grid_block_size, state_block_size))
    allocate(coef_blk_ri(nbf_max, 2*state_block_size), psi_blk_ri(grid_block_size, 2*state_block_size))
    allocate(occ_cache(max(1, nocc_cache)))
    allocate(occ_blk(state_block_size))
    if (n_pw > 0) then
      allocate(pw_tmp_z(grid_block_size, state_block_size))
      allocate(coef_pw_blk(pw_block_size, state_block_size))
    end if
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
    if (system%nspin > 1) then
      allocate(rho_s_bf(rho_s_x_lo:rho_s_x_hi, &
                        rho_s_y_lo:rho_s_y_hi, &
                        rho_s_z_lo:rho_s_z_hi, &
        system%nspin))
    end if
    allocate(rho_send(0:dg_frag%isize-1))
    allocate(rho_recv(0:dg_frag%isize-1))

    need_rhobf_box_cache = (.not. dg_frag%density_rhobf_box_cache_valid)
    need_rhobf_box_cache = need_rhobf_box_cache .or. &
      any(dg_frag%density_rhobf_box_lo /= (/ lbound(rho_bf, 1), lbound(rho_bf, 2), lbound(rho_bf, 3) /))
    need_rhobf_box_cache = need_rhobf_box_cache .or. &
      any(dg_frag%density_rhobf_box_hi /= (/ ubound(rho_bf, 1), ubound(rho_bf, 2), ubound(rho_bf, 3) /))
    need_rhobf_box_cache = need_rhobf_box_cache .or. (.not. allocated(dg_frag%density_grid_bx)) .or. &
                            (.not. allocated(dg_frag%density_grid_by)) .or. (.not. allocated(dg_frag%density_grid_bz))
    if (.not. need_rhobf_box_cache .and. allocated(dg_frag%density_grid_bx)) then
      need_rhobf_box_cache = size(dg_frag%density_grid_bx, 1) /= size(dg_frag%density_grid_points, 1) .or. &
                              size(dg_frag%density_grid_bx, 2) /= size(dg_frag%density_grid_points, 2)
    end if
    if (.not. need_rhobf_box_cache) then
      do irank = 0, dg_frag%isize - 1
        if (dg_frag%density_recv_map(irank)%npts <= 0) cycle
        if ((.not. allocated(dg_frag%density_recv_map(irank)%bx)) .or. &
            (.not. allocated(dg_frag%density_recv_map(irank)%by)) .or. &
            (.not. allocated(dg_frag%density_recv_map(irank)%bz))) then
          need_rhobf_box_cache = .true.
          exit
        end if
      end do
    end if
    if (need_rhobf_box_cache) then
      if (allocated(dg_frag%density_grid_bx)) deallocate(dg_frag%density_grid_bx)
      if (allocated(dg_frag%density_grid_by)) deallocate(dg_frag%density_grid_by)
      if (allocated(dg_frag%density_grid_bz)) deallocate(dg_frag%density_grid_bz)
      allocate(dg_frag%density_grid_bx(size(dg_frag%density_grid_points, 1), size(dg_frag%density_grid_points, 2)))
      allocate(dg_frag%density_grid_by(size(dg_frag%density_grid_points, 1), size(dg_frag%density_grid_points, 2)))
      allocate(dg_frag%density_grid_bz(size(dg_frag%density_grid_points, 1), size(dg_frag%density_grid_points, 2)))
      dg_frag%density_grid_bx(:, :) = 0
      dg_frag%density_grid_by(:, :) = 0
      dg_frag%density_grid_bz(:, :) = 0
      do i_local = 1, ifrag_count
        do point_idx = 1, dg_frag%density_grid_point_count(i_local)
          ixg = dg_frag%density_grid_points(point_idx, i_local)%ixg
          iyg = dg_frag%density_grid_points(point_idx, i_local)%iyg
          izg = dg_frag%density_grid_points(point_idx, i_local)%izg
          dg_frag%density_grid_bx(point_idx, i_local) = &
            map_global_to_phi_box_coord_ham(ixg, lbound(rho_bf, 1), ubound(rho_bf, 1), dg_frag%lgnum_total(1))
          dg_frag%density_grid_by(point_idx, i_local) = &
            map_global_to_phi_box_coord_ham(iyg, lbound(rho_bf, 2), ubound(rho_bf, 2), dg_frag%lgnum_total(2))
          dg_frag%density_grid_bz(point_idx, i_local) = &
            map_global_to_phi_box_coord_ham(izg, lbound(rho_bf, 3), ubound(rho_bf, 3), dg_frag%lgnum_total(3))
        end do
      end do
      do irank = 0, dg_frag%isize - 1
        npts = dg_frag%density_recv_map(irank)%npts
        if (npts <= 0) cycle
        if (allocated(dg_frag%density_recv_map(irank)%bx)) deallocate(dg_frag%density_recv_map(irank)%bx)
        if (allocated(dg_frag%density_recv_map(irank)%by)) deallocate(dg_frag%density_recv_map(irank)%by)
        if (allocated(dg_frag%density_recv_map(irank)%bz)) deallocate(dg_frag%density_recv_map(irank)%bz)
        allocate(dg_frag%density_recv_map(irank)%bx(npts))
        allocate(dg_frag%density_recv_map(irank)%by(npts))
        allocate(dg_frag%density_recv_map(irank)%bz(npts))
        do slot = 1, npts
          dg_frag%density_recv_map(irank)%bx(slot) = map_global_to_phi_box_coord_ham( &
            dg_frag%density_recv_map(irank)%ixg(slot), lbound(rho_bf, 1), ubound(rho_bf, 1), dg_frag%lgnum_total(1))
          dg_frag%density_recv_map(irank)%by(slot) = map_global_to_phi_box_coord_ham( &
            dg_frag%density_recv_map(irank)%iyg(slot), lbound(rho_bf, 2), ubound(rho_bf, 2), dg_frag%lgnum_total(2))
          dg_frag%density_recv_map(irank)%bz(slot) = map_global_to_phi_box_coord_ham( &
            dg_frag%density_recv_map(irank)%izg(slot), lbound(rho_bf, 3), ubound(rho_bf, 3), dg_frag%lgnum_total(3))
        end do
      end do
      dg_frag%density_rhobf_box_lo = (/ lbound(rho_bf, 1), lbound(rho_bf, 2), lbound(rho_bf, 3) /)
      dg_frag%density_rhobf_box_hi = (/ ubound(rho_bf, 1), ubound(rho_bf, 2), ubound(rho_bf, 3) /)
      dg_frag%density_rhobf_box_cache_valid = .true.
    end if
    rho%f = 0.0d0
    rho_bf(:, :, :) = 0.0d0
    if (system%nspin > 1) rho_s_bf(:, :, :, :) = 0.0d0
    rho_blk_reduced(:) = 0.0d0
    if (n_pw > 0 .and. .not. env_phase_block_cache_seen) enable_density_phase_block_cache = .true.
    if (n_pw <= 0) enable_density_phase_block_cache = .false.
    n_frag = dg_frag%n_mat_max
    max_mixed_basis = 0
    if (n_pw > 0 .and. allocated(dg_frag%mixed_basis_dim)) then
      max_mixed_basis = maxval(dg_frag%mixed_basis_dim(1:system%nspin))
    end if
    n_tot = n_frag + n_pw
    if (n_pw > 0) then
      allocate(phase_cache(grid_block_size, n_pw))
      allocate(kpw_hx(n_pw), kpw_hy(n_pw), kpw_hz(n_pw))
      kpw_hx(1:n_pw) = dg_frag%k_pw(1, 1:n_pw) * dg_frag%hgs(1)
      kpw_hy(1:n_pw) = dg_frag%k_pw(2, 1:n_pw) * dg_frag%hgs(2)
      kpw_hz(1:n_pw) = dg_frag%k_pw(3, 1:n_pw) * dg_frag%hgs(3)
    end if

    ! Mixed-basis density reconstruction is only valid when PW channels exist.
    ! For n_pw==0, stay on the pure fragment-basis path.
    use_mixed_density = (n_pw > 0 .and. dg_frag%mixed_basis_ready .and. allocated(dg_frag%mixed_transform) .and. &
      allocated(dg_frag%coef_mix) .and. allocated(dg_frag%mixed_basis_dim))
    ! Density reconstruction uses subgroup-distributed projection and collective reductions on icomm_frag.
    subgroup_root_rank = dg_frag%id - dg_frag%id_frag
    total_send_pts = 0
    do irank = 0, dg_frag%isize - 1
      npts = dg_frag%density_send_count(irank)
      if (npts <= 0) cycle
      allocate(rho_send(irank)%f(density_payload_count * npts, 1, 1))
      rho_send(irank)%f(:, :, :) = 0.0d0
    end do
    do irank = 0, dg_frag%isize - 1
      npts = dg_frag%density_recv_map(irank)%npts
      if (npts <= 0) cycle
    end do
    ! Production density is accumulated from occupied-state projections below.
    ! The raw fragment density-matrix cache is not consumed by current callers,
    ! so avoid rebuilding and zeroing it on every RT density evaluation.
    if (use_mixed_density) then
      use_mixed_density = (max_mixed_basis > 0)
    end if
    if (use_mixed_density) then
      ! The mixed-density projector now forms rho from occupied-state
      ! amplitudes after reducing the fragment-basis contribution across
      ! orbital ranks.  Keep a full occupied-state coefficient view so every
      ! rank participates in the same state batches before the root emits rho.
      need_full_coef_mix_spin = .true.
      allocate(density_mix(max_mixed_basis, max_mixed_basis, system%nspin))
      density_mix(:, :, :) = (0.0d0, 0.0d0)
      allocate(density_mix_partial(max_mixed_basis, max_mixed_basis))
      if (need_full_coef_mix_spin) then
        allocate(coef_mix_spin(max_mixed_basis, max(1, nocc_cache), system%nspin))
        coef_mix_spin(:, :, :) = (0.0d0, 0.0d0)
      end if
      allocate(basis_mix_blk(grid_block_size, max_mixed_basis))
      allocate(density_mix_tmp(grid_block_size, max_mixed_basis))
      allocate(basis_mix_blk_t(max_mixed_basis, grid_block_size))
      allocate(density_mix_tmp_t(max_mixed_basis, grid_block_size))
      allocate(transform_frag_spin(dg_frag%nstate_frag, max_mixed_basis, system%nspin))
      allocate(n_basis_mix_spin(system%nspin))
      n_basis_mix_spin(:) = 0
      transform_frag_spin(:, :, :) = (0.0d0, 0.0d0)
      if (n_pw > 0) then
        allocate(transform_pw_spin(n_pw, max_mixed_basis, system%nspin))
        transform_pw_spin(:, :, :) = (0.0d0, 0.0d0)
      end if
      do ispin = 1, system%nspin
        nocc_spin = dg_frag%nocc_spin(ispin)
        n_basis_mix = min(dg_frag%mixed_basis_dim(ispin), max_mixed_basis, size(dg_frag%coef_mix, 1))
        if (n_basis_mix <= 0 .or. nocc_spin <= 0) cycle
        if (allocated(s_mix)) deallocate(s_mix)
        if (allocated(mix_transform_spin)) deallocate(mix_transform_spin)
        if (allocated(mix_overlap_spin)) deallocate(mix_overlap_spin)
        if (allocated(s_mix_work)) deallocate(s_mix_work)
        if (allocated(coef_mix_eff)) deallocate(coef_mix_eff)
        if (allocated(coef_mix_metric)) deallocate(coef_mix_metric)
        if (allocated(ipiv_mix)) deallocate(ipiv_mix)

        if (rho_mix_mode_kind == 2) then
          allocate(mix_transform_spin(n_tot, n_basis_mix), mix_overlap_spin(n_tot, n_basis_mix), s_mix(n_basis_mix, n_basis_mix))
          mix_transform_spin(:, :) = dg_frag%mixed_transform(1:n_tot, 1:n_basis_mix, ispin)
          call apply_overlap_operator_batch(dg_frag, ispin, mix_transform_spin, mix_overlap_spin, .false.)
          s_mix(:, :) = matmul(conjg(transpose(mix_transform_spin)), mix_overlap_spin)
        end if

        nocc_per_rank_loc = (nocc_spin + max(1, dg_frag%isize_frag) - 1) / max(1, dg_frag%isize_frag)
        io_s_frag = dg_frag%id_frag * nocc_per_rank_loc + 1
        io_e_frag = min((dg_frag%id_frag + 1) * nocc_per_rank_loc, nocc_spin)
        nocc_loc = max(0, io_e_frag - io_s_frag + 1)
        nocc_mix_cols = merge(nocc_spin, nocc_loc, need_full_coef_mix_spin)

        if (rho_mix_mode_kind == 2 .and. allocated(s_mix)) then
          allocate(coef_mix_metric(n_basis_mix, max(1, nocc_mix_cols)), coef_mix_eff(n_basis_mix, max(1, nocc_mix_cols)))
          coef_mix_metric(:, :) = (0.0d0, 0.0d0)
          coef_mix_eff(:, :) = (0.0d0, 0.0d0)
          if (need_full_coef_mix_spin) then
            coef_mix_metric(1:n_basis_mix, 1:nocc_spin) = dg_frag%coef_mix(1:n_basis_mix, 1:nocc_spin, ispin)
          else if (nocc_loc > 0) then
            coef_mix_metric(1:n_basis_mix, 1:nocc_loc) = dg_frag%coef_mix(1:n_basis_mix, io_s_frag:io_e_frag, ispin)
          end if
          if (nocc_mix_cols > 0) then
            allocate(s_mix_work(n_basis_mix, n_basis_mix), ipiv_mix(n_basis_mix))
            s_mix_work(:, :) = s_mix(:, :)
            call zgesv(n_basis_mix, nocc_mix_cols, s_mix_work, n_basis_mix, ipiv_mix, coef_mix_metric, n_basis_mix, info_lapack)
            if (info_lapack /= 0) then
              if (dg_frag%id == 0) then
                write(*,'(1x,a,i0,a,i0)') &
                  ' [FATAL] rho_mix metric_consistent solve failed for ispin=', &
                  ispin, ' zgesv info=', info_lapack
                flush(6)
              end if
              stop "DG-Fragment RT: rho_mix metric_consistent solve failed"
            else
              coef_mix_eff(1:n_basis_mix, 1:nocc_mix_cols) = coef_mix_metric(1:n_basis_mix, 1:nocc_mix_cols)
            end if
          end if
        else
          allocate(coef_mix_eff(n_basis_mix, max(1, nocc_mix_cols)))
          coef_mix_eff(:, :) = (0.0d0, 0.0d0)
          if (need_full_coef_mix_spin) then
            coef_mix_eff(1:n_basis_mix, 1:nocc_spin) = dg_frag%coef_mix(1:n_basis_mix, 1:nocc_spin, ispin)
          else if (nocc_loc > 0) then
            coef_mix_eff(1:n_basis_mix, 1:nocc_loc) = dg_frag%coef_mix(1:n_basis_mix, io_s_frag:io_e_frag, ispin)
          end if
        end if

        if (need_full_coef_mix_spin) then
          coef_mix_spin(1:n_basis_mix, 1:nocc_spin, ispin) = coef_mix_eff(1:n_basis_mix, 1:nocc_spin)
        end if

        ! Hoist occupations once per spin.  The density-matrix and grid
        ! reconstruction loops below then use contiguous cache slices instead of
        ! repeatedly probing system%rocc in state/block hot paths.
        occ_cache(1:nocc_cache) = 0.0d0
        occ_cache(1:nocc_spin) = 1.0d0
        if (allocated(system%rocc)) then
          do io = 1, nocc_spin
            if (io <= size(system%rocc, 1) .and. ispin <= size(system%rocc, 3)) then
              occ_cache(io) = max(0.0d0, system%rocc(io, 1, ispin))
            end if
          end do
        end if

        ! Build the production density matrix from only the local occupied-state
        ! columns.  Diagnostic modes keep a full coefficient copy above, but the
        ! density path still contributes through this local slice and reduces on
        ! icomm_frag.
        density_mix_partial(:, :) = (0.0d0, 0.0d0)
        do io = io_s_frag, io_e_frag
          io_loc = merge(io, io - io_s_frag + 1, need_full_coef_mix_spin)
          occ_factor = occ_cache(io)
          if (occ_factor <= 0.0d0) cycle
          density_mix_partial(1:n_basis_mix, 1:n_basis_mix) = density_mix_partial(1:n_basis_mix, 1:n_basis_mix) + &
            occ_factor * matmul(coef_mix_eff(1:n_basis_mix, io_loc:io_loc), &
                                 conjg(transpose(coef_mix_eff(1:n_basis_mix, io_loc:io_loc))))
        end do
        if (dg_frag%isize_frag > 1 .and. dg_frag%icomm_frag /= COMM_GROUP_NULL) then
          call comm_summation(density_mix_partial(1:n_basis_mix, 1:n_basis_mix), &
                              density_mix(1:n_basis_mix, 1:n_basis_mix, ispin), &
                              n_basis_mix * n_basis_mix, dg_frag%icomm_frag)
        else
          density_mix(1:n_basis_mix, 1:n_basis_mix, ispin) = density_mix_partial(1:n_basis_mix, 1:n_basis_mix)
        end if

      end do
    end if
    ! Plane waves are normalized over the full periodic simulation cell.
    ! mg%num is the MPI-local real-space box size, so using it here makes
    ! PW density depend on the parent/grid or fragment-orbital decomposition.
    boxL(1) = dg_frag%hgs(1) * real(dg_frag%lgnum_total(1), 8)
    boxL(2) = dg_frag%hgs(2) * real(dg_frag%lgnum_total(2), 8)
    boxL(3) = dg_frag%hgs(3) * real(dg_frag%lgnum_total(3), 8)
    inv_sqrt_vol = 1.0d0 / sqrt(max(1.0d-16, boxL(1) * boxL(2) * boxL(3)))
    if (.not. allocated(dg_frag%density_block_nblocks)) then
      allocate(dg_frag%density_block_nblocks(ifrag_count), dg_frag%density_block_first_offset(ifrag_count), &
               dg_frag%density_block_step(ifrag_count))
      block_idx_global = 0
      i_local = 0
      do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
        i_local = i_local + 1
        ngrid = dg_frag%density_grid_point_count(i_local)
        nblocks_ifrag = (ngrid + grid_block_size - 1) / grid_block_size
        dg_frag%density_block_nblocks(i_local) = nblocks_ifrag
        dg_frag%density_block_first_offset(i_local) = 0
        dg_frag%density_block_step(i_local) = 1
        block_idx_global = block_idx_global + nblocks_ifrag
      end do
    end if

    if (n_pw > 0) then
      need_pw_cache_alloc = (.not. allocated(dg_frag%coef_pw_full_cache))
      need_pw_cache_expand = (.not. need_pw_cache_alloc) .and. dg_frag%coef_pw_full_cache_nstate < nocc_cache
      if (need_pw_cache_alloc .or. need_pw_cache_expand) then
        call refresh_pw_coef_cache(dg_frag, nocc_cache)
        rebuilt_pw_cache = .true.
      end if
    end if
    if (n_pw == 0) then
      allocate(D_partial_re(nbf_max, nbf_max))
      allocate(D_frag_re(nbf_max, nbf_max, system%nspin))
      allocate(rho_dmat_q_local(grid_block_size, nbf_max))
      allocate(rho_dmat_q_full(grid_block_size, nbf_max))
      if (enable_wannier_density) then
        if (use_buffer_wannier_density) then
          nkeep_max_density = max(1, maxval(dg_frag%buffer_wannier_nkeep(1:ifrag_count)))
        else
          nkeep_max_density = max(1, maxval(dg_frag%local_wannier_nkeep(1:ifrag_count)))
        end if
        allocate(D_wannier_re(nkeep_max_density, nkeep_max_density, system%nspin))
        allocate(D_wannier_tmp(nbf_max, nkeep_max_density))
        allocate(wannier_coef_owned(nbf_max, nkeep_max_density, system%nspin))
        allocate(wannier_blk_local(grid_block_size, nkeep_max_density))
        allocate(wannier_blk_full(grid_block_size, nkeep_max_density))
        allocate(wannier_dmat_q(grid_block_size, nkeep_max_density))
        allocate(wannier_owned_count_spin(system%nspin))
      end if
    end if
    allocate(coef_re_full(nbf_max, max(1, nocc_cache), system%nspin))
    allocate(coef_im_full(nbf_max, max(1, nocc_cache), system%nspin))
    if (n_pw == 0) allocate(occ_sqrt_cache(max(1, nocc_cache)))
    allocate(coef_c_full(nbf_max, max(1, nocc_cache)), coef_c_frag(nbf_max, max(1, nocc_cache)))
    phi_lb1 = lbound(dg_frag%phi_frag, 1)
    phi_lb2 = lbound(dg_frag%phi_frag, 2)
    phi_lb3 = lbound(dg_frag%phi_frag, 3)
    phi_ub1 = ubound(dg_frag%phi_frag, 1)
    phi_ub2 = ubound(dg_frag%phi_frag, 2)
    phi_ub3 = ubound(dg_frag%phi_frag, 3)
    phi_lg1 = dg_frag%lgnum_total(1)
    phi_lg2 = dg_frag%lgnum_total(2)
    phi_lg3 = dg_frag%lgnum_total(3)
    ! Cache fragment basis values on the density grid block layout.  This avoids
    ! remapping global grid coordinates to buffered phi indices on every density
    ! reconstruction step.
    if (enable_density_phi_block_cache) then
      need_phi_cache_alloc = (.not. allocated(dg_frag%density_phi_block_cache))
      need_phi_count_alloc = (.not. allocated(dg_frag%density_phi_block_count))
      need_phi_cache_invalid = (.not. dg_frag%density_phi_block_cache_valid)
      need_phi_cache_resize = dg_frag%density_phi_block_size /= grid_block_size
      if (.not. need_phi_cache_alloc) then
        need_phi_cache_resize = need_phi_cache_resize .or. &
          size(dg_frag%density_phi_block_cache, 2) /= nbf_frag_cap
      end if
      if (need_phi_cache_alloc .or. need_phi_count_alloc .or. need_phi_cache_invalid .or. need_phi_cache_resize) then
        if (allocated(dg_frag%density_phi_block_cache)) deallocate(dg_frag%density_phi_block_cache)
        if (allocated(dg_frag%density_phi_block_count)) deallocate(dg_frag%density_phi_block_count)
        nblocks_max = max(1, maxval(dg_frag%density_block_nblocks))
        allocate(dg_frag%density_phi_block_cache(grid_block_size, nbf_frag_cap, nblocks_max, ifrag_count))
        allocate(dg_frag%density_phi_block_count(ifrag_count))
        dg_frag%density_phi_block_cache(:, :, :, :) = 0.0d0
        dg_frag%density_phi_block_count(:) = 0
        i_local = 0
        do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
          i_local = i_local + 1
          local_grid_count = dg_frag%density_grid_point_count(i_local)
          dg_frag%density_phi_block_count(i_local) = dg_frag%density_block_nblocks(i_local)
          nbf = maxval(dg_frag%n_basis(ifrag, 1:system%nspin))
          if (dg_frag%parallel_mode_orbital .and. dg_frag%isize_frag > 1) then
            nbf_per_rank_loc = (nbf_max + dg_frag%isize_frag - 1) / dg_frag%isize_frag
            ib_s_frag = dg_frag%id_frag * nbf_per_rank_loc + 1
            ib_e_frag = min((dg_frag%id_frag + 1) * nbf_per_rank_loc, nbf)
          else
            ib_s_frag = 1
            ib_e_frag = nbf
          end if
          nbf_frag_count = max(0, ib_e_frag - ib_s_frag + 1)
          nbf_frag_is = 1
          nbf_frag_ie = nbf_frag_count
          do block_cache_idx = 1, dg_frag%density_phi_block_count(i_local)
            igrid0 = 1 + (block_cache_idx - 1) * grid_block_size
            npt_cache = min(grid_block_size, local_grid_count - igrid0 + 1)
!$omp parallel do private(igrid, ixg, iyg, izg, bx, by, bz, ib_loc, ib_global) schedule(static)
            do igrid = 1, npt_cache
              ixg = dg_frag%density_grid_points(igrid0 + igrid - 1, i_local)%ixg
              iyg = dg_frag%density_grid_points(igrid0 + igrid - 1, i_local)%iyg
              izg = dg_frag%density_grid_points(igrid0 + igrid - 1, i_local)%izg
              bx = map_global_to_phi_box_coord_ham(ixg, phi_lb1, phi_ub1, phi_lg1)
              by = map_global_to_phi_box_coord_ham(iyg, phi_lb2, phi_ub2, phi_lg2)
              bz = map_global_to_phi_box_coord_ham(izg, phi_lb3, phi_ub3, phi_lg3)
              if (bx == 0 .or. by == 0 .or. bz == 0) cycle
              do ib_loc = nbf_frag_is, nbf_frag_ie
                ib_global = ib_s_frag + ib_loc - nbf_frag_is
                dg_frag%density_phi_block_cache(igrid, ib_loc, block_cache_idx, i_local) = &
                  dg_frag%phi_frag(bx, by, bz, ib_global, i_local)
              end do
            end do
!$omp end parallel do
          end do
        end do
        dg_frag%density_phi_block_size = grid_block_size
        dg_frag%density_phi_block_cache_valid = .true.
        rebuilt_phi_block_cache = .true.
      end if
    else
      if (allocated(dg_frag%density_phi_block_cache)) deallocate(dg_frag%density_phi_block_cache)
      if (allocated(dg_frag%density_phi_block_count)) deallocate(dg_frag%density_phi_block_count)
      dg_frag%density_phi_block_size = grid_block_size
      dg_frag%density_phi_block_cache_valid = .false.
    end if
    ! Cache PW phases with the same global-index convention used by operator
    ! construction.  n_pw > 0 revisits these phases for every occupied state,
    ! so keeping them per density block removes repeated trig calls.
    if (enable_density_phase_block_cache .and. n_pw > 0) then
      need_phase_cache_alloc = (.not. allocated(dg_frag%density_phase_block_cache))
      need_phase_cache_invalid = (.not. dg_frag%density_phase_block_cache_valid)
      need_phase_cache_resize = dg_frag%density_phase_block_size /= grid_block_size
      need_phase_cache_npw = dg_frag%density_phase_block_npw /= n_pw
      if (need_phase_cache_alloc .or. need_phase_cache_invalid .or. need_phase_cache_resize .or. need_phase_cache_npw) then
        if (allocated(dg_frag%density_phase_block_cache)) deallocate(dg_frag%density_phase_block_cache)
        nblocks_max = max(1, maxval(dg_frag%density_block_nblocks))
        allocate(dg_frag%density_phase_block_cache(grid_block_size, n_pw, nblocks_max, ifrag_count))
        i_local = 0
        do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
          i_local = i_local + 1
          local_grid_count = dg_frag%density_grid_point_count(i_local)
          do block_cache_idx = 1, dg_frag%density_block_nblocks(i_local)
            igrid0 = 1 + (block_cache_idx - 1) * grid_block_size
            npt_cache = min(grid_block_size, local_grid_count - igrid0 + 1)
!$omp parallel do private(igrid, ixg, iyg, izg, ipw, theta) schedule(static)
            do igrid = 1, npt_cache
              ixg = dg_frag%density_grid_points(igrid0 + igrid - 1, i_local)%ixg
              iyg = dg_frag%density_grid_points(igrid0 + igrid - 1, i_local)%iyg
              izg = dg_frag%density_grid_points(igrid0 + igrid - 1, i_local)%izg
!$omp simd private(theta)
              do ipw = 1, n_pw
                theta = kpw_hx(ipw) * real(ixg, 8) + kpw_hy(ipw) * real(iyg, 8) + kpw_hz(ipw) * real(izg, 8)
                dg_frag%density_phase_block_cache(igrid, ipw, block_cache_idx, i_local) = &
                  cmplx(cos(theta), sin(theta), kind=8) * inv_sqrt_vol
              end do
            end do
!$omp end parallel do
          end do
        end do
        dg_frag%density_phase_block_size = grid_block_size
        dg_frag%density_phase_block_npw = n_pw
        dg_frag%density_phase_block_cache_valid = .true.
        rebuilt_phase_block_cache = .true.
      end if
    else
      if (allocated(dg_frag%density_phase_block_cache)) deallocate(dg_frag%density_phase_block_cache)
      dg_frag%density_phase_block_size = grid_block_size
      dg_frag%density_phase_block_npw = 0
      dg_frag%density_phase_block_cache_valid = .false.
    end if
      i_local = 0
      block_idx_global = 0
      do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
        i_local = i_local + 1
        nxyz(1:3) = dg_frag%nxyz_domain(1:3, ifrag)
        ngrid = dg_frag%density_grid_point_count(i_local)
        nblocks_ifrag = dg_frag%density_block_nblocks(i_local)
        first_block_offset = dg_frag%density_block_first_offset(i_local)
        block_step_blocks = dg_frag%density_block_step(i_local)
        valid_basis_count_spin(:) = 0
        owned_basis_count_spin(:) = 0
        basis_gid_spin(:, :) = 0
        valid_basis_ids_spin(:, :) = 0
        owned_basis_ids_spin(:, :) = 0
        owned_coef_local_ids_spin(:, :) = 0
        do ispin = 1, system%nspin
          nbf = dg_frag%n_basis(ifrag, ispin)
          if (nbf <= 0) cycle
          do istate_frag = 1, nbf
            basis_gid(istate_frag) = dg_frag%index_basis(istate_frag, ifrag, ispin)
            if (basis_gid(istate_frag) < 1 .or. basis_gid(istate_frag) > dg_frag%n_mat_max) cycle
            valid_basis_count_spin(ispin) = valid_basis_count_spin(ispin) + 1
            basis_gid_spin(istate_frag, ispin) = basis_gid(istate_frag)
            valid_basis_ids_spin(valid_basis_count_spin(ispin), ispin) = istate_frag
            if (.not. allocated(dg_frag%coef_global_to_local)) cycle
            if (.not. density_basis_owned_by_rank(basis_gid(istate_frag), ispin)) cycle
            coef_local = dg_frag%coef_global_to_local(basis_gid(istate_frag), ispin)
            if (coef_local < 1 .or. coef_local > size(dg_frag%coef, 1)) cycle
            owned_basis_count_spin(ispin) = owned_basis_count_spin(ispin) + 1
            owned_basis_ids_spin(owned_basis_count_spin(ispin), ispin) = istate_frag
            owned_coef_local_ids_spin(owned_basis_count_spin(ispin), ispin) = coef_local
          end do
        end do
        if (use_mixed_density) then
          n_basis_mix_spin(:) = 0
          do ispin = 1, system%nspin
            nbf = dg_frag%n_basis(ifrag, ispin)
            n_basis_mix = min(dg_frag%mixed_basis_dim(ispin), max_mixed_basis)
            n_basis_mix_spin(ispin) = n_basis_mix
            if (nbf <= 0 .or. n_basis_mix <= 0) cycle
            transform_frag_spin(1:nbf, 1:n_basis_mix, ispin) = (0.0d0, 0.0d0)
            do istate_frag = 1, nbf
              ig_i = dg_frag%index_basis(istate_frag, ifrag, ispin)
              if (ig_i < 1 .or. ig_i > n_frag) cycle
              transform_frag_spin(istate_frag, 1:n_basis_mix, ispin) = dg_frag%mixed_transform(ig_i, 1:n_basis_mix, ispin)
            end do
            if (n_pw > 0) then
              transform_pw_spin(1:n_pw, 1:n_basis_mix, ispin) = dg_frag%mixed_transform(n_frag+1:n_tot, 1:n_basis_mix, ispin)
              ! Phase-fixing test: rotate each mode m so max-abs PW component is real-positive
              if (enable_fp_phase_fix) then
                do im_pf = 1, n_basis_mix
                  ipw_max = 1
                  abs_max_pf = abs(transform_pw_spin(1, im_pf, ispin))
                  do ipw_pf = 2, n_pw
                    if (abs(transform_pw_spin(ipw_pf, im_pf, ispin)) > abs_max_pf) then
                      abs_max_pf = abs(transform_pw_spin(ipw_pf, im_pf, ispin))
                      ipw_max = ipw_pf
                    end if
                  end do
                  if (abs_max_pf > 1.0d-30) then
                    phase_fix = conjg(transform_pw_spin(ipw_max, im_pf, ispin)) / abs_max_pf
                    transform_pw_spin(1:n_pw, im_pf, ispin) = &
                      transform_pw_spin(1:n_pw, im_pf, ispin) * phase_fix
                    transform_frag_spin(1:nbf, im_pf, ispin) = &
                      transform_frag_spin(1:nbf, im_pf, ispin) * phase_fix
                  end if
                end do
              end if
            end if
          end do
        end if
        ! Assemble coefficients once per fragment and form the fragment density
        ! matrix D=C*occ*C^H.  The grid path then evaluates
        ! rho(r)=phi(r)^T D phi(r), avoiding an occupied-state loop per block
        ! where that is cheaper than repeatedly forming psi(r).
        if (n_pw == 0) D_frag_re(:, :, :) = 0.0d0
        if (n_pw == 0 .and. enable_wannier_density) then
          D_wannier_re(:, :, :) = 0.0d0
          wannier_coef_owned(:, :, :) = 0.0d0
          wannier_owned_count_spin(:) = 0
        end if
        if (n_pw == 0) then
          coef_re_full(1:nbf_max, 1:nocc_cache, 1:system%nspin) = 0.0d0
          coef_im_full(1:nbf_max, 1:nocc_cache, 1:system%nspin) = 0.0d0
          do ispin = 1, system%nspin
            nbf = dg_frag%n_basis(ifrag, ispin)
            nocc_spin = dg_frag%nocc_spin(ispin)
            if (nbf <= 0 .or. nocc_spin <= 0) cycle
            occ_cache(1:nocc_spin) = 1.0d0
            if (allocated(system%rocc)) then
              do io = 1, nocc_spin
                if (io <= size(system%rocc, 1) .and. ispin <= size(system%rocc, 3)) then
                  occ_cache(io) = max(0.0d0, system%rocc(io, 1, ispin))
                end if
              end do
            end if
            occ_sqrt_cache(1:nocc_spin) = sqrt(occ_cache(1:nocc_spin))
            call assert_real_vector_finite_density('occ_cache', occ_cache, nocc_spin, ifrag, ispin, 0)
            if (dg_frag%coef_state_block_mode) then
              io_s_frag = max(1, dg_frag%coef_state_start)
              io_e_frag = min(nocc_spin, dg_frag%coef_state_end)
            else
              nocc_per_rank_loc = (nocc_spin + dg_frag%isize_frag - 1) / dg_frag%isize_frag
              io_s_frag = dg_frag%id_frag * nocc_per_rank_loc + 1
              io_e_frag = min((dg_frag%id_frag + 1) * nocc_per_rank_loc, nocc_spin)
            end if
            nocc_loc = max(0, io_e_frag - io_s_frag + 1)
            valid_basis_count = valid_basis_count_spin(ispin)
            coef_c_full(1:nbf_max, 1:nocc_spin) = (0.0d0, 0.0d0)
            if (nocc_loc > 0 .and. valid_basis_count > 0) then
              do idx_local = 1, valid_basis_count
                basis_gid(idx_local) = basis_gid_spin(valid_basis_ids_spin(idx_local, ispin), ispin)
              end do
              coef_c_frag(1:valid_basis_count, 1:nocc_loc) = (0.0d0, 0.0d0)
              if (dg_frag%coef_state_block_mode) then
                do idx_local = 1, valid_basis_count
                  local_coef_idx = dg_frag%coef_global_to_local(basis_gid(idx_local), ispin)
                  if (local_coef_idx < 1 .or. local_coef_idx > size(dg_frag%coef, 1)) cycle
                  do io_loc = 1, nocc_loc
                    local_state_col = io_s_frag + io_loc - 1 - dg_frag%coef_state_start + 1
                    if (local_state_col < 1 .or. local_state_col > size(dg_frag%coef, 2)) cycle
                    coef_c_frag(idx_local, io_loc) = dg_frag%coef(local_coef_idx, local_state_col, ispin)
                  end do
                end do
              else
                call fetch_remote_coef_rows(dg_frag, ispin, basis_gid(1:valid_basis_count), &
                                            coef_c_frag(1:valid_basis_count, 1:nocc_loc), &
                                            io_s_frag, io_e_frag)
              end if
              do io_loc = 1, nocc_loc
                io = io_s_frag + io_loc - 1
                do idx_local = 1, valid_basis_count
                  istate_frag = valid_basis_ids_spin(idx_local, ispin)
                  coef_c_full(istate_frag, io) = coef_c_frag(idx_local, io_loc)
                end do
              end do
            end if
            coef_re_full(1:nbf_max, 1:nocc_spin, ispin) = 0.0d0
            coef_im_full(1:nbf_max, 1:nocc_spin, ispin) = 0.0d0
            if (nocc_loc > 0) then
              coef_re_full(1:nbf, io_s_frag:io_e_frag, ispin) = &
                real(coef_c_full(1:nbf, io_s_frag:io_e_frag), kind=8)
              coef_im_full(1:nbf, io_s_frag:io_e_frag, ispin) = &
                aimag(coef_c_full(1:nbf, io_s_frag:io_e_frag))
            end if
            D_partial_re(1:nbf, 1:nbf) = 0.0d0
            if (nocc_loc > 0 .and. nbf > 0) then
              nbatch = 0
              do io = io_s_frag, io_e_frag
                nbatch = nbatch + 1
                coef_blk_re(1:nbf, nbatch) = occ_sqrt_cache(io) * coef_re_full(1:nbf, io, ispin)
                coef_blk_im(1:nbf, nbatch) = occ_sqrt_cache(io) * coef_im_full(1:nbf, io, ispin)
                if (nbatch == state_block_size .or. io == io_e_frag) then
                  call dsyrk('U', 'N', nbf, nbatch, 1.0d0, coef_blk_re, nbf_max, 1.0d0, D_partial_re, nbf_max)
                  call dsyrk('U', 'N', nbf, nbatch, 1.0d0, coef_blk_im, nbf_max, 1.0d0, D_partial_re, nbf_max)
                  nbatch = 0
                end if
              end do
              call assert_real_matrix_finite_density('D_partial_re', D_partial_re, &
                1, nbf, 1, nbf, ifrag, ispin, 0)
            end if
            if (dg_frag%isize_frag > 1 .and. dg_frag%icomm_frag /= COMM_GROUP_NULL) then
              call comm_summation(D_partial_re(1:nbf, 1:nbf), D_frag_re(1:nbf, 1:nbf, ispin), &
                                  nbf * nbf, dg_frag%icomm_frag)
            else
              D_frag_re(1:nbf, 1:nbf, ispin) = D_partial_re(1:nbf, 1:nbf)
            end if
            call assert_real_matrix_finite_density('D_frag_re-after-reduce', D_frag_re(:, :, ispin), &
              1, nbf, 1, nbf, ifrag, ispin, 0)
            do io = 1, nbf
              do istate_frag = io + 1, nbf
                D_frag_re(istate_frag, io, ispin) = D_frag_re(io, istate_frag, ispin)
              end do
            end do
            call assert_real_matrix_finite_density('D_frag_re-after-sym', D_frag_re(:, :, ispin), &
              1, nbf, 1, nbf, ifrag, ispin, 0)
            if (enable_wannier_density) then
              if (use_buffer_wannier_density) then
                if (i_local >= 1 .and. i_local <= size(dg_frag%buffer_wannier_nkeep)) then
                  nkeep = min(dg_frag%buffer_wannier_nkeep(i_local), size(dg_frag%buffer_wannier_coef, 2))
                else
                  nkeep = 0
                end if
              else if (i_local >= 1 .and. i_local <= size(dg_frag%local_wannier_nkeep)) then
                nkeep = min(dg_frag%local_wannier_nkeep(i_local), size(dg_frag%local_wannier_coef, 2))
              else
                nkeep = 0
              end if
              nw_owned = 0
              if (nkeep > 0) then
                do iw = 1, nkeep
                  if ((.not. use_buffer_wannier_density) .and. allocated(dg_frag%local_wannier_owned)) then
                    if (.not. dg_frag%local_wannier_owned(iw, ispin, i_local)) cycle
                  end if
                  nw_owned = nw_owned + 1
                  if (use_buffer_wannier_density) then
                    wannier_coef_owned(1:nbf, nw_owned, ispin) = &
                      dg_frag%buffer_wannier_coef(1:nbf, iw, ispin, i_local)
                  else
                    wannier_coef_owned(1:nbf, nw_owned, ispin) = &
                      dg_frag%local_wannier_coef(1:nbf, iw, ispin, i_local)
                  end if
                end do
              end if
              wannier_owned_count_spin(ispin) = nw_owned
              if (nw_owned > 0) then
                D_wannier_tmp(1:nbf, 1:nw_owned) = 0.0d0
                call dgemm('N', 'N', nbf, nw_owned, nbf, 1.0d0, D_frag_re(1, 1, ispin), nbf_max, &
                           wannier_coef_owned(1, 1, ispin), nbf_max, 0.0d0, D_wannier_tmp, nbf_max)
                call dgemm('T', 'N', nw_owned, nw_owned, nbf, 1.0d0, wannier_coef_owned(1, 1, ispin), nbf_max, &
                           D_wannier_tmp, nbf_max, 0.0d0, D_wannier_re(1, 1, ispin), nkeep_max_density)
                call assert_real_matrix_finite_density('D_wannier_re', D_wannier_re(:, :, ispin), &
                  1, nw_owned, 1, nw_owned, ifrag, ispin, 0)
              end if
            end if
          end do
        else
          ! n_pw > 0: upfront subgroup assembly of coef_re/im on all icomm_frag ranks
          coef_re_full(1:nbf_max, 1:nocc_cache, 1:system%nspin) = 0.0d0
          coef_im_full(1:nbf_max, 1:nocc_cache, 1:system%nspin) = 0.0d0
          do ispin = 1, system%nspin
            nbf = dg_frag%n_basis(ifrag, ispin)
            nocc_spin = dg_frag%nocc_spin(ispin)
            if (nbf <= 0 .or. nocc_spin <= 0) cycle
            valid_basis_count = owned_basis_count_spin(ispin)
            coef_c_full(1:nbf_max, 1:nocc_spin) = (0.0d0, 0.0d0)
!$omp parallel do private(io, idx_local, istate_frag, io_loc) schedule(static)
            do io = 1, nocc_spin
              do idx_local = 1, valid_basis_count
                istate_frag = owned_basis_ids_spin(idx_local, ispin)
                io_loc = owned_coef_local_ids_spin(idx_local, ispin)
                coef_c_full(istate_frag, io) = dg_frag%coef(io_loc, io, ispin)
              end do
            end do
!$omp end parallel do
            if (dg_frag%isize_frag > 1 .and. dg_frag%icomm_frag /= COMM_GROUP_NULL) then
              coef_c_frag(1:nbf_max, 1:nocc_spin) = (0.0d0, 0.0d0)
              call comm_summation(coef_c_full(1:nbf, 1:nocc_spin), coef_c_frag(1:nbf, 1:nocc_spin), &
                                  nbf * nocc_spin, dg_frag%icomm_frag)
              coef_c_full(1:nbf_max, 1:nocc_spin) = coef_c_frag(1:nbf_max, 1:nocc_spin)
            end if
            coef_re_full(1:nbf_max, 1:nocc_spin, ispin) = real(coef_c_full(1:nbf_max, 1:nocc_spin), kind=8)
            coef_im_full(1:nbf_max, 1:nocc_spin, ispin) = aimag(coef_c_full(1:nbf_max, 1:nocc_spin))
          end do
          ! NOTE: coef_pw_full_cache is already populated on all ranks by refresh_pw_coef_cache
          ! (called before the fragment loop), so no bcast is needed here.

          ! Compute state range for this rank
          nocc_per_rank_loc = (nocc_cache + dg_frag%isize_frag - 1) / dg_frag%isize_frag
          io_s_frag = dg_frag%id_frag * nocc_per_rank_loc + 1
          io_e_frag = min((dg_frag%id_frag + 1) * nocc_per_rank_loc, nocc_cache)
        end if
        ! --- end D pre-pass ---
        do block_offset = first_block_offset, nblocks_ifrag - 1, block_step_blocks
          igrid0 = 1 + block_offset * grid_block_size
          npt_blk = min(grid_block_size, ngrid - igrid0 + 1)
          if (npt_blk <= 0 .or. ngrid <= 0) then
            write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0)') "[FATAL] density block size invalid: rank=", dg_frag%id, &
              " ifrag=", ifrag, " i_local=", i_local, " block_offset=", block_offset, " ngrid=", ngrid, " npt_blk=", npt_blk
            flush(6)
            stop "DG-Fragment RT: density block size invalid"
          end if
          density_on_frag_root = dg_frag%parallel_mode_orbital .and. dg_frag%isize_frag > 1 .and. &
            dg_frag%icomm_frag /= COMM_GROUP_NULL
          local_grid_count = 0
          remote_grid_count = 0
          valid_remote_grid_count = 0
          if (allocated(dg_frag%density_send_slot_map)) then
            call prepare_grid_buffers_owner_map(i_local, igrid0, npt_blk, nxyz, n_pw == 0 .and. .not. dg_frag%parallel_mode_orbital)
          else
            slot_buf(1:npt_blk) = 0
            call prepare_grid_buffers_owner_map_no_slot(i_local, igrid0, npt_blk, nxyz)
          end if
          do igrid = 1, npt_blk
            if (owner_buf(igrid) < 0) cycle
            local_grid_count = local_grid_count + 1
            if (local_grid_count > grid_block_size) then
              write(*,'(1x,a,i0,a,i0,a,i0,a,i0)') "[FATAL] density local_grid_count overflow: rank=", dg_frag%id, &
                " ifrag=", ifrag, " local_grid_count=", local_grid_count, " npt_blk=", npt_blk
              flush(6)
              stop "DG-Fragment RT: density local_grid_count overflow"
            end if
            local_grid_ids(local_grid_count) = igrid
            if (slot_buf(igrid) > 0) then
              remote_grid_count = remote_grid_count + 1
              if (remote_grid_count > grid_block_size) then
                write(*,'(1x,a,i0,a,i0,a,i0,a,i0)') "[FATAL] density remote_grid_count overflow: rank=", dg_frag%id, &
                  " ifrag=", ifrag, " remote_grid_count=", remote_grid_count, " npt_blk=", npt_blk
                flush(6)
                stop "DG-Fragment RT: density remote_grid_count overflow"
              end if
              remote_grid_ids(remote_grid_count) = igrid
              valid_remote_grid_count = valid_remote_grid_count + 1
              if (valid_remote_grid_count > grid_block_size) then
                write(*,'(1x,a,i0,a,i0,a,i0,a,i0)') "[FATAL] density valid_remote_grid_count overflow: rank=", dg_frag%id, &
                  " ifrag=", ifrag, " valid_remote=", valid_remote_grid_count, " npt_blk=", npt_blk
                flush(6)
                stop "DG-Fragment RT: density valid_remote_grid_count overflow"
              end if
              valid_remote_grid_ids(valid_remote_grid_count) = igrid
            end if
          end do

          if (n_pw > 0) then
            if (enable_density_phase_block_cache) then
              phase_cache(1:npt_blk, 1:n_pw) = &
                dg_frag%density_phase_block_cache(1:npt_blk, 1:n_pw, block_offset + 1, i_local)
            else
!$omp parallel do private(ixg, iyg, izg, ipw, theta) schedule(static)
              do igrid = 1, npt_blk
                ixg = ixg_buf(igrid)
                iyg = iyg_buf(igrid)
                izg = izg_buf(igrid)
!$omp simd
                do ipw = 1, n_pw
                  theta = kpw_hx(ipw) * real(ixg, 8) + kpw_hy(ipw) * real(iyg, 8) + kpw_hz(ipw) * real(izg, 8)
                  phase_cache(igrid, ipw) = cmplx(cos(theta), sin(theta), kind=8) * inv_sqrt_vol
                end do
              end do
!$omp end parallel do
            end if
          end if

          do ispin = 1, system%nspin
            nocc_spin = dg_frag%nocc_spin(ispin)
            nbf = dg_frag%n_basis(ifrag, ispin)
            if (nbf <= 0 .or. nocc_spin <= 0) cycle
            valid_basis_count = valid_basis_count_spin(ispin)
            if (dg_frag%parallel_mode_orbital .and. dg_frag%isize_frag > 1) then
              nbf_per_rank_loc = (nbf_max + dg_frag%isize_frag - 1) / dg_frag%isize_frag
              ib_s_frag = dg_frag%id_frag * nbf_per_rank_loc + 1
              ib_e_frag = min((dg_frag%id_frag + 1) * nbf_per_rank_loc, nbf)
            else
              ib_s_frag = 1
              ib_e_frag = nbf
            end if
            nbf_frag_count = max(0, ib_e_frag - ib_s_frag + 1)
            nbf_frag_is = 1
            nbf_frag_ie = nbf_frag_count
          if (enable_density_phi_block_cache) then
            if (nbf_frag_count > 0) then
              phi_blk(1:npt_blk, nbf_frag_is:nbf_frag_ie) = &
                dg_frag%density_phi_block_cache(1:npt_blk, nbf_frag_is:nbf_frag_ie, block_offset + 1, i_local)
            end if
          else
            if (nbf_frag_count > 0) phi_blk(1:npt_blk, nbf_frag_is:nbf_frag_ie) = 0.0d0
!$omp parallel do private(igrid, ixg, iyg, izg, bx, by, bz, ib_loc, ib_global) schedule(static)
            do igrid = 1, npt_blk
              ixg = ixg_buf(igrid)
              iyg = iyg_buf(igrid)
              izg = izg_buf(igrid)
              bx = map_global_to_phi_box_coord_ham(ixg, phi_lb1, phi_ub1, phi_lg1)
              by = map_global_to_phi_box_coord_ham(iyg, phi_lb2, phi_ub2, phi_lg2)
              bz = map_global_to_phi_box_coord_ham(izg, phi_lb3, phi_ub3, phi_lg3)
              if (bx == 0 .or. by == 0 .or. bz == 0) cycle
              do ib_loc = nbf_frag_is, nbf_frag_ie
                ib_global = ib_s_frag + ib_loc - nbf_frag_is
                phi_blk(igrid, ib_loc) = dg_frag%phi_frag(bx, by, bz, ib_global, i_local)
              end do
            end do
!$omp end parallel do
          end if
          if (nbf_frag_count > 0) then
            call assert_real_matrix_finite_density('phi_blk', phi_blk, &
              1, npt_blk, nbf_frag_is, nbf_frag_ie, ifrag, ispin, block_offset)
          end if

            if (use_mixed_density) then
              n_basis_mix = n_basis_mix_spin(ispin)
              if (n_basis_mix <= 0) cycle
              rho_blk(1:npt_blk) = 0.0d0
              if (density_on_frag_root) then
                ! In orbital mode, avoid reducing the full mixed-basis block
                ! (npt_blk x n_basis_mix).  Each rank projects only its local
                ! fragment-basis columns into occupied-state amplitudes, then
                ! the fragment root reduces the smaller psi block and forms rho.
                occ_cache(1:nocc_cache) = 0.0d0
                occ_cache(1:nocc_spin) = 1.0d0
                if (allocated(system%rocc)) then
                  do io = 1, nocc_spin
                    if (io <= size(system%rocc, 1) .and. ispin <= size(system%rocc, 3)) then
                      occ_cache(io) = max(0.0d0, system%rocc(io, 1, ispin))
                    end if
                  end do
                end if
                do io0 = 1, nocc_spin, state_block_size
                  nbatch = min(state_block_size, nocc_spin - io0 + 1)
                  psi_blk_re(1:npt_blk, 1:nbatch) = 0.0d0
                  psi_blk_im(1:npt_blk, 1:nbatch) = 0.0d0
                  if (nbf_frag_count > 0) then
                    coef_c_frag(1:nbf_frag_count, 1:nbatch) = matmul( &
                      transform_frag_spin(ib_s_frag:ib_e_frag, 1:n_basis_mix, ispin), &
                      coef_mix_spin(1:n_basis_mix, io0:io0+nbatch-1, ispin))
                    coef_blk_ri(nbf_frag_is:nbf_frag_ie, 1:nbatch) = real(coef_c_frag(1:nbf_frag_count, 1:nbatch), kind=8)
                    coef_blk_ri(nbf_frag_is:nbf_frag_ie, nbatch+1:2*nbatch) = aimag(coef_c_frag(1:nbf_frag_count, 1:nbatch))
                    call dgemm('N', 'N', npt_blk, 2*nbatch, nbf_frag_count, 1.0d0, phi_blk, grid_block_size, &
                               coef_blk_ri, nbf_max, 0.0d0, psi_blk_ri, grid_block_size)
                    psi_blk_re(1:npt_blk, 1:nbatch) = psi_blk_ri(1:npt_blk, 1:nbatch)
                    psi_blk_im(1:npt_blk, 1:nbatch) = psi_blk_ri(1:npt_blk, nbatch+1:2*nbatch)
                  end if
                  density_mix_tmp(1:npt_blk, 1:nbatch) = cmplx( &
                    psi_blk_re(1:npt_blk, 1:nbatch), psi_blk_im(1:npt_blk, 1:nbatch), kind=8)
                  call comm_summation(density_mix_tmp(1:npt_blk, 1:nbatch), &
                    basis_mix_blk(1:npt_blk, 1:nbatch), npt_blk * nbatch, dg_frag%icomm_frag, 0)
                  if (dg_frag%is_frag_root) then
                    density_mix_tmp(1:npt_blk, 1:nbatch) = basis_mix_blk(1:npt_blk, 1:nbatch)
                    if (n_pw > 0) then
                      do ipw0 = 1, n_pw, pw_block_size
                        npw_blk = min(pw_block_size, n_pw - ipw0 + 1)
                        coef_pw_blk(1:npw_blk, 1:nbatch) = matmul( &
                          transform_pw_spin(ipw0:ipw0+npw_blk-1, 1:n_basis_mix, ispin), &
                          coef_mix_spin(1:n_basis_mix, io0:io0+nbatch-1, ispin))
                        call zgemm('N', 'N', npt_blk, nbatch, npw_blk, zone, phase_cache(1, ipw0), grid_block_size, &
                                   coef_pw_blk, pw_block_size, zone, density_mix_tmp, grid_block_size)
                      end do
                    end if
                    occ_blk(1:nbatch) = occ_cache(io0:io0+nbatch-1)
                    do io = 1, nbatch
                      occ_factor = occ_blk(io)
                      if (occ_factor <= 0.0d0) cycle
!$omp parallel do private(igrid) schedule(static)
                      do igrid = 1, npt_blk
                        rho_blk(igrid) = rho_blk(igrid) + occ_factor * &
                          (real(density_mix_tmp(igrid, io), kind=8)**2 + aimag(density_mix_tmp(igrid, io))**2)
                      end do
!$omp end parallel do
                    end do
                  end if
                end do
              else
                basis_mix_blk(1:npt_blk, 1:n_basis_mix) = (0.0d0, 0.0d0)
                if (nbf_frag_count > 0) then
                  basis_mix_blk(1:npt_blk, 1:n_basis_mix) = matmul(phi_blk(1:npt_blk, nbf_frag_is:nbf_frag_ie), &
                    transform_frag_spin(ib_s_frag:ib_e_frag, 1:n_basis_mix, ispin))
                end if
                if (n_pw > 0) then
                  call zgemm('N', 'N', npt_blk, n_basis_mix, n_pw, zone, phase_cache, grid_block_size, &
                    transform_pw_spin(1, 1, ispin), n_pw, zone, basis_mix_blk, grid_block_size)
                end if
                  occ_cache(1:nocc_cache) = 0.0d0
                  occ_cache(1:nocc_spin) = 1.0d0
                  if (allocated(system%rocc)) then
                    do io = 1, nocc_spin
                      if (io <= size(system%rocc, 1) .and. ispin <= size(system%rocc, 3)) then
                        occ_cache(io) = max(0.0d0, system%rocc(io, 1, ispin))
                      end if
                    end do
                  end if
                  do io0 = 1, nocc_spin, state_block_size
                    nbatch = min(state_block_size, nocc_spin - io0 + 1)
                    call zgemm('N', 'N', npt_blk, nbatch, n_basis_mix, zone, basis_mix_blk, grid_block_size, &
                      coef_mix_spin(1, io0, ispin), max_mixed_basis, zzero, pw_tmp_z, grid_block_size)
                    occ_blk(1:nbatch) = occ_cache(io0:io0+nbatch-1)
                    do io = 1, nbatch
                      occ_factor = occ_blk(io)
                      if (occ_factor <= 0.0d0) cycle
!$omp parallel do private(igrid) schedule(static)
                      do igrid = 1, npt_blk
                        rho_blk(igrid) = rho_blk(igrid) + occ_factor * &
                          (real(pw_tmp_z(igrid, io), kind=8)**2 + aimag(pw_tmp_z(igrid, io))**2)
                      end do
!$omp end parallel do
                    end do
                  end do
              end if
              if (dg_frag%isize_frag > 1 .and. dg_frag%icomm_frag /= COMM_GROUP_NULL .and. &
                  .not. dg_frag%parallel_mode_orbital) then
                call comm_summation(rho_blk(1:npt_blk), rho_blk_reduced(1:npt_blk), npt_blk, dg_frag%icomm_frag)
                rho_blk(1:npt_blk) = rho_blk_reduced(1:npt_blk)
              end if
              call assert_real_vector_finite_density('rho_blk', rho_blk, npt_blk, ifrag, ispin, block_offset)


              ! Mixed orbital block partitioning reduces rho_blk only to the
              ! fragment root; legacy reductions keep the previous all-rank
              ! result.  In both cases only the root emits the final density.
              if (dg_frag%is_frag_root) then
!$omp parallel do private(igrid, ixg, iyg, izg, bx, by, bz, rho_raw_contrib, rho_contrib) schedule(static)
                do igrid = 1, npt_blk
                  ixg = ixg_buf(igrid)
                  iyg = iyg_buf(igrid)
                  izg = izg_buf(igrid)
                  if (slot_buf(igrid) > 0) cycle
                  if (.not. target_rank_owned_by_handler(owner_buf(igrid))) cycle
                    if (ixg < rho_s_x_lo .or. ixg > rho_s_x_hi .or. &
                      iyg < rho_s_y_lo .or. iyg > rho_s_y_hi .or. &
                      izg < rho_s_z_lo .or. izg > rho_s_z_hi) cycle
                  rho_raw_contrib = rho_blk(igrid)
                  rho_contrib = rho_raw_contrib
                  bx = rhobf_bx_buf(igrid)
                  by = rhobf_by_buf(igrid)
                  bz = rhobf_bz_buf(igrid)
                  if (bx == 0 .or. by == 0 .or. bz == 0) then
                    write(*,'(1x,a,i0,a,i0,a,3(i0,1x),a,3(i0,1x),a,3(i0,1x))') &
                      "[FATAL] density local rho_bf map failed: rank=", dg_frag%id, &
                      " id_frag=", dg_frag%id_frag, " idx=", ixg, iyg, izg, &
                      " lb=", lbound(rho_bf, 1), lbound(rho_bf, 2), lbound(rho_bf, 3), &
                      " ub=", ubound(rho_bf, 1), ubound(rho_bf, 2), ubound(rho_bf, 3)
                    flush(6)
                    stop "DG-Fragment RT: density local rho_bf map failed"
                  end if
!$omp atomic update
                  rho_bf(bx, by, bz) = rho_bf(bx, by, bz) + rho_contrib
                  call assert_rho_materialize_value('rho_bf-local-mixed', ifrag, ispin, block_offset, &
                    igrid, bx, by, bz, ixg, iyg, izg, owner_buf(igrid), slot_buf(igrid), &
                    rho_contrib, rho_bf(bx, by, bz))
                  if (system%nspin > 1) then
!$omp atomic update
                    rho_s_bf(ixg, iyg, izg, ispin) = rho_s_bf(ixg, iyg, izg, ispin) + rho_contrib
                  end if
                end do
!$omp end parallel do
              end if
              if (dg_frag%is_frag_root) then
!$omp parallel do private(idx_remote, igrid, owner_rank, slot, rho_contrib, spin_offset) schedule(static)
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
                  if (system%nspin > 1) then
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
                    rho_send(owner_rank)%f(spin_offset + slot, 1, 1) = &
                      rho_send(owner_rank)%f(spin_offset + slot, 1, 1) + rho_contrib
                  end if
                end do
!$omp end parallel do
              end if
            else
              ! D already computed in the pre-pass for n_pw == 0.
              rho_blk_accum(1:npt_blk) = 0.0d0
              if (n_pw == 0) then
                if (.not. allocated(rho_blk_partial)) allocate(rho_blk_partial(grid_block_size))
                rho_blk_partial(1:npt_blk) = 0.0d0
                if (enable_wannier_density) then
                  nw_owned = wannier_owned_count_spin(ispin)
                  rho_blk_accum(1:npt_blk) = 0.0d0
                  if (nw_owned > 0) then
                    wannier_blk_local(1:npt_blk, 1:nw_owned) = 0.0d0
                    if (nbf_frag_count > 0) then
                      call dgemm('N', 'N', npt_blk, nw_owned, nbf_frag_count, 1.0d0, phi_blk, grid_block_size, &
                                 wannier_coef_owned(ib_s_frag, 1, ispin), nbf_max, &
                                 0.0d0, wannier_blk_local, grid_block_size)
                    end if
                    if (density_on_frag_root) then
                      call comm_summation(wannier_blk_local(1:npt_blk, 1:nw_owned), &
                                          wannier_blk_full(1:npt_blk, 1:nw_owned), &
                                          npt_blk * nw_owned, dg_frag%icomm_frag, 0)
                    else if (dg_frag%isize_frag > 1 .and. dg_frag%icomm_frag /= COMM_GROUP_NULL) then
                      call comm_summation(wannier_blk_local(1:npt_blk, 1:nw_owned), &
                                          wannier_blk_full(1:npt_blk, 1:nw_owned), &
                                          npt_blk * nw_owned, dg_frag%icomm_frag)
                    else
                      wannier_blk_full(1:npt_blk, 1:nw_owned) = wannier_blk_local(1:npt_blk, 1:nw_owned)
                    end if
                    if ((.not. density_on_frag_root) .or. dg_frag%is_frag_root) then
                      call assert_real_matrix_finite_density('wannier_blk_full', wannier_blk_full, &
                        1, npt_blk, 1, nw_owned, ifrag, ispin, block_offset)
                      call dgemm('N', 'N', npt_blk, nw_owned, nw_owned, 1.0d0, wannier_blk_full, grid_block_size, &
                                 D_wannier_re(1, 1, ispin), nkeep_max_density, &
                                 0.0d0, wannier_dmat_q, grid_block_size)
!$omp parallel do private(igrid, iw, rho_accum) schedule(static)
                      do igrid = 1, npt_blk
                        rho_accum = 0.0d0
!$omp simd reduction(+:rho_accum)
                        do iw = 1, nw_owned
                          rho_accum = rho_accum + wannier_blk_full(igrid, iw) * wannier_dmat_q(igrid, iw)
                        end do
                        rho_blk_accum(igrid) = rho_accum
                      end do
!$omp end parallel do
                    end if
                  end if
                else if (dg_frag%parallel_mode_orbital) then
                  rho_dmat_q_local(1:npt_blk, 1:nbf) = 0.0d0
                  if (nbf_frag_count > 0) then
                    call dgemm('N', 'N', npt_blk, nbf, nbf_frag_count, 1.0d0, phi_blk, grid_block_size, &
                               D_frag_re(ib_s_frag, 1, ispin), nbf_max, 0.0d0, rho_dmat_q_local, grid_block_size)
                  end if
                  call assert_real_matrix_finite_density('rho_dmat_q_local', rho_dmat_q_local, &
                    1, npt_blk, 1, nbf, ifrag, ispin, block_offset)
                  if (dg_frag%isize_frag > 1 .and. dg_frag%icomm_frag /= COMM_GROUP_NULL) then
                    call comm_summation(rho_dmat_q_local(1:npt_blk, 1:nbf), rho_dmat_q_full(1:npt_blk, 1:nbf), &
                                        npt_blk * nbf, dg_frag%icomm_frag)
                  else
                    rho_dmat_q_full(1:npt_blk, 1:nbf) = rho_dmat_q_local(1:npt_blk, 1:nbf)
                  end if
                  call assert_real_matrix_finite_density('rho_dmat_q_full', rho_dmat_q_full, &
                    1, npt_blk, 1, nbf, ifrag, ispin, block_offset)

                  rho_blk_partial(1:npt_blk) = 0.0d0
                  if (nbf_frag_count > 0) then
!$omp parallel do private(igrid, ib_global, ib_loc, rho_accum) schedule(static)
                    do igrid = 1, npt_blk
                      rho_accum = 0.0d0
!$omp simd reduction(+:rho_accum)
                      do ib_global = ib_s_frag, ib_e_frag
                        ib_loc = ib_global - ib_s_frag + 1
                        rho_accum = rho_accum + phi_blk(igrid, ib_loc) * rho_dmat_q_full(igrid, ib_global)
                      end do
                      rho_blk_partial(igrid) = rho_accum
                    end do
!$omp end parallel do
                  end if
                  call assert_real_vector_finite_density('rho_blk_partial-before-root-reduce', &
                    rho_blk_partial, npt_blk, ifrag, ispin, block_offset)
                  if (density_on_frag_root) then
                    call comm_summation(rho_blk_partial(1:npt_blk), rho_blk_reduced(1:npt_blk), npt_blk, dg_frag%icomm_frag, 0)
                    if (dg_frag%is_frag_root) then
                      call assert_real_vector_finite_density('rho_blk_reduced-root', &
                        rho_blk_reduced, npt_blk, ifrag, ispin, block_offset)
                    end if
                    if (dg_frag%is_frag_root) then
                      rho_blk_accum(1:npt_blk) = rho_blk_reduced(1:npt_blk)
                    else
                      rho_blk_accum(1:npt_blk) = 0.0d0
                    end if
                  else if (dg_frag%isize_frag > 1 .and. dg_frag%icomm_frag /= COMM_GROUP_NULL) then
                    call comm_summation(rho_blk_partial(1:npt_blk), rho_blk_reduced(1:npt_blk), npt_blk, dg_frag%icomm_frag)
                    rho_blk_accum(1:npt_blk) = rho_blk_reduced(1:npt_blk)
                  else
                    rho_blk_accum(1:npt_blk) = rho_blk_partial(1:npt_blk)
                  end if
                else
                ! In orbital mode, each rank contributes its local basis-column
                ! block to every occupied-state batch.  The root forms rho after
                ! reducing the wavefunction amplitude, preserving cross terms.
                if (density_on_frag_root) then
                  io_s_frag = 1
                  io_e_frag = nocc_spin
                else
                  nocc_per_rank_loc = (nocc_spin + dg_frag%isize_frag - 1) / dg_frag%isize_frag
                  io_s_frag = dg_frag%id_frag * nocc_per_rank_loc + 1
                  io_e_frag = min((dg_frag%id_frag + 1) * nocc_per_rank_loc, nocc_spin)
                end if
                do io0 = io_s_frag, io_e_frag, state_block_size
                  nbatch = min(state_block_size, io_e_frag - io0 + 1)
                  psi_blk_re(1:npt_blk, 1:nbatch) = 0.0d0
                  psi_blk_im(1:npt_blk, 1:nbatch) = 0.0d0
                  if (nbf_frag_count > 0) then
                    coef_blk_ri(nbf_frag_is:nbf_frag_ie, 1:nbatch) = &
                      coef_re_full(ib_s_frag:ib_e_frag, io0:io0+nbatch-1, ispin)
                    coef_blk_ri(nbf_frag_is:nbf_frag_ie, nbatch+1:2*nbatch) = &
                      coef_im_full(ib_s_frag:ib_e_frag, io0:io0+nbatch-1, ispin)
                    call dgemm('N', 'N', npt_blk, 2*nbatch, nbf_frag_count, 1.0d0, phi_blk, grid_block_size, &
                               coef_blk_ri, nbf_max, 0.0d0, psi_blk_ri, grid_block_size)
                    psi_blk_re(1:npt_blk, 1:nbatch) = psi_blk_ri(1:npt_blk, 1:nbatch)
                    psi_blk_im(1:npt_blk, 1:nbatch) = psi_blk_ri(1:npt_blk, nbatch+1:2*nbatch)
                  end if
                  if (density_on_frag_root) then
                    call comm_summation(psi_blk_re(1:npt_blk, 1:nbatch), psi_blk_ri(1:npt_blk, 1:nbatch), &
                      npt_blk * nbatch, dg_frag%icomm_frag, 0)
                    call comm_summation(psi_blk_im(1:npt_blk, 1:nbatch), psi_blk_ri(1:npt_blk, nbatch+1:2*nbatch), &
                      npt_blk * nbatch, dg_frag%icomm_frag, 0)
                    if (dg_frag%is_frag_root) then
                      psi_blk_re(1:npt_blk, 1:nbatch) = psi_blk_ri(1:npt_blk, 1:nbatch)
                      psi_blk_im(1:npt_blk, 1:nbatch) = psi_blk_ri(1:npt_blk, nbatch+1:2*nbatch)
                    end if
                  end if
                  if (density_on_frag_root .and. .not. dg_frag%is_frag_root) cycle
                  occ_blk(1:nbatch) = occ_cache(io0:io0+nbatch-1)
!$omp parallel private(io, igrid, occ_factor)
                    do io = 1, nbatch
                      occ_factor = occ_blk(io)
                      if (occ_factor <= 0.0d0) cycle
!$omp do schedule(static)
                      do igrid = 1, npt_blk
                        rho_blk_partial(igrid) = rho_blk_partial(igrid) + occ_factor * &
                          (psi_blk_re(igrid, io) * psi_blk_re(igrid, io) + &
                           psi_blk_im(igrid, io) * psi_blk_im(igrid, io))
                      end do
!$omp end do
                    end do
!$omp end parallel
                end do
                if (density_on_frag_root) then
                  if (dg_frag%is_frag_root) then
                    rho_blk_accum(1:npt_blk) = rho_blk_partial(1:npt_blk)
                  else
                    rho_blk_accum(1:npt_blk) = 0.0d0
                  end if
                else if (dg_frag%isize_frag > 1 .and. dg_frag%icomm_frag /= COMM_GROUP_NULL) then
                  call comm_summation(rho_blk_partial(1:npt_blk), rho_blk_reduced(1:npt_blk), npt_blk, dg_frag%icomm_frag)
                  rho_blk_accum(1:npt_blk) = rho_blk_reduced(1:npt_blk)
                else
                  rho_blk_accum(1:npt_blk) = rho_blk_partial(1:npt_blk)
                end if
                end if
              else
                ! Mixed fragment/PW density path: reduce the fragment-basis
                ! amplitude before adding the replicated PW contribution.
                ! Occupations are applied here because this path builds rho
                ! directly from |psi|**2.
                if (.not. allocated(rho_blk_partial)) allocate(rho_blk_partial(grid_block_size))
                rho_blk_partial(1:npt_blk) = 0.0d0
                occ_cache(1:nocc_cache) = 0.0d0
                occ_cache(1:nocc_spin) = 1.0d0
                if (allocated(system%rocc)) then
                  do io = 1, nocc_spin
                    if (io <= size(system%rocc, 1) .and. ispin <= size(system%rocc, 3)) then
                      occ_cache(io) = max(0.0d0, system%rocc(io, 1, ispin))
                    end if
                  end do
                end if

                if (density_on_frag_root) then
                  io_s_frag = 1
                  io_e_frag = nocc_spin
                else
                  nocc_per_rank_loc = (nocc_spin + dg_frag%isize_frag - 1) / dg_frag%isize_frag
                  io_s_frag = dg_frag%id_frag * nocc_per_rank_loc + 1
                  io_e_frag = min((dg_frag%id_frag + 1) * nocc_per_rank_loc, nocc_spin)
                end if

                do io0 = io_s_frag, min(io_e_frag, nocc_spin), state_block_size
                  nbatch = min(state_block_size, min(io_e_frag, nocc_spin) - io0 + 1)

                  ! Copy fragment coefficients from the full cached view.  The
                  ! PW coefficients are read from the matching full cache below,
                  ! so the fragment and PW pieces describe the same state block.
                  psi_blk_re(1:npt_blk, 1:nbatch) = 0.0d0
                  psi_blk_im(1:npt_blk, 1:nbatch) = 0.0d0
                  if (nbf_frag_count > 0) then
                    coef_blk_re(nbf_frag_is:nbf_frag_ie, 1:nbatch) = &
                      coef_re_full(ib_s_frag:ib_e_frag, io0:io0+nbatch-1, ispin)
                    coef_blk_im(nbf_frag_is:nbf_frag_ie, 1:nbatch) = &
                      coef_im_full(ib_s_frag:ib_e_frag, io0:io0+nbatch-1, ispin)
                    call dgemm('N', 'N', npt_blk, nbatch, nbf_frag_count, 1.0d0, phi_blk, grid_block_size, &
                               coef_blk_re, nbf_max, 0.0d0, psi_blk_re, grid_block_size)
                    call dgemm('N', 'N', npt_blk, nbatch, nbf_frag_count, 1.0d0, phi_blk, grid_block_size, &
                               coef_blk_im, nbf_max, 0.0d0, psi_blk_im, grid_block_size)
                  end if
                  if (density_on_frag_root) then
                    call comm_summation(psi_blk_re(1:npt_blk, 1:nbatch), psi_blk_ri(1:npt_blk, 1:nbatch), &
                      npt_blk * nbatch, dg_frag%icomm_frag, 0)
                    call comm_summation(psi_blk_im(1:npt_blk, 1:nbatch), psi_blk_ri(1:npt_blk, nbatch+1:2*nbatch), &
                      npt_blk * nbatch, dg_frag%icomm_frag, 0)
                    if (dg_frag%is_frag_root) then
                      psi_blk_re(1:npt_blk, 1:nbatch) = psi_blk_ri(1:npt_blk, 1:nbatch)
                      psi_blk_im(1:npt_blk, 1:nbatch) = psi_blk_ri(1:npt_blk, nbatch+1:2*nbatch)
                    else
                      cycle
                    end if
                  end if
                  do ipw0 = 1, n_pw, pw_block_size
                    npw_blk = min(pw_block_size, n_pw - ipw0 + 1)
                    ! Direct access from the replicated PW coefficient cache
                    ! avoids a per-state broadcast inside the hot grid loop.
                    coef_pw_blk(1:npw_blk, 1:nbatch) = &
                      dg_frag%coef_pw_full_cache(ipw0:ipw0+npw_blk-1, io0:io0+nbatch-1, ispin)

                    call zgemm('N', 'N', npt_blk, nbatch, npw_blk, zone, phase_cache(1, ipw0), grid_block_size, &
                               coef_pw_blk, pw_block_size, zzero, pw_tmp_z, grid_block_size)
                    psi_blk_re(1:npt_blk, 1:nbatch) = psi_blk_re(1:npt_blk, 1:nbatch) + real(pw_tmp_z(1:npt_blk, 1:nbatch), kind=8)
                    psi_blk_im(1:npt_blk, 1:nbatch) = psi_blk_im(1:npt_blk, 1:nbatch) + aimag(pw_tmp_z(1:npt_blk, 1:nbatch))
                  end do
                  occ_blk(1:nbatch) = occ_cache(io0:io0+nbatch-1)
                  do io1 = 1, nbatch, rho_state_block_size
                    nstate = min(rho_state_block_size, nbatch - io1 + 1)
                    if (local_grid_count == npt_blk) then
!$omp parallel private(io, igrid, occ_factor)
                      do io = io1, io1 + nstate - 1
                        occ_factor = occ_blk(io)
                        if (occ_factor <= 0.0d0) cycle
!$omp do schedule(static)
                        do igrid = 1, npt_blk
                          rho_blk_partial(igrid) = rho_blk_partial(igrid) + occ_factor * &
                                      (psi_blk_re(igrid, io) * psi_blk_re(igrid, io) + &
                                       psi_blk_im(igrid, io) * psi_blk_im(igrid, io))
                        end do
!$omp end do
                      end do
!$omp end parallel
                    else
!$omp parallel private(io, idx_local, igrid, occ_factor)
                      do io = io1, io1 + nstate - 1
                        occ_factor = occ_blk(io)
                        if (occ_factor <= 0.0d0) cycle
!$omp do schedule(static)
                        do idx_local = 1, local_grid_count
                          igrid = local_grid_ids(idx_local)
                          rho_blk_partial(igrid) = rho_blk_partial(igrid) + occ_factor * &
                                      (psi_blk_re(igrid, io) * psi_blk_re(igrid, io) + &
                                       psi_blk_im(igrid, io) * psi_blk_im(igrid, io))
                        end do
!$omp end do
                      end do
!$omp end parallel
                    end if
                  end do
                end do  ! io0

                ! AllReduce rho_blk_partial across icomm_frag → rho_blk_accum
                ! Note: comm_summation overwrites rho_blk_accum (does not accumulate into it)
                if (density_on_frag_root) then
                  if (dg_frag%is_frag_root) then
                    rho_blk_accum(1:npt_blk) = rho_blk_partial(1:npt_blk)
                  else
                    rho_blk_accum(1:npt_blk) = 0.0d0
                  end if
                else if (dg_frag%isize_frag > 1 .and. dg_frag%icomm_frag /= COMM_GROUP_NULL) then
                  call comm_summation(rho_blk_partial(1:npt_blk), rho_blk_accum(1:npt_blk), npt_blk, dg_frag%icomm_frag)
                else
                  rho_blk_accum(1:npt_blk) = rho_blk_partial(1:npt_blk)
                end if
              end if
              ! rho_blk_accum now holds the fragment contribution for this grid
              ! block.  The first loop materializes points owned by this handler
              ! into rho_bf/rho_s_bf; the second loop packs points whose parent
              ! grid owner lives on another handler rank.
              ! rho_blk_accum: filled by dgemm-path (n_pw==0) or AllReduce (n_pw>0)
              call assert_real_vector_finite_density('rho_blk_accum', rho_blk_accum, npt_blk, &
                ifrag, ispin, block_offset)
                if ((.not. density_on_frag_root) .or. dg_frag%is_frag_root) then
!$omp parallel private(igrid, ixg, iyg, izg, bx, by, bz, rho_contrib, rho_raw_contrib)
!$omp do schedule(static)
                  do igrid = 1, npt_blk
                    ixg = ixg_buf(igrid)
                    iyg = iyg_buf(igrid)
                    izg = izg_buf(igrid)
                    if (slot_buf(igrid) > 0) cycle
                    if ((.not. density_on_frag_root) .and. .not. target_rank_owned_by_handler(owner_buf(igrid))) cycle
                    if (ixg < rho_s_x_lo .or. ixg > rho_s_x_hi .or. &
                      iyg < rho_s_y_lo .or. iyg > rho_s_y_hi .or. &
                      izg < rho_s_z_lo .or. izg > rho_s_z_hi) cycle
                    rho_raw_contrib = rho_blk_accum(igrid)
                    rho_contrib = rho_raw_contrib
                    bx = rhobf_bx_buf(igrid)
                    by = rhobf_by_buf(igrid)
                    bz = rhobf_bz_buf(igrid)
                    if (bx == 0 .or. by == 0 .or. bz == 0) then
                      write(*,'(1x,a,i0,a,i0,a,3(i0,1x),a,3(i0,1x),a,3(i0,1x))') &
                        "[FATAL] density accum rho_bf map failed: rank=", dg_frag%id, &
                        " id_frag=", dg_frag%id_frag, " idx=", ixg, iyg, izg, &
                        " lb=", lbound(rho_bf, 1), lbound(rho_bf, 2), lbound(rho_bf, 3), &
                        " ub=", ubound(rho_bf, 1), ubound(rho_bf, 2), ubound(rho_bf, 3)
                      flush(6)
                      stop "DG-Fragment RT: density accum rho_bf map failed"
                    end if
!$omp atomic update
                    rho_bf(bx, by, bz) = rho_bf(bx, by, bz) + rho_contrib
                    call assert_rho_materialize_value('rho_bf-local-dmat', ifrag, ispin, block_offset, &
                      igrid, bx, by, bz, ixg, iyg, izg, owner_buf(igrid), slot_buf(igrid), &
                      rho_contrib, rho_bf(bx, by, bz))
                    if (system%nspin > 1) then
!$omp atomic update
                      rho_s_bf(ixg, iyg, izg, ispin) = rho_s_bf(ixg, iyg, izg, ispin) + rho_contrib
                    end if
                  end do
!$omp end do
!$omp single
                  if (dg_frag%is_frag_root) then
                    call assert_rho_bf_with_sources('rho_bf-after-local-materialize', &
                      ifrag, ispin, block_offset, npt_blk)
                  end if
!$omp end single
!$omp end parallel
                end if
                if ((.not. density_on_frag_root) .or. dg_frag%is_frag_root) then
!$omp parallel do private(idx_remote, igrid, owner_rank, slot, rho_contrib, spin_offset, ixg, iyg, izg) schedule(static)
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
                    if (dg_frag%density_send_count(owner_rank) <= 0) then
                      if (owner_rank == dg_frag%id) cycle
                      write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0)') &
                        "DG density accum invalid send count: rank=", dg_frag%id, &
                        " id_frag=", dg_frag%id_frag, " i_local=", i_local, &
                        " igrid=", igrid, " owner=", owner_rank, " nsend=", dg_frag%density_send_count(owner_rank)
                      flush(6)
                      stop "DG-Fragment RT: density accum invalid send count"
                    end if
                    if (.not. allocated(rho_send(owner_rank)%f)) then
                      write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0)') &
                        "DG density accum send buffer missing: rank=", dg_frag%id, &
                        " id_frag=", dg_frag%id_frag, " i_local=", i_local, &
                        " igrid=", igrid, " owner=", owner_rank, " nsend=", dg_frag%density_send_count(owner_rank)
                      flush(6)
                      stop "DG-Fragment RT: density accum send buffer missing"
                    end if
                    if (size(rho_send(owner_rank)%f, 1) < dg_frag%density_send_count(owner_rank)) then
                      write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0)') &
                        "DG density accum send buffer size mismatch: rank=", dg_frag%id, &
                        " id_frag=", dg_frag%id_frag, " i_local=", i_local, &
                        " igrid=", igrid, " owner=", owner_rank, &
                        " send_size=", size(rho_send(owner_rank)%f, 1), " nsend=", dg_frag%density_send_count(owner_rank)
                      flush(6)
                      stop "DG-Fragment RT: density accum send buffer size mismatch"
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
                    if (system%nspin > 1) then
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
                      rho_send(owner_rank)%f(spin_offset + slot, 1, 1) = &
                        rho_send(owner_rank)%f(spin_offset + slot, 1, 1) + rho_contrib
                    end if
                  end do
!$omp end parallel do
                end if
                if (dg_frag%is_frag_root) then
                  call assert_rho_bf_with_sources('rho_bf-after-materialize', &
                    ifrag, ispin, block_offset, npt_blk)
                end if
            end if
          end do
        end do
        density_on_frag_root = dg_frag%parallel_mode_orbital .and. dg_frag%isize_frag > 1 .and. &
          dg_frag%icomm_frag /= COMM_GROUP_NULL
        if ((.not. density_on_frag_root) .or. dg_frag%is_frag_root) then
          call assert_rho_bf_with_sources('rho_bf-after-fragment-loop', &
            ifrag, 0, block_idx_global, max(1, min(grid_block_size, ngrid)))
        end if
        block_idx_global = block_idx_global + nblocks_ifrag
      end do
    density_on_frag_root = dg_frag%parallel_mode_orbital .and. dg_frag%isize_frag > 1 .and. &
      dg_frag%icomm_frag /= COMM_GROUP_NULL
    if ((.not. density_on_frag_root) .or. dg_frag%is_frag_root) then
      call assert_rho_bf_with_sources('rho_bf-after-project-before-comm-v2', &
        max(0, dg_frag%ifrag_group), 0, block_idx_global, max(1, min(grid_block_size, ngrid_max)))
    end if
    if (allocated(D_frag_re)) deallocate(D_frag_re)
    if (allocated(D_partial_re)) deallocate(D_partial_re)
    if (allocated(rho_dmat_q_local)) deallocate(rho_dmat_q_local)
    if (allocated(rho_dmat_q_full)) deallocate(rho_dmat_q_full)
    if (allocated(D_wannier_re)) deallocate(D_wannier_re)
    if (allocated(D_wannier_tmp)) deallocate(D_wannier_tmp)
    if (allocated(wannier_coef_owned)) deallocate(wannier_coef_owned)
    if (allocated(wannier_blk_local)) deallocate(wannier_blk_local)
    if (allocated(wannier_blk_full)) deallocate(wannier_blk_full)
    if (allocated(wannier_dmat_q)) deallocate(wannier_dmat_q)
    if (allocated(wannier_owned_count_spin)) deallocate(wannier_owned_count_spin)
    if (allocated(coef_re_full)) deallocate(coef_re_full)
    if (allocated(coef_im_full)) deallocate(coef_im_full)
    if (allocated(coef_c_full)) deallocate(coef_c_full)
    if (allocated(coef_c_frag)) deallocate(coef_c_frag)
    if (allocated(occ_cache)) deallocate(occ_cache)
    if (allocated(occ_blk)) deallocate(occ_blk)
    if (allocated(occ_sqrt_cache)) deallocate(occ_sqrt_cache)
    if (allocated(rho_blk_partial)) deallocate(rho_blk_partial)
    density_on_frag_root = dg_frag%parallel_mode_orbital .and. dg_frag%isize_frag > 1 .and. &
      dg_frag%icomm_frag /= COMM_GROUP_NULL
    density_exchange_active = (.not. density_on_frag_root) .or. dg_frag%is_frag_root
    if (.not. density_exchange_active) then
      send_total_count = 0
      recv_total_count = 0
    else
      allocate(send_counts(1:dg_frag%isize), recv_counts(1:dg_frag%isize))
      allocate(send_displs(1:dg_frag%isize), recv_displs(1:dg_frag%isize))
      send_counts = 0
      recv_counts = 0
      do irank = 0, dg_frag%isize - 1
        if (allocated(rho_send(irank)%f)) send_counts(irank + 1) = size(rho_send(irank)%f)
        npts = dg_frag%density_recv_map(irank)%npts
        if (npts > 0) recv_counts(irank + 1) = density_payload_count * npts
      end do
      send_displs(1) = 0
      recv_displs(1) = 0
      do irank = 2, dg_frag%isize
        send_displs(irank) = send_displs(irank - 1) + send_counts(irank - 1)
        recv_displs(irank) = recv_displs(irank - 1) + recv_counts(irank - 1)
      end do
      send_total_count = sum(send_counts)
      recv_total_count = sum(recv_counts)
      allocate(send_flat(max(1, send_total_count)), recv_flat(max(1, recv_total_count)))
      do irank = 0, dg_frag%isize - 1
        if (.not. allocated(rho_send(irank)%f)) cycle
        if (send_counts(irank + 1) <= 0) cycle
        call assert_real_vector_finite_density('rho_send-peer-before-pack', &
          rho_send(irank)%f(:, 1, 1), send_counts(irank + 1), irank, 0, 0)
        send_flat(send_displs(irank + 1)+1:send_displs(irank + 1)+send_counts(irank + 1)) = rho_send(irank)%f(:, 1, 1)
        deallocate(rho_send(irank)%f)
      end do
      call exchange_density_sparse(send_flat, send_counts, send_displs, recv_flat, recv_counts, recv_displs)
      do irank = 0, dg_frag%isize - 1
        npts = dg_frag%density_recv_map(irank)%npts
        if (npts <= 0) cycle
        if (recv_counts(irank + 1) /= density_payload_count * npts) then
          write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0)') &
            "[FATAL] density unpack recv size mismatch: rank=", dg_frag%id, &
            " id_frag=", dg_frag%id_frag, " peer=", irank, &
            " recv_count=", recv_counts(irank + 1), " npts=", npts
          flush(6)
          stop "DG-Fragment RT: density unpack recv size mismatch"
        end if
        call assert_real_vector_finite_density('recv_flat-peer-after-exchange', &
          recv_flat(recv_displs(irank + 1)+1:recv_displs(irank + 1)+recv_counts(irank + 1)), &
          recv_counts(irank + 1), irank, 0, 0)
        if (system%nspin == 1) then
          do slot = 1, dg_frag%density_recv_map(irank)%npts
            bx = dg_frag%density_recv_map(irank)%bx(slot)
            by = dg_frag%density_recv_map(irank)%by(slot)
            bz = dg_frag%density_recv_map(irank)%bz(slot)
            if (bx == 0 .or. by == 0 .or. bz == 0) then
              ixg = dg_frag%density_recv_map(irank)%ixg(slot)
              iyg = dg_frag%density_recv_map(irank)%iyg(slot)
              izg = dg_frag%density_recv_map(irank)%izg(slot)
              write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,3(i0,1x),a,3(i0,1x),a,3(i0,1x))') &
                "[FATAL] density unpack rho_bf bounds: rank=", dg_frag%id, &
                " id_frag=", dg_frag%id_frag, " peer=", irank, " slot=", slot, &
                " idx=", ixg, iyg, izg, " lb=", lbound(rho_bf, 1), lbound(rho_bf, 2), lbound(rho_bf, 3), &
                " ub=", ubound(rho_bf, 1), ubound(rho_bf, 2), ubound(rho_bf, 3)
              flush(6)
              stop "DG-Fragment RT: density unpack rho_bf bounds"
            end if
            rho_contrib = recv_flat(recv_displs(irank + 1) + slot)
            rho_bf(bx, by, bz) = rho_bf(bx, by, bz) + rho_contrib
            call assert_rho_materialize_value('rho_bf-unpack-fast', irank, 0, 0, &
              slot, bx, by, bz, dg_frag%density_recv_map(irank)%ixg(slot), &
              dg_frag%density_recv_map(irank)%iyg(slot), dg_frag%density_recv_map(irank)%izg(slot), &
              dg_frag%id, slot, rho_contrib, rho_bf(bx, by, bz))
          end do
        else
          do slot = 1, dg_frag%density_recv_map(irank)%npts
            ixg = dg_frag%density_recv_map(irank)%ixg(slot)
            iyg = dg_frag%density_recv_map(irank)%iyg(slot)
            izg = dg_frag%density_recv_map(irank)%izg(slot)
            bx = dg_frag%density_recv_map(irank)%bx(slot)
            by = dg_frag%density_recv_map(irank)%by(slot)
            bz = dg_frag%density_recv_map(irank)%bz(slot)
            if (bx == 0 .or. by == 0 .or. bz == 0) then
              write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,3(i0,1x),a,3(i0,1x),a,3(i0,1x))') &
                "[FATAL] density unpack rho_bf bounds: rank=", dg_frag%id, &
                " id_frag=", dg_frag%id_frag, " peer=", irank, " slot=", slot, &
                " idx=", ixg, iyg, izg, " lb=", lbound(rho_bf, 1), lbound(rho_bf, 2), lbound(rho_bf, 3), &
                " ub=", ubound(rho_bf, 1), ubound(rho_bf, 2), ubound(rho_bf, 3)
              flush(6)
              stop "DG-Fragment RT: density unpack rho_bf bounds"
            end if
            if (ixg < lbound(rho_s_bf, 1) .or. ixg > ubound(rho_s_bf, 1) .or. &
                iyg < lbound(rho_s_bf, 2) .or. iyg > ubound(rho_s_bf, 2) .or. &
                izg < lbound(rho_s_bf, 3) .or. izg > ubound(rho_s_bf, 3)) then
              write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,3(i0,1x),a,3(i0,1x),a,3(i0,1x))') &
                "[FATAL] density unpack rho_s_bf bounds: rank=", dg_frag%id, &
                " id_frag=", dg_frag%id_frag, " peer=", irank, " slot=", slot, &
                " idx=", ixg, iyg, izg, " lb=", lbound(rho_s_bf, 1), lbound(rho_s_bf, 2), lbound(rho_s_bf, 3), &
                " ub=", ubound(rho_s_bf, 1), ubound(rho_s_bf, 2), ubound(rho_s_bf, 3)
              flush(6)
              stop "DG-Fragment RT: density unpack rho_s_bf bounds"
            end if
            rho_raw_contrib = recv_flat(recv_displs(irank + 1) + slot)
            rho_contrib = rho_raw_contrib
            rho_bf(bx, by, bz) = rho_bf(bx, by, bz) + rho_contrib
            call assert_rho_materialize_value('rho_bf-unpack-full', irank, 0, 0, &
              slot, bx, by, bz, ixg, iyg, izg, dg_frag%id, slot, &
              rho_contrib, rho_bf(bx, by, bz))
            do ispin = 1, system%nspin
              spin_offset = ispin * npts
              if (spin_offset + slot < 1 .or. spin_offset + slot > recv_counts(irank + 1)) then
                write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0)') &
                  "[FATAL] density unpack spin slot bounds: rank=", dg_frag%id, &
                  " id_frag=", dg_frag%id_frag, " peer=", irank, " slot=", slot, &
                  " spin_slot=", spin_offset + slot, " recv_size=", recv_counts(irank + 1)
                flush(6)
                stop "DG-Fragment RT: density unpack spin slot bounds"
              end if
              rho_contrib = recv_flat(recv_displs(irank + 1) + spin_offset + slot)
!$omp atomic update
              rho_s_bf(ixg, iyg, izg, ispin) = rho_s_bf(ixg, iyg, izg, ispin) + rho_contrib
            end do
          end do
        end if
        call assert_real_array3_finite_density('rho_bf-after-unpack-peer', rho_bf, &
          lbound(rho_bf, 1), lbound(rho_bf, 2), lbound(rho_bf, 3), &
          rho_x_lo, rho_x_hi, rho_y_lo, rho_y_hi, rho_z_lo, rho_z_hi, 0)
      end do
      call assert_real_array3_finite_density('rho_bf-after-unpack-loop', rho_bf, &
        lbound(rho_bf, 1), lbound(rho_bf, 2), lbound(rho_bf, 3), &
        rho_x_lo, rho_x_hi, rho_y_lo, rho_y_hi, rho_z_lo, rho_z_hi, 0)
      call assert_real_array3_finite_density('rho_bf-after-unpack-timer', rho_bf, &
        lbound(rho_bf, 1), lbound(rho_bf, 2), lbound(rho_bf, 3), &
        rho_x_lo, rho_x_hi, rho_y_lo, rho_y_hi, rho_z_lo, rho_z_hi, 0)
      call assert_real_array3_finite_density('rho_bf-after-exchange', rho_bf, &
        lbound(rho_bf, 1), lbound(rho_bf, 2), lbound(rho_bf, 3), &
        rho_x_lo, rho_x_hi, rho_y_lo, rho_y_hi, rho_z_lo, rho_z_hi, 0)
      call assert_real_array3_finite_density('rho_bf-before-comm-buffer-dealloc', rho_bf, &
        lbound(rho_bf, 1), lbound(rho_bf, 2), lbound(rho_bf, 3), &
        rho_x_lo, rho_x_hi, rho_y_lo, rho_y_hi, rho_z_lo, rho_z_hi, 0)
      deallocate(send_flat, recv_flat)
      deallocate(send_counts, recv_counts, send_displs, recv_displs)
      call assert_real_array3_finite_density('rho_bf-after-comm-buffer-dealloc', rho_bf, &
        lbound(rho_bf, 1), lbound(rho_bf, 2), lbound(rho_bf, 3), &
        rho_x_lo, rho_x_hi, rho_y_lo, rho_y_hi, rho_z_lo, rho_z_hi, 0)
    end if
    if (density_exchange_active) then
    call assert_rho_bf_with_sources('rho_bf-before-copy-v3', &
      max(0, dg_frag%ifrag_group), 0, block_idx_global, max(1, min(grid_block_size, ngrid_max)))

    ! rho_bf -> rho_s boundary: only the authoritative mg-local density is
    ! materialized below for downstream Hartree/XC/reconstruct consumers.
    if (lbound(rho%f, 1) /= rho_x_lo .or. ubound(rho%f, 1) /= rho_x_hi .or. &
        lbound(rho%f, 2) /= rho_y_lo .or. ubound(rho%f, 2) /= rho_y_hi .or. &
        lbound(rho%f, 3) /= rho_z_lo .or. ubound(rho%f, 3) /= rho_z_hi) then
      write(*,'(1x,a,i0,a,3(i0,1x),a,3(i0,1x),a,3(i0,1x),a,3(i0,1x))') &
        "[FATAL] density rho core bounds mismatch: rank=", dg_frag%id, &
        " rho_lb=", lbound(rho%f, 1), lbound(rho%f, 2), lbound(rho%f, 3), &
        " rho_ub=", ubound(rho%f, 1), ubound(rho%f, 2), ubound(rho%f, 3), &
        " core_lo=", rho_x_lo, rho_y_lo, rho_z_lo, &
        " core_hi=", rho_x_hi, rho_y_hi, rho_z_hi
      flush(6)
      stop "DG-Fragment RT: density rho core bounds mismatch"
    end if
    rho%f(rho_x_lo:rho_x_hi, rho_y_lo:rho_y_hi, rho_z_lo:rho_z_hi) = &
      rho_bf(rho_x_lo:rho_x_hi, rho_y_lo:rho_y_hi, rho_z_lo:rho_z_hi)
    if (system%nspin == 1) then
      rho_s(1)%f(rho_x_lo:rho_x_hi, rho_y_lo:rho_y_hi, rho_z_lo:rho_z_hi) = &
        rho%f(rho_x_lo:rho_x_hi, rho_y_lo:rho_y_hi, rho_z_lo:rho_z_hi)
    else
      do ispin = 1, system%nspin
        if (lbound(rho_s(ispin)%f, 1) /= lbound(rho_s_bf, 1) .or. ubound(rho_s(ispin)%f, 1) /= ubound(rho_s_bf, 1) .or. &
            lbound(rho_s(ispin)%f, 2) /= lbound(rho_s_bf, 2) .or. ubound(rho_s(ispin)%f, 2) /= ubound(rho_s_bf, 2) .or. &
            lbound(rho_s(ispin)%f, 3) /= lbound(rho_s_bf, 3) .or. ubound(rho_s(ispin)%f, 3) /= ubound(rho_s_bf, 3)) then
        write(*,'(1x,a,i0,a,i0,a,3(i0,1x),a,3(i0,1x),a,3(i0,1x),a,3(i0,1x))') &
          "[FATAL] density rho_s copy bounds mismatch: rank=", dg_frag%id, " ispin=", ispin, &
          " rho_lb=", lbound(rho_s(ispin)%f, 1), lbound(rho_s(ispin)%f, 2), lbound(rho_s(ispin)%f, 3), &
          " rho_ub=", ubound(rho_s(ispin)%f, 1), ubound(rho_s(ispin)%f, 2), ubound(rho_s(ispin)%f, 3), &
          " bf_lb=", lbound(rho_s_bf, 1), lbound(rho_s_bf, 2), lbound(rho_s_bf, 3), &
          " bf_ub=", ubound(rho_s_bf, 1), ubound(rho_s_bf, 2), ubound(rho_s_bf, 3)
        flush(6)
        stop "DG-Fragment RT: density rho_s copy bounds mismatch"
        end if
        rho_s(ispin)%f(:, :, :) = rho_s_bf(:, :, :, ispin)
      end do
    end if
    call assert_real_array3_finite_density('rho-after-copy', rho%f, &
      lbound(rho%f, 1), lbound(rho%f, 2), lbound(rho%f, 3), &
      rho_x_lo, rho_x_hi, rho_y_lo, rho_y_hi, rho_z_lo, rho_z_hi, 0)
    do ispin = 1, system%nspin
      call assert_real_array3_finite_density('rho_s-after-copy', rho_s(ispin)%f, &
        lbound(rho_s(ispin)%f, 1), lbound(rho_s(ispin)%f, 2), lbound(rho_s(ispin)%f, 3), &
        rho_x_lo, rho_x_hi, rho_y_lo, rho_y_hi, rho_z_lo, rho_z_hi, ispin)
    end do
    total_charge_local = 0.0d0
    do iz = rho_z_lo, rho_z_hi
      do iy = rho_y_lo, rho_y_hi
        do ix = rho_x_lo, rho_x_hi
          total_charge_local = total_charge_local + rho%f(ix, iy, iz)
        end do
      end do
    end do
    total_charge_local = total_charge_local * system%hvol
    call assert_real_scalar_finite_density('total_charge_local-scaled', total_charge_local)
    if (trace_density_charge) then
      rho_trace_min = minval(rho%f(rho_x_lo:rho_x_hi, rho_y_lo:rho_y_hi, rho_z_lo:rho_z_hi))
      rho_trace_max = maxval(rho%f(rho_x_lo:rho_x_hi, rho_y_lo:rho_y_hi, rho_z_lo:rho_z_hi))
      write(*,'(1x,a,i0,a,i0,a,i0,a,l1,a,l1,3(a,es15.6))') &
        '[DG-DENSITY-TRACE] rank=', dg_frag%id, ' itt=', itt_tag, ' frag_group=', dg_frag%ifrag_group, &
        ' sparse_exchange=', density_exchange_active, ' wannier=', enable_wannier_density, &
        ' charge=', total_charge_local, ' rho_min=', rho_trace_min, ' rho_max=', rho_trace_max
      flush(6)
    end if
    ! Keep the raw local electron count.  Avoid a mandatory world-level
    ! collective here: in orbital-fragment RT this diagnostic reduction can be
    ! reached by a different active-rank set than the density construction
    ! itself on Fugaku.  The density grid is already materialized above, and no
    ! normalization is applied from this value.
    total_charge_reduce_in(1) = total_charge_local
    total_charge_reduce_out(1) = total_charge_reduce_in(1)
    total_charge = total_charge_reduce_out(1)
    call assert_real_scalar_finite_density('total_charge_local', total_charge_local)
    call assert_real_scalar_finite_density('total_charge', total_charge)
    dg_frag%rho_global_raw_elec = total_charge
    dg_frag%rho_owned_elec = total_charge
    dg_frag%rho_imported_elec = 0.0d0
    dg_frag%rho_local_all_elec = total_charge
    dg_frag%rho_contract_residual_elec = 0.0d0
    dg_frag%elec_num_raw = total_charge
    dg_frag%rho_scale_factor = 1.0d0
    ! Do not renormalize the reconstructed RT density here.  The raw integrated
    ! charge is the conserved quantity we need to diagnose; rescaling would hide
    ! density reconstruction errors and can spread non-finite scale factors over
    ! the whole real-space grid.
    dg_frag%elec_num_scaled = total_charge
    call assert_real_array3_finite_density('rho-final', rho%f, &
      lbound(rho%f, 1), lbound(rho%f, 2), lbound(rho%f, 3), &
      rho_x_lo, rho_x_hi, rho_y_lo, rho_y_hi, rho_z_lo, rho_z_hi, 0)
    do ispin = 1, system%nspin
      call assert_real_array3_finite_density('rho_s-final', rho_s(ispin)%f, &
        lbound(rho_s(ispin)%f, 1), lbound(rho_s(ispin)%f, 2), lbound(rho_s(ispin)%f, 3), &
        rho_x_lo, rho_x_hi, rho_y_lo, rho_y_hi, rho_z_lo, rho_z_hi, ispin)
    end do
    dg_frag%rho_global_scaled_elec = dg_frag%elec_num_scaled
    end if
    deallocate(ix_buf, iy_buf, iz_buf, owner_buf, ixg_buf, iyg_buf, izg_buf)
    deallocate(rhobf_bx_buf, rhobf_by_buf, rhobf_bz_buf)
    deallocate(slot_buf, local_grid_ids, remote_grid_ids, valid_remote_grid_ids)
    deallocate(basis_gid, valid_basis_ids)
    if (allocated(basis_gid_spin)) deallocate(basis_gid_spin)
    if (allocated(valid_basis_ids_spin)) deallocate(valid_basis_ids_spin)
    if (allocated(valid_basis_count_spin)) deallocate(valid_basis_count_spin)
    if (allocated(owned_basis_ids_spin)) deallocate(owned_basis_ids_spin)
    if (allocated(owned_coef_local_ids_spin)) deallocate(owned_coef_local_ids_spin)
    if (allocated(owned_basis_count_spin)) deallocate(owned_basis_count_spin)
    deallocate(phi_blk, rho_blk, rho_blk_accum, rho_blk_reduced, coef_blk_re, coef_blk_im, psi_blk_re, psi_blk_im)
    if (allocated(coef_blk_ri)) deallocate(coef_blk_ri)
    if (allocated(psi_blk_ri)) deallocate(psi_blk_ri)
    if (allocated(pw_tmp_z)) deallocate(pw_tmp_z)
    if (allocated(coef_pw_blk)) deallocate(coef_pw_blk)
    if (allocated(phase_cache)) deallocate(phase_cache)
    if (allocated(kpw_hx)) deallocate(kpw_hx, kpw_hy, kpw_hz)
    if (allocated(density_mix)) deallocate(density_mix)
    if (allocated(density_mix_partial)) deallocate(density_mix_partial)
    if (allocated(basis_mix_blk)) deallocate(basis_mix_blk)
    if (allocated(density_mix_tmp)) deallocate(density_mix_tmp)
    if (allocated(basis_mix_blk_t)) deallocate(basis_mix_blk_t)
    if (allocated(density_mix_tmp_t)) deallocate(density_mix_tmp_t)
    if (allocated(transform_frag_spin)) deallocate(transform_frag_spin)
    if (allocated(transform_pw_spin)) deallocate(transform_pw_spin)
    if (allocated(mix_transform_spin)) deallocate(mix_transform_spin)
    if (allocated(mix_overlap_spin)) deallocate(mix_overlap_spin)
    if (allocated(s_mix)) deallocate(s_mix)
    if (allocated(s_mix_work)) deallocate(s_mix_work)
    if (allocated(coef_mix_eff)) deallocate(coef_mix_eff)
    if (allocated(coef_mix_metric)) deallocate(coef_mix_metric)
    if (allocated(coef_mix_spin)) deallocate(coef_mix_spin)
    if (allocated(d_raw_ff)) deallocate(d_raw_ff)
    if (allocated(d_raw_fp)) deallocate(d_raw_fp)
    if (allocated(d_raw_pp)) deallocate(d_raw_pp)
    if (allocated(ipiv_mix)) deallocate(ipiv_mix)
    if (allocated(n_basis_mix_spin)) deallocate(n_basis_mix_spin)
    deallocate(rho_bf)
    if (allocated(rho_s_bf)) deallocate(rho_s_bf)
    deallocate(rho_send, rho_recv)
  contains

    logical function density_basis_owned_by_rank(global_basis, ispin_arg) result(is_owned)
      implicit none
      integer, intent(in) :: global_basis, ispin_arg

      is_owned = .true.
      if (.not. dg_frag%parallel_mode_orbital) return
      if (dg_frag%isize_frag <= 1) return
      if (.not. allocated(dg_frag%coef_owner)) return
      if (global_basis < 1 .or. global_basis > size(dg_frag%coef_owner, 1)) then
        is_owned = .false.
        return
      end if
      if (ispin_arg < 1 .or. ispin_arg > size(dg_frag%coef_owner, 2)) then
        is_owned = .false.
        return
      end if
      is_owned = (dg_frag%coef_owner(global_basis, ispin_arg) == dg_frag%id)
    end function density_basis_owned_by_rank

    subroutine assert_real_scalar_finite_density(label, val)
      implicit none
      character(*), intent(in) :: label
      real(8), intent(in) :: val

      if (.not. ieee_is_finite(val)) then
        write(*,'(1x,a,a,a,i0,a,i0,a,1pe20.12)') &
          '[FATAL] invalid DG density scalar: label=', trim(label), &
          ' rank=', dg_frag%id, ' id_frag=', dg_frag%id_frag, ' value=', val
        flush(6)
        stop 'DG-Fragment RT: invalid density scalar'
      end if
    end subroutine assert_real_scalar_finite_density

    subroutine assert_rho_materialize_value(label, ifrag_arg, ispin_arg, block_arg, &
      igrid_arg, bx_arg, by_arg, bz_arg, ixg_arg, iyg_arg, izg_arg, owner_arg, slot_arg, contrib_arg, val_arg)
      implicit none
      character(*), intent(in) :: label
      integer, intent(in) :: ifrag_arg, ispin_arg, block_arg, igrid_arg
      integer, intent(in) :: bx_arg, by_arg, bz_arg, ixg_arg, iyg_arg, izg_arg
      integer, intent(in) :: owner_arg, slot_arg
      real(8), intent(in) :: contrib_arg, val_arg

      if (ieee_is_finite(contrib_arg) .and. ieee_is_finite(val_arg) .and. &
          abs(contrib_arg) <= density_abs_limit .and. abs(val_arg) <= density_abs_limit) return
      write(*,*) '[FATAL] invalid DG rho materialize: label=', trim(label), &
        ' rank=', dg_frag%id, ' id_frag=', dg_frag%id_frag, &
        ' ifrag=', ifrag_arg, ' ispin=', ispin_arg, ' block=', block_arg, &
        ' igrid=', igrid_arg, ' bx=', bx_arg, ' by=', by_arg, ' bz=', bz_arg, &
        ' ixg=', ixg_arg, ' iyg=', iyg_arg, ' izg=', izg_arg, &
        ' owner=', owner_arg, ' slot=', slot_arg, &
        ' contrib=', contrib_arg, ' value=', val_arg, ' limit=', density_abs_limit
      write(0,*) '[FATAL] invalid DG rho materialize: label=', trim(label), &
        ' rank=', dg_frag%id, ' id_frag=', dg_frag%id_frag, &
        ' ifrag=', ifrag_arg, ' ispin=', ispin_arg, ' block=', block_arg, &
        ' igrid=', igrid_arg, ' bx=', bx_arg, ' by=', by_arg, ' bz=', bz_arg, &
        ' ixg=', ixg_arg, ' iyg=', iyg_arg, ' izg=', izg_arg, &
        ' owner=', owner_arg, ' slot=', slot_arg, &
        ' contrib=', contrib_arg, ' value=', val_arg, ' limit=', density_abs_limit
      flush(6)
      flush(0)
      stop 'DG-Fragment RT: invalid rho materialize'
    end subroutine assert_rho_materialize_value

    subroutine assert_rho_bf_with_sources(label, ifrag_arg, ispin_arg, block_arg, npt_arg)
      implicit none
      character(*), intent(in) :: label
      integer, intent(in) :: ifrag_arg, ispin_arg, block_arg, npt_arg
      integer :: ix_a, iy_a, iz_a, nbad
      integer :: first_i, first_j, first_k
      integer :: igrid_src, source_count, first_source_igrid
      integer :: first_source_ixg, first_source_iyg, first_source_izg
      integer :: first_source_owner, first_source_slot
      real(8) :: max_abs, first_value, source_sum, first_source_contrib

      nbad = 0
      max_abs = 0.0d0
      first_i = 0
      first_j = 0
      first_k = 0
      first_value = 0.0d0
      do iz_a = rho_z_lo, rho_z_hi
        do iy_a = rho_y_lo, rho_y_hi
          do ix_a = rho_x_lo, rho_x_hi
            if (ieee_is_finite(rho_bf(ix_a, iy_a, iz_a))) then
              max_abs = max(max_abs, abs(rho_bf(ix_a, iy_a, iz_a)))
              if (abs(rho_bf(ix_a, iy_a, iz_a)) <= density_abs_limit) cycle
            end if
            nbad = nbad + 1
            if (first_i == 0) then
              first_i = ix_a
              first_j = iy_a
              first_k = iz_a
              first_value = rho_bf(ix_a, iy_a, iz_a)
            end if
          end do
        end do
      end do
      if (nbad <= 0) return

      source_count = 0
      first_source_igrid = 0
      first_source_ixg = 0
      first_source_iyg = 0
      first_source_izg = 0
      first_source_owner = -1
      first_source_slot = -1
      first_source_contrib = 0.0d0
      source_sum = 0.0d0
      do igrid_src = 1, npt_arg
        if (rhobf_bx_buf(igrid_src) /= first_i) cycle
        if (rhobf_by_buf(igrid_src) /= first_j) cycle
        if (rhobf_bz_buf(igrid_src) /= first_k) cycle
        source_count = source_count + 1
        source_sum = source_sum + rho_blk_accum(igrid_src)
        if (first_source_igrid == 0) then
          first_source_igrid = igrid_src
          first_source_ixg = ixg_buf(igrid_src)
          first_source_iyg = iyg_buf(igrid_src)
          first_source_izg = izg_buf(igrid_src)
          first_source_owner = owner_buf(igrid_src)
          first_source_slot = slot_buf(igrid_src)
          first_source_contrib = rho_blk_accum(igrid_src)
        end if
      end do

      write(*,*) '[FATAL] invalid DG rho_bf block: label=', trim(label), &
        ' rank=', dg_frag%id, ' id_frag=', dg_frag%id_frag, &
        ' ifrag=', ifrag_arg, ' ispin=', ispin_arg, ' block=', block_arg, &
        ' count=', nbad, ' first_i=', first_i, ' first_j=', first_j, &
        ' first_k=', first_k, ' first_val=', first_value, &
        ' finite_max=', max_abs, ' limit=', density_abs_limit
      write(*,*) '[FATAL] rho_bf source probe: rank=', dg_frag%id, &
        ' source_count=', source_count, ' source_sum=', source_sum, &
        ' igrid=', first_source_igrid, ' ixg=', first_source_ixg, &
        ' iyg=', first_source_iyg, ' izg=', first_source_izg, &
        ' owner=', first_source_owner, ' slot=', first_source_slot, &
        ' contrib=', first_source_contrib
      write(0,*) '[FATAL] invalid DG rho_bf block: label=', trim(label), &
        ' rank=', dg_frag%id, ' id_frag=', dg_frag%id_frag, &
        ' ifrag=', ifrag_arg, ' ispin=', ispin_arg, ' block=', block_arg, &
        ' count=', nbad, ' first_i=', first_i, ' first_j=', first_j, &
        ' first_k=', first_k, ' first_val=', first_value, &
        ' finite_max=', max_abs, ' limit=', density_abs_limit
      write(0,*) '[FATAL] rho_bf source probe: rank=', dg_frag%id, &
        ' source_count=', source_count, ' source_sum=', source_sum, &
        ' igrid=', first_source_igrid, ' ixg=', first_source_ixg, &
        ' iyg=', first_source_iyg, ' izg=', first_source_izg, &
        ' owner=', first_source_owner, ' slot=', first_source_slot, &
        ' contrib=', first_source_contrib
      flush(6)
      flush(0)
      stop 'DG-Fragment RT: invalid rho_bf block'
    end subroutine assert_rho_bf_with_sources

    subroutine assert_real_array3_finite_density(label, vals, vals_ix_lo, vals_iy_lo, vals_iz_lo, &
      ix_lo, ix_hi, iy_lo, iy_hi, iz_lo, iz_hi, ispin_arg)
      implicit none
      character(*), intent(in) :: label
      integer, intent(in) :: vals_ix_lo, vals_iy_lo, vals_iz_lo
      real(8), intent(in) :: vals(vals_ix_lo:, vals_iy_lo:, vals_iz_lo:)
      integer, intent(in) :: ix_lo, ix_hi, iy_lo, iy_hi, iz_lo, iz_hi, ispin_arg
      integer :: ix_a, iy_a, iz_a, nbad
      integer :: first_bad_i, first_bad_j, first_bad_k
      real(8) :: max_abs, first_bad_value

      nbad = 0
      max_abs = 0.0d0
      first_bad_i = 0
      first_bad_j = 0
      first_bad_k = 0
      first_bad_value = 0.0d0
      do iz_a = iz_lo, iz_hi
        do iy_a = iy_lo, iy_hi
          do ix_a = ix_lo, ix_hi
            if (ieee_is_finite(vals(ix_a, iy_a, iz_a))) then
              max_abs = max(max_abs, abs(vals(ix_a, iy_a, iz_a)))
              if (abs(vals(ix_a, iy_a, iz_a)) > density_abs_limit) then
                nbad = nbad + 1
                if (first_bad_i == 0) then
                  first_bad_i = ix_a
                  first_bad_j = iy_a
                  first_bad_k = iz_a
                  first_bad_value = vals(ix_a, iy_a, iz_a)
                end if
              end if
            else
              nbad = nbad + 1
              if (first_bad_i == 0) then
                first_bad_i = ix_a
                first_bad_j = iy_a
                first_bad_k = iz_a
                first_bad_value = vals(ix_a, iy_a, iz_a)
              end if
            end if
          end do
        end do
      end do
      if (nbad > 0) then
        write(*,*) '[FATAL] invalid DG density array: label=', trim(label), &
          ' rank=', dg_frag%id, ' id_frag=', dg_frag%id_frag, ' ispin=', ispin_arg, &
          ' count=', nbad, ' first_i=', first_bad_i, ' first_j=', first_bad_j, &
          ' first_k=', first_bad_k, ' first_val=', first_bad_value, &
          ' finite_max=', max_abs, ' limit=', density_abs_limit
        write(0,*) '[FATAL] invalid DG density array: label=', trim(label), &
          ' rank=', dg_frag%id, ' id_frag=', dg_frag%id_frag, ' ispin=', ispin_arg, &
          ' count=', nbad, ' first_i=', first_bad_i, ' first_j=', first_bad_j, &
          ' first_k=', first_bad_k, ' first_val=', first_bad_value, &
          ' finite_max=', max_abs, ' limit=', density_abs_limit
        flush(6)
        flush(0)
        stop 'DG-Fragment RT: invalid density array'
      end if
    end subroutine assert_real_array3_finite_density

    subroutine assert_real_vector_finite_density(label, vals, nvals, ifrag_arg, ispin_arg, block_arg)
      implicit none
      character(*), intent(in) :: label
      real(8), intent(in) :: vals(:)
      integer, intent(in) :: nvals, ifrag_arg, ispin_arg, block_arg
      integer :: i_a, nbad, first_bad
      real(8) :: max_abs, first_bad_value

      nbad = 0
      first_bad = 0
      max_abs = 0.0d0
      first_bad_value = 0.0d0
      do i_a = 1, nvals
        if (ieee_is_finite(vals(i_a))) then
          max_abs = max(max_abs, abs(vals(i_a)))
          if (abs(vals(i_a)) > density_abs_limit) then
            if (first_bad == 0) then
              first_bad = i_a
              first_bad_value = vals(i_a)
            end if
            nbad = nbad + 1
          end if
        else
          if (first_bad == 0) then
            first_bad = i_a
            first_bad_value = vals(i_a)
          end if
          nbad = nbad + 1
        end if
      end do
      if (nbad > 0) then
        write(*,'(1x,a,a,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,1pe20.12,a,1pe20.12,a,1pe12.4)') &
          '[FATAL] invalid DG density vector: label=', trim(label), &
          ' rank=', dg_frag%id, ' id_frag=', dg_frag%id_frag, &
          ' ifrag=', ifrag_arg, ' ispin=', ispin_arg, &
          ' block=', block_arg, ' count=', nbad, ' first_bad=', first_bad, &
          ' first_val=', first_bad_value, ' finite_max=', max_abs, ' limit=', density_abs_limit
        flush(6)
        write(0,'(1x,a,a,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,1pe20.12,a,1pe20.12,a,1pe12.4)') &
          '[FATAL] invalid DG density vector: label=', trim(label), &
          ' rank=', dg_frag%id, ' id_frag=', dg_frag%id_frag, &
          ' ifrag=', ifrag_arg, ' ispin=', ispin_arg, &
          ' block=', block_arg, ' count=', nbad, ' first_bad=', first_bad, &
          ' first_val=', first_bad_value, ' finite_max=', max_abs, ' limit=', density_abs_limit
        flush(0)
        stop 'DG-Fragment RT: invalid density vector'
      end if
    end subroutine assert_real_vector_finite_density

    subroutine assert_real_matrix_finite_density(label, vals, i_lo, i_hi, j_lo, j_hi, &
      ifrag_arg, ispin_arg, block_arg)
      implicit none
      character(*), intent(in) :: label
      real(8), intent(in) :: vals(:,:)
      integer, intent(in) :: i_lo, i_hi, j_lo, j_hi, ifrag_arg, ispin_arg, block_arg
      integer :: i_a, j_a, nbad, first_bad_i, first_bad_j
      real(8) :: max_abs, first_bad_value

      nbad = 0
      max_abs = 0.0d0
      first_bad_i = 0
      first_bad_j = 0
      first_bad_value = 0.0d0
      do j_a = j_lo, j_hi
        do i_a = i_lo, i_hi
          if (ieee_is_finite(vals(i_a, j_a))) then
            max_abs = max(max_abs, abs(vals(i_a, j_a)))
            if (abs(vals(i_a, j_a)) > density_abs_limit) then
              if (first_bad_i == 0) then
                first_bad_i = i_a
                first_bad_j = j_a
                first_bad_value = vals(i_a, j_a)
              end if
              nbad = nbad + 1
            end if
          else
            if (first_bad_i == 0) then
              first_bad_i = i_a
              first_bad_j = j_a
              first_bad_value = vals(i_a, j_a)
            end if
            nbad = nbad + 1
          end if
        end do
      end do
      if (nbad > 0) then
        write(*,'(1x,a,a,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,1pe20.12,a,1pe20.12,a,1pe12.4)') &
          '[FATAL] invalid DG density matrix: label=', trim(label), &
          ' rank=', dg_frag%id, ' id_frag=', dg_frag%id_frag, &
          ' ifrag=', ifrag_arg, ' ispin=', ispin_arg, &
          ' block=', block_arg, ' count=', nbad, &
          ' first_i=', first_bad_i, ' first_j=', first_bad_j, &
          ' first_val=', first_bad_value, ' finite_max=', max_abs, ' limit=', density_abs_limit
        flush(6)
        write(0,'(1x,a,a,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,1pe20.12,a,1pe20.12,a,1pe12.4)') &
          '[FATAL] invalid DG density matrix: label=', trim(label), &
          ' rank=', dg_frag%id, ' id_frag=', dg_frag%id_frag, &
          ' ifrag=', ifrag_arg, ' ispin=', ispin_arg, &
          ' block=', block_arg, ' count=', nbad, &
          ' first_i=', first_bad_i, ' first_j=', first_bad_j, &
          ' first_val=', first_bad_value, ' finite_max=', max_abs, ' limit=', density_abs_limit
        flush(0)
        stop 'DG-Fragment RT: invalid density matrix'
      end if
    end subroutine assert_real_matrix_finite_density

    subroutine prepare_grid_buffers_owner_map(i_local_grid, igrid0_grid, npt_blk_grid, nxyz_grid, use_subgroup_slot)
      implicit none
      integer, intent(in) :: i_local_grid, igrid0_grid, npt_blk_grid
      integer, intent(in) :: nxyz_grid(3)
      logical, intent(in) :: use_subgroup_slot
      type(density_grid_point_info) :: point

      if (nxyz_grid(1) < 0) return
      if (use_subgroup_slot .and. dg_frag%parallel_mode_orbital) then
        write(*,'(1x,a,i0,a,i0,a)') '[FATAL] orbital mode entered subgroup-slot path: rank=', dg_frag%id, &
          ' id_frag=', dg_frag%id_frag, ' path=prepare_grid_buffers_owner_map'
        flush(6)
        stop 'DG-Fragment RT orbital mode: subgroup real-space slot path is not allowed'
      end if

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
          rhobf_bx_buf(igrid) = 0
          rhobf_by_buf(igrid) = 0
          rhobf_bz_buf(igrid) = 0
          owner_buf(igrid) = -1
          slot_buf(igrid) = 0
          cycle
        end if
        ixg_buf(igrid) = point%ixg
        iyg_buf(igrid) = point%iyg
        izg_buf(igrid) = point%izg
        rhobf_bx_buf(igrid) = dg_frag%density_grid_bx(igrid0_grid + igrid - 1, i_local_grid)
        rhobf_by_buf(igrid) = dg_frag%density_grid_by(igrid0_grid + igrid - 1, i_local_grid)
        rhobf_bz_buf(igrid) = dg_frag%density_grid_bz(igrid0_grid + igrid - 1, i_local_grid)
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

      if (nxyz_grid(1) < 0) return
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
          rhobf_bx_buf(igrid) = 0
          rhobf_by_buf(igrid) = 0
          rhobf_bz_buf(igrid) = 0
          owner_buf(igrid) = -1
          slot_buf(igrid) = 0
          cycle
        end if
        ixg_buf(igrid) = point%ixg
        iyg_buf(igrid) = point%iyg
        izg_buf(igrid) = point%izg
        rhobf_bx_buf(igrid) = dg_frag%density_grid_bx(igrid0_grid + igrid - 1, i_local_grid)
        rhobf_by_buf(igrid) = dg_frag%density_grid_by(igrid0_grid + igrid - 1, i_local_grid)
        rhobf_bz_buf(igrid) = dg_frag%density_grid_bz(igrid0_grid + igrid - 1, i_local_grid)
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
    integer :: frag_size

    frag_size = max(1, dg_frag%isize_frag)
    handler_id_frag = modulo(target_rank, frag_size)
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

  subroutine exchange_density_sparse(sendbuf, scounts, sdispls, recvbuf, rcounts, rdispls)
    use mpi, only: MPI_DOUBLE_PRECISION, MPI_STATUS_SIZE, MPI_REQUEST_NULL
    implicit none
    real(8), intent(in) :: sendbuf(:)
    integer, intent(in) :: scounts(:), sdispls(:)
    real(8), intent(inout) :: recvbuf(:)
    integer, intent(in) :: rcounts(:), rdispls(:)
    integer, parameter :: density_exchange_tag = 26017
    integer, parameter :: density_exchange_chunk_size = 524288
    integer :: dest, src, scount, rcount, sidx, ridx
    integer :: soff, roff, send_chunk, recv_chunk
    integer :: ierr, peer_idx, chunk_idx, max_chunks, req_count, max_req
    integer :: send_peer_count, recv_peer_count
    integer, allocatable :: send_peers(:), recv_peers(:)
    integer, allocatable :: requests(:), statuses(:,:)

    if (dg_frag%isize <= 1) then
      if (size(recvbuf) >= 1 .and. size(sendbuf) >= 1 .and. scounts(1) > 0) then
        recvbuf(1:scounts(1)) = sendbuf(1:scounts(1))
      end if
      return
    end if

    scount = scounts(dg_frag%id + 1)
    rcount = rcounts(dg_frag%id + 1)
    if (scount /= rcount) then
      write(*,'(1x,a,i0,a,i0,a,i0)') &
        "[FATAL] DG density self-exchange size mismatch: rank=", dg_frag%id, &
        " send=", scount, " recv=", rcount
      flush(6)
      stop "DG-Fragment RT: density self-exchange size mismatch"
    end if
    if (scount > 0) then
      sidx = sdispls(dg_frag%id + 1) + 1
      ridx = rdispls(dg_frag%id + 1) + 1
      recvbuf(ridx:ridx + scount - 1) = sendbuf(sidx:sidx + scount - 1)
    end if

    allocate(send_peers(dg_frag%isize), recv_peers(dg_frag%isize))
    send_peer_count = 0
    recv_peer_count = 0
    max_chunks = 0
    do peer_idx = 0, dg_frag%isize - 1
      if (peer_idx == dg_frag%id) cycle
      scount = scounts(peer_idx + 1)
      rcount = rcounts(peer_idx + 1)
      if (scount > 0) then
        send_peer_count = send_peer_count + 1
        send_peers(send_peer_count) = peer_idx
        max_chunks = max(max_chunks, (scount + density_exchange_chunk_size - 1) / density_exchange_chunk_size)
      end if
      if (rcount > 0) then
        recv_peer_count = recv_peer_count + 1
        recv_peers(recv_peer_count) = peer_idx
        max_chunks = max(max_chunks, (rcount + density_exchange_chunk_size - 1) / density_exchange_chunk_size)
      end if
    end do

    if (send_peer_count == 0 .and. recv_peer_count == 0) then
      deallocate(send_peers, recv_peers)
      return
    end if

    ! Only communicate with ranks that actually own sparse density payload.
    ! Receives for the current chunk are posted before sends, so each phase is
    ! deadlock-free without forcing every rank through all MPI ranks as in the
    ! old ring schedule.  Chunking keeps each Tofu-registered buffer bounded.
    max_req = max(1, send_peer_count + recv_peer_count)
    allocate(requests(max_req), statuses(MPI_STATUS_SIZE, max_req))
    do chunk_idx = 0, max_chunks - 1
      requests(:) = MPI_REQUEST_NULL
      req_count = 0

      do peer_idx = 1, recv_peer_count
        src = recv_peers(peer_idx)
        rcount = rcounts(src + 1)
        roff = chunk_idx * density_exchange_chunk_size
        if (roff >= rcount) cycle
        recv_chunk = min(density_exchange_chunk_size, rcount - roff)
        ridx = rdispls(src + 1) + 1
        req_count = req_count + 1
        call MPI_Irecv(recvbuf(ridx + roff), recv_chunk, MPI_DOUBLE_PRECISION, src, density_exchange_tag, &
                       dg_frag%icomm, requests(req_count), ierr)
        if (ierr /= 0) then
          write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0)') &
            "[FATAL] DG density sparse irecv failed: rank=", dg_frag%id, &
            " chunk=", chunk_idx, " src=", src, " count=", recv_chunk, " ierr=", ierr
          flush(6)
          stop "DG-Fragment RT: density sparse irecv failed"
        end if
      end do

      do peer_idx = 1, send_peer_count
        dest = send_peers(peer_idx)
        scount = scounts(dest + 1)
        soff = chunk_idx * density_exchange_chunk_size
        if (soff >= scount) cycle
        send_chunk = min(density_exchange_chunk_size, scount - soff)
        sidx = sdispls(dest + 1) + 1
        req_count = req_count + 1
        call MPI_Isend(sendbuf(sidx + soff), send_chunk, MPI_DOUBLE_PRECISION, dest, density_exchange_tag, &
                       dg_frag%icomm, requests(req_count), ierr)
        if (ierr /= 0) then
          write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0)') &
            "[FATAL] DG density sparse isend failed: rank=", dg_frag%id, &
            " chunk=", chunk_idx, " dest=", dest, " count=", send_chunk, " ierr=", ierr
          flush(6)
          stop "DG-Fragment RT: density sparse isend failed"
        end if
      end do

      if (req_count > 0) then
        call MPI_Waitall(req_count, requests, statuses, ierr)
        if (ierr /= 0) then
          write(*,'(1x,a,i0,a,i0,a,i0)') &
            "[FATAL] DG density sparse waitall failed: rank=", dg_frag%id, &
            " chunk=", chunk_idx, " ierr=", ierr
          flush(6)
          stop "DG-Fragment RT: density sparse waitall failed"
        end if
      end if
    end do
    deallocate(requests, statuses, send_peers, recv_peers)
  end subroutine exchange_density_sparse

  end subroutine calculate_density_from_fragments
