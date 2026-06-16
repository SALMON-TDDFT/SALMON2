  subroutine update_density_hamiltonian_stage(dg_frag, system, info, rt, itt, Ac_tot, &
                                              lg, mg, stencil, xc_func, srg, srg_scalar, fg, poisson, pp, ppg, ppn, &
                                              rho, rho_s, Vh, Vxc, Vpsl, energy, update_energy)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use structures
    use sendrecv_grid, only: s_sendrecv_grid
    use salmon_xc, only: s_xc_functional, exchange_correlation, &
                         salmon_xctype_none, salmon_xctype_pz
    use poisson_dg_distributed, only: hartree_dg_distributed
    use density_matrix_and_energy_plusU_sub, only: calc_density_matrix_and_energy_plusU, PLUS_U_ON
    use salmon_global, only: yn_spinorbit, ae_shape1, ae_shape2, theory
    use misc_routines, only: get_wtime
    use rt_dg_plane_wave, only: compute_fragment_pw_overlap, &
      compute_fragment_pw_hamiltonian, build_mixed_hamiltonian
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system),     intent(inout) :: system
    type(s_parallel_info),  intent(in)    :: info
    type(s_rt),             intent(inout) :: rt
    integer,                intent(in)    :: itt
    real(8),                intent(in)    :: Ac_tot(3)
    type(s_rgrid),          intent(in)    :: lg, mg
    type(s_stencil),        intent(in)    :: stencil
    type(s_xc_functional),  intent(in)    :: xc_func
    type(s_sendrecv_grid),  intent(inout) :: srg, srg_scalar
    type(s_reciprocal_grid),intent(in)    :: fg
    type(s_poisson),        intent(inout) :: poisson
    type(s_pp_info),        intent(in)    :: pp
    type(s_pp_grid),        intent(in)    :: ppg
    type(s_pp_nlcc),        intent(in)    :: ppn
    type(s_scalar),         intent(inout) :: rho, Vh, Vpsl
    type(s_scalar),         intent(inout) :: rho_s(system%nspin), Vxc(system%nspin)
    type(s_dft_energy),     intent(inout) :: energy
    logical, optional,      intent(in)    :: update_energy
    integer :: n_frag, n_pw, ispin
    integer :: n_floor_local
    integer :: n_nan_rho_local
    integer :: n_nan_vh_local
    integer :: n_nan_vhbuf_local
    integer :: n_bad_vxcbuf_local
    logical :: use_rank_buffered_potential, use_fragment_xc
    logical :: need_energy_update
    logical :: need_stage_buffer_alloc
    real(8) :: t0, t1
    real(8) :: time_density, time_hartree, time_xc, time_reconstruct, time_pw_mix
    integer :: trace_call_id
    integer, save :: trace_stage_call_count = 0
    logical :: trace_stage
    logical, save :: timing_initialized = .false.
    logical, save :: enable_stage_timing = .false.
    character(16) :: env_timing
    integer :: env_status

    if (.not. timing_initialized) then
      env_timing = ''
      call get_environment_variable('SALMON_DG_STAGE_TIMING', env_timing, status=env_status)
      if (env_status == 0) then
        select case(trim(adjustl(env_timing)))
        case('1','y','Y','yes','YES','true','TRUE','on','ON')
          enable_stage_timing = .true.
        end select
      end if
      timing_initialized = .true.
    end if
    n_pw = max(0, dg_frag%n_plane_waves)
    need_energy_update = .true.
    if (present(update_energy)) need_energy_update = update_energy
    time_density = 0.0d0
    time_hartree = 0.0d0
    time_xc = 0.0d0
    time_reconstruct = 0.0d0
    time_pw_mix = 0.0d0
    trace_call_id = 0
    trace_stage = .false.
    if (enable_stage_timing .and. itt == 1 .and. dg_frag%id == 0 .and. trace_stage_call_count < 5) then
      trace_stage_call_count = trace_stage_call_count + 1
      trace_call_id = trace_stage_call_count
      trace_stage = .true.
      write(*,'(1x,a,i0,a,i0)') '[DG-STAGE] enter itt=', itt, ' call=', trace_call_id
      flush(6)
    end if

    if (trace_stage) then
      write(*,'(1x,a)') '[DG-STAGE] density start'
      flush(6)
    end if
    t0 = get_wtime()
    call calculate_density_from_fragments(dg_frag, system, mg, rho, rho_s, itt)
    t1 = get_wtime()
    time_density = time_density + (t1 - t0)
    n_nan_rho_local = count(.not. ieee_is_finite(rho%f(mg%is(1):mg%ie(1), &
                                                       mg%is(2):mg%ie(2), mg%is(3):mg%ie(3))))
    do ispin = 1, system%nspin
      n_nan_rho_local = n_nan_rho_local + &
        count(.not. ieee_is_finite(rho_s(ispin)%f(mg%is(1):mg%ie(1), &
                                                  mg%is(2):mg%ie(2), mg%is(3):mg%ie(3))))
    end do
    if (n_nan_rho_local > 0) then
      write(*,'(1x,a,i0,a,i0,a,i0,a,i0)') '[FATAL] NaN/Inf in DG stage density: rank=', &
        dg_frag%id, ' itt=', itt, ' call=', trace_call_id, ' count=', n_nan_rho_local
      flush(6)
      stop "NaN/Inf in DG stage density"
    end if
    if (trace_stage) then
      write(*,'(1x,a,1pe12.4)') '[DG-STAGE] density done time=', t1 - t0
      flush(6)
    end if

    use_rank_buffered_potential = all(dg_frag%rank_buf_hi(:) >= dg_frag%rank_buf_lo(:))
    use_fragment_xc = use_rank_buffered_potential .and. dg_frag%parallel_mode_orbital .and. n_pw == 0 .and. &
                      system%nspin == 1 .and. xc_func%ispin == 0 .and. .not. xc_func%use_gradient .and. &
                      .not. xc_func%use_laplacian .and. .not. xc_func%use_kinetic_energy .and. &
                      .not. xc_func%use_current .and. xc_func%xctype(1) == salmon_xctype_pz .and. &
                      xc_func%xctype(2) == salmon_xctype_none .and. &
                      xc_func%xctype(3) == salmon_xctype_none .and. .not. pp%flag_nlcc .and. &
                      yn_spinorbit /= 'y'
    if (use_rank_buffered_potential) then
      need_stage_buffer_alloc = (.not. allocated(dg_frag%stage_Vh_buffer)) .or. &
                                (.not. allocated(dg_frag%stage_Vpsl_buffer)) .or. &
                                (.not. allocated(dg_frag%stage_Vxc_buffer))
      if (.not. need_stage_buffer_alloc) then
        need_stage_buffer_alloc = any(lbound(dg_frag%stage_Vh_buffer) /= dg_frag%rank_buf_lo) .or. &
                                  any(ubound(dg_frag%stage_Vh_buffer) /= dg_frag%rank_buf_hi) .or. &
                                  any(lbound(dg_frag%stage_Vpsl_buffer) /= dg_frag%rank_buf_lo) .or. &
                                  any(ubound(dg_frag%stage_Vpsl_buffer) /= dg_frag%rank_buf_hi) .or. &
                                  lbound(dg_frag%stage_Vxc_buffer, 1) /= dg_frag%rank_buf_lo(1) .or. &
                                  lbound(dg_frag%stage_Vxc_buffer, 2) /= dg_frag%rank_buf_lo(2) .or. &
                                  lbound(dg_frag%stage_Vxc_buffer, 3) /= dg_frag%rank_buf_lo(3) .or. &
                                  ubound(dg_frag%stage_Vxc_buffer, 1) /= dg_frag%rank_buf_hi(1) .or. &
                                  ubound(dg_frag%stage_Vxc_buffer, 2) /= dg_frag%rank_buf_hi(2) .or. &
                                  ubound(dg_frag%stage_Vxc_buffer, 3) /= dg_frag%rank_buf_hi(3) .or. &
                                  size(dg_frag%stage_Vxc_buffer, 4) /= system%nspin
      end if
      if (need_stage_buffer_alloc) then
        if (allocated(dg_frag%stage_Vh_buffer)) deallocate(dg_frag%stage_Vh_buffer)
        if (allocated(dg_frag%stage_Vpsl_buffer)) deallocate(dg_frag%stage_Vpsl_buffer)
        if (allocated(dg_frag%stage_Vxc_buffer)) deallocate(dg_frag%stage_Vxc_buffer)
        if (allocated(dg_frag%stage_gx_map)) deallocate(dg_frag%stage_gx_map)
        if (allocated(dg_frag%stage_gy_map)) deallocate(dg_frag%stage_gy_map)
        if (allocated(dg_frag%stage_gz_map)) deallocate(dg_frag%stage_gz_map)
        allocate(dg_frag%stage_Vh_buffer(dg_frag%rank_buf_lo(1):dg_frag%rank_buf_hi(1), &
                                         dg_frag%rank_buf_lo(2):dg_frag%rank_buf_hi(2), &
                                         dg_frag%rank_buf_lo(3):dg_frag%rank_buf_hi(3)))
        allocate(dg_frag%stage_Vpsl_buffer(dg_frag%rank_buf_lo(1):dg_frag%rank_buf_hi(1), &
                                           dg_frag%rank_buf_lo(2):dg_frag%rank_buf_hi(2), &
                                           dg_frag%rank_buf_lo(3):dg_frag%rank_buf_hi(3)))
        allocate(dg_frag%stage_Vxc_buffer(dg_frag%rank_buf_lo(1):dg_frag%rank_buf_hi(1), &
                                          dg_frag%rank_buf_lo(2):dg_frag%rank_buf_hi(2), &
                                          dg_frag%rank_buf_lo(3):dg_frag%rank_buf_hi(3), system%nspin))
        allocate(dg_frag%stage_gx_map(dg_frag%rank_buf_lo(1):dg_frag%rank_buf_hi(1)))
        allocate(dg_frag%stage_gy_map(dg_frag%rank_buf_lo(2):dg_frag%rank_buf_hi(2)))
        allocate(dg_frag%stage_gz_map(dg_frag%rank_buf_lo(3):dg_frag%rank_buf_hi(3)))
        dg_frag%stage_vpsl_buffer_valid = .false.
        dg_frag%stage_map_valid = .false.
      end if
      if (.not. dg_frag%stage_vpsl_buffer_valid) then
        ! Vpsl is static for DG-RT; keep the buffer copy out of the RK hot path.
        call copy_periodic_global_scalar_to_rank_buffer(dg_frag, mg, Vpsl, dg_frag%stage_Vpsl_buffer)
        dg_frag%stage_vpsl_buffer_valid = .true.
      end if
    end if
    if (trace_stage) then
      write(*,'(1x,a)') '[DG-STAGE] Hartree start'
      flush(6)
    end if
    t0 = get_wtime()
    call hartree_dg_distributed(info, lg, mg, fg, poisson, dg_frag, rho, Vh)
    t1 = get_wtime()
    time_hartree = time_hartree + (t1 - t0)
    n_nan_vh_local = count(.not. ieee_is_finite(Vh%f(mg%is(1):mg%ie(1), &
                                                     mg%is(2):mg%ie(2), mg%is(3):mg%ie(3))))
    if (n_nan_vh_local > 0) then
      write(*,'(1x,a,i0,a,i0,a,i0,a,i0)') '[FATAL] NaN/Inf in DG Hartree Vh: rank=', &
        dg_frag%id, ' itt=', itt, ' call=', trace_call_id, ' count=', n_nan_vh_local
      flush(6)
      stop "NaN/Inf in DG Hartree Vh"
    end if
    if (trace_stage) then
      write(*,'(1x,a,1pe12.4)') '[DG-STAGE] Hartree done time=', t1 - t0
      flush(6)
    end if
    if (use_rank_buffered_potential) then
      call copy_periodic_global_scalar_to_rank_buffer(dg_frag, mg, Vh, dg_frag%stage_Vh_buffer)
      n_nan_vhbuf_local = count(.not. ieee_is_finite(dg_frag%stage_Vh_buffer( &
        dg_frag%rank_buf_lo(1):dg_frag%rank_buf_hi(1), &
        dg_frag%rank_buf_lo(2):dg_frag%rank_buf_hi(2), &
        dg_frag%rank_buf_lo(3):dg_frag%rank_buf_hi(3))))
      if (n_nan_vhbuf_local > 0) then
        write(*,'(1x,a,i0,a,i0,a,i0,a,i0)') '[FATAL] NaN/Inf in DG Hartree Vh buffer: rank=', &
          dg_frag%id, ' itt=', itt, ' call=', trace_call_id, ' count=', n_nan_vhbuf_local
        flush(6)
        stop "NaN/Inf in DG Hartree Vh buffer"
      end if
    end if
    ! Guard against tiny negative rho before XC; this prevents non-physical Vxc NaN under aggressive FP reassociation.
    ! Validation note: O3 and O3+no-unsafe-math runs matched key observables (J_para/J_total/Ne_raw) in this case.
    n_floor_local = 0
    do ispin = 1, system%nspin
      n_floor_local = n_floor_local + count(rho_s(ispin)%f(:, :, :) < 0.0d0)
      where (rho_s(ispin)%f(:, :, :) < 0.0d0) rho_s(ispin)%f(:, :, :) = 0.0d0
    end do
    if (trace_stage) then
      write(*,'(1x,a)') '[DG-STAGE] XC start'
      if (use_fragment_xc) then
        write(*,'(1x,a)') '[DG-STAGE] XC route=fragment-pz-buffer'
      else
        write(*,'(1x,a)') '[DG-STAGE] XC route=generic'
      end if
      flush(6)
    end if
    t0 = get_wtime()
    if (use_fragment_xc) then
      call compute_fragment_pz_xc_buffer(dg_frag, system, info, mg, rho_s(1), &
                                         dg_frag%stage_Vxc_buffer(:, :, :, 1), energy%E_xc, need_energy_update)
    else
      call exchange_correlation(system, xc_func, mg, srg_scalar, srg, rho_s, pp, ppn, &
                   info, rt%tpsi0, stencil, Vxc, energy%E_xc)
    end if
    t1 = get_wtime()
    time_xc = time_xc + (t1 - t0)
    if (trace_stage) then
      write(*,'(1x,a,1pe12.4)') '[DG-STAGE] XC done time=', t1 - t0
      flush(6)
    end if
    if (system%nspin > 0 .and. .not. use_fragment_xc .and. use_rank_buffered_potential) then
      do n_frag = 1, system%nspin
        call copy_periodic_global_scalar_to_rank_buffer(dg_frag, mg, Vxc(n_frag), &
                                                        dg_frag%stage_Vxc_buffer(:, :, :, n_frag))
      end do
    end if
    if (use_rank_buffered_potential) then
      n_bad_vxcbuf_local = count(.not. ieee_is_finite(dg_frag%stage_Vxc_buffer(:, :, :, 1:system%nspin)))
      if (n_bad_vxcbuf_local > 0) then
        write(*,'(1x,a,i0,a,i0,a,i0,a,i0)') '[FATAL] NaN/Inf in DG Vxc buffer: rank=', &
          dg_frag%id, ' itt=', itt, ' call=', trace_call_id, ' count=', n_bad_vxcbuf_local
        flush(6)
        stop "NaN/Inf in DG Vxc buffer"
      end if
    end if

    if (dg_frag%use_plusu .and. PLUS_U_ON) then
      call calc_density_matrix_and_energy_plusU(rt%tpsi0, ppg, info, system, energy%E_U)
    else
      energy%E_U = 0.0d0
    end if
    if (use_rank_buffered_potential) then
      if (trace_stage) then
        write(*,'(1x,a)') '[DG-STAGE] reconstruct start'
        flush(6)
      end if
      t0 = get_wtime()
      call reconstruct_hamiltonian_matrix(dg_frag, system, stencil, Vh, Vxc, Vpsl, Ac_tot, &
                                          dg_frag%stage_Vh_buffer, dg_frag%stage_Vxc_buffer, &
                                          dg_frag%stage_Vpsl_buffer)
      t1 = get_wtime()
      time_reconstruct = time_reconstruct + (t1 - t0)
    else
      if (trace_stage) then
        write(*,'(1x,a)') '[DG-STAGE] reconstruct start'
        flush(6)
      end if
      t0 = get_wtime()
      call reconstruct_hamiltonian_matrix(dg_frag, system, stencil, Vh, Vxc, Vpsl, Ac_tot)
      t1 = get_wtime()
      time_reconstruct = time_reconstruct + (t1 - t0)
    end if
    if (trace_stage) then
      write(*,'(1x,a,1pe12.4)') '[DG-STAGE] reconstruct done time=', time_reconstruct
      flush(6)
    end if

    ! Nonlocal PP is applied directly by the derivative/observable projector
    ! route.  Building H_nl_blocks here would force a global block reduction at
    ! every RT stage and breaks the rank-distributed DG-RT design.

    if (dg_frag%use_plane_wave_basis .and. dg_frag%n_plane_waves > 0) then
      if (trace_stage) then
        write(*,'(1x,a)') '[DG-STAGE] PW mixed update start'
        flush(6)
      end if
      t0 = get_wtime()
      n_frag = size(dg_frag%coef, 1)
      n_pw = dg_frag%n_plane_waves

      if (.not. allocated(dg_frag%S_mat_frag_pw)) then
        allocate(dg_frag%S_mat_frag_pw(n_frag, n_pw, dg_frag%nspin))
        call compute_fragment_pw_overlap(dg_frag, dg_frag%S_mat_frag_pw)
      else if (size(dg_frag%S_mat_frag_pw,1) /= n_frag .or. size(dg_frag%S_mat_frag_pw,2) /= n_pw .or. &
               size(dg_frag%S_mat_frag_pw,3) /= dg_frag%nspin) then
        deallocate(dg_frag%S_mat_frag_pw)
        allocate(dg_frag%S_mat_frag_pw(n_frag, n_pw, dg_frag%nspin))
        call compute_fragment_pw_overlap(dg_frag, dg_frag%S_mat_frag_pw)
      end if

      if (.not. allocated(dg_frag%H_mat_frag_pw)) then
        allocate(dg_frag%H_mat_frag_pw(n_frag, n_pw, dg_frag%nspin))
      else if (size(dg_frag%H_mat_frag_pw,1) /= n_frag .or. size(dg_frag%H_mat_frag_pw,2) /= n_pw .or. &
               size(dg_frag%H_mat_frag_pw,3) /= dg_frag%nspin) then
        deallocate(dg_frag%H_mat_frag_pw)
        allocate(dg_frag%H_mat_frag_pw(n_frag, n_pw, dg_frag%nspin))
      end if

      call compute_fragment_pw_hamiltonian(dg_frag, Vh, Vxc, Vpsl, dg_frag%H_mat_frag_pw)
      call build_mixed_hamiltonian(dg_frag, lg, Vh, Vxc, Vpsl, Ac_tot, dg_frag%S_mat_frag_pw, dg_frag%H_mat_frag_pw)
      t1 = get_wtime()
      time_pw_mix = time_pw_mix + (t1 - t0)
      if (trace_stage) then
        write(*,'(1x,a,1pe12.4)') '[DG-STAGE] PW mixed update done time=', t1 - t0
        flush(6)
      end if
    end if

    if (enable_stage_timing .and. dg_frag%id == 0) then
      write(*,'(1x,a,i0,5(a,1pe12.4))') '        stage timing: itt=', itt, &
        ' density=', time_density, ' hartree=', time_hartree, ' xc=', time_xc, &
        ' reconstruct=', time_reconstruct, ' pw_mix=', time_pw_mix
      flush(6)
    end if

  end subroutine update_density_hamiltonian_stage

  subroutine compute_fragment_pz_xc_buffer(dg_frag, system, info, mg, rho_s_in, Vxc_buffer, E_xc, update_energy)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use structures, only: s_dft_system, s_parallel_info, s_rgrid, s_scalar
    use communication, only: comm_summation
    use builtin_pz, only: PZxc
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system), intent(in) :: system
    type(s_parallel_info), intent(in) :: info
    type(s_rgrid), intent(in) :: mg
    type(s_scalar), intent(in) :: rho_s_in
    real(8), intent(out) :: Vxc_buffer(dg_frag%rank_buf_lo(1):dg_frag%rank_buf_hi(1), &
                                       dg_frag%rank_buf_lo(2):dg_frag%rank_buf_hi(2), &
                                       dg_frag%rank_buf_lo(3):dg_frag%rank_buf_hi(3))
    real(8), intent(inout) :: E_xc
    logical, intent(in) :: update_energy
    integer :: ix, iy, iz
    integer :: gx, gy, gz
    integer :: rho_nan_count
    real(8) :: rho_val, exc, dexc_drho
    real(8) :: e_xc_local, e_xc_global

    if (.not. allocated(dg_frag%stage_gx_map)) then
      allocate(dg_frag%stage_gx_map(dg_frag%rank_buf_lo(1):dg_frag%rank_buf_hi(1)))
    end if
    if (.not. allocated(dg_frag%stage_gy_map)) then
      allocate(dg_frag%stage_gy_map(dg_frag%rank_buf_lo(2):dg_frag%rank_buf_hi(2)))
    end if
    if (.not. allocated(dg_frag%stage_gz_map)) then
      allocate(dg_frag%stage_gz_map(dg_frag%rank_buf_lo(3):dg_frag%rank_buf_hi(3)))
    end if
    if (.not. dg_frag%stage_map_valid) then
      do ix = dg_frag%rank_buf_lo(1), dg_frag%rank_buf_hi(1)
        dg_frag%stage_gx_map(ix) = map_global_to_stage_box_coord(ix, mg%is(1), mg%ie(1), dg_frag%lgnum_total(1))
      end do
      do iy = dg_frag%rank_buf_lo(2), dg_frag%rank_buf_hi(2)
        dg_frag%stage_gy_map(iy) = map_global_to_stage_box_coord(iy, mg%is(2), mg%ie(2), dg_frag%lgnum_total(2))
      end do
      do iz = dg_frag%rank_buf_lo(3), dg_frag%rank_buf_hi(3)
        dg_frag%stage_gz_map(iz) = map_global_to_stage_box_coord(iz, mg%is(3), mg%ie(3), dg_frag%lgnum_total(3))
      end do
      dg_frag%stage_map_valid = .true.
    end if

    rho_nan_count = count(.not. ieee_is_finite(rho_s_in%f(mg%is(1):mg%ie(1), &
                                                          mg%is(2):mg%ie(2), mg%is(3):mg%ie(3))))
    if (rho_nan_count > 0) then
      write(*,'(1x,a,i0,a,i0)') "[FATAL] NaN/Inf density before fragment PZ XC buffer: rank=", &
        dg_frag%id, " count=", rho_nan_count
      flush(6)
      stop "NaN/Inf density before fragment PZ XC buffer"
    end if

    ! The DG Hamiltonian rebuild only needs Vxc on the fragment buffer.
    ! The scalar XC energy is integrated only when the caller requests an
    ! observable energy update; RK intermediate stages only need the potential.
!$omp parallel do private(ix,iy,gx,gy,gz,rho_val,exc,dexc_drho) schedule(static)
    do iz = dg_frag%rank_buf_lo(3), dg_frag%rank_buf_hi(3)
      gz = dg_frag%stage_gz_map(iz)
      do iy = dg_frag%rank_buf_lo(2), dg_frag%rank_buf_hi(2)
        gy = dg_frag%stage_gy_map(iy)
!$omp simd private(gx,rho_val,exc,dexc_drho)
        do ix = dg_frag%rank_buf_lo(1), dg_frag%rank_buf_hi(1)
          gx = dg_frag%stage_gx_map(ix)
          if (gx == 0 .or. gy == 0 .or. gz == 0) then
            Vxc_buffer(ix, iy, iz) = 0.0d0
          else
            rho_val = max(rho_s_in%f(gx, gy, gz), 0.0d0)
            if (rho_val <= 0.0d0) then
              Vxc_buffer(ix, iy, iz) = 0.0d0
            else
              call PZxc(rho_val, exc, dexc_drho)
              Vxc_buffer(ix, iy, iz) = exc + rho_val * dexc_drho
            end if
          end if
        end do
      end do
    end do
!$omp end parallel do

    if (.not. update_energy) return

    e_xc_local = 0.0d0
!$omp parallel do private(ix,iy,rho_val,exc,dexc_drho) reduction(+:e_xc_local) schedule(static)
    do iz = mg%is(3), mg%ie(3)
      do iy = mg%is(2), mg%ie(2)
!$omp simd private(rho_val,exc,dexc_drho) reduction(+:e_xc_local)
        do ix = mg%is(1), mg%ie(1)
          rho_val = max(rho_s_in%f(ix, iy, iz), 0.0d0)
          if (rho_val > 0.0d0) then
            call PZxc(rho_val, exc, dexc_drho)
            e_xc_local = e_xc_local + exc * rho_val
          end if
        end do
      end do
    end do
!$omp end parallel do
    e_xc_local = e_xc_local * system%hvol
    call comm_summation(e_xc_local, e_xc_global, info%icomm_r)
    E_xc = e_xc_global
  end subroutine compute_fragment_pz_xc_buffer

  integer function map_global_to_stage_box_coord(ig, lb, ub, lgtot) result(iloc)
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
  end function map_global_to_stage_box_coord
