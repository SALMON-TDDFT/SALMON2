  subroutine update_density_hamiltonian_stage(dg_frag, system, info, rt, itt, Ac_tot, &
                                              lg, mg, stencil, xc_func, srg, srg_scalar, fg, poisson, pp, ppg, ppn, &
                                              rho, rho_s, Vh, Vxc, Vpsl, energy)
    use structures
    use sendrecv_grid, only: s_sendrecv_grid
    use communication, only: comm_summation, comm_get_min
    use salmon_xc, only: s_xc_functional, exchange_correlation
    use poisson_dg_distributed, only: hartree_dg_distributed
    use density_matrix_and_energy_plusU_sub, only: calc_density_matrix_and_energy_plusU, PLUS_U_ON
    use rt_dg_fragment_ops, only: ensure_nonlocal_pp_matrix_A
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
    integer :: n_frag, n_pw, ispin
    integer :: n_floor_local, n_floor_global
    real(8) :: rho_floor_min_local, rho_floor_min_global
    real(8) :: rho_floor_min_buf(1)
    real(8), allocatable :: rho_buffer(:,:,:), Vh_buffer(:,:,:)
    logical :: use_rank_buffered_potential
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
    time_density = 0.0d0
    time_hartree = 0.0d0
    time_xc = 0.0d0
    time_reconstruct = 0.0d0
    time_pw_mix = 0.0d0
    trace_call_id = 0
    trace_stage = .false.
    if (itt == 1 .and. dg_frag%id == 0 .and. trace_stage_call_count < 5) then
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
    call cpu_time(t0)
    call calculate_density_from_fragments(dg_frag, system, mg, rho, rho_s, itt)
    call cpu_time(t1)
    time_density = time_density + (t1 - t0)
    if (trace_stage) then
      write(*,'(1x,a,1pe12.4)') '[DG-STAGE] density done time=', t1 - t0
      flush(6)
    end if

    dg_frag%rho_frag(:, :, :) = rho%f(:, :, :)
    if (system%nspin > 0) then
      dg_frag%rho_s_frag(:, :, :, 1:system%nspin) = 0.0d0
      do n_frag = 1, system%nspin
        dg_frag%rho_s_frag(:, :, :, n_frag) = rho_s(n_frag)%f(:, :, :)
      end do
    end if

    use_rank_buffered_potential = all(dg_frag%rank_buf_hi(:) >= dg_frag%rank_buf_lo(:))
    if (use_rank_buffered_potential) then
      allocate(rho_buffer(dg_frag%rank_buf_lo(1):dg_frag%rank_buf_hi(1), &
                          dg_frag%rank_buf_lo(2):dg_frag%rank_buf_hi(2), &
                          dg_frag%rank_buf_lo(3):dg_frag%rank_buf_hi(3)))
      allocate(Vh_buffer(dg_frag%rank_buf_lo(1):dg_frag%rank_buf_hi(1), &
                         dg_frag%rank_buf_lo(2):dg_frag%rank_buf_hi(2), &
                         dg_frag%rank_buf_lo(3):dg_frag%rank_buf_hi(3)))
      call copy_periodic_global_scalar_to_rank_buffer(dg_frag, mg, rho, rho_buffer)
    end if
    if (trace_stage) then
      write(*,'(1x,a)') '[DG-STAGE] Hartree start'
      flush(6)
    end if
    call cpu_time(t0)
    call hartree_dg_distributed(info, lg, mg, fg, poisson, dg_frag, rho, Vh)
    call cpu_time(t1)
    time_hartree = time_hartree + (t1 - t0)
    if (trace_stage) then
      write(*,'(1x,a,1pe12.4)') '[DG-STAGE] Hartree done time=', t1 - t0
      flush(6)
    end if
    dg_frag%Vh_frag(:, :, :) = Vh%f(:, :, :)
    if (use_rank_buffered_potential) then
      call copy_periodic_global_scalar_to_rank_buffer(dg_frag, mg, Vh, Vh_buffer)
    end if
    ! Guard against tiny negative rho before XC; this prevents non-physical Vxc NaN under aggressive FP reassociation.
    ! Validation note: O3 and O3+no-unsafe-math runs matched key observables (J_para/J_total/Ne_raw) in this case.
    n_floor_local = 0
    rho_floor_min_local = huge(1.0d0)
    do ispin = 1, system%nspin
      n_floor_local = n_floor_local + count(rho_s(ispin)%f(:, :, :) < 0.0d0)
      if (any(rho_s(ispin)%f(:, :, :) < 0.0d0)) then
        rho_floor_min_local = min(rho_floor_min_local, minval(rho_s(ispin)%f(:, :, :), &
          mask=rho_s(ispin)%f(:, :, :) < 0.0d0))
      end if
      where (rho_s(ispin)%f(:, :, :) < 0.0d0) rho_s(ispin)%f(:, :, :) = 0.0d0
    end do
    call comm_summation(n_floor_local, n_floor_global, dg_frag%icomm)
    if (n_floor_global > 0) then
      rho_floor_min_buf(1) = rho_floor_min_local
      call comm_get_min(rho_floor_min_buf, rho_floor_min_buf, 1, dg_frag%icomm)
      rho_floor_min_global = rho_floor_min_buf(1)
      if (dg_frag%id == 0) then
        write(*,'(1x,a,i0,a,i0,a,1pe14.6)') '[RHO-FLOOR] itt=', itt, ' clamped_points=', n_floor_global, &
          ' min_rho_before=', rho_floor_min_global
        flush(6)
      end if
    end if
    if (trace_stage) then
      write(*,'(1x,a)') '[DG-STAGE] XC start'
      flush(6)
    end if
    call cpu_time(t0)
    call exchange_correlation(system, xc_func, mg, srg_scalar, srg, rho_s, pp, ppn, &
                 info, rt%tpsi0, stencil, Vxc, energy%E_xc)
    call cpu_time(t1)
    time_xc = time_xc + (t1 - t0)
    if (trace_stage) then
      write(*,'(1x,a,1pe12.4)') '[DG-STAGE] XC done time=', t1 - t0
      flush(6)
    end if
    if (system%nspin > 0) then
      dg_frag%Vxc_frag(:, :, :, 1:system%nspin) = 0.0d0
      do n_frag = 1, system%nspin
        dg_frag%Vxc_frag(:, :, :, n_frag) = Vxc(n_frag)%f(:, :, :)
      end do
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
      call cpu_time(t0)
      call reconstruct_hamiltonian_matrix(dg_frag, system, stencil, Vh, Vxc, Vpsl, Ac_tot, Vh_buffer)
      call cpu_time(t1)
      time_reconstruct = time_reconstruct + (t1 - t0)
    else
      if (trace_stage) then
        write(*,'(1x,a)') '[DG-STAGE] reconstruct start'
        flush(6)
      end if
      call cpu_time(t0)
      call reconstruct_hamiltonian_matrix(dg_frag, system, stencil, Vh, Vxc, Vpsl, Ac_tot)
      call cpu_time(t1)
      time_reconstruct = time_reconstruct + (t1 - t0)
    end if
    if (trace_stage) then
      write(*,'(1x,a,1pe12.4)') '[DG-STAGE] reconstruct done time=', time_reconstruct
      flush(6)
    end if

    if (dg_frag%use_plane_wave_basis .and. dg_frag%n_plane_waves > 0) then
      if (trace_stage) then
        write(*,'(1x,a)') '[DG-STAGE] PW mixed update start'
        flush(6)
      end if
      call cpu_time(t0)
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

      ! Mixed FP/PP nonlocal terms use the subgroup-reduced block cache.
      ! Requesting the dense cache here scales as n_mat_max**2 and is not
      ! viable for large DC->DG runs.
      call ensure_nonlocal_pp_matrix_A(dg_frag, mg, ppg, system, Ac_tot)

      call compute_fragment_pw_hamiltonian(dg_frag, Vh, Vxc, Vpsl, dg_frag%H_mat_frag_pw)
      call build_mixed_hamiltonian(dg_frag, lg, Vh, Vxc, Vpsl, Ac_tot, dg_frag%S_mat_frag_pw, dg_frag%H_mat_frag_pw)
      call cpu_time(t1)
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

    if (allocated(rho_buffer)) deallocate(rho_buffer)
    if (allocated(Vh_buffer)) deallocate(Vh_buffer)

  end subroutine update_density_hamiltonian_stage
