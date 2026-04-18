  subroutine time_evolution_rk(dg_frag, system, info, rt, itt, dt, &
                             lg, mg, stencil, xc_func, srg, srg_scalar, fg, poisson, pp, ppg, ppn, &
                             rho, rho_s, Vh, Vxc, Vpsl, energy)
    use structures
    use communication, only: comm_summation
    use salmon_global, only: yn_fix_func
    use sendrecv_grid, only: s_sendrecv_grid
    use salmon_xc, only: s_xc_functional
    use rt_dg_fragment_ops, only: sync_mixed_coef_from_raw, apply_overlap_operator, gather_full_coef_view, &
                    zero_nonowned_coefficients
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system),     intent(inout) :: system
    type(s_parallel_info),  intent(in)    :: info
    type(s_rt),             intent(inout) :: rt
    integer,                intent(in)    :: itt
    real(8),                intent(in)    :: dt
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
    
    complex(8), allocatable, save :: k(:,:,:,:)
    complex(8), allocatable, save :: k_pw(:,:,:,:)
    complex(8), allocatable :: coef_ref(:,:,:)
    complex(8), allocatable :: coef_pw_ref(:,:,:)
    integer :: istage, io, jo, ispin
    integer :: jo_local, jo_global, owned_s, owned_e, n_owned
    real(8) :: Ac_tot(3), t_stage
    integer :: n, n_pw
    logical :: use_mixed_rt
    logical :: found_nan
    integer :: nan_jo, nan_io, nan_ispin
    complex(8) :: nan_val
    logical, parameter :: enable_rk_trace = .false.
    logical, parameter :: enable_rk_nan_check = .false.
    logical :: enable_rk_stage_metric, zero_nonowned_after_stage, renormalize_occ_after_rk
    character(len=64) :: env_stage_metric, env_stage_metric_maxitt
    character(len=64) :: env_zero_nonowned, env_renorm_occ
    integer :: env_stage_metric_len, env_stage_metric_status
    integer :: stage_metric_max_itt
    logical :: enable_hmat_s_check
    character(len=64) :: env_hmat_s_check
    
    ! Use n_mat_max (global basis size) instead of nstate_frag (local basis size)
    n = dg_frag%n_mat_max
    n_pw = 0
    if (dg_frag%use_plane_wave_basis .and. allocated(dg_frag%coef_pw)) n_pw = dg_frag%n_plane_waves
    use_mixed_rt = (n_pw > 0 .and. dg_frag%mixed_basis_ready .and. allocated(dg_frag%coef_mix))
    if (enable_rk_trace .and. itt == 1 .and. dg_frag%id == 0) then
      write(*,'(1x,a,l1,a,i0,a,l1,a,i0,a,i0,a,i0)') "        rk basis probe: use_pw=", &
        dg_frag%use_plane_wave_basis, " n_plane_waves=", dg_frag%n_plane_waves, " coef_pw_alloc=", &
        allocated(dg_frag%coef_pw), " coef_pw_dim1=", merge(size(dg_frag%coef_pw,1), 0, allocated(dg_frag%coef_pw)), &
        " n_mat_max=", n, " nstate_tot=", dg_frag%nstate_tot
      flush(6)
    end if
    enable_rk_stage_metric = .false.
    zero_nonowned_after_stage = .false.
    renormalize_occ_after_rk = (n_pw == 0)
    stage_metric_max_itt = 5
    env_stage_metric = ''
    call get_environment_variable("SALMON_DG_RK_STAGE_METRIC", env_stage_metric, length=env_stage_metric_len, status=env_stage_metric_status)
    if (env_stage_metric_status == 0 .and. env_stage_metric_len > 0) then
      if (env_stage_metric(1:1) == '1' .or. env_stage_metric(1:1) == 'y' .or. env_stage_metric(1:1) == 'Y' .or. &
          env_stage_metric(1:1) == 't' .or. env_stage_metric(1:1) == 'T') then
        enable_rk_stage_metric = .true.
      end if
    end if
    env_stage_metric_maxitt = ''
    call get_environment_variable("SALMON_DG_RK_STAGE_METRIC_MAXITT", env_stage_metric_maxitt, &
                                  length=env_stage_metric_len, status=env_stage_metric_status)
    if (env_stage_metric_status == 0 .and. env_stage_metric_len > 0) then
      read(env_stage_metric_maxitt(1:env_stage_metric_len), *, err=100, end=100) stage_metric_max_itt
100   continue
      stage_metric_max_itt = max(1, stage_metric_max_itt)
    end if
    env_zero_nonowned = ''
    call get_environment_variable("SALMON_DG_ZERO_NONOWNED_AFTER_RK_STAGE", env_zero_nonowned, &
                                  length=env_stage_metric_len, status=env_stage_metric_status)
    if (env_stage_metric_status == 0 .and. env_stage_metric_len > 0) then
      if (env_zero_nonowned(1:1) == '1' .or. env_zero_nonowned(1:1) == 'y' .or. env_zero_nonowned(1:1) == 'Y' .or. &
          env_zero_nonowned(1:1) == 't' .or. env_zero_nonowned(1:1) == 'T') then
        zero_nonowned_after_stage = .true.
      end if
    end if
    env_renorm_occ = ''
    call get_environment_variable("SALMON_DG_RENORMALIZE_OCC_AFTER_RK", env_renorm_occ, &
                                  length=env_stage_metric_len, status=env_stage_metric_status)
    if (env_stage_metric_status == 0 .and. env_stage_metric_len > 0) then
      if (env_renorm_occ(1:1) == '1' .or. env_renorm_occ(1:1) == 'y' .or. env_renorm_occ(1:1) == 'Y' .or. &
          env_renorm_occ(1:1) == 't' .or. env_renorm_occ(1:1) == 'T') then
        renormalize_occ_after_rk = .true.
      else if (env_renorm_occ(1:1) == '0' .or. env_renorm_occ(1:1) == 'n' .or. env_renorm_occ(1:1) == 'N' .or. &
               env_renorm_occ(1:1) == 'f' .or. env_renorm_occ(1:1) == 'F') then
        renormalize_occ_after_rk = .false.
      end if
    end if
    
    enable_hmat_s_check = .false.
    env_hmat_s_check = ''
    call get_environment_variable("SALMON_DG_RK_HMAT_S_CHECK", env_hmat_s_check, length=env_stage_metric_len, status=env_stage_metric_status)
    if (env_stage_metric_status == 0 .and. env_stage_metric_len > 0) then
      if (env_hmat_s_check(1:1) == '1' .or. env_hmat_s_check(1:1) == 'y' .or. env_hmat_s_check(1:1) == 'Y' .or. &
          env_hmat_s_check(1:1) == 't' .or. env_hmat_s_check(1:1) == 'T') then
        enable_hmat_s_check = .true.
      end if
    end if
    if (itt == 1 .and. dg_frag%id == 0) then
      write(*,'(1x,a,l1,a,l1)') "        DEBUG: RK init enable_rk_stage_metric=", enable_rk_stage_metric, &
        " enable_hmat_s_check=", enable_hmat_s_check
      flush(6)
    end if
    
    ! Reuse RK stage arrays across calls to reduce per-step allocation overhead.
    if (.not. allocated(k)) then
      allocate(k(n, dg_frag%nstate_tot, dg_frag%nspin, dg_frag%rk_stages))
    else if (size(k, 1) /= n .or. size(k, 2) /= dg_frag%nstate_tot .or. &
             size(k, 3) /= dg_frag%nspin .or. size(k, 4) /= dg_frag%rk_stages) then
      deallocate(k)
      allocate(k(n, dg_frag%nstate_tot, dg_frag%nspin, dg_frag%rk_stages))
    end if
    if (n_pw > 0) then
      if (.not. allocated(k_pw)) then
        allocate(k_pw(n_pw, dg_frag%nstate_tot, dg_frag%nspin, dg_frag%rk_stages))
      else if (size(k_pw, 1) /= n_pw .or. size(k_pw, 2) /= dg_frag%nstate_tot .or. &
               size(k_pw, 3) /= dg_frag%nspin .or. size(k_pw, 4) /= dg_frag%rk_stages) then
        deallocate(k_pw)
        allocate(k_pw(n_pw, dg_frag%nstate_tot, dg_frag%nspin, dg_frag%rk_stages))
      end if
    end if
    
    if (dg_frag%time_integrator == 3) then
      ! Classical RK4 (paper-aligned): k1@t, k2@t+dt/2, k3@t+dt/2, k4@t+dt
      owned_s = max(1, dg_frag%owned_coef_start)
      owned_e = min(n, dg_frag%owned_coef_end)
      if (owned_e >= owned_s) then
        n_owned = owned_e - owned_s + 1
      else
        n_owned = 0
      end if
      allocate(coef_ref(max(1, n_owned), dg_frag%nstate_tot, dg_frag%nspin))
      coef_ref = (0.0d0, 0.0d0)
      if (n_owned > 0) then
        coef_ref(1:n_owned, :, :) = dg_frag%coef(owned_s:owned_e, :, :)
      end if
      if (n_pw > 0) then
        allocate(coef_pw_ref(n_pw, dg_frag%nstate_tot, dg_frag%nspin))
        coef_pw_ref = dg_frag%coef_pw
      end if

      ! Stage 1
      Ac_tot = rt%Ac_tot(:, itt)
      if (n_owned > 0) then
        dg_frag%coef(owned_s:owned_e, :, :) = coef_ref(1:n_owned, :, :)
      end if
      if (use_mixed_rt) then
        do ispin = 1, dg_frag%nspin
          call sync_mixed_coef_from_raw(dg_frag, ispin)
        end do
      end if
      call debug_rk_stage_metric("rk4-stage1-state")
      if (enable_rk_trace) then
        write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,a)') "        rk trace: rank=", dg_frag%id, &
          " id_frag=", dg_frag%id_frag, " ifrag_group=", dg_frag%ifrag_group, " itt=", itt, &
          " stage=", "rk4-stage1-entry"
        flush(6)
      end if
      if (yn_fix_func == 'n') then
        if (enable_rk_trace) then
          write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,a)') "        rk trace: rank=", dg_frag%id, &
            " id_frag=", dg_frag%id_frag, " ifrag_group=", dg_frag%ifrag_group, " itt=", itt, &
            " stage=", "rk4-stage1-before-update"
          flush(6)
        end if
        call update_density_hamiltonian_stage(dg_frag, system, info, rt, itt, Ac_tot, &
                                              lg, mg, stencil, xc_func, srg, srg_scalar, fg, poisson, pp, ppg, ppn, &
                                              rho, rho_s, Vh, Vxc, Vpsl, energy)
        if (enable_rk_trace) then
          write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,a)') "        rk trace: rank=", dg_frag%id, &
            " id_frag=", dg_frag%id_frag, " ifrag_group=", dg_frag%ifrag_group, " itt=", itt, &
            " stage=", "rk4-stage1-after-update"
          flush(6)
        end if
      end if
      if (enable_rk_trace) then
        write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,a)') "        rk trace: rank=", dg_frag%id, &
          " id_frag=", dg_frag%id_frag, " ifrag_group=", dg_frag%ifrag_group, " itt=", itt, &
          " stage=", "rk4-stage1-before-derivative"
        flush(6)
      end if
      if (n_pw > 0) then
        call calculate_time_derivative(dg_frag, system, mg, stencil, ppg, Ac_tot, itt, k(:,:,:,1), k_pw(:,:,:,1))
      else
        call calculate_time_derivative(dg_frag, system, mg, stencil, ppg, Ac_tot, itt, k(:,:,:,1))
      end if
      if (enable_rk_trace) then
        write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,a)') "        rk trace: rank=", dg_frag%id, &
          " id_frag=", dg_frag%id_frag, " ifrag_group=", dg_frag%ifrag_group, " itt=", itt, &
          " stage=", "rk4-stage1-after-derivative"
        flush(6)
      end if

      ! Stage 2
      Ac_tot = 0.5d0 * (rt%Ac_tot(:, itt) + rt%Ac_tot(:, itt+1))
      if (enable_rk_trace) then
        write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,a)') "        rk trace: rank=", dg_frag%id, &
          " id_frag=", dg_frag%id_frag, " ifrag_group=", dg_frag%ifrag_group, " itt=", itt, &
          " stage=", "rk4-stage2-entry"
        flush(6)
      end if
      if (n_pw > 0) then
!$omp parallel private(jo,jo_local,jo_global)
!$omp do collapse(2) schedule(static)
        do ispin = 1, dg_frag%nspin
          do io = 1, dg_frag%nstate_tot
!$omp simd
            do jo_local = 1, n_owned
              jo_global = owned_s + jo_local - 1
              dg_frag%coef(jo_global, io, ispin) = coef_ref(jo_local, io, ispin) + 0.5d0 * dt * k(jo_global, io, ispin, 1)
            end do
          end do
        end do
!$omp end do
!$omp do collapse(2) schedule(static)
        do ispin = 1, dg_frag%nspin
          do io = 1, dg_frag%nstate_tot
!$omp simd
            do jo = 1, n_pw
              dg_frag%coef_pw(jo, io, ispin) = coef_pw_ref(jo, io, ispin) + 0.5d0 * dt * k_pw(jo, io, ispin, 1)
            end do
          end do
        end do
!$omp end do
!$omp end parallel
      else
!$omp parallel do collapse(2) private(jo,jo_local,jo_global) schedule(static)
        do ispin = 1, dg_frag%nspin
          do io = 1, dg_frag%nstate_tot
!$omp simd
            do jo_local = 1, n_owned
              jo_global = owned_s + jo_local - 1
              dg_frag%coef(jo_global, io, ispin) = coef_ref(jo_local, io, ispin) + 0.5d0 * dt * k(jo_global, io, ispin, 1)
            end do
          end do
        end do
!$omp end parallel do
      end if
      if (use_mixed_rt) then
        do ispin = 1, dg_frag%nspin
          call sync_mixed_coef_from_raw(dg_frag, ispin)
        end do
      end if
      if (zero_nonowned_after_stage) call zero_nonowned_coefficients(dg_frag)
      call debug_rk_stage_metric("rk4-stage2-state")
      if (yn_fix_func == 'n') then
        if (enable_rk_trace) then
          write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,a)') "        rk trace: rank=", dg_frag%id, &
            " id_frag=", dg_frag%id_frag, " ifrag_group=", dg_frag%ifrag_group, " itt=", itt, &
            " stage=", "rk4-stage2-before-update"
          flush(6)
        end if
        call update_density_hamiltonian_stage(dg_frag, system, info, rt, itt, Ac_tot, &
                                              lg, mg, stencil, xc_func, srg, srg_scalar, fg, poisson, pp, ppg, ppn, &
                                              rho, rho_s, Vh, Vxc, Vpsl, energy)
        if (enable_rk_trace) then
          write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,a)') "        rk trace: rank=", dg_frag%id, &
            " id_frag=", dg_frag%id_frag, " ifrag_group=", dg_frag%ifrag_group, " itt=", itt, &
            " stage=", "rk4-stage2-after-update"
          flush(6)
        end if
      end if
      if (enable_rk_trace) then
        write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,a)') "        rk trace: rank=", dg_frag%id, &
          " id_frag=", dg_frag%id_frag, " ifrag_group=", dg_frag%ifrag_group, " itt=", itt, &
          " stage=", "rk4-stage2-before-derivative"
        flush(6)
      end if
      if (n_pw > 0) then
        call calculate_time_derivative(dg_frag, system, mg, stencil, ppg, Ac_tot, itt, k(:,:,:,2), k_pw(:,:,:,2))
      else
        call calculate_time_derivative(dg_frag, system, mg, stencil, ppg, Ac_tot, itt, k(:,:,:,2))
      end if
      if (enable_rk_trace) then
        write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,a)') "        rk trace: rank=", dg_frag%id, &
          " id_frag=", dg_frag%id_frag, " ifrag_group=", dg_frag%ifrag_group, " itt=", itt, &
          " stage=", "rk4-stage2-after-derivative"
        flush(6)
      end if

      ! Stage 3
      if (n_pw > 0) then
!$omp parallel private(jo,jo_local,jo_global)
!$omp do collapse(2) schedule(static)
        do ispin = 1, dg_frag%nspin
          do io = 1, dg_frag%nstate_tot
!$omp simd
            do jo_local = 1, n_owned
              jo_global = owned_s + jo_local - 1
              dg_frag%coef(jo_global, io, ispin) = coef_ref(jo_local, io, ispin) + 0.5d0 * dt * k(jo_global, io, ispin, 2)
            end do
          end do
        end do
!$omp end do
!$omp do collapse(2) schedule(static)
        do ispin = 1, dg_frag%nspin
          do io = 1, dg_frag%nstate_tot
!$omp simd
            do jo = 1, n_pw
              dg_frag%coef_pw(jo, io, ispin) = coef_pw_ref(jo, io, ispin) + 0.5d0 * dt * k_pw(jo, io, ispin, 2)
            end do
          end do
        end do
!$omp end do
!$omp end parallel
      else
!$omp parallel do collapse(2) private(jo,jo_local,jo_global) schedule(static)
        do ispin = 1, dg_frag%nspin
          do io = 1, dg_frag%nstate_tot
!$omp simd
            do jo_local = 1, n_owned
              jo_global = owned_s + jo_local - 1
              dg_frag%coef(jo_global, io, ispin) = coef_ref(jo_local, io, ispin) + 0.5d0 * dt * k(jo_global, io, ispin, 2)
            end do
          end do
        end do
!$omp end parallel do
      end if
      if (use_mixed_rt) then
        do ispin = 1, dg_frag%nspin
          call sync_mixed_coef_from_raw(dg_frag, ispin)
        end do
      end if
      if (zero_nonowned_after_stage) call zero_nonowned_coefficients(dg_frag)
      call debug_rk_stage_metric("rk4-stage3-state")
      if (yn_fix_func == 'n') then
        call update_density_hamiltonian_stage(dg_frag, system, info, rt, itt, Ac_tot, &
                                              lg, mg, stencil, xc_func, srg, srg_scalar, fg, poisson, pp, ppg, ppn, &
                                              rho, rho_s, Vh, Vxc, Vpsl, energy)
      end if
      if (n_pw > 0) then
        call calculate_time_derivative(dg_frag, system, mg, stencil, ppg, Ac_tot, itt, k(:,:,:,3), k_pw(:,:,:,3))
      else
        call calculate_time_derivative(dg_frag, system, mg, stencil, ppg, Ac_tot, itt, k(:,:,:,3))
      end if

      ! Stage 4
      Ac_tot = rt%Ac_tot(:, itt+1)
      if (n_pw > 0) then
!$omp parallel private(jo,jo_local,jo_global)
!$omp do collapse(2) schedule(static)
        do ispin = 1, dg_frag%nspin
          do io = 1, dg_frag%nstate_tot
!$omp simd
            do jo_local = 1, n_owned
              jo_global = owned_s + jo_local - 1
              dg_frag%coef(jo_global, io, ispin) = coef_ref(jo_local, io, ispin) + dt * k(jo_global, io, ispin, 3)
            end do
          end do
        end do
!$omp end do
!$omp do collapse(2) schedule(static)
        do ispin = 1, dg_frag%nspin
          do io = 1, dg_frag%nstate_tot
!$omp simd
            do jo = 1, n_pw
              dg_frag%coef_pw(jo, io, ispin) = coef_pw_ref(jo, io, ispin) + dt * k_pw(jo, io, ispin, 3)
            end do
          end do
        end do
!$omp end do
!$omp end parallel
      else
!$omp parallel do collapse(2) private(jo,jo_local,jo_global) schedule(static)
        do ispin = 1, dg_frag%nspin
          do io = 1, dg_frag%nstate_tot
!$omp simd
            do jo_local = 1, n_owned
              jo_global = owned_s + jo_local - 1
              dg_frag%coef(jo_global, io, ispin) = coef_ref(jo_local, io, ispin) + dt * k(jo_global, io, ispin, 3)
            end do
          end do
        end do
!$omp end parallel do
      end if
      if (use_mixed_rt) then
        do ispin = 1, dg_frag%nspin
          call sync_mixed_coef_from_raw(dg_frag, ispin)
        end do
      end if
      if (zero_nonowned_after_stage) call zero_nonowned_coefficients(dg_frag)
      call debug_rk_stage_metric("rk4-stage4-state")
      call debug_rk_stage_metric("rk4-final-state")
      if (enable_hmat_s_check .and. itt <= stage_metric_max_itt) then
        call debug_hmat_s_consistency("rk4-stage4-before-update")
      end if
      if (yn_fix_func == 'n') then
        call update_density_hamiltonian_stage(dg_frag, system, info, rt, itt, Ac_tot, &
                                              lg, mg, stencil, xc_func, srg, srg_scalar, fg, poisson, pp, ppg, ppn, &
                                              rho, rho_s, Vh, Vxc, Vpsl, energy)
      end if
      if (enable_hmat_s_check .and. itt <= stage_metric_max_itt) then
        call debug_hmat_s_consistency("rk4-stage4-after-update")
      end if
      if (n_pw > 0) then
        call calculate_time_derivative(dg_frag, system, mg, stencil, ppg, Ac_tot, itt, k(:,:,:,4), k_pw(:,:,:,4))
      else
        call calculate_time_derivative(dg_frag, system, mg, stencil, ppg, Ac_tot, itt, k(:,:,:,4))
      end if

      ! Final RK4 combination
      if (n_pw > 0) then
!$omp parallel private(jo,jo_local,jo_global)
!$omp do collapse(2) schedule(static)
        do ispin = 1, dg_frag%nspin
          do io = 1, dg_frag%nstate_tot
!$omp simd
            do jo_local = 1, n_owned
              jo_global = owned_s + jo_local - 1
              dg_frag%coef(jo_global, io, ispin) = coef_ref(jo_local, io, ispin) + dt * ( &
                                             k(jo_global, io, ispin, 1) + 2.0d0 * k(jo_global, io, ispin, 2) + &
                                             2.0d0 * k(jo_global, io, ispin, 3) + k(jo_global, io, ispin, 4)) / 6.0d0
            end do
          end do
        end do
!$omp end do
!$omp do collapse(2) schedule(static)
        do ispin = 1, dg_frag%nspin
          do io = 1, dg_frag%nstate_tot
!$omp simd
            do jo = 1, n_pw
              dg_frag%coef_pw(jo, io, ispin) = coef_pw_ref(jo, io, ispin) + dt * ( &
                                                k_pw(jo, io, ispin, 1) + 2.0d0 * k_pw(jo, io, ispin, 2) + &
                                                2.0d0 * k_pw(jo, io, ispin, 3) + k_pw(jo, io, ispin, 4)) / 6.0d0
            end do
          end do
        end do
!$omp end do
!$omp end parallel
      else
!$omp parallel do collapse(2) private(jo,jo_local,jo_global) schedule(static)
        do ispin = 1, dg_frag%nspin
          do io = 1, dg_frag%nstate_tot
!$omp simd
            do jo_local = 1, n_owned
              jo_global = owned_s + jo_local - 1
              dg_frag%coef(jo_global, io, ispin) = coef_ref(jo_local, io, ispin) + dt * ( &
                                             k(jo_global, io, ispin, 1) + 2.0d0 * k(jo_global, io, ispin, 2) + &
                                             2.0d0 * k(jo_global, io, ispin, 3) + k(jo_global, io, ispin, 4)) / 6.0d0
            end do
          end do
        end do
!$omp end parallel do
      end if
      if (use_mixed_rt) then
        do ispin = 1, dg_frag%nspin
          call sync_mixed_coef_from_raw(dg_frag, ispin)
        end do
      end if
      if (zero_nonowned_after_stage) call zero_nonowned_coefficients(dg_frag)
      if (renormalize_occ_after_rk) call renormalize_occupied_norms_after_rk()
      if (yn_fix_func == 'n') then
        ! Rebuild H at the final RK4 state/time so next step starts from consistent rho/H.
        Ac_tot = rt%Ac_tot(:, itt+1)
        call update_density_hamiltonian_stage(dg_frag, system, info, rt, itt, Ac_tot, &
                                              lg, mg, stencil, xc_func, srg, srg_scalar, fg, poisson, pp, ppg, ppn, &
                                              rho, rho_s, Vh, Vxc, Vpsl, energy)
      end if
      deallocate(coef_ref)
      if (allocated(coef_pw_ref)) deallocate(coef_pw_ref)

    else
      ! SSPRK3 stages.
      ! Store initial coefficients for Shu-Osher blending.
      if (.not. allocated(dg_frag%coef_work)) then
        allocate(dg_frag%coef_work(n, dg_frag%nstate_tot, dg_frag%nspin))
      end if
      dg_frag%coef_work = dg_frag%coef
      ! Save initial PW coefficients for the Shu-Osher alpha*coef_work term
      ! (analogous to dg_frag%coef_work which is already set above for fragment coef).
      if (n_pw > 0) then
        allocate(coef_pw_ref(n_pw, dg_frag%nstate_tot, dg_frag%nspin))
        coef_pw_ref = dg_frag%coef_pw
      end if

      if (n_pw > 0) then
        do istage = 1, dg_frag%rk_stages
          ! Get vector potential at this time (velocity gauge)
          ! For RK stages, interpolate between itt and itt+1
          if (istage == 1) then
            Ac_tot = rt%Ac_tot(:, itt)
          else
            t_stage = dble(istage-1) / dble(dg_frag%rk_stages)
            Ac_tot = (1.0d0 - t_stage) * rt%Ac_tot(:, itt) + t_stage * rt%Ac_tot(:, itt+1)
          end if

          if (enable_rk_trace) then
            write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0,a,a)') "        rk trace: rank=", dg_frag%id, &
              " id_frag=", dg_frag%id_frag, " ifrag_group=", dg_frag%ifrag_group, " itt=", itt, &
              " stage=", istage, " phase=", "ssprk-before-derivative(mixed)"
            flush(6)
          end if

          ! Calculate time derivative: d/dt coef = -i*(H_0 + A^2/2)*coef + A·<∇>*coef
          ! In velocity gauge: H(t) = H_0 - i*A(t)·∇ + A(t)^2/2
          call calculate_time_derivative(dg_frag, system, mg, stencil, ppg, Ac_tot, itt, k(:,:,:,istage), k_pw(:,:,:,istage))

          if (enable_rk_trace) then
            write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0,a,a)') "        rk trace: rank=", dg_frag%id, &
              " id_frag=", dg_frag%id_frag, " ifrag_group=", dg_frag%ifrag_group, " itt=", itt, &
              " stage=", istage, " phase=", "ssprk-after-derivative(mixed)"
            flush(6)
          end if

          ! Update coefficients for next stage
          if (istage < dg_frag%rk_stages) then
            ! OpenMP parallelization for coefficient update
!$omp parallel private(jo,jo_local,jo_global)
!$omp do collapse(2) schedule(static)
            do ispin = 1, dg_frag%nspin
              do io = 1, dg_frag%nstate_tot
!$omp simd
                do jo = 1, n
                  dg_frag%coef(jo, io, ispin) = dg_frag%rk_alpha(istage) * dg_frag%coef_work(jo, io, ispin) + &
                                                 dg_frag%rk_beta(istage) * dg_frag%coef(jo, io, ispin) + &
                                                 dg_frag%rk_gamma(istage) * dt * k(jo, io, ispin, istage)
                end do
              end do
            end do
!$omp end do
!$omp do collapse(2) schedule(static)
            do ispin = 1, dg_frag%nspin
              do io = 1, dg_frag%nstate_tot
!$omp simd
                do jo = 1, n_pw
                  ! BUG FIX: was coef_pw += gamma*dt*k_pw (missing alpha/beta terms).
                  ! Apply full Shu-Osher formula, consistent with fragment update above.
                  dg_frag%coef_pw(jo, io, ispin) = dg_frag%rk_alpha(istage) * coef_pw_ref(jo, io, ispin) + &
                                                   dg_frag%rk_beta(istage)  * dg_frag%coef_pw(jo, io, ispin) + &
                                                   dg_frag%rk_gamma(istage) * dt * k_pw(jo, io, ispin, istage)
                end do
              end do
            end do
!$omp end do
!$omp end parallel
            if (use_mixed_rt) then
              do ispin = 1, dg_frag%nspin
                call sync_mixed_coef_from_raw(dg_frag, ispin)
              end do
            end if
            if (enable_rk_trace) then
              write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0,a,a)') "        rk trace: rank=", dg_frag%id, &
                " id_frag=", dg_frag%id_frag, " ifrag_group=", dg_frag%ifrag_group, " itt=", itt, &
                " stage=", istage, " phase=", "ssprk-after-update(mixed)"
              flush(6)
            end if
          end if
        end do
      else
        do istage = 1, dg_frag%rk_stages
          ! Get vector potential at this time (velocity gauge)
          ! For RK stages, interpolate between itt and itt+1
          if (istage == 1) then
            Ac_tot = rt%Ac_tot(:, itt)
          else
            t_stage = dble(istage-1) / dble(dg_frag%rk_stages)
            Ac_tot = (1.0d0 - t_stage) * rt%Ac_tot(:, itt) + t_stage * rt%Ac_tot(:, itt+1)
          end if

          if (enable_rk_trace) then
            write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0,a,a)') "        rk trace: rank=", dg_frag%id, &
              " id_frag=", dg_frag%id_frag, " ifrag_group=", dg_frag%ifrag_group, " itt=", itt, &
              " stage=", istage, " phase=", "ssprk-before-derivative"
            flush(6)
          end if

          ! Calculate time derivative: d/dt coef = -i*(H_0 + A^2/2)*coef + A·<∇>*coef
          ! In velocity gauge: H(t) = H_0 - i*A(t)·∇ + A(t)^2/2
          call calculate_time_derivative(dg_frag, system, mg, stencil, ppg, Ac_tot, itt, k(:,:,:,istage))

          if (enable_rk_trace) then
            write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0,a,a)') "        rk trace: rank=", dg_frag%id, &
              " id_frag=", dg_frag%id_frag, " ifrag_group=", dg_frag%ifrag_group, " itt=", itt, &
              " stage=", istage, " phase=", "ssprk-after-derivative"
            flush(6)
          end if

          ! Update coefficients for next stage
          if (istage < dg_frag%rk_stages) then
            ! OpenMP parallelization for coefficient update
!$omp parallel do collapse(2) private(jo,jo_local,jo_global) schedule(static)
            do ispin = 1, dg_frag%nspin
              do io = 1, dg_frag%nstate_tot
!$omp simd
                do jo = 1, n
                  dg_frag%coef(jo, io, ispin) = dg_frag%rk_alpha(istage) * dg_frag%coef_work(jo, io, ispin) + &
                                                 dg_frag%rk_beta(istage) * dg_frag%coef(jo, io, ispin) + &
                                                 dg_frag%rk_gamma(istage) * dt * k(jo, io, ispin, istage)
                end do
              end do
            end do
!$omp end parallel do
            if (use_mixed_rt) then
              do ispin = 1, dg_frag%nspin
                call sync_mixed_coef_from_raw(dg_frag, ispin)
              end do
            end if
            if (enable_rk_trace) then
              write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0,a,a)') "        rk trace: rank=", dg_frag%id, &
                " id_frag=", dg_frag%id_frag, " ifrag_group=", dg_frag%ifrag_group, " itt=", itt, &
                " stage=", istage, " phase=", "ssprk-after-update"
              flush(6)
            end if
          end if
        end do
      end if
      
      ! Final update: apply Shu-Osher formula for stage rk_stages.
      !
      ! BUG FIX (fragment): the previous code reset coef = coef_work then accumulated
      !   coef += Σ_s gamma_s * dt * k_s  →  coef_work + γ1*dt*k1 + γ2*dt*k2 + γ3*dt*k3
      ! which is WRONG.  Expanding SSPRK3 (α=[1,3/4,1/3], β=[0,1/4,2/3], γ=[1,1/4,2/3])
      ! the correct final result is
      !   u^{n+1} = 1/3*u^0 + 2/3*u^{(2)} + 2/3*dt*k3
      ! which is exactly the Shu-Osher formula for stage rk_stages applied to the
      ! current coef (= u^{(rk_stages-1)} after the last intermediate update).
      !
      ! BUG FIX (PW): final step for coef_pw was missing entirely; added here.
      associate(rs => dg_frag%rk_stages)
!$omp parallel private(jo,jo_local,jo_global)
!$omp do collapse(2) schedule(static)
      do ispin = 1, dg_frag%nspin
        do io = 1, dg_frag%nstate_tot
!$omp simd
          do jo = 1, n
            dg_frag%coef(jo, io, ispin) = dg_frag%rk_alpha(rs) * dg_frag%coef_work(jo, io, ispin) + &
                                          dg_frag%rk_beta(rs)  * dg_frag%coef(jo, io, ispin) + &
                                          dg_frag%rk_gamma(rs) * dt * k(jo, io, ispin, rs)
          end do
        end do
      end do
!$omp end do
      if (n_pw > 0) then
!$omp do collapse(2) schedule(static)
        do ispin = 1, dg_frag%nspin
          do io = 1, dg_frag%nstate_tot
!$omp simd
            do jo = 1, n_pw
              dg_frag%coef_pw(jo, io, ispin) = dg_frag%rk_alpha(rs) * coef_pw_ref(jo, io, ispin) + &
                                               dg_frag%rk_beta(rs)  * dg_frag%coef_pw(jo, io, ispin) + &
                                               dg_frag%rk_gamma(rs) * dt * k_pw(jo, io, ispin, rs)
            end do
          end do
        end do
!$omp end do
      end if
!$omp end parallel
      end associate
      if (use_mixed_rt) then
        do ispin = 1, dg_frag%nspin
          call sync_mixed_coef_from_raw(dg_frag, ispin)
        end do
      end if
      if (allocated(coef_pw_ref)) deallocate(coef_pw_ref)
    end if

    if (enable_rk_nan_check) then
      found_nan = .false.
      nan_jo = 0
      nan_io = 0
      nan_ispin = 0
      nan_val = (0.0d0, 0.0d0)
      do ispin = 1, dg_frag%nspin
        do io = 1, dg_frag%nstate_tot
          do jo = 1, n
            if (real(dg_frag%coef(jo, io, ispin)) /= real(dg_frag%coef(jo, io, ispin)) .or. &
                aimag(dg_frag%coef(jo, io, ispin)) /= aimag(dg_frag%coef(jo, io, ispin))) then
              found_nan = .true.
              nan_jo = jo
              nan_io = io
              nan_ispin = ispin
              nan_val = dg_frag%coef(jo, io, ispin)
              exit
            end if
          end do
          if (found_nan) exit
        end do
        if (found_nan) exit
      end do
      if (found_nan) then
        write(*,'(a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,2es12.4)') "[NaN] coef detected (final): rank=", dg_frag%id, &
          " itt=", itt, " stage=", dg_frag%rk_stages, " ispin=", nan_ispin, " io=", nan_io, " jo=", nan_jo, &
          " val=", real(nan_val), aimag(nan_val)
        stop "NaN in coef"
      end if
    end if

  contains

  subroutine renormalize_occupied_norms_after_rk()
    integer :: ispin_norm, io_norm, jo_norm, n_frag_norm, n_pw_norm, n_tot_norm, nocc_norm
    real(8) :: norm2_norm, norm2_local, norm_scale
    logical :: use_s_norm
    complex(8), allocatable :: vec_norm(:), svec_norm(:)
    complex(8), allocatable :: coef_frag_norm(:,:), coef_pw_norm(:,:)

    do ispin_norm = 1, dg_frag%nspin
      n_frag_norm = dg_frag%n_mat(ispin_norm)
      n_pw_norm = 0
      if (dg_frag%use_plane_wave_basis .and. allocated(dg_frag%coef_pw)) n_pw_norm = dg_frag%n_plane_waves
      use_s_norm = allocated(dg_frag%S_mat_frag_pw) .or. allocated(dg_frag%S_mat_prop_c) .or. &
                   allocated(dg_frag%S_mat_prop) .or. allocated(dg_frag%S_mat_c) .or. allocated(dg_frag%S_mat) .or. &
                   allocated(dg_frag%S_mat_prop_blocks) .or. allocated(dg_frag%S_mat_blocks)
      if (n_pw_norm > 0 .or. use_s_norm) n_frag_norm = dg_frag%n_mat_max
      n_tot_norm = n_frag_norm + n_pw_norm
      if (n_tot_norm <= 0) cycle
      nocc_norm = dg_frag%nstate_tot
      if (allocated(dg_frag%nocc_spin) .and. ispin_norm <= size(dg_frag%nocc_spin)) then
        nocc_norm = min(nocc_norm, max(0, dg_frag%nocc_spin(ispin_norm)))
      end if
      if (nocc_norm <= 0) cycle

      call gather_full_coef_view(dg_frag, ispin_norm, n_frag_norm, nocc_norm, coef_frag_norm, coef_pw_norm, 1, nocc_norm)
      allocate(vec_norm(n_tot_norm), svec_norm(n_tot_norm))
      do io_norm = 1, nocc_norm
        vec_norm(:) = (0.0d0, 0.0d0)
        vec_norm(1:n_frag_norm) = coef_frag_norm(1:n_frag_norm, io_norm)
        if (n_pw_norm > 0) vec_norm(n_frag_norm+1:n_tot_norm) = coef_pw_norm(1:n_pw_norm, io_norm)
        if (use_s_norm) then
          call apply_overlap_operator(dg_frag, ispin_norm, vec_norm, svec_norm, .true.)
          norm2_local = real(sum(conjg(vec_norm) * svec_norm), kind=8)
          call comm_summation(norm2_local, norm2_norm, dg_frag%icomm)
        else
          norm2_norm = sum(abs(vec_norm)**2)
        end if
        if (norm2_norm <= 1.0d-24) cycle
        norm_scale = 1.0d0 / sqrt(norm2_norm)
        coef_frag_norm(1:n_frag_norm, io_norm) = coef_frag_norm(1:n_frag_norm, io_norm) * norm_scale
        if (n_pw_norm > 0) coef_pw_norm(1:n_pw_norm, io_norm) = coef_pw_norm(1:n_pw_norm, io_norm) * norm_scale
      end do

      do io_norm = 1, nocc_norm
        do jo_norm = 1, n_frag_norm
          if (allocated(dg_frag%coef_owner)) then
            if (dg_frag%coef_owner(jo_norm, ispin_norm) /= dg_frag%id) cycle
          end if
          dg_frag%coef(jo_norm, io_norm, ispin_norm) = coef_frag_norm(jo_norm, io_norm)
        end do
        do jo_norm = 1, n_pw_norm
          if (allocated(dg_frag%coef_pw_owner)) then
            if (dg_frag%coef_pw_owner(jo_norm) /= dg_frag%id) cycle
          end if
          dg_frag%coef_pw(jo_norm, io_norm, ispin_norm) = coef_pw_norm(jo_norm, io_norm)
        end do
      end do
      deallocate(vec_norm, svec_norm)
      if (allocated(coef_frag_norm)) deallocate(coef_frag_norm)
      if (allocated(coef_pw_norm)) deallocate(coef_pw_norm)
    end do
    call zero_nonowned_coefficients(dg_frag)
  end subroutine renormalize_occupied_norms_after_rk

  subroutine debug_rk_stage_metric(stage_label)
    character(len=*), intent(in) :: stage_label
    integer :: ispin_probe, io, n_frag_rows, n_pw_probe, n_tot_probe, nocc_probe
    real(8) :: cs2_sum, cs2_min, cs2_max, sval
    complex(8), allocatable :: vec(:), svec(:)
    complex(8), allocatable :: coef_frag_probe(:,:), coef_pw_probe(:,:)

    if (.not. enable_rk_stage_metric) return
    if (itt > stage_metric_max_itt) return
    if (dg_frag%nspin <= 0 .or. dg_frag%nstate_tot <= 0) return

    do ispin_probe = 1, dg_frag%nspin
      nocc_probe = min(dg_frag%nstate_tot, dg_frag%nocc_spin(ispin_probe))
      if (nocc_probe <= 0) cycle
      n_frag_rows = dg_frag%n_mat_max
      n_pw_probe = 0
      if (dg_frag%use_plane_wave_basis .and. allocated(dg_frag%coef_pw)) n_pw_probe = dg_frag%n_plane_waves
      n_tot_probe = n_frag_rows + n_pw_probe
      if (n_tot_probe <= 0) cycle

      call gather_full_coef_view(dg_frag, ispin_probe, n_frag_rows, nocc_probe, coef_frag_probe, coef_pw_probe, 1, nocc_probe)
      allocate(vec(n_tot_probe), svec(n_tot_probe))
      cs2_sum = 0.0d0
      cs2_min = huge(1.0d0)
      cs2_max = -huge(1.0d0)
      do io = 1, nocc_probe
        vec(:) = (0.0d0, 0.0d0)
        vec(1:n_frag_rows) = coef_frag_probe(1:n_frag_rows, io)
        if (n_pw_probe > 0) vec(n_frag_rows+1:n_tot_probe) = coef_pw_probe(1:n_pw_probe, io)
        call apply_overlap_operator(dg_frag, ispin_probe, vec, svec, .true.)
        sval = real(sum(conjg(vec) * svec), kind=8)
        cs2_sum = cs2_sum + sval
        cs2_min = min(cs2_min, sval)
        cs2_max = max(cs2_max, sval)
      end do
      if (dg_frag%id == 0) then
        write(*,'(1x,a,a,a,i0,a,i0,a,3(1x,es12.4))') "        rk stage metric: stage=", trim(stage_label), &
          " ispin=", ispin_probe, " nocc=", nocc_probe, " cs2(sum,min,max)=", cs2_sum, cs2_min, cs2_max
        flush(6)
      end if
      deallocate(vec, svec)
      if (allocated(coef_frag_probe)) deallocate(coef_frag_probe)
      if (allocated(coef_pw_probe)) deallocate(coef_pw_probe)
    end do
  end subroutine debug_rk_stage_metric

  subroutine debug_hmat_s_consistency(label)
    character(len=*), intent(in) :: label
    integer :: ispin, io, jo, nmat, max_inspect
    real(8) :: h_diag_min, h_diag_max, s_diag_min, s_diag_max, h_off_max, s_off_max
    real(8) :: h_trace, s_trace, ratio_hs
    logical :: hmat_ok, smat_ok
    
    if (dg_frag%id /= 0) return  ! Only rank 0 prints
    
    ! Defensive checks
    if (dg_frag%nspin <= 0) return
    if (.not. allocated(dg_frag%H_mat)) return
    if (.not. allocated(dg_frag%S_mat)) return
    
    nmat = size(dg_frag%H_mat, 1)
    if (nmat <= 0) return
    if (size(dg_frag%H_mat, 3) /= dg_frag%nspin) return
    if (size(dg_frag%S_mat, 3) /= dg_frag%nspin) return
    if (size(dg_frag%S_mat, 1) /= nmat) return
    if (size(dg_frag%S_mat, 2) /= nmat) return
    
    max_inspect = min(nmat, 200)  ! Limit inspection to first 200 states
    
    do ispin = 1, dg_frag%nspin
      ! Compute diagnostics for this spin
      h_diag_min = huge(1.0d0)
      h_diag_max = -huge(1.0d0)
      s_diag_min = huge(1.0d0)
      s_diag_max = -huge(1.0d0)
      h_off_max = 0.0d0
      s_off_max = 0.0d0
      h_trace = 0.0d0
      s_trace = 0.0d0
      
      do io = 1, nmat
        h_trace = h_trace + dg_frag%H_mat(io, io, ispin)
        s_trace = s_trace + dg_frag%S_mat(io, io, ispin)
        h_diag_min = min(h_diag_min, dg_frag%H_mat(io, io, ispin))
        h_diag_max = max(h_diag_max, dg_frag%H_mat(io, io, ispin))
        s_diag_min = min(s_diag_min, dg_frag%S_mat(io, io, ispin))
        s_diag_max = max(s_diag_max, dg_frag%S_mat(io, io, ispin))
      end do
      
      do io = 1, max_inspect
        do jo = io+1, max_inspect
          h_off_max = max(h_off_max, abs(dg_frag%H_mat(io, jo, ispin)))
          s_off_max = max(s_off_max, abs(dg_frag%S_mat(io, jo, ispin)))
        end do
      end do
      
      hmat_ok = (h_diag_max > 0.0d0)
      smat_ok = (s_diag_min > 1.0d-10 .and. s_diag_max < 1.0d2)
      ratio_hs = merge(h_trace / s_trace, 0.0d0, s_trace /= 0.0d0)
      
      write(*,'(1x,a,a,a,i0,a,i0)') "        hmat-s-check: label=", trim(label), &
        " ispin=", ispin, " nmat=", nmat
      write(*,'(1x,a,4(1x,es12.4))') "          H_diag(min,max) S_diag(min,max):", &
        h_diag_min, h_diag_max, s_diag_min, s_diag_max
      write(*,'(1x,a,2(1x,es12.4))') "          H_offmax(first200) S_offmax(first200):", h_off_max, s_off_max
      write(*,'(1x,a,2(1x,es12.4),1x,l1,1x,l1)') "          H_trace S_trace H/S_ratio hmat_ok smat_ok:", &
        h_trace, s_trace, ratio_hs, hmat_ok, smat_ok
      flush(6)
    end do
  end subroutine debug_hmat_s_consistency
    
  end subroutine time_evolution_rk
