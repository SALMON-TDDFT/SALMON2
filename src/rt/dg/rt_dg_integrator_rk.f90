  subroutine time_evolution_rk(dg_frag, system, info, rt, itt, dt, &
                             lg, mg, stencil, xc_func, srg, srg_scalar, fg, poisson, pp, ppg, ppn, &
                             rho, rho_s, Vh, Vxc, Vpsl, energy)
    use structures
    use salmon_global, only: yn_fix_func
    use sendrecv_grid, only: s_sendrecv_grid
    use salmon_xc, only: s_xc_functional
    use rt_dg_fragment_ops, only: sync_mixed_coef_from_raw, sync_raw_coef_from_mixed
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
    
    complex(8), allocatable, save :: k(:,:,:), k_accum(:,:,:)
    complex(8), allocatable, save :: k_pw(:,:,:), k_pw_accum(:,:,:)
    complex(8), allocatable :: coef_ref(:,:,:)
    complex(8), allocatable :: coef_pw_ref(:,:,:)
    integer :: istage, io, jo, ispin
    real(8) :: Ac_tot(3), t_stage
    real(8) :: t0, t1
    real(8) :: time_sync, time_stage_update, time_deriv, time_coef_update
    integer :: n, n_pw
    logical :: use_mixed_rt
    logical :: trace_first_step
    logical, save :: timing_initialized = .false.
    logical, save :: enable_rk_timing = .false.
    character(16) :: env_timing
    integer :: env_status
    
    ! Coefficients are stored only for the fragment rows local to this rank.
    n = size(dg_frag%coef, 1)
    n_pw = 0
    if (dg_frag%use_plane_wave_basis .and. allocated(dg_frag%coef_pw)) n_pw = dg_frag%n_plane_waves
    use_mixed_rt = (n_pw > 0 .and. dg_frag%mixed_basis_ready .and. allocated(dg_frag%coef_mix))
    trace_first_step = (itt == 1 .and. dg_frag%id == 0)
    if (trace_first_step) then
      write(*,'(1x,a,i0,4(a,i0),a,l1)') '[DG-RK] enter itt=', itt, &
        ' local_basis_rows=', n, ' nstate_tot=', dg_frag%nstate_tot, &
        ' n_pw=', n_pw, ' nspin=', dg_frag%nspin, ' use_mixed=', use_mixed_rt
      flush(6)
    end if
    if (.not. timing_initialized) then
      env_timing = ''
      call get_environment_variable('SALMON_DG_RK_TIMING', env_timing, status=env_status)
      if (env_status == 0) then
        select case(trim(adjustl(env_timing)))
        case('1','y','Y','yes','YES','true','TRUE','on','ON')
          enable_rk_timing = .true.
        end select
      end if
      timing_initialized = .true.
    end if
    time_sync = 0.0d0
    time_stage_update = 0.0d0
    time_deriv = 0.0d0
    time_coef_update = 0.0d0
    
    ! Reuse RK work arrays across calls.  Classical RK4 only needs the current
    ! stage derivative and its weighted sum, not all four stage derivatives.
    if (.not. allocated(k)) then
      allocate(k(n, dg_frag%nstate_tot, dg_frag%nspin))
      allocate(k_accum(n, dg_frag%nstate_tot, dg_frag%nspin))
    else if (size(k, 1) /= n .or. size(k, 2) /= dg_frag%nstate_tot .or. &
             size(k, 3) /= dg_frag%nspin) then
      deallocate(k, k_accum)
      allocate(k(n, dg_frag%nstate_tot, dg_frag%nspin))
      allocate(k_accum(n, dg_frag%nstate_tot, dg_frag%nspin))
    end if
    if (n_pw > 0) then
      if (.not. allocated(k_pw)) then
        allocate(k_pw(n_pw, dg_frag%nstate_tot, dg_frag%nspin))
        allocate(k_pw_accum(n_pw, dg_frag%nstate_tot, dg_frag%nspin))
      else if (size(k_pw, 1) /= n_pw .or. size(k_pw, 2) /= dg_frag%nstate_tot .or. &
               size(k_pw, 3) /= dg_frag%nspin) then
        deallocate(k_pw, k_pw_accum)
        allocate(k_pw(n_pw, dg_frag%nstate_tot, dg_frag%nspin))
        allocate(k_pw_accum(n_pw, dg_frag%nstate_tot, dg_frag%nspin))
      end if
    end if
    if (trace_first_step) then
      write(*,'(1x,a)') '[DG-RK] work arrays ready'
      flush(6)
    end if
    
    if (dg_frag%time_integrator == 3) then
      ! Classical RK4 (paper-aligned): k1@t, k2@t+dt/2, k3@t+dt/2, k4@t+dt
      if (use_mixed_rt) then
        ! Orbital-parallel matrix/density construction splits basis/state work
        ! across ranks, so the RK reference state must be the canonical mixed
        ! view on every rank before any stage-local column work starts.
        call cpu_time(t0)
        do ispin = 1, dg_frag%nspin
          call sync_mixed_coef_from_raw(dg_frag, ispin)
        end do
        do ispin = 1, dg_frag%nspin
          call sync_raw_coef_from_mixed(dg_frag, ispin)
        end do
        call cpu_time(t1)
        time_sync = time_sync + (t1 - t0)
      end if
      allocate(coef_ref(n, dg_frag%nstate_tot, dg_frag%nspin))
      coef_ref = dg_frag%coef
      if (n_pw > 0) then
        allocate(coef_pw_ref(n_pw, dg_frag%nstate_tot, dg_frag%nspin))
        coef_pw_ref = dg_frag%coef_pw
      end if
      if (trace_first_step) then
        write(*,'(1x,a)') '[DG-RK] RK4 reference copied'
        flush(6)
      end if

      ! Stage 1
      Ac_tot = rt%Ac_tot(:, itt)
      dg_frag%coef = coef_ref
      ! The mixed/raw views were canonicalized immediately before coef_ref was
      ! captured, so Stage 1 can reuse that state without another sync pair.
      if (yn_fix_func == 'n') then
        if (trace_first_step) then
          write(*,'(1x,a)') '[DG-RK] stage 1 density/H update start'
          flush(6)
        end if
        call cpu_time(t0)
        call update_density_hamiltonian_stage(dg_frag, system, info, rt, itt, Ac_tot, &
                                              lg, mg, stencil, xc_func, srg, srg_scalar, fg, poisson, pp, ppg, ppn, &
                                              rho, rho_s, Vh, Vxc, Vpsl, energy)
        call cpu_time(t1)
        time_stage_update = time_stage_update + (t1 - t0)
        if (trace_first_step) then
          write(*,'(1x,a,1pe12.4)') '[DG-RK] stage 1 density/H update done time=', t1 - t0
          flush(6)
        end if
      end if
      if (trace_first_step) then
        write(*,'(1x,a)') '[DG-RK] stage 1 derivative start'
        flush(6)
      end if
      call cpu_time(t0)
      if (n_pw > 0) then
        call calculate_time_derivative(dg_frag, system, mg, stencil, ppg, Ac_tot, itt, k, k_pw)
      else
        call calculate_time_derivative(dg_frag, system, mg, stencil, ppg, Ac_tot, itt, k)
      end if
      call cpu_time(t1)
      time_deriv = time_deriv + (t1 - t0)
      if (trace_first_step) then
        write(*,'(1x,a,1pe12.4)') '[DG-RK] stage 1 derivative done time=', t1 - t0
        flush(6)
      end if
      k_accum(:, :, :) = k(:, :, :)
      if (n_pw > 0) k_pw_accum(:, :, :) = k_pw(:, :, :)

      ! Stage 2
      Ac_tot = 0.5d0 * (rt%Ac_tot(:, itt) + rt%Ac_tot(:, itt+1))
      call cpu_time(t0)
      if (n_pw > 0) then
!$omp parallel private(jo)
!$omp do collapse(2) schedule(static)
        do ispin = 1, dg_frag%nspin
          do io = 1, dg_frag%nstate_tot
!$omp simd
            do jo = 1, n
              dg_frag%coef(jo, io, ispin) = coef_ref(jo, io, ispin) + 0.5d0 * dt * k(jo, io, ispin)
            end do
          end do
        end do
!$omp end do
!$omp do collapse(2) schedule(static)
        do ispin = 1, dg_frag%nspin
          do io = 1, dg_frag%nstate_tot
!$omp simd
            do jo = 1, n_pw
              dg_frag%coef_pw(jo, io, ispin) = coef_pw_ref(jo, io, ispin) + 0.5d0 * dt * k_pw(jo, io, ispin)
            end do
          end do
        end do
!$omp end do
!$omp end parallel
      else
!$omp parallel do collapse(2) private(jo) schedule(static)
        do ispin = 1, dg_frag%nspin
          do io = 1, dg_frag%nstate_tot
!$omp simd
            do jo = 1, n
              dg_frag%coef(jo, io, ispin) = coef_ref(jo, io, ispin) + 0.5d0 * dt * k(jo, io, ispin)
            end do
          end do
        end do
!$omp end parallel do
      end if
      call cpu_time(t1)
      time_coef_update = time_coef_update + (t1 - t0)
      if (use_mixed_rt) then
        call cpu_time(t0)
        ! The provisional RK update changed raw coefficients; rebuild coef_mix
        ! and then refresh raw coefficients from the canonical mixed view.
        do ispin = 1, dg_frag%nspin
          call sync_mixed_coef_from_raw(dg_frag, ispin)
        end do
        do ispin = 1, dg_frag%nspin
          call sync_raw_coef_from_mixed(dg_frag, ispin)
        end do
        call cpu_time(t1)
        time_sync = time_sync + (t1 - t0)
      end if
      if (yn_fix_func == 'n') then
        call cpu_time(t0)
        call update_density_hamiltonian_stage(dg_frag, system, info, rt, itt, Ac_tot, &
                                              lg, mg, stencil, xc_func, srg, srg_scalar, fg, poisson, pp, ppg, ppn, &
                                              rho, rho_s, Vh, Vxc, Vpsl, energy)
        call cpu_time(t1)
        time_stage_update = time_stage_update + (t1 - t0)
      end if
      call cpu_time(t0)
      if (n_pw > 0) then
        call calculate_time_derivative(dg_frag, system, mg, stencil, ppg, Ac_tot, itt, k, k_pw)
      else
        call calculate_time_derivative(dg_frag, system, mg, stencil, ppg, Ac_tot, itt, k)
      end if
      call cpu_time(t1)
      time_deriv = time_deriv + (t1 - t0)
!$omp parallel do collapse(2) private(jo) schedule(static)
      do ispin = 1, dg_frag%nspin
        do io = 1, dg_frag%nstate_tot
!$omp simd
          do jo = 1, n
            k_accum(jo, io, ispin) = k_accum(jo, io, ispin) + 2.0d0 * k(jo, io, ispin)
          end do
        end do
      end do
!$omp end parallel do
      if (n_pw > 0) then
!$omp parallel do collapse(2) private(jo) schedule(static)
        do ispin = 1, dg_frag%nspin
          do io = 1, dg_frag%nstate_tot
!$omp simd
            do jo = 1, n_pw
              k_pw_accum(jo, io, ispin) = k_pw_accum(jo, io, ispin) + 2.0d0 * k_pw(jo, io, ispin)
            end do
          end do
        end do
!$omp end parallel do
      end if

      ! Stage 3
      call cpu_time(t0)
      if (n_pw > 0) then
!$omp parallel private(jo)
!$omp do collapse(2) schedule(static)
        do ispin = 1, dg_frag%nspin
          do io = 1, dg_frag%nstate_tot
!$omp simd
            do jo = 1, n
              dg_frag%coef(jo, io, ispin) = coef_ref(jo, io, ispin) + 0.5d0 * dt * k(jo, io, ispin)
            end do
          end do
        end do
!$omp end do
!$omp do collapse(2) schedule(static)
        do ispin = 1, dg_frag%nspin
          do io = 1, dg_frag%nstate_tot
!$omp simd
            do jo = 1, n_pw
              dg_frag%coef_pw(jo, io, ispin) = coef_pw_ref(jo, io, ispin) + 0.5d0 * dt * k_pw(jo, io, ispin)
            end do
          end do
        end do
!$omp end do
!$omp end parallel
      else
!$omp parallel do collapse(2) private(jo) schedule(static)
        do ispin = 1, dg_frag%nspin
          do io = 1, dg_frag%nstate_tot
!$omp simd
            do jo = 1, n
              dg_frag%coef(jo, io, ispin) = coef_ref(jo, io, ispin) + 0.5d0 * dt * k(jo, io, ispin)
            end do
          end do
        end do
!$omp end parallel do
      end if
      call cpu_time(t1)
      time_coef_update = time_coef_update + (t1 - t0)
      if (use_mixed_rt) then
        call cpu_time(t0)
        do ispin = 1, dg_frag%nspin
          call sync_mixed_coef_from_raw(dg_frag, ispin)
        end do
        do ispin = 1, dg_frag%nspin
          call sync_raw_coef_from_mixed(dg_frag, ispin)
        end do
        call cpu_time(t1)
        time_sync = time_sync + (t1 - t0)
      end if
      if (yn_fix_func == 'n') then
        call cpu_time(t0)
        call update_density_hamiltonian_stage(dg_frag, system, info, rt, itt, Ac_tot, &
                                              lg, mg, stencil, xc_func, srg, srg_scalar, fg, poisson, pp, ppg, ppn, &
                                              rho, rho_s, Vh, Vxc, Vpsl, energy)
        call cpu_time(t1)
        time_stage_update = time_stage_update + (t1 - t0)
      end if
      call cpu_time(t0)
      if (n_pw > 0) then
        call calculate_time_derivative(dg_frag, system, mg, stencil, ppg, Ac_tot, itt, k, k_pw)
      else
        call calculate_time_derivative(dg_frag, system, mg, stencil, ppg, Ac_tot, itt, k)
      end if
      call cpu_time(t1)
      time_deriv = time_deriv + (t1 - t0)
!$omp parallel do collapse(2) private(jo) schedule(static)
      do ispin = 1, dg_frag%nspin
        do io = 1, dg_frag%nstate_tot
!$omp simd
          do jo = 1, n
            k_accum(jo, io, ispin) = k_accum(jo, io, ispin) + 2.0d0 * k(jo, io, ispin)
          end do
        end do
      end do
!$omp end parallel do
      if (n_pw > 0) then
!$omp parallel do collapse(2) private(jo) schedule(static)
        do ispin = 1, dg_frag%nspin
          do io = 1, dg_frag%nstate_tot
!$omp simd
            do jo = 1, n_pw
              k_pw_accum(jo, io, ispin) = k_pw_accum(jo, io, ispin) + 2.0d0 * k_pw(jo, io, ispin)
            end do
          end do
        end do
!$omp end parallel do
      end if

      ! Stage 4
      Ac_tot = rt%Ac_tot(:, itt+1)
      call cpu_time(t0)
      if (n_pw > 0) then
!$omp parallel private(jo)
!$omp do collapse(2) schedule(static)
        do ispin = 1, dg_frag%nspin
          do io = 1, dg_frag%nstate_tot
!$omp simd
            do jo = 1, n
              dg_frag%coef(jo, io, ispin) = coef_ref(jo, io, ispin) + dt * k(jo, io, ispin)
            end do
          end do
        end do
!$omp end do
!$omp do collapse(2) schedule(static)
        do ispin = 1, dg_frag%nspin
          do io = 1, dg_frag%nstate_tot
!$omp simd
            do jo = 1, n_pw
              dg_frag%coef_pw(jo, io, ispin) = coef_pw_ref(jo, io, ispin) + dt * k_pw(jo, io, ispin)
            end do
          end do
        end do
!$omp end do
!$omp end parallel
      else
!$omp parallel do collapse(2) private(jo) schedule(static)
        do ispin = 1, dg_frag%nspin
          do io = 1, dg_frag%nstate_tot
!$omp simd
            do jo = 1, n
              dg_frag%coef(jo, io, ispin) = coef_ref(jo, io, ispin) + dt * k(jo, io, ispin)
            end do
          end do
        end do
!$omp end parallel do
      end if
      call cpu_time(t1)
      time_coef_update = time_coef_update + (t1 - t0)
      if (use_mixed_rt) then
        call cpu_time(t0)
        do ispin = 1, dg_frag%nspin
          call sync_mixed_coef_from_raw(dg_frag, ispin)
        end do
        do ispin = 1, dg_frag%nspin
          call sync_raw_coef_from_mixed(dg_frag, ispin)
        end do
        call cpu_time(t1)
        time_sync = time_sync + (t1 - t0)
      end if
      if (yn_fix_func == 'n') then
        call cpu_time(t0)
        call update_density_hamiltonian_stage(dg_frag, system, info, rt, itt, Ac_tot, &
                                              lg, mg, stencil, xc_func, srg, srg_scalar, fg, poisson, pp, ppg, ppn, &
                                              rho, rho_s, Vh, Vxc, Vpsl, energy)
        call cpu_time(t1)
        time_stage_update = time_stage_update + (t1 - t0)
      end if
      call cpu_time(t0)
      if (n_pw > 0) then
        call calculate_time_derivative(dg_frag, system, mg, stencil, ppg, Ac_tot, itt, k, k_pw)
      else
        call calculate_time_derivative(dg_frag, system, mg, stencil, ppg, Ac_tot, itt, k)
      end if
      call cpu_time(t1)
      time_deriv = time_deriv + (t1 - t0)

      ! Final RK4 combination
      call cpu_time(t0)
      if (n_pw > 0) then
!$omp parallel private(jo)
!$omp do collapse(2) schedule(static)
        do ispin = 1, dg_frag%nspin
          do io = 1, dg_frag%nstate_tot
!$omp simd
            do jo = 1, n
              dg_frag%coef(jo, io, ispin) = coef_ref(jo, io, ispin) + &
                                             (dt / 6.0d0) * (k_accum(jo, io, ispin) + k(jo, io, ispin))
            end do
          end do
        end do
!$omp end do
!$omp do collapse(2) schedule(static)
        do ispin = 1, dg_frag%nspin
          do io = 1, dg_frag%nstate_tot
!$omp simd
            do jo = 1, n_pw
              dg_frag%coef_pw(jo, io, ispin) = coef_pw_ref(jo, io, ispin) + &
                                                (dt / 6.0d0) * (k_pw_accum(jo, io, ispin) + k_pw(jo, io, ispin))
            end do
          end do
        end do
!$omp end do
!$omp end parallel
      else
!$omp parallel do collapse(2) private(jo) schedule(static)
        do ispin = 1, dg_frag%nspin
          do io = 1, dg_frag%nstate_tot
!$omp simd
            do jo = 1, n
              dg_frag%coef(jo, io, ispin) = coef_ref(jo, io, ispin) + &
                                             (dt / 6.0d0) * (k_accum(jo, io, ispin) + k(jo, io, ispin))
            end do
          end do
        end do
!$omp end parallel do
      end if
      call cpu_time(t1)
      time_coef_update = time_coef_update + (t1 - t0)
      if (use_mixed_rt) then
        call cpu_time(t0)
        do ispin = 1, dg_frag%nspin
          call sync_mixed_coef_from_raw(dg_frag, ispin)
        end do
        do ispin = 1, dg_frag%nspin
          call sync_raw_coef_from_mixed(dg_frag, ispin)
        end do
        call cpu_time(t1)
        time_sync = time_sync + (t1 - t0)
      end if
      if (yn_fix_func == 'n') then
        ! Rebuild H at the final RK4 state/time so next step starts from consistent rho/H.
        Ac_tot = rt%Ac_tot(:, itt+1)
        call cpu_time(t0)
        call update_density_hamiltonian_stage(dg_frag, system, info, rt, itt, Ac_tot, &
                                              lg, mg, stencil, xc_func, srg, srg_scalar, fg, poisson, pp, ppg, ppn, &
                                              rho, rho_s, Vh, Vxc, Vpsl, energy)
        call cpu_time(t1)
        time_stage_update = time_stage_update + (t1 - t0)
      end if
      if (enable_rk_timing .and. dg_frag%id == 0) then
        write(*,'(1x,a,i0,4(a,1pe12.4))') '        rk timing: itt=', itt, &
          ' sync=', time_sync, ' coef=', time_coef_update, ' stage_update=', time_stage_update, ' deriv=', time_deriv
        flush(6)
      end if
      deallocate(coef_ref)
      if (allocated(coef_pw_ref)) deallocate(coef_pw_ref)

    else
      ! SSPRK3 stages.
      ! Store initial coefficients for Shu-Osher blending.
      if (use_mixed_rt) then
        ! Keep the saved Shu-Osher reference in the same canonical coefficient
        ! view used by the orbital-parallel matrix and density builders.
        do ispin = 1, dg_frag%nspin
          call sync_mixed_coef_from_raw(dg_frag, ispin)
        end do
        do ispin = 1, dg_frag%nspin
          call sync_raw_coef_from_mixed(dg_frag, ispin)
        end do
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

          ! Calculate time derivative: d/dt coef = -i*(H_0 + A^2/2)*coef + A·<∇>*coef
          ! In velocity gauge: H(t) = H_0 - i*A(t)·∇ + A(t)^2/2
          call calculate_time_derivative(dg_frag, system, mg, stencil, ppg, Ac_tot, itt, k, k_pw)

          ! Update coefficients for next stage
          if (istage < dg_frag%rk_stages) then
            ! OpenMP parallelization for coefficient update
!$omp parallel private(jo)
!$omp do collapse(2) schedule(static)
            do ispin = 1, dg_frag%nspin
              do io = 1, dg_frag%nstate_tot
!$omp simd
                do jo = 1, n
                  dg_frag%coef(jo, io, ispin) = dg_frag%rk_alpha(istage) * dg_frag%coef_work(jo, io, ispin) + &
                                                 dg_frag%rk_beta(istage) * dg_frag%coef(jo, io, ispin) + &
                                                 dg_frag%rk_gamma(istage) * dt * k(jo, io, ispin)
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
                                                   dg_frag%rk_gamma(istage) * dt * k_pw(jo, io, ispin)
                end do
              end do
            end do
!$omp end do
!$omp end parallel
            if (use_mixed_rt) then
              do ispin = 1, dg_frag%nspin
                call sync_mixed_coef_from_raw(dg_frag, ispin)
              end do
              do ispin = 1, dg_frag%nspin
                call sync_raw_coef_from_mixed(dg_frag, ispin)
              end do
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

          ! Calculate time derivative: d/dt coef = -i*(H_0 + A^2/2)*coef + A·<∇>*coef
          ! In velocity gauge: H(t) = H_0 - i*A(t)·∇ + A(t)^2/2
          call calculate_time_derivative(dg_frag, system, mg, stencil, ppg, Ac_tot, itt, k)

          ! Update coefficients for next stage
          if (istage < dg_frag%rk_stages) then
            ! OpenMP parallelization for coefficient update
!$omp parallel do collapse(2) private(jo) schedule(static)
            do ispin = 1, dg_frag%nspin
              do io = 1, dg_frag%nstate_tot
!$omp simd
                do jo = 1, n
                  dg_frag%coef(jo, io, ispin) = dg_frag%rk_alpha(istage) * dg_frag%coef_work(jo, io, ispin) + &
                                                 dg_frag%rk_beta(istage) * dg_frag%coef(jo, io, ispin) + &
                                                 dg_frag%rk_gamma(istage) * dt * k(jo, io, ispin)
                end do
              end do
            end do
!$omp end parallel do
            if (use_mixed_rt) then
              do ispin = 1, dg_frag%nspin
                call sync_mixed_coef_from_raw(dg_frag, ispin)
              end do
              do ispin = 1, dg_frag%nspin
                call sync_raw_coef_from_mixed(dg_frag, ispin)
              end do
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
!$omp parallel private(jo)
!$omp do collapse(2) schedule(static)
      do ispin = 1, dg_frag%nspin
        do io = 1, dg_frag%nstate_tot
!$omp simd
          do jo = 1, n
            dg_frag%coef(jo, io, ispin) = dg_frag%rk_alpha(rs) * dg_frag%coef_work(jo, io, ispin) + &
                                          dg_frag%rk_beta(rs)  * dg_frag%coef(jo, io, ispin) + &
                                          dg_frag%rk_gamma(rs) * dt * k(jo, io, ispin)
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
                                               dg_frag%rk_gamma(rs) * dt * k_pw(jo, io, ispin)
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
        do ispin = 1, dg_frag%nspin
          call sync_raw_coef_from_mixed(dg_frag, ispin)
        end do
      end if
      if (allocated(coef_pw_ref)) deallocate(coef_pw_ref)
    end if
  end subroutine time_evolution_rk
