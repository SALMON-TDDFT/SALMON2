  subroutine time_evolution_rk(dg_frag, system, info, rt, itt, dt, &
                             lg, mg, stencil, xc_func, srg, srg_scalar, fg, poisson, pp, ppg, ppn, &
                             rho, rho_s, Vh, Vxc, Vpsl, energy)
    use structures
    use salmon_global, only: yn_fix_func
    use sendrecv_grid, only: s_sendrecv_grid
    use salmon_xc, only: s_xc_functional
    use rt_dg_fragment_types, only: s_dg_fragment_rt
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
    
    complex(8), allocatable, save :: k_blk(:,:,:)
    complex(8), allocatable, save :: next_coef(:,:,:)
    complex(8), allocatable, save :: acc(:,:,:)
    integer :: istage
    integer :: state0, state_s, state_e, nstate_blk, nstate_work, nstate_prop
    real(8) :: Ac_tot(3), t_stage
    integer :: n
    integer, parameter :: state_work_target_mb = 512
    integer, parameter :: state_work_vectors = 5
    integer(8) :: target_bytes, bytes_per_state
    logical, parameter :: trace_first_step = .false.
    
    ! Coefficients are rank-distributed in orbital mode.  The derivative
    ! gathers global rows in state blocks internally and scatters back here.
    n = size(dg_frag%coef, 1)
    if (dg_frag%use_plane_wave_basis .or. allocated(dg_frag%coef_pw)) then
      stop "DG RK now supports the pure fragment block-sparse route only"
    end if

    nstate_prop = dg_frag%nstate_tot
    if (allocated(dg_frag%nocc_spin)) then
      nstate_prop = min(dg_frag%nstate_tot, max(1, maxval(dg_frag%nocc_spin(1:dg_frag%nspin))))
    end if

    if (.not. allocated(dg_frag%coef_work)) then
      allocate(dg_frag%coef_work(n, dg_frag%nstate_tot, dg_frag%nspin))
    else if (size(dg_frag%coef_work, 1) /= n .or. size(dg_frag%coef_work, 2) /= dg_frag%nstate_tot .or. &
             size(dg_frag%coef_work, 3) /= dg_frag%nspin) then
      deallocate(dg_frag%coef_work)
      allocate(dg_frag%coef_work(n, dg_frag%nstate_tot, dg_frag%nspin))
    end if

    target_bytes = int(state_work_target_mb, kind=8) * 1024_8 * 1024_8
    bytes_per_state = 16_8 * int(max(1, dg_frag%n_mat_max), kind=8) * int(max(1, dg_frag%nspin), kind=8) * &
                      int(state_work_vectors, kind=8)
    nstate_work = max(1, min(nstate_prop, int(max(1_8, target_bytes / max(1_8, bytes_per_state)))))

    if (.not. allocated(k_blk)) then
      allocate(k_blk(n, nstate_work, dg_frag%nspin))
    else if (size(k_blk, 1) /= n .or. size(k_blk, 2) /= nstate_work .or. size(k_blk, 3) /= dg_frag%nspin) then
      deallocate(k_blk)
      allocate(k_blk(n, nstate_work, dg_frag%nspin))
    end if
    if (.not. allocated(next_coef)) then
      allocate(next_coef(n, nstate_prop, dg_frag%nspin))
    else if (size(next_coef, 1) /= n .or. size(next_coef, 2) /= nstate_prop .or. &
             size(next_coef, 3) /= dg_frag%nspin) then
      deallocate(next_coef)
      allocate(next_coef(n, nstate_prop, dg_frag%nspin))
    end if

    if (dg_frag%time_integrator == 3) then
      ! Classical RK4 (paper-aligned): k1@t, k2@t+dt/2, k3@t+dt/2, k4@t+dt
      dg_frag%coef_work(:, 1:nstate_prop, :) = dg_frag%coef(:, 1:nstate_prop, :)
      if (.not. allocated(acc)) then
        allocate(acc(n, nstate_prop, dg_frag%nspin))
      else if (size(acc, 1) /= n .or. size(acc, 2) /= nstate_prop .or. size(acc, 3) /= dg_frag%nspin) then
        deallocate(acc)
        allocate(acc(n, nstate_prop, dg_frag%nspin))
      end if
      acc(:, :, :) = (0.0d0, 0.0d0)

      ! Stage 1
      Ac_tot = rt%Ac_tot(:, itt)
      dg_frag%coef(:, 1:nstate_prop, :) = dg_frag%coef_work(:, 1:nstate_prop, :)
      if (trace_first_step) then
        write(*,'(1x,a,i0,a,i0)') '[DG-RT-FIRST] enter RK4 step=', itt, ' nstate_work=', nstate_work
        flush(6)
      end if
      if (yn_fix_func == 'n') then
        call update_density_hamiltonian_stage(dg_frag, system, info, rt, itt, Ac_tot, &
                                              lg, mg, stencil, xc_func, srg, srg_scalar, fg, poisson, pp, ppg, ppn, &
                                              rho, rho_s, Vh, Vxc, Vpsl, energy)
      end if
      next_coef(:, :, :) = dg_frag%coef_work(:, 1:nstate_prop, :)
      do state0 = 1, nstate_prop, nstate_work
        nstate_blk = min(nstate_work, nstate_prop - state0 + 1)
        state_s = state0
        state_e = state0 + nstate_blk - 1
        if (trace_first_step .and. state0 == 1) then
          write(*,'(1x,a,i0,a,i0)') '[DG-RT-FIRST] stage=1 derivative begin state_e=', state_e, ' nlocal=', n
          flush(6)
        end if
        call calculate_time_derivative(dg_frag, system, mg, ppg, Ac_tot, k_blk(:, 1:nstate_blk, :), state_s, state_e)
        call check_finite_derivative_block(1, state_s, state_e, k_blk(:, 1:nstate_blk, :))
        if (trace_first_step .and. state0 == 1) then
          write(*,'(1x,a)') '[DG-RT-FIRST] stage=1 derivative done'
          flush(6)
        end if
        acc(:, state_s:state_e, :) = acc(:, state_s:state_e, :) + k_blk(:, 1:nstate_blk, :) / 6.0d0
        next_coef(:, state_s:state_e, :) = dg_frag%coef_work(:, state_s:state_e, :) + &
                                           0.5d0 * dt * k_blk(:, 1:nstate_blk, :)
      end do
      dg_frag%coef(:, 1:nstate_prop, :) = next_coef(:, :, :)

      ! Stage 2
      Ac_tot = 0.5d0 * (rt%Ac_tot(:, itt) + rt%Ac_tot(:, itt+1))
      if (yn_fix_func == 'n') then
        call update_density_hamiltonian_stage(dg_frag, system, info, rt, itt, Ac_tot, &
                                              lg, mg, stencil, xc_func, srg, srg_scalar, fg, poisson, pp, ppg, ppn, &
                                              rho, rho_s, Vh, Vxc, Vpsl, energy)
      end if
      next_coef(:, :, :) = dg_frag%coef_work(:, 1:nstate_prop, :)
      do state0 = 1, nstate_prop, nstate_work
        nstate_blk = min(nstate_work, nstate_prop - state0 + 1)
        state_s = state0
        state_e = state0 + nstate_blk - 1
        if (trace_first_step .and. state0 == 1) then
          write(*,'(1x,a)') '[DG-RT-FIRST] stage=2 derivative begin'
          flush(6)
        end if
        call calculate_time_derivative(dg_frag, system, mg, ppg, Ac_tot, k_blk(:, 1:nstate_blk, :), state_s, state_e)
        call check_finite_derivative_block(2, state_s, state_e, k_blk(:, 1:nstate_blk, :))
        if (trace_first_step .and. state0 == 1) then
          write(*,'(1x,a)') '[DG-RT-FIRST] stage=2 derivative done'
          flush(6)
        end if
        acc(:, state_s:state_e, :) = acc(:, state_s:state_e, :) + k_blk(:, 1:nstate_blk, :) / 3.0d0
        next_coef(:, state_s:state_e, :) = dg_frag%coef_work(:, state_s:state_e, :) + &
                                           0.5d0 * dt * k_blk(:, 1:nstate_blk, :)
      end do
      dg_frag%coef(:, 1:nstate_prop, :) = next_coef(:, :, :)

      ! Stage 3
      if (yn_fix_func == 'n') then
        call update_density_hamiltonian_stage(dg_frag, system, info, rt, itt, Ac_tot, &
                                              lg, mg, stencil, xc_func, srg, srg_scalar, fg, poisson, pp, ppg, ppn, &
                                              rho, rho_s, Vh, Vxc, Vpsl, energy)
      end if
      next_coef(:, :, :) = dg_frag%coef_work(:, 1:nstate_prop, :)
      do state0 = 1, nstate_prop, nstate_work
        nstate_blk = min(nstate_work, nstate_prop - state0 + 1)
        state_s = state0
        state_e = state0 + nstate_blk - 1
        call calculate_time_derivative(dg_frag, system, mg, ppg, Ac_tot, k_blk(:, 1:nstate_blk, :), state_s, state_e)
        call check_finite_derivative_block(3, state_s, state_e, k_blk(:, 1:nstate_blk, :))
        acc(:, state_s:state_e, :) = acc(:, state_s:state_e, :) + k_blk(:, 1:nstate_blk, :) / 3.0d0
        next_coef(:, state_s:state_e, :) = dg_frag%coef_work(:, state_s:state_e, :) + &
                                           dt * k_blk(:, 1:nstate_blk, :)
      end do
      dg_frag%coef(:, 1:nstate_prop, :) = next_coef(:, :, :)

      ! Stage 4
      Ac_tot = rt%Ac_tot(:, itt+1)
      if (yn_fix_func == 'n') then
        call update_density_hamiltonian_stage(dg_frag, system, info, rt, itt, Ac_tot, &
                                              lg, mg, stencil, xc_func, srg, srg_scalar, fg, poisson, pp, ppg, ppn, &
                                              rho, rho_s, Vh, Vxc, Vpsl, energy)
      end if
      next_coef(:, :, :) = dg_frag%coef_work(:, 1:nstate_prop, :)
      do state0 = 1, nstate_prop, nstate_work
        nstate_blk = min(nstate_work, nstate_prop - state0 + 1)
        state_s = state0
        state_e = state0 + nstate_blk - 1
        call calculate_time_derivative(dg_frag, system, mg, ppg, Ac_tot, k_blk(:, 1:nstate_blk, :), state_s, state_e)
        call check_finite_derivative_block(4, state_s, state_e, k_blk(:, 1:nstate_blk, :))
        acc(:, state_s:state_e, :) = acc(:, state_s:state_e, :) + k_blk(:, 1:nstate_blk, :) / 6.0d0
        next_coef(:, state_s:state_e, :) = dg_frag%coef_work(:, state_s:state_e, :) + &
                                           dt * acc(:, state_s:state_e, :)
      end do
      dg_frag%coef(:, 1:nstate_prop, :) = next_coef(:, :, :)
      if (yn_fix_func == 'n') then
        ! Rebuild H at the final RK4 state/time so next step starts from consistent rho/H.
        Ac_tot = rt%Ac_tot(:, itt+1)
        call update_density_hamiltonian_stage(dg_frag, system, info, rt, itt, Ac_tot, &
                                              lg, mg, stencil, xc_func, srg, srg_scalar, fg, poisson, pp, ppg, ppn, &
                                              rho, rho_s, Vh, Vxc, Vpsl, energy)
      end if

    else
      ! SSPRK3 stages.
      ! Store initial coefficients for Shu-Osher blending.
      dg_frag%coef_work(:, 1:nstate_prop, :) = dg_frag%coef(:, 1:nstate_prop, :)

      do istage = 1, dg_frag%rk_stages
        ! Get vector potential at this time (velocity gauge)
        ! For RK stages, interpolate between itt and itt+1
        if (istage == 1) then
          Ac_tot = rt%Ac_tot(:, itt)
        else
          t_stage = dble(istage-1) / dble(dg_frag%rk_stages)
          Ac_tot = (1.0d0 - t_stage) * rt%Ac_tot(:, itt) + t_stage * rt%Ac_tot(:, itt+1)
        end if

        next_coef(:, :, :) = dg_frag%coef_work(:, 1:nstate_prop, :)
        do state0 = 1, nstate_prop, nstate_work
          nstate_blk = min(nstate_work, nstate_prop - state0 + 1)
          state_s = state0
          state_e = state0 + nstate_blk - 1
          call calculate_time_derivative(dg_frag, system, mg, ppg, Ac_tot, k_blk(:, 1:nstate_blk, :), state_s, state_e)
          call check_finite_derivative_block(istage, state_s, state_e, k_blk(:, 1:nstate_blk, :))
          next_coef(:, state_s:state_e, :) = dg_frag%rk_alpha(istage) * dg_frag%coef_work(:, state_s:state_e, :) + &
                                             dg_frag%rk_beta(istage)  * dg_frag%coef(:, state_s:state_e, :) + &
                                             dg_frag%rk_gamma(istage) * dt * k_blk(:, 1:nstate_blk, :)
        end do
        dg_frag%coef(:, 1:nstate_prop, :) = next_coef(:, :, :)
      end do
    end if

  contains

    subroutine check_finite_derivative_block(stage_id, state_s_in, state_e_in, block)
      integer, intent(in) :: stage_id, state_s_in, state_e_in
      complex(8), intent(in) :: block(:, :, :)

      if (any(real(block, kind=8) /= real(block, kind=8)) .or. &
          any(aimag(block) /= aimag(block))) then
        write(*,'(1x,a,i0,a,i0,a,i0,a,i0)') '[FATAL] NaN in DG RK derivative: rank=', &
          dg_frag%id, ' stage=', stage_id, ' state_s=', state_s_in, ' state_e=', state_e_in
        stop 'DG-Fragment RT: NaN derivative block'
      end if
    end subroutine check_finite_derivative_block

  end subroutine time_evolution_rk
