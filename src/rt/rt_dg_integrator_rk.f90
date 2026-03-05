  subroutine time_evolution_rk(dg_frag, system, info, rt, itt, dt, &
                             lg, mg, stencil, xc_func, srg, srg_scalar, fg, poisson, pp, ppg, ppn, &
                             rho, rho_s, Vh, Vxc, Vpsl, energy)
    use structures
    use salmon_global, only: yn_fix_func
    use sendrecv_grid, only: s_sendrecv_grid
    use salmon_xc, only: s_xc_functional
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
    
    complex(8), allocatable :: k(:,:,:,:)
    complex(8), allocatable :: k_pw(:,:,:,:)
    complex(8), allocatable :: coef_ref(:,:,:)
    complex(8), allocatable :: coef_pw_ref(:,:,:)
    integer :: istage, io, jo, ispin
    real(8) :: Ac_tot(3), t_stage
    integer :: n, n_pw
    logical :: found_nan
    integer :: nan_jo, nan_io, nan_ispin
    complex(8) :: nan_val
    
    ! Use n_mat_max (global basis size) instead of nstate_frag (local basis size)
    n = dg_frag%n_mat_max
    n_pw = 0
    if (dg_frag%use_plane_wave_basis .and. allocated(dg_frag%coef_pw)) n_pw = dg_frag%n_plane_waves
    
    ! Allocate RK stage arrays (complex for proper phase evolution)
    allocate(k(n, dg_frag%nstate_tot, dg_frag%nspin, dg_frag%rk_stages))
    if (n_pw > 0) allocate(k_pw(n_pw, dg_frag%nstate_tot, dg_frag%nspin, dg_frag%rk_stages))
    
    ! Store initial coefficients
    dg_frag%coef_work = dg_frag%coef
    
    if (dg_frag%time_integrator == 3) then
      ! Classical RK4 (paper-aligned): k1@t, k2@t+dt/2, k3@t+dt/2, k4@t+dt
      allocate(coef_ref(n, dg_frag%nstate_tot, dg_frag%nspin))
      coef_ref = dg_frag%coef
      if (n_pw > 0) then
        allocate(coef_pw_ref(n_pw, dg_frag%nstate_tot, dg_frag%nspin))
        coef_pw_ref = dg_frag%coef_pw
      end if

      ! Stage 1
      Ac_tot = rt%Ac_tot(:, itt)
      dg_frag%coef = coef_ref
      if (yn_fix_func == 'n') then
        call update_density_hamiltonian_stage(dg_frag, system, info, rt, itt, Ac_tot, &
                                              lg, mg, stencil, xc_func, srg, srg_scalar, fg, poisson, pp, ppg, ppn, &
                                              rho, rho_s, Vh, Vxc, Vpsl, energy)
      end if
      if (n_pw > 0) then
        call calculate_time_derivative(dg_frag, system, mg, stencil, ppg, Ac_tot, itt, k(:,:,:,1), k_pw(:,:,:,1))
      else
        call calculate_time_derivative(dg_frag, system, mg, stencil, ppg, Ac_tot, itt, k(:,:,:,1))
      end if

      ! Stage 2
      Ac_tot = 0.5d0 * (rt%Ac_tot(:, itt) + rt%Ac_tot(:, itt+1))
      if (n_pw > 0) then
!$omp parallel do private(io,jo) schedule(static)
        do ispin = 1, dg_frag%nspin
          do io = 1, dg_frag%nstate_tot
            do jo = 1, n
              dg_frag%coef(jo, io, ispin) = coef_ref(jo, io, ispin) + 0.5d0 * dt * k(jo, io, ispin, 1)
            end do
          end do
          do io = 1, dg_frag%nstate_tot
            do jo = 1, n_pw
              dg_frag%coef_pw(jo, io, ispin) = coef_pw_ref(jo, io, ispin) + 0.5d0 * dt * k_pw(jo, io, ispin, 1)
            end do
          end do
        end do
!$omp end parallel do
      else
!$omp parallel do private(io,jo) schedule(static)
        do ispin = 1, dg_frag%nspin
          do io = 1, dg_frag%nstate_tot
            do jo = 1, n
              dg_frag%coef(jo, io, ispin) = coef_ref(jo, io, ispin) + 0.5d0 * dt * k(jo, io, ispin, 1)
            end do
          end do
        end do
!$omp end parallel do
      end if
      if (yn_fix_func == 'n') then
        call update_density_hamiltonian_stage(dg_frag, system, info, rt, itt, Ac_tot, &
                                              lg, mg, stencil, xc_func, srg, srg_scalar, fg, poisson, pp, ppg, ppn, &
                                              rho, rho_s, Vh, Vxc, Vpsl, energy)
      end if
      if (n_pw > 0) then
        call calculate_time_derivative(dg_frag, system, mg, stencil, ppg, Ac_tot, itt, k(:,:,:,2), k_pw(:,:,:,2))
      else
        call calculate_time_derivative(dg_frag, system, mg, stencil, ppg, Ac_tot, itt, k(:,:,:,2))
      end if

      ! Stage 3
      if (n_pw > 0) then
!$omp parallel do private(io,jo) schedule(static)
        do ispin = 1, dg_frag%nspin
          do io = 1, dg_frag%nstate_tot
            do jo = 1, n
              dg_frag%coef(jo, io, ispin) = coef_ref(jo, io, ispin) + 0.5d0 * dt * k(jo, io, ispin, 2)
            end do
          end do
          do io = 1, dg_frag%nstate_tot
            do jo = 1, n_pw
              dg_frag%coef_pw(jo, io, ispin) = coef_pw_ref(jo, io, ispin) + 0.5d0 * dt * k_pw(jo, io, ispin, 2)
            end do
          end do
        end do
!$omp end parallel do
      else
!$omp parallel do private(io,jo) schedule(static)
        do ispin = 1, dg_frag%nspin
          do io = 1, dg_frag%nstate_tot
            do jo = 1, n
              dg_frag%coef(jo, io, ispin) = coef_ref(jo, io, ispin) + 0.5d0 * dt * k(jo, io, ispin, 2)
            end do
          end do
        end do
!$omp end parallel do
      end if
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
!$omp parallel do private(io,jo) schedule(static)
        do ispin = 1, dg_frag%nspin
          do io = 1, dg_frag%nstate_tot
            do jo = 1, n
              dg_frag%coef(jo, io, ispin) = coef_ref(jo, io, ispin) + dt * k(jo, io, ispin, 3)
            end do
          end do
          do io = 1, dg_frag%nstate_tot
            do jo = 1, n_pw
              dg_frag%coef_pw(jo, io, ispin) = coef_pw_ref(jo, io, ispin) + dt * k_pw(jo, io, ispin, 3)
            end do
          end do
        end do
!$omp end parallel do
      else
!$omp parallel do private(io,jo) schedule(static)
        do ispin = 1, dg_frag%nspin
          do io = 1, dg_frag%nstate_tot
            do jo = 1, n
              dg_frag%coef(jo, io, ispin) = coef_ref(jo, io, ispin) + dt * k(jo, io, ispin, 3)
            end do
          end do
        end do
!$omp end parallel do
      end if
      if (yn_fix_func == 'n') then
        call update_density_hamiltonian_stage(dg_frag, system, info, rt, itt, Ac_tot, &
                                              lg, mg, stencil, xc_func, srg, srg_scalar, fg, poisson, pp, ppg, ppn, &
                                              rho, rho_s, Vh, Vxc, Vpsl, energy)
      end if
      if (n_pw > 0) then
        call calculate_time_derivative(dg_frag, system, mg, stencil, ppg, Ac_tot, itt, k(:,:,:,4), k_pw(:,:,:,4))
      else
        call calculate_time_derivative(dg_frag, system, mg, stencil, ppg, Ac_tot, itt, k(:,:,:,4))
      end if

      ! Final RK4 combination
      if (n_pw > 0) then
!$omp parallel do private(io,jo) schedule(static)
        do ispin = 1, dg_frag%nspin
          do io = 1, dg_frag%nstate_tot
            do jo = 1, n
              dg_frag%coef(jo, io, ispin) = coef_ref(jo, io, ispin) + dt * ( &
                                             k(jo, io, ispin, 1) + 2.0d0 * k(jo, io, ispin, 2) + &
                                             2.0d0 * k(jo, io, ispin, 3) + k(jo, io, ispin, 4)) / 6.0d0
            end do
          end do
          do io = 1, dg_frag%nstate_tot
            do jo = 1, n_pw
              dg_frag%coef_pw(jo, io, ispin) = coef_pw_ref(jo, io, ispin) + dt * ( &
                                                k_pw(jo, io, ispin, 1) + 2.0d0 * k_pw(jo, io, ispin, 2) + &
                                                2.0d0 * k_pw(jo, io, ispin, 3) + k_pw(jo, io, ispin, 4)) / 6.0d0
            end do
          end do
        end do
!$omp end parallel do
      else
!$omp parallel do private(io,jo) schedule(static)
        do ispin = 1, dg_frag%nspin
          do io = 1, dg_frag%nstate_tot
            do jo = 1, n
              dg_frag%coef(jo, io, ispin) = coef_ref(jo, io, ispin) + dt * ( &
                                             k(jo, io, ispin, 1) + 2.0d0 * k(jo, io, ispin, 2) + &
                                             2.0d0 * k(jo, io, ispin, 3) + k(jo, io, ispin, 4)) / 6.0d0
            end do
          end do
        end do
!$omp end parallel do
      end if
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
      ! SSPRK3 stages
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
        if (n_pw > 0) then
          call calculate_time_derivative(dg_frag, system, mg, stencil, ppg, Ac_tot, itt, k(:,:,:,istage), k_pw(:,:,:,istage))
        else
          call calculate_time_derivative(dg_frag, system, mg, stencil, ppg, Ac_tot, itt, k(:,:,:,istage))
        end if
        
        ! Update coefficients for next stage
        if (istage < dg_frag%rk_stages) then
          ! OpenMP parallelization for coefficient update
          if (n_pw > 0) then
!$omp parallel do private(io,jo) schedule(static)
            do ispin = 1, dg_frag%nspin
              do io = 1, dg_frag%nstate_tot
                do jo = 1, n
                  dg_frag%coef(jo, io, ispin) = dg_frag%rk_alpha(istage) * dg_frag%coef_work(jo, io, ispin) + &
                                                 dg_frag%rk_beta(istage) * dg_frag%coef(jo, io, ispin) + &
                                                 dg_frag%rk_gamma(istage) * dt * k(jo, io, ispin, istage)
                end do
              end do
              do io = 1, dg_frag%nstate_tot
                do jo = 1, n_pw
                  dg_frag%coef_pw(jo, io, ispin) = dg_frag%coef_pw(jo, io, ispin) + &
                                                   dg_frag%rk_gamma(istage) * dt * k_pw(jo, io, ispin, istage)
                end do
              end do
            end do
!$omp end parallel do
          else
!$omp parallel do private(io,jo) schedule(static)
            do ispin = 1, dg_frag%nspin
              do io = 1, dg_frag%nstate_tot
                do jo = 1, n
                  dg_frag%coef(jo, io, ispin) = dg_frag%rk_alpha(istage) * dg_frag%coef_work(jo, io, ispin) + &
                                                 dg_frag%rk_beta(istage) * dg_frag%coef(jo, io, ispin) + &
                                                 dg_frag%rk_gamma(istage) * dt * k(jo, io, ispin, istage)
                end do
              end do
            end do
!$omp end parallel do
          end if
        end if
      end do
      
      ! Final update
      dg_frag%coef = dg_frag%coef_work
      ! OpenMP parallelization for final coefficient update
!$omp parallel do collapse(2) private(io,jo,istage) schedule(static)
      do ispin = 1, dg_frag%nspin
        do io = 1, dg_frag%nstate_tot
          do jo = 1, n
            do istage = 1, dg_frag%rk_stages
              dg_frag%coef(jo, io, ispin) = dg_frag%coef(jo, io, ispin) + &
                                            dg_frag%rk_gamma(istage) * dt * k(jo, io, ispin, istage)
            end do
          end do
        end do
      end do
!$omp end parallel do
    end if

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
    
    deallocate(k)
    if (allocated(k_pw)) deallocate(k_pw)
    
  end subroutine time_evolution_rk
