  subroutine time_evolution_aetrs(dg_frag, system, mg, stencil, ppg, rt, itt, dt)
    use structures
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system),     intent(in)    :: system
    type(s_rgrid),          intent(in)    :: mg
    type(s_stencil),        intent(in)    :: stencil
    type(s_pp_grid),        intent(in)    :: ppg
    type(s_rt),             intent(inout) :: rt
    integer,                intent(in)    :: itt
    real(8),                intent(in)    :: dt

    ! Approximated ETRS (AETRS) for the fragment basis.
    !
    ! Algorithm (frozen-H midpoint RK2, 2nd-order accurate):
    !   1. k1 = f(coef(t),  A(t))      — derivative at current state with A(t)
    !   2. coef_mid = coef(t) + dt/2 * k1    — Euler half-step (predictor)
    !   3. k_mid = f(coef_mid, A(t+dt/2))   — derivative at midpoint with A_mid
    !   4. coef(t+dt) = coef(t) + dt * k_mid — corrector
    !
    ! This uses A(t+dt/2) = (A(t)+A(t+dt))/2 for the external field and the
    ! current Hamiltonian H_mat (frozen; H is updated by the self-consistent loop
    ! outside this subroutine).  The scheme is 2nd-order and approximately
    ! time-reversal symmetric in the frozen-H limit.

    integer :: n, n_pw
    real(8) :: Ac_t(3), Ac_mid(3)
    complex(8), allocatable :: coef_save(:,:,:)
    complex(8), allocatable :: k1(:,:,:), k_mid(:,:,:)
    complex(8), allocatable :: coef_pw_save(:,:,:)
    complex(8), allocatable :: k1_pw(:,:,:), k_mid_pw(:,:,:)

    n = size(dg_frag%coef, 1)
    if (n <= 0) return
    n_pw = 0
    if (dg_frag%use_plane_wave_basis .and. allocated(dg_frag%coef_pw)) n_pw = dg_frag%n_plane_waves

    allocate(coef_save(n, dg_frag%nstate_tot, dg_frag%nspin))
    allocate(k1      (n, dg_frag%nstate_tot, dg_frag%nspin))
    allocate(k_mid   (n, dg_frag%nstate_tot, dg_frag%nspin))
    if (n_pw > 0) then
      allocate(coef_pw_save(n_pw, dg_frag%nstate_tot, dg_frag%nspin))
      allocate(k1_pw      (n_pw, dg_frag%nstate_tot, dg_frag%nspin))
      allocate(k_mid_pw   (n_pw, dg_frag%nstate_tot, dg_frag%nspin))
    end if

    ! Vector potentials
    Ac_t   = rt%Ac_tot(:, itt)
    Ac_mid = 0.5d0 * (rt%Ac_tot(:, itt) + rt%Ac_tot(:, itt+1))

    ! Save current state
    coef_save(1:n, :, :) = dg_frag%coef(1:n, :, :)
    if (n_pw > 0) coef_pw_save(1:n_pw, :, :) = dg_frag%coef_pw(1:n_pw, :, :)

    ! Step 1: derivative at coef(t) with A(t)
    if (n_pw > 0) then
      call calculate_time_derivative(dg_frag, system, mg, stencil, ppg, Ac_t, itt, k1, k1_pw)
    else
      call calculate_time_derivative(dg_frag, system, mg, stencil, ppg, Ac_t, itt, k1)
    end if

    ! Step 2: advance to midpoint (predictor half-step)
    dg_frag%coef(1:n, :, :) = coef_save(1:n, :, :) + 0.5d0 * dt * k1(1:n, :, :)
    if (n_pw > 0) then
      dg_frag%coef_pw(1:n_pw, :, :) = coef_pw_save(1:n_pw, :, :) + 0.5d0 * dt * k1_pw(1:n_pw, :, :)
    end if

    ! Step 3: derivative at coef_mid with A(t+dt/2)
    if (n_pw > 0) then
      call calculate_time_derivative(dg_frag, system, mg, stencil, ppg, Ac_mid, itt, k_mid, k_mid_pw)
    else
      call calculate_time_derivative(dg_frag, system, mg, stencil, ppg, Ac_mid, itt, k_mid)
    end if

    ! Step 4: apply midpoint propagation from original coef(t)
    dg_frag%coef(1:n, :, :) = coef_save(1:n, :, :) + dt * k_mid(1:n, :, :)
    if (n_pw > 0) then
      dg_frag%coef_pw(1:n_pw, :, :) = coef_pw_save(1:n_pw, :, :) + dt * k_mid_pw(1:n_pw, :, :)
    end if

    deallocate(coef_save, k1, k_mid)
    if (allocated(coef_pw_save)) deallocate(coef_pw_save)
    if (allocated(k1_pw)) deallocate(k1_pw)
    if (allocated(k_mid_pw)) deallocate(k_mid_pw)

  end subroutine time_evolution_aetrs
