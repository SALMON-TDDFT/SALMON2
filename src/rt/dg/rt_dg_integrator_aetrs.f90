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

    integer :: n
    real(8) :: Ac_t(3), Ac_mid(3)
    complex(8), allocatable :: coef_save(:,:,:)
    complex(8), allocatable :: k1(:,:,:), k_mid(:,:,:)

    n = dg_frag%n_mat_max
    if (n <= 0) return
    if (dg_frag%use_plane_wave_basis .or. allocated(dg_frag%coef_pw)) then
      stop "DG AETRS now supports the pure fragment block-sparse route only"
    end if
    if (allocated(dg_frag%local_coef_global_ids) .and. size(dg_frag%coef, 1) < n) then
      stop "DG AETRS is disabled for row-split coefficients; use rk4 or ssprk3"
    end if

    allocate(coef_save(n, dg_frag%nstate_tot, dg_frag%nspin))
    allocate(k1      (n, dg_frag%nstate_tot, dg_frag%nspin))
    allocate(k_mid   (n, dg_frag%nstate_tot, dg_frag%nspin))

    ! Vector potentials
    Ac_t   = rt%Ac_tot(:, itt)
    Ac_mid = 0.5d0 * (rt%Ac_tot(:, itt) + rt%Ac_tot(:, itt+1))

    ! Save current state
    coef_save(1:n, :, :) = dg_frag%coef(1:n, :, :)

    ! Step 1: derivative at coef(t) with A(t)
    call calculate_time_derivative(dg_frag, system, mg, ppg, Ac_t, k1)

    ! Step 2: advance to midpoint (predictor half-step)
    dg_frag%coef(1:n, :, :) = coef_save(1:n, :, :) + 0.5d0 * dt * k1(1:n, :, :)

    ! Step 3: derivative at coef_mid with A(t+dt/2)
    call calculate_time_derivative(dg_frag, system, mg, ppg, Ac_mid, k_mid)

    ! Step 4: apply midpoint propagation from original coef(t)
    dg_frag%coef(1:n, :, :) = coef_save(1:n, :, :) + dt * k_mid(1:n, :, :)

    deallocate(coef_save, k1, k_mid)

  end subroutine time_evolution_aetrs
