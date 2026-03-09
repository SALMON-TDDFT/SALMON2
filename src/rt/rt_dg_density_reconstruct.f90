  subroutine calculate_density_from_fragments(dg_frag, system, mg, rho, rho_s)
    use structures
    use salmon_global, only: nelec, nelec_spin
    use communication, only: comm_summation
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system),     intent(in)    :: system
    type(s_rgrid),          intent(in)    :: mg
    type(s_scalar),         intent(inout) :: rho
    type(s_scalar),         intent(inout) :: rho_s(system%nspin)

    integer :: ifrag, io, i_local, ispin
    integer :: istate_frag
    integer :: ix, iy, iz, ixg, iyg, izg
    integer :: ig_i, nbf, ipw, n_pw
    integer :: nxyz(3), ixyz0(3)
    integer :: nocc_per_spin, nocc_spin
    real(8) :: occ_factor
    complex(8) :: coef_i, psi_val, phase_pw, ci
    real(8) :: phi_i, rho_contrib
    real(8) :: total_charge, scale_rho
    real(8) :: rx, ry, rz, boxL(3), inv_sqrt_vol, theta
    real(8), allocatable :: rho_local(:,:,:), w_local(:,:,:), w_global(:,:,:)

    rho%f = 0.0d0
    do ispin = 1, system%nspin
      rho_s(ispin)%f = 0.0d0
    end do

    if (.not. allocated(dg_frag%phi_frag)) return

    allocate(rho_local(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))
    allocate(w_local(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))
    allocate(w_global(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))
    rho_local = 0.0d0
    w_local = 0.0d0
    w_global = 0.0d0

    ! Closed-shell fallback: nelec = 2 * nocc_per_spin.
    nocc_per_spin = min(dg_frag%nstate_tot, int(nelec / 2.0d0 + 1.0d-12))
    occ_factor = 2.0d0 / real(system%nspin, 8)
    n_pw = 0
    if (dg_frag%use_plane_wave_basis .and. allocated(dg_frag%coef_pw) .and. allocated(dg_frag%k_pw)) then
      n_pw = dg_frag%n_plane_waves
    end if
    boxL(1) = dg_frag%hgs(1) * real(mg%num(1), 8)
    boxL(2) = dg_frag%hgs(2) * real(mg%num(2), 8)
    boxL(3) = dg_frag%hgs(3) * real(mg%num(3), 8)
    inv_sqrt_vol = 1.0d0 / sqrt(max(1.0d-16, boxL(1) * boxL(2) * boxL(3)))

    i_local = 0
    do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
      i_local = i_local + 1

      nxyz(1:3) = dg_frag%nxyz_domain(1:3, ifrag)
      ixyz0(1:3) = dg_frag%ixyz_frag(1:3, ifrag)

      ! Count overlap multiplicity on the global grid for later averaging.
      !$omp parallel do collapse(3) private(ix, iy, iz, ixg, iyg, izg) schedule(static)
      do iz = 1, nxyz(3)
        do iy = 1, nxyz(2)
          do ix = 1, nxyz(1)
            ixg = ixyz0(1) + ix - 1
            iyg = ixyz0(2) + iy - 1
            izg = ixyz0(3) + iz - 1

            ixg = mod(ixg - 1, mg%num(1)) + 1
            iyg = mod(iyg - 1, mg%num(2)) + 1
            izg = mod(izg - 1, mg%num(3)) + 1

            if (ixg < mg%is(1) .or. ixg > mg%ie(1)) cycle
            if (iyg < mg%is(2) .or. iyg > mg%ie(2)) cycle
            if (izg < mg%is(3) .or. izg > mg%ie(3)) cycle

            !$omp atomic update
            w_local(ixg, iyg, izg) = w_local(ixg, iyg, izg) + 1.0d0
          end do
        end do
      end do
      !$omp end parallel do

      do ispin = 1, system%nspin
        nocc_spin = nocc_per_spin
        if (system%nspin == 2 .and. sum(nelec_spin(:)) > 0) then
          nocc_spin = min(dg_frag%nstate_tot, nelec_spin(ispin))
        end if
        nbf = dg_frag%n_basis(ifrag, ispin)
        do io = 1, nocc_spin
          !$omp parallel do collapse(3) private(ix, iy, iz, ixg, iyg, izg, istate_frag, ig_i, coef_i, phi_i, psi_val, rho_contrib) schedule(static)
          do iz = 1, nxyz(3)
            do iy = 1, nxyz(2)
              do ix = 1, nxyz(1)
                ixg = ixyz0(1) + ix - 1
                iyg = ixyz0(2) + iy - 1
                izg = ixyz0(3) + iz - 1

                ixg = mod(ixg - 1, mg%num(1)) + 1
                iyg = mod(iyg - 1, mg%num(2)) + 1
                izg = mod(izg - 1, mg%num(3)) + 1

                if (ixg < mg%is(1) .or. ixg > mg%ie(1)) cycle
                if (iyg < mg%is(2) .or. iyg > mg%ie(2)) cycle
                if (izg < mg%is(3) .or. izg > mg%ie(3)) cycle

                psi_val = (0.0d0, 0.0d0)
                do istate_frag = 1, nbf
                  ig_i = dg_frag%index_basis(istate_frag, ifrag, ispin)
                  if (ig_i < 1 .or. ig_i > dg_frag%n_mat_max) cycle
                  coef_i = dg_frag%coef(ig_i, io, ispin)
                  phi_i = dg_frag%phi_frag(ix, iy, iz, istate_frag, i_local)
                  psi_val = psi_val + coef_i * phi_i
                end do
                if (n_pw > 0) then
                  rx = real(ixg - 1, 8) * dg_frag%hgs(1)
                  ry = real(iyg - 1, 8) * dg_frag%hgs(2)
                  rz = real(izg - 1, 8) * dg_frag%hgs(3)
                  do ipw = 1, n_pw
                    theta = dg_frag%k_pw(1, ipw) * rx + dg_frag%k_pw(2, ipw) * ry + dg_frag%k_pw(3, ipw) * rz
                    phase_pw = cmplx(cos(theta), sin(theta), kind=8)
                    ci = dg_frag%coef_pw(ipw, io, ispin)
                    psi_val = psi_val + ci * phase_pw * inv_sqrt_vol
                  end do
                end if

                rho_contrib = occ_factor * real(conjg(psi_val) * psi_val, kind=8)

                !$omp atomic update
                rho_local(ixg, iyg, izg) = rho_local(ixg, iyg, izg) + rho_contrib
              end do
            end do
          end do
          !$omp end parallel do
        end do
      end do
    end do

    call comm_summation(rho_local, rho%f, &
                       (mg%ie(1)-mg%is(1)+1)*(mg%ie(2)-mg%is(2)+1)*(mg%ie(3)-mg%is(3)+1), &
                       dg_frag%icomm)

    call comm_summation(w_local, w_global, &
                       (mg%ie(1)-mg%is(1)+1)*(mg%ie(2)-mg%is(2)+1)*(mg%ie(3)-mg%is(3)+1), &
                       dg_frag%icomm)

    where (w_global > 0.5d0)
      rho%f = rho%f / w_global
    end where

    total_charge = sum(rho%f) * system%hvol
    if (total_charge > 1.0d-14 .and. total_charge == total_charge) then
      scale_rho = nelec / total_charge
      rho%f = rho%f * scale_rho
    end if

    do ispin = 1, system%nspin
      rho_s(ispin)%f = rho%f / real(system%nspin, 8)
    end do

    deallocate(rho_local, w_local, w_global)

  end subroutine calculate_density_from_fragments
