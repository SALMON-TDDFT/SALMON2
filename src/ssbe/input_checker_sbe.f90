!
!  Copyright 2020 SALMON developers
!
!  Licensed under the Apache License, Version 2.0 (the "License");
!  you may not use this file except in compliance with the License.
!  You may obtain a copy of the License at
!
!      http://www.apache.org/licenses/LICENSE-2.0
!
!  Unless required by applicable law or agreed to in writing, software
!  distributed under the License is distributed on an "AS IS" BASIS,
!  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!  See the License for the specific language governing permissions and
!  limitations under the License.
!

!=======================================================================


module input_checker_sbe
    use salmon_global
    use communication, only: comm_is_root
    use parallelization, only: nproc_id_global, nproc_size_global
    implicit none
contains

function check_input_variables_sbe() result(flag)
    implicit none
    logical :: flag
    integer :: i
    logical :: flag_spinor
    integer :: nvb_abs   ! occupied bands of the ABSOLUTE (unwindowed) manifold

    flag = .true.

    ! SO-SBE spinor mode: spin='noncollinear' (the global input checker
    ! already enforces yn_spinorbit='y' with it).  nstate/nstate_sbe/nelec are
    ! then SPINOR counts: 1 electron per occupied band, Kramers partners are
    ! separate band indices, so the occupied manifold is bands 1..nelec
    ! (spinless: 1..nelec/2).
    flag_spinor = (trim(spin) == 'noncollinear')
    if (flag_spinor) then
        nvb_abs = nelec
    else
        nvb_abs = nelec / 2
    end if

    if (flag_spinor) then
        if (trim(theory) == "maxwell_sbe") then
            call raise("ERROR! spin='noncollinear' (SO-SBE) is not supported for theory='maxwell_sbe' yet!")
        end if
        ! Kramers theorem (SOC with time-reversal symmetry) makes EVERY band
        ! doubly degenerate at EVERY k: the legacy length-gauge machinery
        ! (un-gated scalar -q/t_2 dephasing + |dnm| dipole phases) is
        ! basis-dependent / divergent inside those degenerate pairs.  Only the
        ! gauge-covariant path ('gicov': Wilson-line transport + the
        ! delta-omega-GATED T2, which correctly never dephases WITHIN a
        ! Kramers pair) is well-defined there.  The velocity gauge carries no
        ! scalar T2 and needs no gate.
        if (trim(gauge_sbe) == "length_gauge" .and. trim(sbe_lg_degen) /= "gicov") then
            call raise("ERROR! spin='noncollinear' (SO-SBE) with gauge_sbe='length_gauge' requires " // &
                & "sbe_lg_degen='gicov' (Kramers degeneracy needs the gauge-covariant transport " // &
                & "and the delta-omega-gated T2)!")
        end if
        if (comm_is_root(nproc_id_global)) then
            write(*, '(a)') "# SO-SBE (spin='noncollinear'): spinor bands, occupation 1 per band; " // &
                & "occupied manifold = bands 1..nelec."
        end if
    end if

    if(num_sbe>=2)then
        if (nstate_sbe(2) <= 0) call raise("ERROR! nstate_sbe must be specified when num_sbe>=2!")
        if (nk_sbe(2) <= 0) call raise("ERROR! nk_sbe must be specified when num_sbe>=2!")
        if (nelec_sbe(2) <= 0) call raise("ERROR! nelec_sbe must be specified when num_sbe>=2!")
        if ((norm2(al_vec1_sbe(:,2)) <= 1d-9) .and. &
            (norm2(al_vec2_sbe(:,2)) <= 1d-9) .and. &
            (norm2(al_vec3_sbe(:,2)) <= 1d-9) .and. &
            (norm2(al_sbe(:,2)) <= 1d-9)) then
            call raise("ERROR! 'al_sbe' or 'al_vec1..3_sbe' must be specified when num_sbe>=2!")
        end if
        if(nx_m*ny_m*nz_m < nproc_size_global) then
            call raise("ERROR! Number of macro points must be larger than or equal to nproc_size_global when num_sbe>=2!")
        end if
    end if

    if ((norm2(al_vec1) <= 1d-9) .and. (norm2(al_vec2) <= 1d-9) .and. (norm2(al_vec3) <= 1d-9)) then
        if (norm2(al) > 1d-9) then
            al_vec1(1:3) = (/ al(1), 0.0d0, 0.0d0 /)
            al_vec2(1:3) = (/ 0.0d0, al(2), 0.0d0 /)
            al_vec3(1:3) = (/ 0.0d0, 0.0d0, al(3) /)
        else
            if ((norm2(al_vec1_sbe(:,1)) <= 1d-9) .and. &
                (norm2(al_vec2_sbe(:,1)) <= 1d-9) .and. &
                (norm2(al_vec3_sbe(:,1)) <= 1d-9) .and. &
                (norm2(al_sbe(:,1)) <= 1d-9)) then
                call raise("ERROR! Either one of 'al', 'al_vec1..3', 'al_sbe' or 'al_vec1..3_sbe' must be specified!")
            end if
        end if
    end if

    if ((norm2(al_vec1_sbe(:,1)) <= 1d-9) .and. &
        (norm2(al_vec2_sbe(:,1)) <= 1d-9) .and. &
        (norm2(al_vec3_sbe(:,1)) <= 1d-9)) then
        if(norm2(al_sbe(:,1)) > 1d-9) then
            al_vec1_sbe(1:3,1) = (/ al_sbe(1,1), 0.0d0, 0.0d0 /)
            al_vec2_sbe(1:3,1) = (/ 0.0d0, al_sbe(2,1), 0.0d0 /)
            al_vec3_sbe(1:3,1) = (/ 0.0d0, 0.0d0, al_sbe(3,1) /)
        else
            al_vec1_sbe(1:3,1) = al_vec1(1:3)
            al_vec2_sbe(1:3,1) = al_vec2(1:3)
            al_vec3_sbe(1:3,1) = al_vec3(1:3)
        end if
    end if

    do i=2,num_sbe
      if ((norm2(al_vec1_sbe(:,i)) <= 1d-9) .and. &
          (norm2(al_vec2_sbe(:,i)) <= 1d-9) .and. &
          (norm2(al_vec3_sbe(:,i)) <= 1d-9)) then
          al_vec1_sbe(1:3,i) = (/ al_sbe(1,i), 0.0d0, 0.0d0 /)
          al_vec2_sbe(1:3,i) = (/ 0.0d0, al_sbe(2,i), 0.0d0 /)
          al_vec3_sbe(1:3,i) = (/ 0.0d0, 0.0d0, al_sbe(3,i) /)
      end if
    end do

    if (nstate_sbe(1) <= 0) nstate_sbe(1) = nstate
    if (nstate_sbe(1) > nstate) then
        call raise("ERROR! 'nstate_sbe' must be smaller than 'nstate'!")
    end if

    if (sysname_sbe(1) == "default") sysname_sbe(1) = sysname
    if (nk_sbe(1) <= 0) nk_sbe(1) = num_kgrid(1)*num_kgrid(2)*num_kgrid(3)
    if (nelec_sbe(1) <= 0) nelec_sbe(1) = nelec

    ! Kramers-pair alignment (spinor SO-SBE): with time-reversal symmetry the
    ! spinor export pairs adjacent bands (2m-1, 2m) into exactly degenerate
    ! Kramers doublets.  A half-occupied doublet (odd nelec) or a frozen cut
    ! through a doublet (even nband_sbe_min-1 violated) makes the basis
    ! choice INSIDE an exactly degenerate pair observable -- the same
    ! pathology the gicov guard exists for -- so both are rejected.  An odd
    ! window TOP merely truncates the partner of the highest conduction
    ! doublet (a basis-truncation artifact, not an occupied-manifold one):
    ! warn, do not fail.
    if (flag_spinor) then
        if (mod(nelec, 2) /= 0) then
            call raise("ERROR! spin='noncollinear' (SO-SBE) requires even 'nelec' " // &
                & "(a half-occupied Kramers doublet is basis-dependent)!")
        end if
        if (mod(nband_sbe_min - 1, 2) /= 0) then
            call raise("ERROR! spin='noncollinear' (SO-SBE): 'nband_sbe_min'-1 must be even " // &
                & "(the frozen cut must not split a Kramers doublet)!")
        end if
        if (mod(nstate_sbe(1), 2) /= 0) then
            if (comm_is_root(nproc_id_global)) then
                write(*, '(a)') "WARNING: odd 'nstate_sbe' splits the highest Kramers doublet " // &
                    & "across the window top (basis-truncation artifact risk)."
            end if
        end if
    end if

    ! Band-window lower-cut: SBE propagates the contiguous window
    ! [nband_sbe_min, nstate_sbe(1)]; bands 1..nband_sbe_min-1 are frozen as
    ! inert fully-occupied (spinless: 2 electrons each; spinor: 1 each), so
    ! every frozen band must lie inside the occupied manifold 1..nvb_abs
    ! (nvb_abs = nelec/2 spinless, = nelec spinor).
    if (trim(theory) /= "maxwell_sbe" .and. nstate_sbe(1) < nvb_abs) then
        call raise("ERROR! 'nstate_sbe' must be >= the occupied band count " // &
            & "(nelec/2 spinless, nelec for spin='noncollinear') " // &
            & "(the SBE band window must contain the occupied manifold)!")
    end if
    if (nband_sbe_min < 1) then
        call raise("ERROR! 'nband_sbe_min' must be >= 1!")
    end if
    if (nband_sbe_min > nstate_sbe(1)) then
        call raise("ERROR! 'nband_sbe_min' must not exceed 'nstate_sbe' (empty band window)!")
    end if
    if (nband_sbe_min - 1 > nvb_abs) then
        call raise("ERROR! 'nband_sbe_min'-1 must be <= the occupied band count " // &
            & "(nelec/2 spinless, nelec for spin='noncollinear') " // &
            & "(frozen bands must be fully occupied)!")
    end if
    if (nband_sbe_min > 1) then
        if (trim(theory) == "maxwell_sbe") then
            call raise("ERROR! 'nband_sbe_min' > 1 is not supported for theory='maxwell_sbe'!")
        end if
        if (comm_is_root(nproc_id_global)) then
            write(*, '(a,i0,a)') "WARNING: 'nband_sbe_min' freezes bands 1..", &
                & nband_sbe_min - 1, " as inert fully-occupied (no dynamics, zero current)."
            if (flag_spinor) then
                write(*, '(a,i0,a,i0)') "  This is a frozen-core approximation: valid only if the " // &
                    & "frozen bands lie far below the band gap. nelec_eff = ", &
                    & nelec - (nband_sbe_min - 1), " of nelec = ", nelec
            else
                write(*, '(a,i0,a,i0)') "  This is a frozen-core approximation: valid only if the " // &
                    & "frozen bands lie far below the band gap. nelec_eff = ", &
                    & nelec - 2 * (nband_sbe_min - 1), " of nelec = ", nelec
            end if
            if (nband_sbe_min - 1 == nvb_abs) then
                write(*, '(a)') "  WARNING: the band window contains NO occupied bands " // &
                    & "(the entire valence manifold is frozen)."
            end if
        end if
    end if

    if (trim(gauge_sbe) == "length_gauge") then
        if(t_2 <= 0.d0) call raise("ERROR! 't_2' must be positive.")
        if(am_s /= 4 .and. am_s /= 8) then
            call raise("ERROR! Supported only when 'am_s' is 4 or 8.")
        end if
    end if

    select case(trim(sbe_lg_degen))
    case("off")
    case("gi", "gifix", "gicov")
        if (sbe_lg_diag == 2 .or. sbe_lg_diag == 3) then
            call raise("ERROR! 'sbe_lg_degen'='"//trim(sbe_lg_degen)//"' is incompatible with 'sbe_lg_diag'=2 or 3.")
        end if
    case default
        call raise("ERROR! 'sbe_lg_degen' must be 'off', 'gi', 'gifix', or 'gicov'.")
    end select

    if (sbe_lg_degen_floor <= 0d0) call raise("ERROR! 'sbe_lg_degen_floor' must be positive.")

    select case(yn_sbe_gs_current_subtract)
    case('y', 'n')
    case default
        call raise("ERROR! 'yn_sbe_gs_current_subtract' must be 'y' or 'n'!")
    end select
    if (yn_sbe_gs_current_subtract == 'y' .and. trim(gauge_sbe) == "length_gauge") then
        ! WARN, not error: the f-sum-deficiency subtraction acts only on the
        ! velocity-gauge current readout (calc_current_bloch); the length-gauge
        ! current has no A(t)-proportional truncation artifact to subtract, so
        ! the flag is a silent no-op there.  Follow the nband_sbe_min style:
        ! tell the user rather than fail.
        if (comm_is_root(nproc_id_global)) then
            write(*, '(a)') "WARNING: 'yn_sbe_gs_current_subtract'='y' has no effect for " // &
                & "gauge_sbe='length_gauge' (velocity-gauge current readout only)."
        end if
    end if
    ! --- VG completion (yn_sbe_vnl_exact): all-order nonlocal V_nl(k+A) ---
    select case(yn_sbe_vnl_exact)
    case('y', 'n')
    case default
        call raise("ERROR! 'yn_sbe_vnl_exact' must be 'y' or 'n'!")
    end select
    if (yn_sbe_vnl_exact == 'y') then
        if (trim(gauge_sbe) /= "velocity_gauge") then
            call raise("ERROR! 'yn_sbe_vnl_exact'='y' requires gauge_sbe='velocity_gauge' " // &
                & "(the kappa-stencil completion is a velocity-gauge construction)!")
        end if
        if (yn_vnl_correction == 'y') then
            ! HARD error, not a priority rule: the first-order A.rvnl term and
            ! the all-order DeltaV(k,A) are the SAME physics -- combining them
            ! double-counts the nonlocal coupling (a physics bug, not a
            ! preference).  The exact mode supersedes yn_vnl_correction.
            call raise("ERROR! 'yn_sbe_vnl_exact'='y' is mutually exclusive with " // &
                & "'yn_vnl_correction'='y' (the all-order DeltaV supersedes the first-order " // &
                & "A.rvnl term; combining them double-counts the nonlocal coupling)!")
        end if
        if (norder_correction /= 0) then
            call raise("ERROR! 'yn_sbe_vnl_exact'='y' requires 'norder_correction'=0 " // &
                & "(the norder readout corrections are a perturbative variant of the same " // &
                & "H(k+A) completion)!")
        end if
        if (trim(theory) == "maxwell_sbe") then
            call raise("ERROR! 'yn_sbe_vnl_exact'='y' is not supported for theory='maxwell_sbe' " // &
                & "(macroscopic A(t) cannot be pre-validated against the 1D stencil)!")
        end if
        if (flag_spinor) then
            call raise("ERROR! 'yn_sbe_vnl_exact'='y' does not support spin='noncollinear' yet " // &
                & "(scalar first stage; the file contract is spinor-ready)!")
        end if
        if (len_trim(file_sbe_vnl_kappa) == 0) then
            call raise("ERROR! 'file_sbe_vnl_kappa' must be specified when 'yn_sbe_vnl_exact'='y'!")
        end if
    end if

    if (yn_sbe_gs_current_subtract == 'y' .and. norder_correction >= 1) then
        ! WARN, not error: the deficiency tensor D is derived against the
        ! norder_correction=0 readout (the production default).  At norder>=1
        ! the core's own A-linear velocity-operator correction supplies (at
        ! rho ~ rho0) a bare-p, 1d-3-floored variant of the captured-window
        ! paramagnetic term that D is built against, so combining both changes
        ! the A-linear readout beyond the derived correction.  Tell the user;
        ! do not fail.
        if (comm_is_root(nproc_id_global)) then
            write(*, '(a)') "WARNING: 'yn_sbe_gs_current_subtract'='y' is derived for " // &
                & "'norder_correction'=0; with norder_correction>=1 the readout's own " // &
                & "A-linear correction overlaps the subtracted D*A(t) term."
        end if
    end if

    select case(yn_sbe_out_occ)
    case('y', 'n')
    case default
        call raise("ERROR! 'yn_sbe_out_occ' must be 'y' or 'n'!")
    end select

    if (trim(theory) /= "maxwell_sbe") return

    if (nx_m < 1) &
        call raise("ERROR! 'nx_m' must be larger than 1!")
    if (ny_m < 1) &
        call raise("ERROR! 'ny_m' must be larger than 1!")
    if (nz_m < 1) &
        call raise("ERROR! 'nz_m' must be larger than 1!")

    if ((nxvac_m(1) > 0) .and. (nxvacl_m > 0)) &
        call raise("ERROR! 'nxvac_m' and 'nxvacl_m' must not be specified simultaneously!")
    if ((nxvac_m(2) > 0) .and. (nxvacr_m > 0)) &
        call raise("ERROR! 'nxvac_m' and 'nxvacr_m' must not be specified simultaneously!")

    if ((hx_m > 0) .and. (dl_em(1) > 0)) &
        call raise("ERROR! 'hx_m' and 'dl_em' must not be specified simultaneously!")
    if ((hy_m > 0) .and. (dl_em(2) > 0)) &
        call raise("ERROR! 'hy_m' and 'dl_em' must not be specified simultaneously!")
    if ((hz_m > 0) .and. (dl_em(3) > 0)) &
        call raise("ERROR! 'hz_m' and 'dl_em' must not be specified simultaneously!")
    if ((hx_m <= 0) .and. (dl_em(1) <= 0)) &
        call raise("ERROR! 'hx_m' or 'dl_em(1)' must be specified!")
    if ((hy_m <= 0) .and. (dl_em(2) <= 0)) &
        call raise("ERROR! 'hy_m' or 'dl_em(2)' must be specified!")
    if ((hz_m <= 0) .and. (dl_em(3) <= 0)) &
        call raise("ERROR! 'hz_m' or 'dl_em(1)' must be specified!")
        
    select case(trim(fdtddim))
    case("1d", "1D")
    case("3d", "3D")
    case default
        call raise("ERROR! 'fdtddim' must be specified ('1d' or '3d')")
    end select

    if (abs(nxvacl_m) >= 1) nxvac_m(1) = abs(nxvacl_m)
    if (abs(nxvacr_m) >= 1) nxvac_m(2) = abs(nxvacr_m)

    return

    contains

    subroutine raise(msg)
        implicit none
        character(*), intent(in) :: msg
        if (comm_is_root(nproc_id_global)) then
            write(*, '(a)') msg
        end if
        flag = .false.
    end subroutine raise

end function check_input_variables_sbe
end module input_checker_sbe


