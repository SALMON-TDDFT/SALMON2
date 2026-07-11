module realtime_ssbe
    implicit none
contains

subroutine main_realtime_ssbe(icomm)
    use salmon_global
    use communication
    use gs_info_ssbe
    use bloch_solver_ssbe
    use sbe_collision_gw, only: load_gw_rate
    use em_field
    use datafile_ssbe
    use input_checker_sbe
    use filesystem, only: get_filehandle
    use sbe_lg_mode_ssbe, only: uses_integral_gicov
    use gicov_integral_ssbe, only: gicov_int_axis_single, gicov_int_jmax
    implicit none
    integer, intent(in) :: icomm

    type(s_sbe_bloch_solver) :: sbe
    type(s_sbe_gs_info) :: gs
    real(8) :: t,  E(3), jmat(3)
    real(8), allocatable :: Ac_ext_t(:, :)
    integer :: it, i
    real(8) :: energy, tr_all, tr_vb
    integer :: nproc, irank, ierr
    integer :: fh_sbe_rt, fh_sbe_rt_energy, fh_sbe_nex
    integer :: fh_sbe_diag
    integer :: fh_sbe_occ, fh_sbe_edist
    integer :: ib, jb, ik
    real(8) :: herm_norm(1), herm_norm_l(1), trace_re, trace_re_l
    integer :: nk
    integer :: nb_sbe_eff, nelec_eff
    real(8), allocatable :: gamma_abs(:, :)
    real(8) :: bj_am(8,8)
    ! Band-resolved population + energy-resolved distribution (yn_sbe_out_occ='y')
    real(8), allocatable :: nex_b(:), edist(:)
    real(8) :: ed_lo, ed_de, ed_sigma
    integer :: ed_nbin
    ! integral (covariant-Houston) transport mode (sbe_lg_degen='gicov_int')
    logical :: gi_mode, gi_ok_l
    integer :: gi_axis_l, gi_jmax_l, jt
    real(8) :: q_mid, q_now, a_sh(3), qmx, tol_ax, twopi_l
    real(8), allocatable :: q_all(:, :)

    call comm_get_groupinfo(icomm, irank, nproc)

    if (.not. check_input_variables_sbe()) return

    ! Band-window lower-cut (nband_sbe_min, default 1 = full window):
    ! the SBE propagates the contiguous window [nband_sbe_min, nstate_sbe(1)]
    ! (window index w <-> absolute band w + nband_sbe_min - 1); the frozen
    ! bands 1..nband_sbe_min-1 are inert fully-occupied and enter the
    ! trace / n_ex bookkeeping only through nelec_eff.
    nb_sbe_eff = nstate_sbe(1) - (nband_sbe_min - 1)

    ! Read ground state electronic system: the exports carry nstate bands
    ! (all consumed by the readers); gs%* stores the window
    ! [nband_sbe_min : nstate_sbe(1)], so gs%nb == nb_sbe_eff == sbe%nb.
    nk = num_kgrid(1)*num_kgrid(2)*num_kgrid(3)
    call init_sbe_gs_info(gs, sysname, base_directory, &
        & nk, nstate, nband_sbe_min, nstate_sbe(1), nelec, &
        & al_vec1, al_vec2, al_vec3, &
        & .false., icomm)

    ! Electrons inside the window, per the occupation convention resolved by
    ! init_sbe_gs_info: spinless nelec - 2*(nband_sbe_min-1); spinor
    ! (noncollinear/SOC, 1 electron per band) nelec - (nband_sbe_min-1).
    nelec_eff = gs%ne

    ! Initialization of SBE solver and density matrix:
    call init_sbe_bloch_solver(sbe, gs, nb_sbe_eff, icomm)
    sbe%flag_vnl_correction = (yn_vnl_correction == 'y')

    ! One-time diagnostic (flag-on only: the flag-off stdout is byte-unchanged):
    ! the window f-sum deficiency tensor D applied as J(t) -= D*A(t).
    if (yn_sbe_gs_current_subtract == 'y' .and. irank == 0) then
        write(*, '(a)') "# yn_sbe_gs_current_subtract: window f-sum deficiency tensor D [a.u.]"
        do i = 1, 3
            write(*, '(a,3es24.15e3)') "#   D(i,1:3) = ", sbe%fsum_D(i, 1:3)
        end do
        write(*, '(a,es14.6,a)') "#   (frozen-window limit Ne/V = ", &
            & dble(nelec_eff) / gs%volume, ")"
    end if

    if (trim(gauge_sbe) == "length_gauge") then
        ! Prepare qnm
        call prepare_qnm(sbe, gs, icomm)
        call adams_moulton_coefs(bj_am)
    end if

    ! GW collision-term setup (Phase 2): load Gamma(n,k) and the cold reference.
    ! The rate file is indexed by ABSOLUTE band, so load the full band range
    ! 1..nstate_sbe(1) and keep the window slice [nband_sbe_min : nstate_sbe(1)]
    ! (reader-side windowing, same as the gs exports).
    if (yn_sbe_gw_collision == 'y') then
        allocate(gs%gamma_gw(1:nb_sbe_eff, 1:nk))
        allocate(gs%f0_ref  (1:nb_sbe_eff, 1:nk))
        gs%gamma_gw(:, :) = 0d0
        gs%f0_ref(1:nb_sbe_eff, 1:nk) = gs%occup(1:nb_sbe_eff, 1:nk)
        allocate(gamma_abs(1:nstate_sbe(1), 1:nk))
        call load_gw_rate(trim(file_sbe_gw_rate), nstate_sbe(1), nk, gamma_abs)
        gs%gamma_gw(1:nb_sbe_eff, 1:nk) = gamma_abs(nband_sbe_min:nstate_sbe(1), 1:nk)
        deallocate(gamma_abs)
        if (irank == 0) write(*,'(a,a)') "# SBE GW collision ON, mode = ", trim(sbe_deph_mode)
    end if

    ! Prepare external pulse
    allocate(Ac_ext_t(1:3, -1:nt+1))
    call calc_Ac_ext_t(0.0d0, dt, -1, nt+1, Ac_ext_t)
    ! VG completion: validate the WHOLE trajectory against the kappa-stencil
    ! up front (axis collinearity + range incl. the 4-point support), so an
    ! under-sized stencil fails here and not mid-run.
    if (yn_sbe_vnl_exact == 'y') then
        call sbe_vnl_validate_trajectory(gs, Ac_ext_t, -1, nt+1, irank)
    end if

    ! Integral (covariant-Houston) transport setup (sbe_lg_degen='gicov_int'):
    ! analyse the WHOLE trajectory up front.  The reduced-mesh displacement is
    !   q_i(t) = num_kgrid(i) * (a(t) . a_i) / 2pi,   a(t) = -[A_ext(t)-A_ext(0)],
    ! (SALMON field convention E = -dA/dt, so da/dt = E; the mesh shift a(t) is
    ! the negative shifted external vector potential -- getting this sign wrong
    ! is a pi phase flip).  Require exactly one reduced axis to move over the
    ! entire pulse (runtime guard, not a deck-vector test), fix j_max from its
    ! peak span, and build the transport cache once.
    gi_mode = uses_integral_gicov(sbe_lg_degen)
    gi_axis_l = 0;  gi_jmax_l = 0;  twopi_l = 2d0 * acos(-1d0)
    if (gi_mode) then
        allocate(q_all(3, 0:nt))
        do jt = 0, nt
            a_sh(1:3) = -(Ac_ext_t(1:3, jt) - Ac_ext_t(1:3, 0))
            do i = 1, 3
                q_all(i, jt) = dble(num_kgrid(i)) &
                    & * dot_product(a_sh(1:3), gs%a_matrix(1:3, i)) / twopi_l
            end do
        end do
        qmx = maxval(abs(q_all))
        tol_ax = max(1d-12, 1d-8 * qmx)
        call gicov_int_axis_single(q_all, 3, nt + 1, tol_ax, gi_ok_l, gi_axis_l)
        if (.not. gi_ok_l) then
            if (irank == 0) write(*, '(a)') "ERROR(realtime_ssbe): sbe_lg_degen=" // &
                & "'gicov_int' requires single-axis linear polarization on a reciprocal-" // &
                & "mesh axis; the field trajectory drives 0 or >=2 reduced axes " // &
                & "([110]/elliptic/circular are v1-unsupported)."
            error stop 1
        end if
        gi_jmax_l = gicov_int_jmax(maxval(abs(q_all(gi_axis_l, :))), 1d0)
        call build_gicov_integral_cache(sbe, gs, gi_axis_l, gi_jmax_l, icomm)
        if (irank == 0) write(*, '(a,i0,a,i0)') &
            & "# gicov_int: driven reduced axis = ", gi_axis_l, &
            & ", transport span j_max = ", gi_jmax_l
        deallocate(q_all)
    end if

    ! Initial energy
    energy = 0.0d0
    E(:) = 0.0d0; Jmat(:) = 0.0d0;

    ! Population-output buffers + energy grid (yn_sbe_out_occ='y').  Allocated on
    ! EVERY rank: calc_band_population / calc_energy_distribution are collective
    ! (comm_summation), so every rank fills a local partial and reduces.  The
    ! grid comes from gs%eigen (replicated) => identical on all ranks.
    if (yn_sbe_out_occ == 'y') then
        allocate(nex_b(1:nb_sbe_eff))
        ed_nbin = 300
        call sbe_edist_grid(gs, nb_sbe_eff, ed_nbin, ed_lo, ed_de, ed_sigma)
        allocate(edist(1:ed_nbin))
    end if

    if (irank == 0) then
        ! SYSNAME_sbe_rt.data
        fh_sbe_rt = get_filehandle()
        open(unit=fh_sbe_rt, file=trim(base_directory)//trim(sysname)//"_sbe_rt.data", action="write")
        call write_sbe_rt_header(fh_sbe_rt)
        ! SYSNAME_sbe_rt_energy.data
        fh_sbe_rt_energy = get_filehandle()
        open(unit=fh_sbe_rt_energy, file=trim(base_directory)//trim(sysname)//"_sbe_rt_energy.data", action="write")
        call write_sbe_rt_energy_header(fh_sbe_rt_energy)
        ! SYSNAME_sbe_nex.data
        fh_sbe_nex = get_filehandle()
        open(unit=fh_sbe_nex, file=trim(base_directory)//trim(sysname)//"_sbe_nex.data", action="write")
        call write_sbe_nex_header(fh_sbe_nex)
        if (gs%is_metal) then
            ! Diagnostic-only caveat (does not affect propagation/current):
            ! n_ex/n_hole below are read out against the rigid valence split
            ! 1..gs%nvb, which is not the real per-k occupied-band set for a
            ! metal-like ground state -- expect a nonzero, non-physical
            ! n_ex/n_hole already at t=0 (before any pulse). Occupation-
            ! weighted generalization of this diagnostic is a separate
            ! follow-up; calc_band_population (SYSNAME_sbe_occ.data) already
            ! reads out against the real gs%occup and is unaffected.
            write(*, '(a)') "# WARNING: metal-like occupation detected -- " // &
                & "_sbe_nex.data n_ex/n_hole use the insulator-only 1..nvb " // &
                & "valence split and are not physically meaningful here " // &
                & "(current/dynamics are unaffected; see calc_band_population)."
        end if
        ! SYSNAME_sbe_occ.data / _sbe_edist.data (band- and energy-resolved population)
        if (yn_sbe_out_occ == 'y') then
            fh_sbe_occ = get_filehandle()
            open(unit=fh_sbe_occ, file=trim(base_directory)//trim(sysname)//"_sbe_occ.data", action="write")
            call write_sbe_occ_header(fh_sbe_occ, nb_sbe_eff)
            fh_sbe_edist = get_filehandle()
            open(unit=fh_sbe_edist, file=trim(base_directory)//trim(sysname)//"_sbe_edist.data", action="write")
            call write_sbe_edist_header(fh_sbe_edist)
        end if
        ! SYSNAME_sbe_diag.data (LG-SBE mechanism diagnostics; length_gauge only)
        if (trim(gauge_sbe) == "length_gauge") then
            fh_sbe_diag = get_filehandle()
            open(unit=fh_sbe_diag, file=trim(base_directory)//trim(sysname)//"_sbe_diag.data", action="write")
            write(fh_sbe_diag, '(a)') "# LG-SBE mechanism diagnostics"
            write(fh_sbe_diag, '(a)') "# 1:time[a.u.]  2:herm_norm  3:trace_re"
        end if
        ! Stdout logs:
        write(*, "(a)") "  time-step  time[fs] Current(xyz)[a.u.]                     electrons   Total energy[au]"
        write(*, "(a)") "-----------------------------------------------------------------------------------------"
    end if

    call comm_sync_all(icomm)

    ! Realtime calculation
    do it = 1, nt
        t = dt * it
        E(:) = -(Ac_ext_t(:, it + 1) - Ac_ext_t(:, it - 1)) / (2 * dt)
        if (trim(gauge_sbe) == "velocity_gauge") then
            call dt_evolve_bloch(sbe, gs, Ac_ext_t(:, it), dt)
            call calc_current_bloch(sbe, gs, Ac_ext_t(:, it), Jmat, icomm)
        else if (gi_mode) then ! length_gauge, integral covariant-Houston transport
            ! reduced-mesh shift at the STEP MIDPOINT t_it - dt/2 (propagation)
            a_sh(1:3) = -(0.5d0 * (Ac_ext_t(1:3, it - 1) + Ac_ext_t(1:3, it)) &
                &         - Ac_ext_t(1:3, 0))
            q_mid = dble(num_kgrid(gi_axis_l)) &
                & * dot_product(a_sh(1:3), gs%a_matrix(1:3, gi_axis_l)) / twopi_l
            call dt_evolve_bloch_lg_integral(sbe, gs, q_mid, dt)
            ! reduced-mesh shift at t_it (current readout)
            a_sh(1:3) = -(Ac_ext_t(1:3, it) - Ac_ext_t(1:3, 0))
            q_now = dble(num_kgrid(gi_axis_l)) &
                & * dot_product(a_sh(1:3), gs%a_matrix(1:3, gi_axis_l)) / twopi_l
            call calc_current_bloch_lg_integral(sbe, gs, q_now, Jmat, icomm)
        else ! trim(gauge_sbe) == "length_gauge") (FD / gi / gifix / gicov)
            call dt_evolve_bloch_lg(sbe, gs, E(:), bj_am, dt, icomm)
            call calc_current_bloch_lg(sbe, gs, Jmat, icomm)
        end if
        energy = energy + dot_product(E(1:3), -Jmat(1:3)) * gs%volume * dt
        
        if (irank == 0) then
            call write_sbe_rt_line(fh_sbe_rt, &
                & t, Ac_ext_t(1:3, it), E(1:3), Ac_ext_t(1:3, it), E(1:3), Jmat(1:3))
        end if

        if (mod(it, 10) == 0) then
            tr_all = calc_trace(sbe, gs, nb_sbe_eff, icomm)
            if (irank == 0) then
                call write_sbe_rt_energy_line(fh_sbe_rt_energy, t, energy, energy)
                write(*, "(i8,f12.3,3es12.3,2f12.3)") it, t, Jmat(1:3), tr_all, energy
            end if
            ! LG-SBE mechanism diagnostics: Hermiticity norm & trace of qnm
            if (trim(gauge_sbe) == "length_gauge") then
                herm_norm_l(1) = 0.0d0
                trace_re_l = 0.0d0
                do ik = sbe%ik_min, sbe%ik_max
                    do ib = 1, sbe%nb
                        trace_re_l = trace_re_l + gs%kweight(ik) * real(sbe%qnm(ib, ib, ik))
                        do jb = 1, sbe%nb
                            herm_norm_l(1) = max(herm_norm_l(1), &
                                & abs(sbe%qnm(ib, jb, ik) - conjg(sbe%qnm(jb, ib, ik))))
                        end do
                    end do
                end do
                call comm_get_max(herm_norm_l, herm_norm, 1, icomm)
                call comm_summation(trace_re_l, trace_re, icomm)
                if (irank == 0) then
                    write(fh_sbe_diag, '(3es24.15e3)') t, herm_norm(1), trace_re
                end if
            end if
        end if
        
        if (mod(it, out_projection_step) == 0) then
            ! window bookkeeping: the occupied bands inside the window are
            ! 1..gs%nvb (window indices; nvb = nelec_eff/2 spinless, =
            ! nelec_eff spinor); the frozen bands are inert, so n_ex and
            ! n_hole are exact against nelec_eff (not nelec).
            tr_all = calc_trace(sbe, gs, nb_sbe_eff, icomm)
            tr_vb = calc_trace(sbe, gs, gs%nvb, icomm)
            if (irank == 0) then
                call write_sbe_nex_line(fh_sbe_nex, t, tr_all - tr_vb, nelec_eff - tr_vb)
            end if
            ! Band-resolved population + energy-resolved distribution.  Both
            ! calc_* are collective (comm_summation) => every rank calls them;
            ! only rank 0 writes.
            if (yn_sbe_out_occ == 'y') then
                if (gi_mode) then
                    ! gicov_int: occupations by projection onto the INSTANTANEOUS
                    ! eigenbasis at x = kappa - a(t) (diag(rho~) is basis-dependent
                    ! under transport).  The static-eps energy distribution below is
                    ! left as an approximate diagnostic in the co-moving frame
                    ! (moving-frame energy binning is a follow-up).
                    call calc_band_population_integral(sbe, gs, q_now, nex_b, icomm)
                else
                    call calc_band_population(sbe, gs, nex_b, icomm)
                end if
                if (irank == 0) call write_sbe_occ_line(fh_sbe_occ, t, nb_sbe_eff, nex_b)
                call calc_energy_distribution(sbe, gs, ed_lo, ed_de, ed_nbin, ed_sigma, edist, icomm)
                if (irank == 0) call write_sbe_edist_block(fh_sbe_edist, t, ed_nbin, ed_lo, ed_de, edist)
            end if
        end if

        if (mod(it, 500) == 0) then
            if (irank == 0) then
                flush(fh_sbe_rt)
                flush(fh_sbe_rt_energy)
                flush(fh_sbe_nex)
                if (yn_sbe_out_occ == 'y') then
                    flush(fh_sbe_occ)
                    flush(fh_sbe_edist)
                end if
                if (trim(gauge_sbe) == "length_gauge") flush(fh_sbe_diag)
            end if
        end if
    end do
    
    call comm_sync_all(icomm)

    if (trim(gauge_sbe) == "length_gauge") call sbe_print_timers(icomm)

    if (irank == 0) then
        close(fh_sbe_rt)
        close(fh_sbe_rt_energy)
        close(fh_sbe_nex)
        if (yn_sbe_out_occ == 'y') then
            close(fh_sbe_occ)
            close(fh_sbe_edist)
        end if
        if (trim(gauge_sbe) == "length_gauge") close(fh_sbe_diag)
    end if

    if (yn_sbe_out_occ == 'y') then
        deallocate(nex_b, edist)
    end if

    return
end subroutine main_realtime_ssbe

end module realtime_ssbe
