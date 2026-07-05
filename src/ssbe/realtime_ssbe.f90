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
    integer :: ib, jb, ik
    real(8) :: herm_norm(1), herm_norm_l(1), trace_re, trace_re_l
    integer :: nk
    real(8) :: bj_am(8,8)

    call comm_get_groupinfo(icomm, irank, nproc)

    if (.not. check_input_variables_sbe()) return

    ! Read ground state electronic system:
    nk = num_kgrid(1)*num_kgrid(2)*num_kgrid(3)
    call init_sbe_gs_info(gs, sysname, base_directory, &
        & nk, nstate, nelec, &
        & al_vec1, al_vec2, al_vec3, &
        & .false., icomm)        
    
    ! Initialization of SBE solver and density matrix:
    call init_sbe_bloch_solver(sbe, gs, nstate_sbe(1), icomm)
    sbe%flag_vnl_correction = (yn_vnl_correction == 'y')

    if (trim(gauge_sbe) == "length_gauge") then
        ! Prepare qnm
        call prepare_qnm(sbe, gs, icomm)
        call adams_moulton_coefs(bj_am)
    end if

    ! GW collision-term setup (Phase 2): load Gamma(n,k) and the cold reference.
    if (yn_sbe_gw_collision == 'y') then
        allocate(gs%gamma_gw(1:nstate_sbe(1), 1:nk))
        allocate(gs%f0_ref  (1:nstate_sbe(1), 1:nk))
        gs%gamma_gw(:, :) = 0d0
        gs%f0_ref(1:nstate_sbe(1), 1:nk) = gs%occup(1:nstate_sbe(1), 1:nk)
        call load_gw_rate(trim(file_sbe_gw_rate), nstate_sbe(1), nk, gs%gamma_gw)
        if (irank == 0) write(*,'(a,a)') "# SBE GW collision ON, mode = ", trim(sbe_deph_mode)
    end if

    ! Prepare external pulse
    allocate(Ac_ext_t(1:3, -1:nt+1))
    call calc_Ac_ext_t(0.0d0, dt, -1, nt+1, Ac_ext_t)
    ! Initial energy
    energy = 0.0d0
    E(:) = 0.0d0; Jmat(:) = 0.0d0;

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
        else ! trim(gauge_sbe) == "length_gauge")
            call dt_evolve_bloch_lg(sbe, gs, E(:), bj_am, dt, icomm)
            call calc_current_bloch_lg(sbe, gs, Jmat, icomm)
        end if
        energy = energy + dot_product(E(1:3), -Jmat(1:3)) * gs%volume * dt
        
        if (irank == 0) then
            call write_sbe_rt_line(fh_sbe_rt, &
                & t, Ac_ext_t(1:3, it), E(1:3), Ac_ext_t(1:3, it), E(1:3), Jmat(1:3))
        end if

        if (mod(it, 10) == 0) then
            tr_all = calc_trace(sbe, gs, nstate_sbe(1), icomm)
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
            tr_all = calc_trace(sbe, gs, nstate_sbe(1), icomm)
            tr_vb = calc_trace(sbe, gs, nelec / 2, icomm)    
            if (irank == 0) then
                call write_sbe_nex_line(fh_sbe_nex, t, tr_all - tr_vb, nelec - tr_vb)
            end if
        end if

        if (mod(it, 500) == 0) then
            if (irank == 0) then
                flush(fh_sbe_rt)
                flush(fh_sbe_rt_energy)
                flush(fh_sbe_nex)
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
        if (trim(gauge_sbe) == "length_gauge") close(fh_sbe_diag)
    end if

    return
end subroutine main_realtime_ssbe

end module realtime_ssbe
