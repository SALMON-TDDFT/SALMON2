module realtime_ssbe
    implicit none
contains

subroutine main_realtime_ssbe(icomm)
    use salmon_global
    use communication
    use gs_info_ssbe
    use bloch_solver_ssbe
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

    call comm_get_groupinfo(icomm, irank, nproc)

    if (.not. check_input_variables_sbe()) return

    ! Read ground state electronic system:
    call init_sbe_gs_info(gs, sysname, base_directory, &
        & num_kgrid, nstate, nelec, &
        & al_vec1, al_vec2, al_vec3, &
        & .false., icomm)        
    
    ! Initialization of SBE solver and density matrix:
    call init_sbe_bloch_solver(sbe, gs, nstate_sbe, icomm)
    sbe%flag_vnl_correction = (yn_vnl_correction == 'y')

    ! Prepare external pulse
    allocate(Ac_ext_t(1:3, -1:nt+1))
    call calc_Ac_ext_t(0.0d0, dt, 0, nt+1, Ac_ext_t)
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
        ! Stdout logs:
        write(*, "(a)") " time-step time[fs] Current(xyz)[a.u.]                     electrons   Total energy[au]"
        write(*, "(a)") "---------------------------------------------------------------------------------------"
    end if

    call comm_sync_all(icomm)

    ! Realtime calculation
    do it = 1, nt
        t = dt * it
        call dt_evolve_bloch(sbe, gs, Ac_ext_t(:, it), dt)
        call calc_current_bloch(sbe, gs, Ac_ext_t(:, it), Jmat, icomm)
        E(:) = -(Ac_ext_t(:, it + 1) - Ac_ext_t(:, it - 1)) / (2 * dt)
        energy = energy + dot_product(E(1:3), -Jmat(1:3)) * gs%volume * dt
        
        if (irank == 0) then
            call write_sbe_rt_line(fh_sbe_rt, &
                & t, Ac_ext_t(1:3, it), E(1:3), Ac_ext_t(1:3, it), E(1:3), Jmat(1:3))
        end if

        if (mod(it, 10) == 0) then
            tr_all = calc_trace(sbe, gs, nstate_sbe, icomm)
            if (irank == 0) then
                call write_sbe_rt_energy_line(fh_sbe_rt_energy, t, energy, energy)
                write(*, "(i6,f12.3,3es12.3,2f12.3)") it, t, Jmat(1:3), tr_all, energy
            end if
        end if
        
        if (mod(it, out_projection_step) == 0) then
            tr_all = calc_trace(sbe, gs, nstate_sbe, icomm)
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
            end if
        end if
    end do
    
    call comm_sync_all(icomm)

    if (irank == 0) then
        close(fh_sbe_rt)
        close(fh_sbe_rt_energy)
        close(fh_sbe_nex)
    end if

    return
end subroutine main_realtime_ssbe

end module realtime_ssbe
