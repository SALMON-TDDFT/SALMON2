module realtime_ssbe
    implicit none
contains
subroutine realtime_main_ssbe(icomm)
    use salmon_global
    use communication
    use gs_info_ssbe
    use bloch_solver_ssbe
    use em_field
    implicit none
    integer, intent(in) :: icomm

    type(s_sbe_bloch_solver) :: sbe
    type(s_sbe_gs_info) :: gs
    real(8) :: t,  E(3), jmat(3)
    real(8), allocatable :: Ac_ext_t(:, :)
    integer :: it, i
    real(8) :: energy0, energy
    real(8) :: tr_all, tr_vb
    integer :: nproc, irank, ierr
    integer :: nthread
    integer :: nstate_sbe

    call comm_get_groupinfo(icomm, irank, nproc)

    if (0.0d0 < al(1)) al_vec1(1:3) = (/ al(1), 0.0d0, 0.0d0 /)
    if (0.0d0 < al(2)) al_vec2(1:3) = (/ 0.0d0, al(2), 0.0d0 /)
    if (0.0d0 < al(3)) al_vec3(1:3) = (/ 0.0d0, 0.0d0, al(3) /)
    ! if (nstate_sbe < 1) nstate_sbe = nstate
    nstate_sbe = nstate

    ! Read ground state electronic system:
    call init_sbe_gs_info(gs, sysname, base_directory, &
        & num_kgrid, nstate, nelec, &
        & al_vec1, al_vec2, al_vec3, &
        & .false., icomm)        
    
    ! Initialization of SBE solver and density matrix:
    call init_sbe_bloch_solver(sbe, gs, nstate_sbe, icomm)

    ! Prepare external pulse
    allocate(Ac_ext_t(1:3, -1:nt+1))
    call calc_Ac_ext_t(0.0d0, dt, 0, nt, Ac_ext_t)

    energy0 = calc_energy(sbe, gs, Ac_ext_t(:, 0), icomm)

    ! Realtime calculation
    if (irank == 0) then
        open(unit=100, file=trim(base_directory)//trim(sysname)//"_sbe_rt.data")

        write(100, '(a)') "# Real time calculation:"
        write(100, '(a)') "# Ac_ext: External vector potential field"
        write(100, '(a)') "# E_ext: External electric field"
        write(100, '(a)') "# Ac_tot: Total vector potential field"
        write(100, '(a)') "# E_tot: Total electric field"
        write(100, '(a)') "# Jm: Matter current density (electrons)"
        write(100, '(4a)') "# 1:Time[a.u.] 2:Ac_ext_x[a.u.] 3:Ac_ext_y[a.u.] 4:Ac_ext_z[a.u.] ", &
            & "5:E_ext_x[a.u.] 6:E_ext_y[a.u.] 7:E_ext_z[a.u.] 8:Ac_tot_x[a.u.] ", &
            & "9:Ac_tot_y[a.u.] 10:Ac_tot_z[a.u.] 11:E_tot_x[a.u.] 12:E_tot_y[a.u.] ", &
            & "13:E_tot_z[a.u.]  14:Jm_x[a.u.] 15:Jm_y[a.u.] 16:Jm_z[a.u.]"

        open(unit=101, file=trim(base_directory)//trim(sysname)//"_sbe_rt_energy.data")
        write(101, '(a)') "# 1:Time[a.u.] 2:Eall[a.u.] 3:Eall-Eall0[a.u.]"
        write(101, '(f12.6,2(es24.15e3))') 0.0d0, energy0, 0.0d0

        open(unit=102, file=trim(base_directory)//trim(sysname)//"_sbe_nex.data")
        write(102, '(a)') "# 1:Time[a.u.] 2:nelec[a.u.] 3:nhole[a.u.]"

    end if
    

    do it = 1, nt
        t = dt * it
        call dt_evolve_bloch(sbe, gs, Ac_ext_t(:, it), dt)

        ! if (mod(it, 1) == 0) then
            E(:) = -(Ac_ext_t(:, it + 1) - Ac_ext_t(:, it - 1)) / (2 * dt)
            call calc_current_bloch(sbe, gs, Ac_ext_t(:, it), Jmat, icomm)
            energy = calc_energy(sbe, gs, Ac_ext_t(:, it), icomm)
            tr_all = calc_trace(sbe, gs, nstate_sbe, icomm)
            tr_vb = calc_trace(sbe, gs, nelec / 2, icomm)
            
            if (irank == 0) then
                write(100, '(f12.6,15(es24.15e3))') t, Ac_ext_t(:, it), E(:), Ac_ext_t(:, it), E(:), Jmat(:)
            end if
        ! end if

        if (mod(it, 100) == 0) then
            if (irank == 0) then
                write(101, '(f12.6,2(es24.15e3))') t, energy, energy-energy0
                write(102, '(f12.6,2(es24.15e3))') t, tr_all - tr_vb, nelec - tr_vb 
                write(*, '(i6,f12.6,es24.15e3)') it, t, tr_all
            end if
        end if

    end do

    if (irank == 0) then
        close(100)
        close(101)
        close(102)
    end if

    return
end subroutine realtime_main_ssbe

end module realtime_ssbe
