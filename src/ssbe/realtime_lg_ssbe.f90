module realtime_lg_ssbe
    implicit none
contains

subroutine main_realtime_lg_ssbe(icomm)
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

    ! Prepare qnm
    call prepare_qnm(sbe, gs, icomm)

    call adams_moulton_coefs(bj_am)

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
        write(*, "(a)") " time-step time[fs] Current(xyz)[a.u.]                     Electric Field(z) [au]"
        write(*, "(a)") "---------------------------------------------------------------------------------------"
    end if

    call comm_sync_all(icomm)


    ! Realtime calculation
    do it = 1, nt
        t = dt * it
        E(:) = -(Ac_ext_t(:, it + 1) - Ac_ext_t(:, it - 1)) / (2 * dt)
        call dt_evolve_bloch_lg(sbe, gs, E(:), bj_am, dt)
        call calc_current_bloch_lg(sbe, gs, jmat)

        if (mod(it, 10) == 0) then
            if (irank == 0) then
                write(*, "(i6,f12.3,3es12.3,es12.3)") it, t, Jmat(1:3), E(3)
            end if
        end if
    end do

    return
end subroutine main_realtime_lg_ssbe


subroutine adams_moulton_coefs(bj_am)
  implicit none
  real(8) :: bj_am(8,8)

  bj_am(1,1) = 1.d0

  !bj_am(1,2) = 3.d0/2.d0
  !bj_am(2,2) = -1.d0/2.d0

  !bj_am(1,3) = 23.d0/12.d0
  !bj_am(2,3) = -4.d0/3.d0
  !bj_am(3,3) = 5.d0/12.d0

  bj_am(1,4) = 55.d0/24.d0
  bj_am(2,4) = -59.d0/24.d0
  bj_am(3,4) = 37.d0/24.d0
  bj_am(4,4) = -3.d0/8.d0

  !bj_am(1,5) = 1901.d0/720.d0
  !bj_am(2,5) = -1387.d0/360.d0
  !bj_am(3,5) = 109.d0/30.d0
  !bj_am(4,5) = -637.d0/360.d0
  !bj_am(5,5) = 251.d0/720.d0

  !bj_am(1,6) = 4277.d0/1440.d0
  !bj_am(2,6) = -2641.d0/480.d0
  !bj_am(3,6) = 4991.d0/720.d0
  !bj_am(4,6) = -3649.d0/720.d0
  !bj_am(5,6) = 959.d0/480.d0
  !bj_am(6,6) = -95.d0/288.d0

  !bj_am(1,7) = 198721.d0/60480.d0
  !bj_am(2,7) = -18637.d0/2520.d0
  !bj_am(3,7) = 235183.d0/20160.d0
  !bj_am(4,7) = -10754.d0/945.d0
  !bj_am(5,7) = 135713.d0/20160.d0
  !bj_am(6,7) = -5603.d0/2520.d0
  !bj_am(7,7) = 19087.d0/60480.d0

  bj_am(1,8) = 16083.d0/4480.d0
  bj_am(2,8) = -1152169.d0/120960.d0
  bj_am(3,8) = 242653.d0/13440.d0
  bj_am(4,8) = -296053.d0/13440.d0
  bj_am(5,8) = 2102243.d0/120960.d0
  bj_am(6,8) = -115747.d0/13440.d0
  bj_am(7,8) = 32863.d0/13440.d0
  bj_am(8,8) = -5257.d0/17280.d0

end subroutine adams_moulton_coefs


end module realtime_lg_ssbe
