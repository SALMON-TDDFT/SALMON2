module multiscale_ssbe
    implicit none
contains

subroutine main_multiscale_ssbe(icomm)
    use omp_lib
    use communication
    use gs_info_ssbe
    use bloch_solver_ssbe
    use util_ssbe
    use em_field
    use fdtd_weyl
    use math_constants, only: pi
    use phys_constants, only: cspeed_au
    use salmon_global
    use filesystem, only: get_filehandle
    use datafile_ssbe
    implicit none
    integer, intent(in) :: icomm

    type(s_sbe_gs_info) :: gs
    type(s_sbe_bloch_solver), allocatable :: sbe(:)
    real(8) :: t
    real(8), allocatable :: Ac_ext_t(:, :)
    integer :: it, i
    integer :: nproc, irank, icomm_macro, nproc_macro, irank_macro
    real(8), allocatable :: Ac_macro(:, :), E_macro(:, :)
    real(8), allocatable :: Jmat_macro_tmp(:, :), Jmat_macro(:, :)
    integer, allocatable :: itbl_macro_coord(:, :)
    integer :: nmacro, nmacro_max
    integer :: imacro_min, imacro_max
    integer :: ix, iy, iz, mt, imacro, iobs, ii, jj
    real(8) :: jmat(3)

    type(s_fdtd_system) :: fs
    type(ls_fdtd_weyl) :: fw
    character(256) :: tmp

    logical :: flag_1d_model

    integer :: nstate_sbe

    integer :: fh_sbe_wave
    integer, allocatable :: fh_sbe_obs(:)
    integer, allocatable :: fh_sbe_rt(:)

    call comm_get_groupinfo(icomm, irank, nproc)

    ! FDTD setup
    fs%mg%nd = 1
    fs%mg%is(1) = 1 - abs(nxvac_m(1))
    fs%mg%ie(1) = nx_m + abs(nxvac_m(2))
    fs%mg%is(2) = 1 - abs(nyvac_m(1))
    fs%mg%ie(2) = ny_m + abs(nyvac_m(2))
    fs%mg%is(3) = 1 - abs(nzvac_m(1))
    fs%mg%ie(3) = nz_m + abs(nzvac_m(2))
    fs%mg%is_array(1:3) = fs%mg%is(1:3) - fs%mg%nd
    fs%mg%ie_array(1:3) = fs%mg%ie(1:3) + fs%mg%nd
    fs%hgs(1:3) = (/ hx_m, hy_m, hz_m /)
    fw%dt = dt
    fw%fdtddim = trim(fdtddim)
    fs%rlsize(1) = 1d0
    fs%rlsize(2) = 1d0
    fs%rlsize(3) = 1d0
    fs%origin(1) = 0d0
    fs%origin(2) = 0d0
    fs%origin(3) = 0d0
    do ii = 1, 3
        do jj = 1, 2
            fs%a_bc(ii,jj) = trim(boundary_em(ii,jj))
        end do
    end do
    call weyl_init(fs, fw)

    flag_1d_model = ((ny_m == 1) .and. (nz_m == 1))
    nstate_sbe = nstate

    ! Prepare external pulse
    mt = max(nt, int(abs(nxvac_m(1)) * fs%hgs(1) / cspeed_au / dt))
    allocate(Ac_ext_t(1:3, -1:mt+1))
    call calc_Ac_ext_t(0.0d0, dt, 0, mt, Ac_ext_t)    
    call set_incident_field(mt, Ac_ext_t, fs, fw)

    ! Macropoint and media setup
    nmacro_max = nx_m * ny_m * nz_m
    allocate(itbl_macro_coord(3, nmacro_max))
    if (irank == 0) then
        call read_media_info(nmacro_max, itbl_macro_coord, nmacro, fw)
    end if
    call comm_bcast(itbl_macro_coord, icomm, 0)
    call comm_bcast(nmacro, icomm, 0)

    if (nmacro > 0) then
        if (0.0d0 < al(1)) al_vec1(1:3) = (/ al(1), 0.0d0, 0.0d0 /)
        if (0.0d0 < al(2)) al_vec2(1:3) = (/ 0.0d0, al(2), 0.0d0 /)
        if (0.0d0 < al(3)) al_vec3(1:3) = (/ 0.0d0, 0.0d0, al(3) /)
        if (nstate_sbe < 1) nstate_sbe = nstate

        ! Read ground state electronic system:
        call init_sbe_gs_info(gs, sysname, base_directory, &
            & num_kgrid, nstate, nelec, &
            & al_vec1, al_vec2, al_vec3, &
            & .false., icomm)
        
        ! Distribute
        call distribute_macropoints(irank, nmacro, nproc, imacro_min, imacro_max)
        icomm_macro = comm_create_group(icomm, imacro_min, 0)
        call comm_get_groupinfo(icomm_macro, irank_macro, nproc_macro) 
        allocate(Ac_macro(1:3, nmacro))
        allocate(E_macro(1:3, nmacro))
        allocate(Jmat_macro_tmp(1:3, nmacro))
        allocate(Jmat_macro(1:3, nmacro))
        allocate(sbe(imacro_min:imacro_max))    

        ! Initialization of SBE solver and density matrix:
        do i = imacro_min, imacro_max
            call init_sbe_bloch_solver(sbe(i), gs, nstate_sbe, icomm_macro)
        end do
    end if

    call comm_sync_all(icomm)

    if (irank == 0) then
        write(tmp, "(a,a,a,a)") "mkdir ", trim(base_directory), trim(sysname), "_sbe_RT_Ac"
        call system(trim(tmp))
        write(tmp, "(a,a,a,a)") "mkdir ", trim(base_directory), trim(sysname), "_sbe_m"
        call system(trim(tmp))
        do imacro = 1, nmacro
            write(tmp, "(a,a,a,a,i6.6)") "mkdir ", trim(base_directory), trim(sysname), "_sbe_m/m", imacro
            call system(trim(tmp))
        end do
    end if

    call comm_sync_all(icomm)

    if (irank == 0) then
        ! _sbe_wave_rt.data
        write(tmp, "(a,a,a)") trim(base_directory), trim(sysname), "_sbe_wave.data"
        fh_sbe_wave = get_filehandle()
        open(unit=fh_sbe_wave, file=trim(tmp), action="write")
        call write_sbe_wave_header(fh_sbe_wave)
        ! _sbe_obs_?????.data
        allocate(fh_sbe_obs(1:obs_num_em))
        do iobs = 1, obs_num_em
            write(tmp, "(a,a,a,i3.3,a)") trim(base_directory), trim(sysname), "_sbe_obs_", iobs, "_at_point_rt.data"
            fh_sbe_obs(iobs) = get_filehandle()
            open(unit=fh_sbe_obs(iobs), file=trim(tmp), action="write")
            call write_sbe_obs_header(fh_sbe_obs(iobs))
        end do
    end if

    if (irank_macro == 0) then
        ! _sbe_rt.data
        allocate(fh_sbe_rt(imacro_min:imacro_max))
        do imacro = imacro_min, imacro_max
            write(tmp, "(a,a,a,i6.6,a,a,a)") trim(base_directory), &
                & trim(sysname), "_sbe_m/m", imacro, "/", trim(sysname), "_sbe_rt.data"
            fh_sbe_rt(imacro) = get_filehandle()
            open(unit=fh_sbe_rt(imacro), file=trim(tmp), action="write")
            call write_sbe_rt_header(fh_sbe_rt(imacro))
        end do
    end if

    call comm_sync_all(icomm)

    do it = 1, nt
        t = dt * it

        if (nmacro > 0) then
            ! Get vector potential field
            if (irank == 0) then
                do imacro = 1, nmacro
                    ix = itbl_macro_coord(1, imacro)
                    iy = itbl_macro_coord(2, imacro)
                    iz = itbl_macro_coord(3, imacro)
                    Ac_macro(1:3, imacro) = fw%vec_Ac_new%v(1:3, ix, iy, iz)
                    E_macro(1:3, imacro) = fw%vec_E%v(1:3, ix, iy, iz)
                end do
            end if
            call comm_bcast(Ac_macro, icomm, 0)
            call comm_bcast(E_macro, icomm, 0)
            
            Jmat_macro_tmp = 0.0d0
            do imacro = imacro_min, imacro_max
                call dt_evolve_bloch(sbe(imacro), gs, Ac_macro(1:3, imacro), dt)
                call calc_current_bloch(sbe(imacro), gs, Ac_macro(1:3, imacro), jmat, icomm_macro)
                if (irank_macro == 0) then
                    Jmat_macro_tmp(1:3, imacro) = jmat(1:3)
                end if
            end do
            call comm_summation(Jmat_macro_tmp, Jmat_macro, 3*nmacro, icomm)

            do imacro = 1, nmacro
                ix = itbl_macro_coord(1, imacro)
                iy = itbl_macro_coord(2, imacro)
                iz = itbl_macro_coord(3, imacro)
                fw%vec_j_em%v(:, ix, iy, iz) = -Jmat_macro(:, imacro)
            end do
        end if

        call weyl_calc(fs, fw)

        if (irank == 0) then
            if (mod(it, out_ms_step) == 0) then
                call write_Ac_field(it, fs, fw)
            end if
            if (mod(it, 10) == 0) then
                write(*, "(a,i6)") "Time step = ", it
                if (flag_1d_model) call write_wave_data_file(fh_sbe_wave, it, fs, fw)
                do iobs = 1, obs_num_em
                    ix = int(obs_loc_em(iobs, 1) / fs%hgs(1))
                    iy = int(obs_loc_em(iobs, 2) / fs%hgs(2))
                    iz = int(obs_loc_em(iobs, 3) / fs%hgs(3))
                    call write_sbe_obs_line(fh_sbe_obs(iobs), it * dt, &
                        &  -(fw%vec_Ac_new%v(:, ix, iy, iz) - fw%vec_Ac_old%v(:, ix, iy, iz)) / (2 * dt))
                end do            
            end if
            if (mod(it, 500) == 0) then
                flush(fh_sbe_wave)
                do iobs = 1, obs_num_em
                    flush(fh_sbe_obs(iobs))
                end do
            end if
        end if

        if (irank_macro == 0) then
            if (nmacro > 0) then
                do imacro = imacro_min, imacro_max
                    call write_sbe_rt_line(fh_sbe_rt(imacro), t, &
                        & Ac_macro(1:3, imacro), E_macro(1:3, imacro), &
                        & Ac_macro(1:3, imacro), E_macro(1:3, imacro), & 
                        & Jmat_macro(1:3, imacro))
                end do
                if (mod(it, 500) == 0) then
                    do imacro = imacro_min, imacro_max
                        flush(fh_sbe_rt(imacro))
                    end do
                end if
            end if
        end if
    end do

    call comm_sync_all(icomm)

    if (irank == 0) then
        close(fh_sbe_wave)
        do iobs = 1, obs_num_em
            close(fh_sbe_obs(iobs))
        end do
    end if

    if (irank_macro == 0) then
        if (nmacro > 0) then
            do imacro = imacro_min, imacro_max
                close(fh_sbe_rt(imacro))
            end do
        end if
    end if

    call comm_sync_all(icomm)

    return
end subroutine main_multiscale_ssbe


subroutine distribute_macropoints(irank, nmacro, nproc, imacro_min, imacro_max)
    use util_ssbe, only: split_range
    use communication, only: comm_create_group, comm_get_groupinfo
    use salmon_global
    implicit none
    integer, intent(in) :: irank
    integer, intent(in) :: nmacro
    integer, intent(in) :: nproc
    integer, intent(out) :: imacro_min
    integer, intent(out) :: imacro_max

    integer :: i
    integer :: itbl_macro_min(0:nproc-1)
    integer :: itbl_macro_max(0:nproc-1)
    integer :: itbl_rank_min(1:nmacro)
    integer :: itbl_rank_max(1:nmacro)

    if (nproc <= nmacro) then
        call split_range(1, nmacro, nproc, itbl_macro_min, itbl_macro_max)
        imacro_min = itbl_macro_min(irank)
        imacro_max = itbl_macro_max(irank)    
    else
        call split_range(0, nproc-1, nmacro, itbl_rank_min, itbl_rank_max)
        do i = 1, nmacro
            if ((itbl_rank_min(i) <= irank) .and. (irank <= itbl_rank_max(i))) then
                imacro_min = i
                imacro_max = i
            end if
        end do
    end if

    return
end subroutine distribute_macropoints


subroutine read_media_info(nmacro_max, itbl_macro_coord, nmacro, fw)
    use fdtd_weyl, only: ls_fdtd_weyl
    use salmon_global
    implicit none
    integer, intent(in) :: nmacro_max
    integer, intent(out) :: itbl_macro_coord(1:3, nmacro_max)
    integer, intent(out) :: nmacro
    type(ls_fdtd_weyl), intent(inout) :: fw

    integer :: i, imacro, n, ix, iy, iz, itype

    fw%epsilon%f(:, :, :) = 1.0d0
    imacro = 0
    if (len_trim(file_macropoint) > 0) then
        open(99, file=trim(file_macropoint), action="read")
        read(99, *) n
        do i = 1, n
            read(99, *) ix, iy, iz, itype
            if (ix < 1) stop "Error: invalid range!"
            if (iy < 1) stop "Error: invalid range!"
            if (iz < 1) stop "Error: invalid range!"
            if (ix > nx_m) stop "Error: invalid range!"
            if (iy > ny_m) stop "Error: invalid range!"
            if (iz > nz_m) stop "Error: invalid range!"
            if (itype > 0) then
                fw%epsilon%f(ix, iy, iz) = epsilon_em(itype)
            else if (itype == 0) then
                fw%epsilon%f(ix, iy, iz) = 1.0d0
            else
                imacro = imacro + 1
                if (imacro > nmacro_max) stop "Error: number of macropoints is too large!" 
                itbl_macro_coord(1:3, imacro) = (/ ix, iy, iz /)
            end if
        end do
        close(99)
    else
        do iz = 1, nz_m
        do iy = 1, ny_m
        do ix = 1, nx_m
            imacro = imacro + 1
            if (imacro > nmacro_max) stop "Error: number of macropoints is too large!" 
            itbl_macro_coord(1:3, imacro) = (/ ix, iy, iz /)
        end do
        end do
        end do
    end if
    nmacro = imacro
end subroutine 


! subroutine read_shape_info(fs, fw)
!     use input_parameter, only: file_macropoint, epsilon_em, nx_m, ny_m, nz_m
!     use fdtd_weyl, only: ls_fdtd_weyl
!     implicit none
!     integer, intent(in) :: nmacro_max
!     integer, intent(out) :: itbl_macro_coord(1:3, nmacro_max)
!     integer, intent(out) :: nmacro
!     type(s_fdtd_system), intent(in) :: fs
!     type(ls_fdtd_weyl), intent(inout) :: fw

!     integer :: itbl_media( &
!         & fw%mg%is(1):fw%mg%ie(1), &
!         & fw%mg%is(2):fw%mg%ie(2), &
!         & fw%mg%is(3):fw%mg%ie(3), &
!         & )

!     do i = 1, n_s
!         select case trim(typ_s(i))
!         case "ellipsoid"
!         end select
!     end do
! end subroutine read_shape_info
        

subroutine set_incident_field(mt, Ac, fs, fw)
    use fdtd_weyl, only: s_fdtd_system, ls_fdtd_weyl
    use phys_constants, only: cspeed_au
    use salmon_global
    implicit none
    integer, intent(in) :: mt
    real(8), intent(in) :: Ac(1:3, 0:mt)
    type(s_fdtd_system), intent(in) :: fs
    type(ls_fdtd_weyl), intent(inout) :: fw

    integer :: ix, iy, iz, it
    real(8) :: r_it, w, x

    do ix = fs%mg%is_array(1), 0
        do iy = fs%mg%is_array(2), fs%mg%ie_array(2)
            do iz = fs%mg%is_array(3), fs%mg%ie_array(3)
                x = ix * fs%hgs(1)
                r_it = (-x / cspeed_au / fw%dt)
                it = int(r_it)
                if ((0 <= it) .and. (it < mt)) then
                    w = r_it - it
                    fw%vec_Ac_new%v(:, ix, iy, iz) = w * Ac(:, it+1) + (1.0d0-w) * Ac(:, it)
                end if
                r_it = r_it - 1.0d0
                it = int(r_it)
                if ((0 <= it) .and. (it < mt)) then
                    w = r_it - it
                    fw%vec_Ac%v(:, ix, iy, iz) = w * Ac(:, it+1) + (1.0d0-w) * Ac(:, it)
                end if
                r_it = r_it - 1.0d0
                it = int(r_it)
                if ((0 <= it) .and. (it < mt)) then
                    w = r_it - it
                    fw%vec_Ac_old%v(:, ix, iy, iz) = w * Ac(:, it+1) + (1.0d0-w) * Ac(:, it)
                end if
            end do
        end do
    end do
end subroutine set_incident_field

subroutine write_Ac_field(iit, fs, fw)
    use fdtd_weyl, only: s_fdtd_system, ls_fdtd_weyl
    use phys_constants, only: cspeed_au
    use filesystem
    use salmon_global
    implicit none
    integer, intent(in) :: iit
    character(256) :: file_Ac_data
    type(s_fdtd_system), intent(in) :: fs
    type(ls_fdtd_weyl), intent(in) :: fw
    integer :: ix, iy, iz
    integer :: fh

    write(file_Ac_data, "(a,a,a,a,i6.6,a)") trim(sysname), "_sbe_RT_Ac/", trim(sysname), "_Ac_", iit, ".data"
    fh = get_filehandle()
    open(unit=fh, file=trim(file_Ac_data), action="write")
    write(fh, '(a)') "# Multiscale TDDFT calculation"
    write(fh, '(a)') "# IX, IY, IZ: FDTD Grid index"
    write(fh, '(a)') "# x, y, z: Coordinates"
    write(fh, '(a)') "# Ac: Vector potential field"
    write(fh, '(a)') "# E: Electric field"
    write(fh, '(a)') "# J_em: Electromagnetic current density"
    write(fh, '("#",99(1X,I0,":",A,"[",A,"]"))') &
        & 1, "IX", "none", &
        & 2, "IY", "none", &
        & 3, "IZ", "none", &
        & 4, "Ac_x", "[au]", &
        & 5, "Ac_y", "[au]", &
        & 6, "Ac_z", "[au]", &
        & 7, "E_x", "[au]", &
        & 8, "E_y", "[au]", &
        & 9, "E_z", "[au]", &
        & 10, "B_x", "a.u.", &
        & 11, "B_y", "a.u.", &
        & 12, "B_z", "a.u.", &
        & 13, "Jem_x", "[au]", &
        & 14, "Jem_y", "[au]", &
        & 15, "Jem_z", "[au]", &
        & 16, "E_em", "[au]" // "/vol", &
        & 17, "E_abs", "[au]" //  "/vol"
    do iz = max(fs%mg%is(3), out_ms_region_iz_m(1)), min(fs%mg%ie(3), out_ms_region_iz_m(2))
    do iy = max(fs%mg%is(2), out_ms_region_iy_m(1)), min(fs%mg%ie(2), out_ms_region_iy_m(2))
    do ix = max(fs%mg%is(1), out_ms_region_ix_m(1)), min(fs%mg%ie(1), out_ms_region_ix_m(2))
        write(fh, "(3i9,99es25.15e4)") &
        ix, iy, iz, &
        fw%vec_Ac%v(:, ix, iy, iz), &
        fw%vec_e%v(:, ix, iy, iz), &
        fw%vec_h%v(:, ix, iy, iz), &
        fw%vec_j_em%v(:, ix, iy, iz)
    end do
    end do
    end do
    close(fh)
    return
end subroutine write_Ac_field

 ! Experimetal Implementation of Incident/Reflection/Transmit field output
subroutine write_wave_data_file(fh, iit, fs, fw)
    use fdtd_weyl, only: s_fdtd_system, ls_fdtd_weyl
    use phys_constants, only: cspeed_au
    use salmon_global
    use datafile_ssbe
    implicit none
    integer, intent(in) :: fh
    integer, intent(in) :: iit
    type(s_fdtd_system), intent(in) :: fs
    type(ls_fdtd_weyl), intent(in) :: fw
    ! real(8) :: dt
    real(8) :: e_inc(3)
    real(8) :: e_ref(3)
    real(8) :: e_tra(3)
    real(8) :: dt_Ac(3)
    real(8) :: dx_Ac(3)
    integer :: iiy, iiz
    iiy = fs%mg%is(2)
    iiz = fs%mg%is(3)

    ! Left side boundary:
    dx_Ac(:) = (fw%vec_Ac%v(:,0,iiy,iiz) - fw%vec_Ac%v(:,-1,iiy,iiz)) / fs%hgs(1)
    dt_Ac(:) = (0.5d0 * (fw%vec_Ac_new%v(:,0,iiy,iiz) + fw%vec_Ac_new%v(:,-1,iiy,iiz)) & 
        & - 0.5d0 * (fw%vec_Ac_old%v(:,0,iiy,iiz) + fw%vec_Ac_old%v(:,-1,iiy,iiz))) / (2 * dt)
    e_inc(:) = -0.5d0 * (dt_Ac - cspeed_au * dx_Ac)
    e_ref(:) = -0.5d0 * (dt_Ac + cspeed_au * dx_Ac)
    ! Right side boundary:
    dx_Ac(:) = (fw%vec_Ac%v(:,nx_m+2,iiy,iiz) - fw%vec_Ac%v(:,nx_m+1,iiy,iiz)) / fs%hgs(1)
    dt_Ac(:) = (0.5d0 * (fw%vec_Ac_new%v(:,nx_m+2,iiy,iiz) + fw%vec_Ac_new%v(:,nx_m+1,iiy,iiz)) & 
        & - 0.5d0 * (fw%vec_Ac_old%v(:,nx_m+2,iiy,iiz) + fw%vec_Ac_old%v(:,nx_m+1,iiy,iiz))) / (2 * dt)
    e_tra(:) = -0.5d0 * (dt_Ac - cspeed_au * dx_Ac)
    call write_sbe_wave_line(fh, iit * dt, e_inc, e_ref, e_tra)
    return
end subroutine write_wave_data_file



end module multiscale_ssbe
