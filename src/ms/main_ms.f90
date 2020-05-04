!
!  Copyright 2019-2020 SALMON developers
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

#include "config.h"

subroutine main_ms
use math_constants, only: pi
use salmon_global
use structures
use inputoutput, only: nx_m, ny_m, nz_m, dt
use communication, only: comm_is_root, comm_sync_all, comm_create_group_byid, comm_get_groupinfo, comm_sync_all, comm_summation, comm_bcast
use salmon_xc, only: finalize_xc
use timer
use write_sub, only: write_response_0d,write_response_3d,write_pulse_0d,write_pulse_3d
use initialization_rt_sub
use checkpoint_restart_sub
use fdtd_weyl, only: ls_fdtd_weyl, weyl_init, weyl_calc, weyl_finalize
use parallelization, only: nproc_id_global, nproc_size_global, nproc_group_global
use filesystem, only: create_directory, get_filehandle
use phys_constants, only: cspeed_au
use em_field, only: calc_Ac_ext
implicit none

type(s_rgrid) :: lg
type(s_rgrid) :: mg
type(s_dft_system)  :: system
type(s_rt) :: rt
type(s_parallel_info) :: info
type(s_poisson) :: poisson
type(s_stencil) :: stencil
type(s_xc_functional) :: xc_func
type(s_reciprocal_grid) :: fg
type(s_ewald_ion_ion) :: ewald
type(s_dft_energy) :: energy
type(s_md) :: md
type(s_ofile) :: ofl
type(s_scalar) :: sVpsl
type(s_scalar) :: srho,sVh,sVh_stock1,sVh_stock2,Vbox
type(s_scalar),allocatable :: srho_s(:),V_local(:),sVxc(:)
type(s_dmatrix) :: dmat
type(s_orbital) :: spsi_in,spsi_out
type(s_orbital) :: tpsi ! temporary wavefunctions
type(s_sendrecv_grid) :: srg,srg_scalar
type(s_pp_info) :: pp
type(s_pp_grid) :: ppg
type(s_pp_nlcc) :: ppn
type(s_singlescale) :: singlescale

type(s_fdtd_system) :: fs
type(ls_fdtd_weyl) :: fw
type(s_multiscale) :: ms


integer :: Mit, itt, itotNtime
integer :: nntime

integer :: i, ix, iy, iz

integer, allocatable :: iranklists(:)

real(8), allocatable :: Ac_inc(:, :)

character(256) :: file_debug_log
integer :: nmacro_mygroup, isize_mygroup

if (.not. check_input_variables()) return

!! Open logfile for debugging
! write(file_debug_log, "('ms_debug', i3.3, '.log')") nproc_id_global
! open(unit=9999, file=file_debug_log)
! write(9999, *) 'logging start'; flush(9999)

call timer_begin(LOG_TOTAL)
! Initialization

call initialization_ms()
! if (comm_is_root(ms%id_ms_world)) call print_header()

call comm_sync_all
call timer_enable_sub
call timer_begin(LOG_RT_ITERATION)

TE : do itt=Mit+1,itotNtime
    call time_evolution_step_ms()
    call write_RT_Ac_file()

   if((checkpoint_interval >= 1) .and. (mod(itt,checkpoint_interval) == 0)) &
        & call checkpoint_ms()
end do TE

call timer_end(LOG_RT_ITERATION)
call timer_disable_sub

call timer_end(LOG_TOTAL)

if(write_rt_wfn_k=='y')then
    call checkpoint_ms(ofl%dir_out_restart)
end if

call finalize_xc(xc_func)

! close(gp)

contains





function macropoint_in_mygroup(imacro) result(r)
    implicit none
    integer, intent(in) :: imacro
    logical :: r
    r = (ms%imacro_mygroup_s <= imacro) &
        & .and. (imacro <= ms%imacro_mygroup_e)
    return
end function macropoint_in_mygroup





! Create the base_directory name for each macropoints:
function base_directory_macro(imacro) result(r)
    implicit none
    integer, intent(in) :: imacro
    character(256) :: r
    write(r, '(a, a, a, i6.6, a)') trim(ms%base_directory), trim(sysname), '_m/m', imacro, '/'
    return
end function base_directory_macro






function directory_read_data_macro(basedir, imacro) result(r)
    implicit none
    character(*), intent(in) :: basedir
    integer, intent(in) :: imacro
    character(256) :: r
    write(r, '(a,a,i6.6,a)') trim(basedir), 'm', imacro, '/'
    return
end function directory_read_data_macro





function restart_directory_macro(iter, imacro) result(r)
    implicit none
    integer, intent(in) :: iter
    integer, intent(in) :: imacro
    character(256) :: r
    write(r,'(a,i6.6,a)') "checkpoint_ms_", iter, "/"
    r = trim(directory_read_data_macro(trim(r), imacro))
    return
end function restart_directory_macro


subroutine initialization_ms()
    implicit none
    integer :: ii, jj

    ! Store global information
    ms%nmacro = nx_m * ny_m * nz_m
    ms%base_directory = trim(base_directory)
    ms%base_directory_RT_Ac = trim(ms%base_directory) // trim(sysname) // "_RT_Ac/"
    ms%icomm_ms_world = nproc_group_global
    ms%isize_ms_world = nproc_size_global
    ms%id_ms_world = nproc_id_global
    ms%directory_read_data = trim(directory_read_data)
    

    if (ms%nmacro <= ms%isize_ms_world) then
        if (mod(ms%isize_ms_world, ms%nmacro) == 0) then
            nmacro_mygroup = 1
            isize_mygroup = ms%isize_ms_world / ms%nmacro
            ms%imacro_mygroup_s = ms%id_ms_world / isize_mygroup + 1
            ms%imacro_mygroup_e = ms%imacro_mygroup_s
            ms%id_mygroup_s = (ms%imacro_mygroup_s - 1) * isize_mygroup
            ms%id_mygroup_e = ms%id_mygroup_s + isize_mygroup - 1
        else
            stop "mod(number_of_procs, number_of_points) must be 0!"
        end if
    else
        stop "number of procs must be larger than number of points!"
    end if
    
    if (ms%nmacro < 1) stop "Invalid macropoint number"
    
    allocate(iranklists(isize_mygroup))
    do i = 1, isize_mygroup
        iranklists(i) = ms%id_mygroup_s + (i - 1)
    end do
    
    ms%icomm_macropoint = comm_create_group_byid(ms%icomm_ms_world, iranklists)
    call comm_get_groupinfo(ms%icomm_macropoint, ms%id_macropoint, ms%isize_macropoint)

    fs%mg%ndir = 3
    fs%mg%nd = 1
    
    fs%hgs(1:3) = dl_em(1:3)
    if (0d0 < hx_m) fs%hgs(1) = hx_m
    if (0d0 < hy_m) fs%hgs(2) = hy_m
    if (0d0 < hz_m) fs%hgs(3) = hz_m

    fw%dt = dt
    fw%fdtddim = '1d'
    fs%mg%is(1) = -nxvacl_m
    fs%mg%ie(1) = nx_m+nxvacr_m
    fs%mg%is(2) = 1
    fs%mg%ie(2) = ny_m
    fs%mg%is(3) = 1
    fs%mg%ie(3) = nz_m
    fs%mg%is_overlap(1:3) = fs%mg%is(1:3) - fs%mg%nd
    fs%mg%ie_overlap(1:3) = fs%mg%ie(1:3) + fs%mg%nd
    fs%mg%is_array(1:3) = fs%mg%is_overlap(1:3)
    fs%mg%ie_array(1:3) = fs%mg%ie_overlap(1:3)
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
    fs%imedia(:,:,:) = 0

    call Weyl_init(fs, fw)

    allocate(ms%curr(1:3, 1:ms%nmacro))
    allocate(ms%vec_Ac(1:3, 1:ms%nmacro))
    allocate(ms%vec_Ac_old(1:3, 1:ms%nmacro))
    allocate(ms%ixyz_tbl(1:3, 1:ms%nmacro))
    allocate(ms%imacro_tbl( &
        & fs%mg%is(1):fs%mg%ie(1), &
        & fs%mg%is(2):fs%mg%ie(2), &
        & fs%mg%is(3):fs%mg%ie(3)))

    i = 1
    do iz = 1, nz_m
        do iy = 1, ny_m
            do ix = 1, nx_m
                ms%imacro_tbl(ix, iy, iz) = i
                ms%ixyz_tbl(1, i) = ix
                ms%ixyz_tbl(2, i) = iy
                ms%ixyz_tbl(3, i) = iz
                i = i + 1
            end do
        end do
    end do


    if (comm_is_root(ms%id_ms_world)) then
        call create_directory(trim(ms%base_directory_RT_Ac))
        do i = 1, ms%nmacro
            call create_directory(trim(base_directory_macro(i)))
        end do
    end if

    do i = 1, ms%nmacro

        call comm_sync_all()

        if (macropoint_in_mygroup(i)) then
            ! Override global variables
            base_directory = trim(base_directory_macro(i))
            nproc_group_global = ms%icomm_macropoint
            nproc_id_global = ms%id_macropoint
            nproc_size_global = ms%isize_macropoint

            if (yn_restart == 'y') then
                directory_read_data = trim(directory_read_data_macro(trim(ms%directory_read_data), i))
            end if
                    
            ! Initializa TDKS system
            call initialization_rt( Mit, itotNtime, system, energy, ewald, rt, md, &
                                    singlescale,  &
                                    stencil, fg, poisson,  &
                                    lg, mg,   &
                                    info,  &
                                    xc_func, dmat, ofl,  &
                                    srg, srg_scalar,  &
                                    spsi_in, spsi_out, tpsi, srho, srho_s,  &
                                    V_local, Vbox, sVh, sVh_stock1, sVh_stock2, sVxc, sVpsl,&
                                    pp, ppg, ppn )

            ! Override global variables (restore)
            base_directory = trim(ms%base_directory)
            nproc_group_global = ms%icomm_ms_world
            nproc_id_global = ms%id_ms_world
            nproc_size_global = ms%isize_ms_world        
            directory_read_data = trim(ms%directory_read_data)
        end if
    end do

    call incident()

    return    
end subroutine initialization_ms





subroutine print_header()
    implicit none
    !(header of standard output)
    ! if(comm_is_root(nproc_id_global))then
        write(*,*)
        select case(iperiodic)
        case(0)
        write(*,'(1x,a10,a11,a11,a10,a48,a15,a18,a10)') &
                    "time-step ", "time[fs]", "macropoint", &
                    "Dipole moment(xyz)[A]"     &
                    ,"electrons", "Total energy[eV]", "iterVh"
        case(3)
        write(*,'(1x,a10,a11,a10,a48,a15,a18)')   &
                    "time-step", "time[fs] ", "macropoint", &
                    "Current(xyz)[a.u.]",     &
                    "electrons", "Total energy[eV] "
        end select
        write(*,'("#",7("----------"))')
    ! endif
end subroutine print_header





subroutine time_evolution_step_ms
    implicit none
    integer :: ii, iimacro, iix, iiy, iiz
    real(8) :: curr_tmp(3, ms%nmacro)

    ! Override Global Variables
    nproc_group_global = ms%icomm_macropoint
    nproc_id_global = ms%id_macropoint
    nproc_size_global = ms%isize_macropoint

    fw%Ac_inc(:) = Ac_inc(:, itt)
    fw%Ac_inc_old(:) = Ac_inc(:, itt-1)
    call weyl_calc(fs, fw)

    do ii = 1, ms%nmacro
        iix = ms%ixyz_tbl(1, ii)
        iiy = ms%ixyz_tbl(2, ii)
        iiz = ms%ixyz_tbl(3, ii)
        ms%vec_Ac(1:3, ii) = fw%vec_Ac%v(1:3, iix, iiy, iiz)
        ms%vec_Ac_old(1:3, ii) = fw%vec_Ac_old%v(1:3, iix, iiy, iiz)
    end do

    if (ms%imacro_mygroup_e - ms%imacro_mygroup_s + 1 > 1) then
        stop "ERROR! Unsupported paralization scheme!"
    else
        iimacro = ms%imacro_mygroup_s

        rt%Ac_ext(:, itt-1) = ms%vec_Ac_old(1:3, iimacro)
        rt%Ac_ext(:, itt) = ms%vec_Ac(1:3, iimacro)
        rt%Ac_ext(:, itt+1) = ms%vec_Ac(1:3, iimacro) &
            & + (ms%vec_Ac(1:3, iimacro) -  ms%vec_Ac_old(1:3, iimacro))

        if(mod(itt,2)==1)then
            call time_evolution_step(Mit,itotNtime,itt,lg,mg,system,rt,info,stencil,xc_func &
            & ,srg,srg_scalar,pp,ppg,ppn,spsi_in,spsi_out,tpsi,srho,srho_s,V_local,Vbox,sVh,sVh_stock1,sVh_stock2,sVxc &
            & ,sVpsl,dmat,fg,energy,ewald,md,ofl,poisson,singlescale)
        else
            call time_evolution_step(Mit,itotNtime,itt,lg,mg,system,rt,info,stencil,xc_func &
            & ,srg,srg_scalar,pp,ppg,ppn,spsi_out,spsi_in,tpsi,srho,srho_s,V_local,Vbox,sVh,sVh_stock1,sVh_stock2,sVxc &
            & ,sVpsl,dmat,fg,energy,ewald,md,ofl,poisson,singlescale)
        end if

        curr_tmp(:, :) = 0d0
        if (comm_is_root(ms%id_macropoint)) &
            & curr_tmp(1:3, iimacro) = rt%curr(1:3, itt)
        call comm_summation(curr_tmp, ms%curr, 3 * ms%nmacro, ms%icomm_ms_world)

    end if

    do ii = 1, ms%nmacro
        iix = ms%ixyz_tbl(1, ii)
        iiy = ms%ixyz_tbl(2, ii)
        iiz = ms%ixyz_tbl(3, ii)
        fw%vec_j_em%v(1:3, iix, iiy, iiz) = - ms%curr(1:3, ii)
    end do

    ! Override Global Variables (Repair)
    nproc_group_global = ms%icomm_ms_world
    nproc_id_global = ms%id_ms_world
    nproc_size_global = ms%isize_ms_world

    return
end subroutine time_evolution_step_ms





subroutine checkpoint_ms(odir)
    implicit none
    character(*),optional, intent(in) :: odir
    character(256) :: idir
    integer :: fh_bin
    
    ! Override Global Variables
    base_directory = trim(base_directory_macro(i))
    nproc_group_global = ms%icomm_macropoint
    nproc_id_global = ms%id_macropoint
    nproc_size_global = ms%isize_macropoint
    
    call timer_begin(LOG_CHECKPOINT_SYNC)
    call timer_begin(LOG_CHECKPOINT_SELF)
    
    do i = 1, ms%nmacro
        if (present(odir)) then
            idir = trim(directory_read_data_macro(trim(odir), i))
        else
            idir = trim(restart_directory_macro(itt, i))
        end if

        if (comm_is_root(ms%id_ms_world)) call create_directory(trim(idir))

        call comm_sync_all()

        if (macropoint_in_mygroup(i)) then
            if (mod(itt,2)==1) then
                call checkpoint_rt(lg,mg,system,info,spsi_out,itt,sVh_stock1,sVh_stock2,singlescale,idir)
            else
                call checkpoint_rt(lg,mg,system,info,spsi_in, itt,sVh_stock1,sVh_stock2,singlescale,idir)
            endif
        end if
    end do

    if (comm_is_root(ms%id_ms_world)) then
        fh_bin = get_filehandle()

        open(fh_bin,file=trim(idir) // '../vec_Ac.bin',form='unformatted')
        write(fh_bin) fw%vec_Ac%v
        close(fh_bin)

        open(fh_bin,file=trim(idir) // '../vec_Ac_old.bin',form='unformatted')
        write(fh_bin) fw%vec_Ac_old%v
        close(fh_bin)

        open(fh_bin,file=trim(idir) // '../vec_j_em.bin',form='unformatted')
        write(fh_bin) fw%vec_j_em%v
        close(fh_bin)
    end if


    call timer_end(LOG_CHECKPOINT_SELF)
    call comm_sync_all
    call timer_end(LOG_CHECKPOINT_SYNC)

    ! Override Global Variables (Repair)
    base_directory = trim(ms%base_directory)
    nproc_group_global = ms%icomm_ms_world
    nproc_id_global = ms%id_ms_world
    nproc_size_global = ms%isize_ms_world        

    return
end subroutine checkpoint_ms




subroutine write_RT_Ac_file()
    use inputoutput, only: t_unit_ac, t_unit_current, t_unit_elec, t_unit_length
    implicit none
    integer ::  iix, iiy, iiz
    character(256) :: filename


    if (ms%id_ms_world == 0) then
        if (mod(itt, 100) == 0) then
            write(filename, '(a, a, a, i6.6, a)') trim(ms%base_directory_RT_Ac), trim(sysname), "_Ac_",  itt, '.data'
            open(8888, file=trim(filename))
            write(8888, '(a)') "# Multiscale TDDFT calculation"
            write(8888, '(a)') "# IX, IY, IZ: FDTD Grid index"
            write(8888, '(a)') "# x, y, z: Coordinates"
            write(8888, '(a)') "# Ac: Vector potential field"
            write(8888, '(a)') "# E: Electric field"
            write(8888, '(a)') "# J_em: Electromagnetic current density"
            write(8888, '("#",99(1X,I0,":",A,"[",A,"]"))') &
              & 1, "IX", "none", &
              & 2, "IY", "none", &
              & 3, "IZ", "none", &
            !   & 4, "x", trim(t_unit_length%name), &
            !   & 5, "y", trim(t_unit_length%name), &
            !   & 6, "z", trim(t_unit_length%name), &
              & 4, "Ac_x", trim(t_unit_ac%name), &
              & 5, "Ac_y", trim(t_unit_ac%name), &
              & 6, "Ac_z", trim(t_unit_ac%name), &
              & 7, "E_x", trim(t_unit_elec%name), &
              & 8, "E_y", trim(t_unit_elec%name), &
              & 9, "E_z", trim(t_unit_elec%name), &
              & 10, "B_x", "a.u.", &
              & 11, "B_y", "a.u.", &
              & 12, "B_z", "a.u.", &
              & 13, "Jm_x", trim(t_unit_current%name), &
              & 14, "Jm_y", trim(t_unit_current%name), &
              & 15, "Jm_z", trim(t_unit_current%name)

            do iiz = fs%mg%is(3), fs%mg%ie(3)
            do iiy = fs%mg%is(2), fs%mg%ie(2)
            do iix = fs%mg%is(1), fs%mg%ie(1)
                write(8888, '(3(i6, 1x), 12(e23.15e3, 1x))')  &
                    & iix, & 
                    & iiy, & 
                    & iiz, & 
                    ! & iix * fs%hgs(1) * t_unit_length%conv, &
                    ! & iiy * fs%hgs(2) * t_unit_length%conv, &
                    ! & iiz * fs%hgs(3) * t_unit_length%conv, &
                    & fw%vec_Ac%v(1, iix, iiy, iiz) * t_unit_ac%conv, &
                    & fw%vec_Ac%v(2, iix, iiy, iiz) * t_unit_ac%conv, &
                    & fw%vec_Ac%v(3, iix, iiy, iiz) * t_unit_ac%conv, &
                    & fw%vec_E%v(1, iix, iiy, iiz) * t_unit_elec%conv, &
                    & fw%vec_E%v(2, iix, iiy, iiz) * t_unit_elec%conv, &
                    & fw%vec_E%v(3, iix, iiy, iiz) * t_unit_elec%conv, &
                    & fw%vec_H%v(1, iix, iiy, iiz), &
                    & fw%vec_H%v(2, iix, iiy, iiz), &
                    & fw%vec_H%v(3, iix, iiy, iiz), &
                    & fw%vec_j_em%v(1, iix, iiy, iiz) * t_unit_current%conv, &
                    & fw%vec_j_em%v(2, iix, iiy, iiz) * t_unit_current%conv, &
                    & fw%vec_j_em%v(3, iix, iiy, iiz) * t_unit_current%conv
            end do
            end do
            end do
            close(8888)
        end if
    end if


end subroutine write_RT_Ac_file


subroutine incident()
    use em_field, only: calc_Ac_ext_t
    implicit none
    integer :: iix, iiy, iiz
    real(8), allocatable :: Ac(:, :)
    real(8), allocatable :: Ac_old(:, :)
    integer :: fh_bin

    fw%vec_Ac%v = 0d0
    fw%vec_Ac_old%v = 0d0

    ! x-directed incident
    allocate(Ac_inc(1:3, -1:itotNtime))

    call calc_Ac_ext_t(-(fs%mg%is(1)-0.5d0)*fs%hgs(1) / cspeed_au, fw%dt, &
        & -1, itotNtime, Ac_inc)

    if (yn_restart == 'y') then

        if (comm_is_root(ms%id_ms_world)) then
            fh_bin = get_filehandle()

            open(fh_bin,file=trim(ms%directory_read_data) // 'vec_Ac.bin',form='unformatted')
            read(fh_bin) fw%vec_Ac%v
            close(fh_bin)

            open(fh_bin,file=trim(ms%directory_read_data) // 'vec_Ac_old.bin',form='unformatted')
            read(fh_bin) fw%vec_Ac_old%v
            close(fh_bin)

            open(fh_bin,file=trim(ms%directory_read_data) // 'vec_j_em.bin',form='unformatted')
            read(fh_bin) fw%vec_j_em%v
            close(fh_bin)
        end if

        call comm_bcast(fw%vec_Ac%v, ms%icomm_ms_world)
        call comm_bcast(fw%vec_Ac_old%v, ms%icomm_ms_world)

    else

        allocate(Ac(1:3, fs%mg%is_overlap(1):0))
        allocate(Ac_old(1:3, fs%mg%is_overlap(1):0))    

        call calc_Ac_ext_t(-fw%dt*1, -fs%hgs(1) / cspeed_au, &
            & fs%mg%is_overlap(1), 0, Ac)
        call calc_Ac_ext_t(-fw%dt*2, -fs%hgs(1) / cspeed_au, &
            & fs%mg%is_overlap(1), 0, Ac_old)

        do iiz = fs%mg%is_overlap(3), fs%mg%ie_overlap(3)
            do iiy = fs%mg%is_overlap(2), fs%mg%ie_overlap(2)
                fw%vec_Ac%v(1:3,  fs%mg%is_overlap(1):0, iiy, iiz) = Ac(1:3,  fs%mg%is_overlap(1):0)
                fw%vec_Ac_old%v(1:3,  fs%mg%is_overlap(1):0, iiy, iiz) = Ac_old(1:3,  fs%mg%is_overlap(1):0)
            end do
        end do

        ! Evolve the first (0-th timestep)
        fw%Ac_inc_old(1:3) = Ac_inc(1:3,-1)
        fw%Ac_inc(1:3) = Ac_inc(1:3,0)
        call weyl_calc(fs, fw)

        if (comm_is_root(ms%id_ms_world)) &
            call write_RT_Ac_file()

        deallocate(Ac, Ac_old)
    end if

    return
 end subroutine incident




function check_input_variables() result(r)
    implicit none
    logical :: r
    r = .true.
    if (nx_m < 1) then
        if (comm_is_root(nproc_id_global)) &
            & write(*, *) "ERROR! 'nx_m' must be larger than 1!"
        r = .false.
    end if
    if (ny_m /= 1) then
        if (comm_is_root(nproc_id_global)) &
            & write(*, *) "ERROR! 'ny_m' must be 1!"
        r = .false.
    end if
    if (nz_m /= 1) then
        if (comm_is_root(nproc_id_global)) &
            & write(*, *) "ERROR! 'nz_m' must be 1!"
        r = .false.
    end if
    if (nproc_size_global < nx_m * ny_m * nz_m) then
        if (comm_is_root(nproc_id_global)) &
            & write(*, *) "ERROR! MPI procs is too small!"
        r = .false.
    end if
    if (mod(nproc_size_global, nx_m * ny_m * nz_m) > 0) then
        if (comm_is_root(nproc_id_global)) &
            & write(*, *) "ERROR! MPI procs number of processes is an integer multiple of macropoints!"
        r = .false.
    end if
    if (hx_m < 1d-6 .and. dl_em(1) < 1d-6) then
        if (comm_is_root(nproc_id_global)) &
            & write(*, *) "ERROR! 'hx_m' or 'dl_em(1)' must be specified!"
        r = .false.
    end if
    if (hy_m < 1d-6 .and. dl_em(2) < 1d-6) then
        if (comm_is_root(nproc_id_global)) &
            & write(*, *) "ERROR! 'hy_m' or 'dl_em(2)' must be specified!"
        r = .false.
    end if
    if (hz_m < 1d-6 .and. dl_em(3) < 1d-6) then
        if (comm_is_root(nproc_id_global)) &
            & write(*, *) "ERROR! 'hz_m' or 'dl_em(3)' must be specified!"
        r = .false.
    end if
    if (dt < 1d-6) then
        if (comm_is_root(nproc_id_global)) &
            & write(*, *) "ERROR! 'dt' must be specified!"
        r = .false.
    end if
    if (dt_em > 1d-6) then
        if (comm_is_root(nproc_id_global)) &
            & write(*, *) "ERROR! 'dt_em' must not be specified!"
        r = .false.
    end if
    if (nxvacl_m < 1) then
        if (comm_is_root(nproc_id_global)) &
            & write(*, *) "ERROR! 'nxvacl_m' must not larger than 1!"
        r = .false.
    end if
    if (nxvacr_m < 1) then
        if (comm_is_root(nproc_id_global)) &
            & write(*, *) "ERROR! 'nxvacr_m' must not larger than 1!"
        r = .false.
    end if
    if (trim(boundary_em(1,1)) .ne. 'periodic' &
        & .and. trim(boundary_em(1,1)) .ne. 'pec' &
        & .and. trim(boundary_em(1,1)) .ne. 'abc') then
        if (comm_is_root(nproc_id_global)) &
            & write(*, *) "ERROR! 'boundary_em(1,1)' unknown boundary condition!"
        r = .false.
    end if
    if (trim(boundary_em(1,2)) .ne. 'periodic' &
        & .and. trim(boundary_em(1,2)) .ne. 'pec' &
        & .and. trim(boundary_em(1,2)) .ne. 'abc') then
        if (comm_is_root(nproc_id_global)) &
            & write(*, *) "ERROR! 'boundary_em(1,2)' unknown boundary condition!"
        r = .false.
    end if
    return
end function check_input_variables

end subroutine main_ms

