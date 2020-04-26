!
!  Copyright 2019 SALMON developers
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
use communication, only: comm_is_root, comm_sync_all, comm_create_group_byid, comm_get_groupinfo, comm_sync_all, comm_summation
use salmon_xc, only: finalize_xc
use timer
use write_sub, only: write_response_0d,write_response_3d,write_pulse_0d,write_pulse_3d
use initialization_rt_sub
use checkpoint_restart_sub
use fdtd_weyl, only: ls_fdtd_weyl, weyl_init, weyl_calc, weyl_finalize
use parallelization, only: nproc_id_global, nproc_size_global, nproc_group_global
use filesystem, only: create_directory
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

character(256) :: file_debug_log
integer :: nmacro_mygroup, isize_mygroup

if (.not. check_input_variables()) return

! Open logfile for debugging
write(file_debug_log, "('ms_debug', i3.3, '.log')") nproc_id_global
open(unit=9999, file=file_debug_log)
write(9999, *) 'logging start'; flush(9999)

call timer_begin(LOG_TOTAL)
! Initialization
call initialization_ms()
! if (comm_is_root(ms%id_ms_world)) call print_header()


call comm_sync_all
call timer_enable_sub
call timer_begin(LOG_RT_ITERATION)

TE : do itt=Mit+1,itotNtime

    write(9999, *) 'Step', itt; flush(9999)

    write(9999, *) "StartTime evolution"; flush(9999)
    call time_evolution_step_ms_macro()
    write(9999, *) "EndTime evolution"; flush(9999)
    call write_RT_Ac_file()


!   if((checkpoint_interval >= 1) .and. (mod(itt,checkpoint_interval) == 0)) then
!     call timer_begin(LOG_CHECKPOINT_SYNC)
!     call timer_begin(LOG_CHECKPOINT_SELF)
!     if (mod(itt,2)==1) then
!       call checkpoint_rt(lg,mg,system,info,spsi_out,itt,sVh_stock1,sVh_stock2,singlescale)
!     else
!       call checkpoint_rt(lg,mg,system,info,spsi_in, itt,sVh_stock1,sVh_stock2,singlescale)
!     endif
!     call timer_end(LOG_CHECKPOINT_SELF)
!     call comm_sync_all
!     call timer_end(LOG_CHECKPOINT_SYNC)
!   endif


end do TE

write(9999, *) 'Loop Complete'
flush(9999)

call timer_end(LOG_RT_ITERATION)
call timer_disable_sub

! #ifdef __FUJITSU
! call fapp_stop('time_evol',1,0) ! performance profiling
! #endif

! close(030) ! laser


! !--------------------------------- end of time-evolution


! !------------ Writing part -----------

! ! write RT: 
! call timer_begin(LOG_WRITE_RT_RESULTS)

! !
! select case(iperiodic)
! case(0)
!   if(theory=="tddft_response")then
!     call write_response_0d(ofl,rt)
!   else
!     call write_pulse_0d(ofl,rt)
!   end if
! case(3)
!   if(theory=="tddft_response")then
!     call write_response_3d(ofl,rt)
!   else
!     call write_pulse_3d(ofl,rt)
!   end if
! end select

! call timer_end(LOG_WRITE_RT_RESULTS)
call timer_end(LOG_TOTAL)

! if(write_rt_wfn_k=='y')then
!   call checkpoint_rt(lg,mg,system,info,spsi_out,Mit,sVh_stock1,sVh_stock2,singlescale,ofl%dir_out_restart)
! end if

call finalize_xc(xc_func)

close(9999)

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
    ! character(256) :: tmp
    character(256) :: r
    ! write(9999, *) trim(ms%base_directory)
    ! write(9999, *) trim(sysname)
    ! write(9999, *) imacro
    ! FLUSH(9999)
    ! write(9999, '(a, a, a, i6.6, a)') trim(ms%base_directory), trim(sysname), '_m/', imacro, '/'
    ! flush(9999)
    ! write(9999, *) trim(tmp)
    ! flush(9999)
    write(r, '(a, a, a, i6.6, a)') trim(ms%base_directory), trim(sysname), '_m/', imacro, '/'
    return
end function base_directory_macro


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
    return
end function check_input_variables




subroutine initialization_ms()
    implicit none

    ! Store global information
    ms%nmacro = nx_m * ny_m * nz_m
    ms%base_directory = trim(base_directory)
    ms%base_directory_RT_Ac = trim(ms%base_directory) // trim(sysname) // "_RT_Ac/"
    ms%icomm_ms_world = nproc_group_global
    ms%isize_ms_world = nproc_size_global
    ms%id_ms_world = nproc_id_global
    

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

    ! allocate(iranklists())
    write(9999, *) 'nmacro_mygroup:', nmacro_mygroup
    write(9999, *) 'isize_mygroup:', isize_mygroup
    write(9999, *) 'ms%imacro_mygroup_s:', ms%imacro_mygroup_s
    write(9999, *) 'ms%imacro_mygroup_e:', ms%imacro_mygroup_e
    write(9999, *) 'ms%id_mygroup_s:', ms%id_mygroup_s
    write(9999, *) 'ms%id_mygroup_e:', ms%id_mygroup_e
    write(9999, *) 'ms%base_directory:', trim(ms%base_directory)
    flush(9999)
    
    if (ms%nmacro < 1) stop "Invalid macropoint number"
    
    allocate(iranklists(isize_mygroup))
    do i = 1, isize_mygroup
        iranklists(i) = ms%id_mygroup_s + (i - 1)
    end do
    
    write(9999, *) 'iranklists:', iranklists(:)
    flush(9999)
    
    ms%icomm_macropoint = comm_create_group_byid(ms%icomm_ms_world, iranklists)
    call comm_get_groupinfo(ms%icomm_macropoint, ms%id_macropoint, ms%isize_macropoint)
    
    write(9999, *) 'ms%icomm_ms_world:', ms%icomm_ms_world
    write(9999, *) 'ms%id_ms_world:', ms%id_ms_world
    write(9999, *) 'ms%isize_ms_world:', ms%isize_ms_world
    write(9999, *) 'ms%icomm_macropoint:', ms%icomm_macropoint
    write(9999, *) 'ms%id_macropoint:', ms%id_macropoint
    write(9999, *) 'ms%isize_macropoint:', ms%isize_macropoint
    flush(9999)
    
    fs%mg%ndir = 3
    fs%mg%nd = 1
    
    fs%hgs(1:3) = dl_em(1:3)
    if (0d0 < hx_m) fs%hgs(1) = hx_m
    if (0d0 < hy_m) fs%hgs(2) = hy_m
    if (0d0 < hz_m) fs%hgs(3) = hz_m
        
    write(9999, *) 'fs%mg%ndir:', fs%mg%ndir
    write(9999, *) 'fs%mg%nd:', fs%mg%nd
    write(9999, *) 'fs%hgs(1:3):', fs%hgs(1:3)
    flush(9999)

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
    fs%a_bc(1,1:2) = 'periodic'
    fs%a_bc(2,1:2) = 'periodic'
    fs%a_bc(3,1:2) = 'periodic'
    fs%imedia(:,:,:) = 0

    write(9999, *) 'Initialization FDTD Weyl start';flush(9999)
    call Weyl_init(fs, fw)
    write(9999, *) 'Initialization FDTD Weyl end';flush(9999)

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

    write(9999, *) trim(ms%base_directory_RT_Ac); flush(9999)
    if (comm_is_root(ms%id_ms_world)) &
        & call create_directory(trim(ms%base_directory_RT_Ac))

    write(9999, *) 'Macrpoint initialization start';flush(9999)
    do i = 1, ms%nmacro
        if (comm_is_root(ms%id_ms_world)) then
            write(9999, *) trim(base_directory_macro(i)); flush(9999)
            call create_directory(trim(base_directory_macro(i)))
            write(*, '(a)') "################################"
            write(*, '(a, i6)') "# Initialization of macropoint #", i
        end if

        call comm_sync_all()

        if (macropoint_in_mygroup(i)) then
            ! Override global variables
            base_directory = trim(base_directory_macro(i))
            nproc_group_global = ms%icomm_macropoint
            nproc_id_global = ms%id_macropoint
            nproc_size_global = ms%isize_macropoint
            
            write(9999, *) "trim(base_directory)", trim(base_directory); flush(9999)
            write(9999, *) "nproc_group_global", nproc_group_global; flush(9999)
            write(9999, *) "nproc_id_global", nproc_id_global; flush(9999)
            write(9999, *) "nproc_size_global", nproc_size_global; flush(9999)
                    
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
        end if
    end do
    write(9999, *) 'Macrpoint initialization end';flush(9999)

    write(9999, *) 'Incident field setup start';flush(9999)
    call incident()
    write(9999, *) 'Incident field setup end';flush(9999)
    call write_RT_Ac_file()

    write(9999, *) 'Initialization Complete'
    flush(9999)
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


subroutine time_evolution_step_ms_macro
    implicit none
    integer :: ii, iimacro, iix, iiy, iiz
    real(8) :: curr_tmp(3, ms%nmacro)

    ! Override Global Variables
    nproc_group_global = ms%icomm_macropoint
    nproc_id_global = ms%id_macropoint
    nproc_size_global = ms%isize_macropoint

    write(9999, *) "FDTD iteration start"; flush(9999)
    call weyl_calc(fs, fw)

    do ii = 1, ms%nmacro
        iix = ms%ixyz_tbl(1, ii)
        iiy = ms%ixyz_tbl(2, ii)
        iiz = ms%ixyz_tbl(3, ii)
        ms%vec_Ac(1:3, ii) = fw%vec_Ac%v(1:3, iix, iiy, iiz)
        ms%vec_Ac_old(1:3, ii) = fw%vec_Ac_old%v(1:3, iix, iiy, iiz)
    end do
    write(9999, *) "FDTD iteration end"; flush(9999)

    if (ms%imacro_mygroup_e - ms%imacro_mygroup_s + 1 > 1) then
        stop "ERROR! Unsupported paralization scheme!"
    else
        iimacro = ms%imacro_mygroup_s

        rt%Ac_ext(:, itt) = ms%vec_Ac(1:3, iimacro)
        rt%Ac_ext(:, itt+1) = ms%vec_Ac(1:3, iimacro) &
            & + (ms%vec_Ac(1:3, iimacro) -  ms%vec_Ac_old(1:3, iimacro))

        write(9999, *) "RT iteration start"; flush(9999)
        if(mod(itt,2)==1)then
            call time_evolution_step(Mit,itotNtime,itt,lg,mg,system,rt,info,stencil,xc_func &
            & ,srg,srg_scalar,pp,ppg,ppn,spsi_in,spsi_out,tpsi,srho,srho_s,V_local,Vbox,sVh,sVh_stock1,sVh_stock2,sVxc &
            & ,sVpsl,dmat,fg,energy,ewald,md,ofl,poisson,singlescale)
        else
            call time_evolution_step(Mit,itotNtime,itt,lg,mg,system,rt,info,stencil,xc_func &
            & ,srg,srg_scalar,pp,ppg,ppn,spsi_out,spsi_in,tpsi,srho,srho_s,V_local,Vbox,sVh,sVh_stock1,sVh_stock2,sVxc &
            & ,sVpsl,dmat,fg,energy,ewald,md,ofl,poisson,singlescale)
        end if
        write(9999, *) "RT iteration end"; flush(9999)

        write(9999, *) "Current summation start"; flush(9999)
        curr_tmp(:, :) = 0d0
        if (comm_is_root(ms%id_macropoint)) &
            & curr_tmp(1:3, iimacro) = rt%curr(1:3, itt)
        write(9999, *) "iimacro, curr_tmp, curr", iimacro, curr_tmp(3, iimacro), rt%curr(3, itt)
        call comm_summation(curr_tmp, ms%curr, 3 * ms%nmacro, ms%icomm_ms_world)
        write(9999, *) "Current summation end"; flush(9999)

    end if

    do ii = 1, ms%nmacro
        iix = ms%ixyz_tbl(1, ii)
        iiy = ms%ixyz_tbl(2, ii)
        iiz = ms%ixyz_tbl(3, ii)
        fw%vec_j_em%v(1:3, iix, iiy, iiz) = - ms%curr(1:3, ii)
        write(9999, *) "macromove", iix, iiy, iiz, ms%curr(1, ii), ms%curr(2, ii), ms%curr(3, ii)
    end do

    ! write(9999, *) "Current write start"; flush(9999)
    ! do ii = 1, ms%nmacro
    !     write(9999, "(4(e23.15e3, 1x), a, i3.3)") itt * dt, ms%curr(1, ii),  ms%curr(2, ii),  ms%curr(3, ii), " #CURRM", ii
    ! end do
    ! write(9999, *) "Current write end"; flush(9999)

    ! Override Global Variables (Repair)
    nproc_group_global = ms%icomm_ms_world
    nproc_id_global = ms%id_ms_world
    nproc_size_global = ms%isize_ms_world

    return
end subroutine time_evolution_step_ms_macro





subroutine write_RT_Ac_file()
    use inputoutput, only: t_unit_ac, t_unit_current, t_unit_elec, t_unit_length
    implicit none
    integer ::  iix, iiy, iiz
    character(256) :: filename

    write(9999, *) "Start write_RT_Ac_file"; flush(9999)

    if (ms%id_ms_world == 0) then
        if (mod(itt, 100) == 0) then
            write(filename, '(a, a, a, i6.6, a)') trim(ms%base_directory_RT_Ac), trim(sysname), "_Ac_",  itt, '.data'
            write(9999, *) trim(filename); flush(9999)
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

    write(9999, *) "End write_RT_Ac_file"; flush(9999)

end subroutine write_RT_Ac_file


subroutine incident()
    use em_field, only: calc_Ac_ext_t
    implicit none
    integer :: iix, iiy, iiz
    real(8), allocatable :: Ac(:, :)
    real(8), allocatable :: Ac_old(:, :)

    fw%vec_Ac%v = 0d0
    fw%vec_Ac_old%v = 0d0

    ! x-directed incident
    allocate(Ac(1:3, fs%mg%is_overlap(1):0))
    allocate(Ac_old(1:3, fs%mg%is_overlap(1):0))

    call calc_Ac_ext_t(0d0, -fs%hgs(1) / cspeed_au, &
        & fs%mg%is_overlap(1), 0, Ac)
    call calc_Ac_ext_t(-fw%dt, -fs%hgs(1) / cspeed_au, &
        & fs%mg%is_overlap(1), 0, Ac_old)

    do iiz = fs%mg%is_overlap(3), fs%mg%ie_overlap(3)
        do iiy = fs%mg%is_overlap(2), fs%mg%ie_overlap(2)
            fw%vec_Ac%v(1:3,  fs%mg%is_overlap(1):0, iiy, iiz) = Ac(1:3,  fs%mg%is_overlap(1):0)
            fw%vec_Ac_old%v(1:3,  fs%mg%is_overlap(1):0, iiy, iiz) = Ac_old(1:3,  fs%mg%is_overlap(1):0)
        end do
    end do

    deallocate(Ac, Ac_old)

    return
 end subroutine incident

end subroutine main_ms

