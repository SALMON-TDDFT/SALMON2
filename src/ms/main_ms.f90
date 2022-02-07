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
use communication, only: comm_is_root, comm_sync_all, comm_create_group_byid, &
                         comm_get_groupinfo, comm_sync_all, comm_summation, &
                         comm_bcast
use salmon_xc, only: finalize_xc
use timer
use write_sub, only: write_response_0d,write_response_3d,write_pulse_0d,write_pulse_3d
use initialization_rt_sub
use checkpoint_restart_sub
use fdtd_weyl, only: ls_fdtd_weyl, weyl_init, weyl_calc, weyl_finalize
use parallelization, only: nproc_id_global, nproc_size_global, nproc_group_global, &
                           adjust_elapse_time
use filesystem, only: create_directory, get_filehandle
use phys_constants, only: cspeed_au
use em_field, only: calc_Ac_ext
use input_checker_ms, only: check_input_variables_ms

use inputoutput, only: t_unit_ac, t_unit_elec, t_unit_current, t_unit_time, t_unit_length, t_unit_energy
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
type(s_scalar) :: Vpsl
type(s_scalar) :: rho,rho_jm,Vh,Vh_stock1,Vh_stock2,Vbox
type(s_scalar),allocatable :: rho_s(:),V_local(:),Vxc(:)
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


integer :: Mit, itt
integer :: nntime

integer :: i, ix, iy, iz

integer, allocatable :: iranklists(:)

real(8), allocatable :: Ac_inc(:, :), rNe(:)

integer :: nmacro_mygroup, isize_mygroup

logical :: is_checkpoint_iter, is_shutdown_time
! Only for 1D calculation outputs:
integer :: fh_wave


! character(256) :: file_debug_log
!! Open logfile for debugging
! write(file_debug_log, "('ms_debug', i3.3, '.log')") nproc_id_global
! open(unit=9999, file=file_debug_log)
! write(9999, *) 'logging start'; flush(9999)

if (.not. check_input_variables_ms()) return

call timer_begin(LOG_TOTAL)
! Initialization

call initialization_ms()

call comm_sync_all
call timer_enable_sub
call timer_begin(LOG_RT_ITERATION)

call print_header()

TE : do itt=Mit+1,nt
    call time_evolution_step_ms()
    

    is_checkpoint_iter = (checkpoint_interval >= 1) .and. (mod(itt,checkpoint_interval) == 0)
    is_shutdown_time   = (time_shutdown > 0d0) .and. (adjust_elapse_time(timer_now(LOG_TOTAL)) > time_shutdown)

    if (is_checkpoint_iter .or. is_shutdown_time) then
      if (is_shutdown_time .and. comm_is_root(info%id_rko)) then
        print *, 'shutdown the calculation, iter =', itt
      end if

      call checkpoint_ms()

      if (is_shutdown_time) then
        exit TE
      end if
    end if
end do TE

call timer_end(LOG_RT_ITERATION)
call timer_disable_sub

call timer_end(LOG_TOTAL)

if(write_rt_wfn_k=='y')then
    call checkpoint_ms(ofl%dir_out_restart)
end if

call finalize_xc(xc_func)

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
    integer :: iimacro_s, iimacro_e, iimacro

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
    fw%fdtddim = trim(fdtddim)

    ! Experimental implementation for oblique incidence
    fw%theta = theta_oblique_deg / 180d0 * pi
    fw%nsmooth = nsmooth_oblique
    
    ! For compatibility with previous versions 
    ! (nxvacl_m and nxvacr_m will be removed in future)
    if ((nxvac1_m == 0) .and. (nxvacl_m /= 0)) nxvac1_m = abs(nxvacl_m) + 1
    if ((nxvac2_m == 0) .and. (nxvacr_m /= 0)) nxvac2_m = abs(nxvacr_m)

    fs%mg%is(1) = 1 - abs(nxvac1_m)
    fs%mg%ie(1) = nx_m + abs(nxvac2_m)
    fs%mg%is(2) = 1 - abs(nyvac1_m)
    fs%mg%ie(2) = ny_m + abs(nyvac2_m)
    fs%mg%is(3) = 1 - abs(nzvac1_m)
    fs%mg%ie(3) = nz_m + abs(nzvac2_m)
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

    ! incident field
    call Weyl_init(fs, fw)

    allocate(rNe(1:ms%nmacro))
    allocate(fs%imedia(fs%mg%is_array(1):fs%mg%ie_array(1), &
                       fs%mg%is_array(2):fs%mg%ie_array(2), &
                       fs%mg%is_array(3):fs%mg%ie_array(3)))
    fs%imedia(:,:,:) = 0

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

    call comm_sync_all(ms%icomm_ms_world)

    do iimacro_s = 1, ms%nmacro, nmacro_chunk
        iimacro_e = min(iimacro_s + nmacro_chunk - 1, ms%nmacro)
        if (comm_is_root(ms%id_ms_world)) &
            write(*, '(a, i6, "-", i6)') "Initializing macropoint:", iimacro_s, iimacro_e
        
        do iimacro = iimacro_s, iimacro_e

            if (macropoint_in_mygroup(iimacro)) then
                ! Override global variables
                base_directory = trim(base_directory_macro(iimacro))
                nproc_group_global = ms%icomm_macropoint
                nproc_id_global = ms%id_macropoint
                nproc_size_global = ms%isize_macropoint
                quiet = (1 < iimacro)

                if (yn_restart == 'y') then
                    directory_read_data = trim(directory_read_data_macro(trim(ms%directory_read_data), iimacro))
                end if
                
                ! Initializa TDKS system
                call initialization_rt( Mit, system, energy, ewald, rt, md, &
                                        singlescale,  &
                                        stencil, fg, poisson,  &
                                        lg, mg,   &
                                        info,  &
                                        xc_func, ofl,  &
                                        srg, srg_scalar,  &
                                        spsi_in, spsi_out, tpsi, rho, rho_jm, rho_s,  &
                                        V_local, Vbox, Vh, Vh_stock1, Vh_stock2, Vxc, Vpsl,&
                                        pp, ppg, ppn )

                ! Override global variables (restore)
                base_directory = trim(ms%base_directory)
                nproc_group_global = ms%icomm_ms_world
                nproc_id_global = ms%id_ms_world
                nproc_size_global = ms%isize_ms_world        
                directory_read_data = trim(ms%directory_read_data)
                quiet = .false.
            end if
        end do
        call comm_sync_all()
    end do

    itt = mit
    allocate(Ac_inc(1:3, -1:nt+2))
    call incident()

    ! Experimental implementation
    if (comm_is_root(ms%id_ms_world)) then
        if (trim(fdtddim) == '1d' .or. trim(fdtddim) == '1D') then
            call open_wave_data_file()
        end if
    end if

    call calc_weight(fs%mg%is(1), fs%mg%ie(1), &
        fw%weight(fs%mg%is(1):fs%mg%ie(1)), &
        nx_m, nsmooth_oblique)

        
    


    return    
end subroutine initialization_ms



subroutine print_header()
    implicit none
    if (comm_is_root(ms%id_ms_world)) then
        write(*,*)
        write(*,'(a7,a7,a9,a33,a15,a11)') &
            & "Step", "Macro", "Time", "Current","Electrons", "Eabs/cell"
        write(*,'(7x,7x,a9,a33,a15,a11)') &
            & trim(t_unit_time%name), &
            & trim(t_unit_current%name), &
            & "", &
            & trim(t_unit_energy%name)
        write(*,'("#",9("----------"))')
    endif
end subroutine print_header


subroutine print_linelog()
    implicit none
    integer :: iimacro, iix, iiy, iiz
    if (comm_is_root(ms%id_ms_world)) then
        do iimacro = 1, ms%nmacro
            iix = ms%ixyz_tbl(1, iimacro)
            iiy = ms%ixyz_tbl(2, iimacro)
            iiz = ms%ixyz_tbl(3, iimacro)
            write(*,'(i7,i7,f9.3,3(es11.2e3),f15.8,es11.2e3)') &
                & itt, iimacro, &
                & itt * dt * t_unit_time%conv, &
                & fw%vec_j_em_new%v(1:3, iix, iiy, iiz) * t_unit_current%conv, &
                & rNe(iimacro), &
                & fw%edensity_absorb%f(iix, iiy, iiz) * system%det_a * t_unit_energy%conv
        end do
    end if
end subroutine





subroutine time_evolution_step_ms
    implicit none
    integer :: iimacro, iix, iiy, iiz
    real(8) :: curr_tmp(3, ms%nmacro), curr(3, ms%nmacro), rNe_tmp(ms%nmacro)

    real(8) :: c1, c2, c3

    ! ----------------------------------------
    ! Time Evolution of FDTD System
    ! ----------------------------------------
    fw%Ac_inc_new(:) = Ac_inc(:, itt)
    fw%Ac_inc(:) = Ac_inc(:, itt-1)
    call weyl_calc(fs, fw)
    if (mod(itt, out_ms_step) == 0) call write_RT_Ac_file()
    
    if(yn_ms_ld .ne. 'y') then

    ! ----------------------------------------
    ! Time Evolution of TDDFT System
    ! ----------------------------------------
    ! Override Global Variables
    nproc_group_global = ms%icomm_macropoint
    nproc_id_global = ms%id_macropoint
    nproc_size_global = ms%isize_macropoint
    quiet = .true.

    curr_tmp(:,:) = 0d0
    rNe_tmp(:) = 0d0
    do iimacro = ms%imacro_mygroup_s, ms%imacro_mygroup_s
        iix = ms%ixyz_tbl(1, iimacro)
        iiy = ms%ixyz_tbl(2, iimacro)
        iiz = ms%ixyz_tbl(3, iimacro)
        rt%Ac_tot(1:3, itt-1) = fw%vec_Ac%v(1:3, iix, iiy, iiz)
        rt%Ac_tot(1:3, itt)   = fw%vec_Ac_new%v(1:3, iix, iiy, iiz)
        rt%Ac_ext(1:3, itt-1) = rt%Ac_tot(1:3, itt-1)
        rt%Ac_ext(1:3, itt)   = rt%Ac_tot(1:3, itt)

        if(mod(itt,2)==1)then
            call time_evolution_step(Mit,nt,itt,lg,mg,system,rt,info,stencil,xc_func &
            & ,srg,srg_scalar,pp,ppg,ppn,spsi_in,spsi_out,tpsi,rho,rho_jm,rho_s,V_local,Vbox,Vh,Vh_stock1,Vh_stock2,Vxc &
            & ,Vpsl,fg,energy,ewald,md,ofl,poisson,singlescale)
        else
            call time_evolution_step(Mit,nt,itt,lg,mg,system,rt,info,stencil,xc_func &
            & ,srg,srg_scalar,pp,ppg,ppn,spsi_out,spsi_in,tpsi,rho,rho_jm,rho_s,V_local,Vbox,Vh,Vh_stock1,Vh_stock2,Vxc &
            & ,Vpsl,fg,energy,ewald,md,ofl,poisson,singlescale)
        end if
            
        if (comm_is_root(ms%id_macropoint)) then
            curr_tmp(1:3, iimacro) = rt%curr(1:3, itt)
            rNe_tmp(iimacro) = rt%rIe(itt)
        endif
    end do
    ! Override Global Variables (Repair)
    nproc_group_global = ms%icomm_ms_world
    nproc_id_global = ms%id_ms_world
    nproc_size_global = ms%isize_ms_world
    quiet = .false.
    
    call comm_summation(curr_tmp, curr, 3 * ms%nmacro, ms%icomm_ms_world)
    call comm_summation(rNe_tmp, rNe, ms%nmacro, ms%icomm_ms_world)
    

    do iimacro = 1, ms%nmacro
        iix = ms%ixyz_tbl(1, iimacro)
        iiy = ms%ixyz_tbl(2, iimacro)
        iiz = ms%ixyz_tbl(3, iimacro)
        fw%vec_j_em_new%v(1:3, iix, iiy, iiz) = -1.0d0 * curr(1:3, iimacro)
    end do

    if (mod(itt, out_ms_step) == 0) call print_linelog()

    else
    ! Lorentz-Drude oscilator
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    c1 = (1.0d0 - 0.5d0 * (ms_omega_ld * dt) ** 2) &
        & / (1.0d0 + 0.5d0 * ms_gamma_ld * dt)
    c2 = (1.0d0 - 0.5d0 * (ms_gamma_ld * dt)) &
        & / (1.0d0 + 0.5d0 * (ms_gamma_ld * dt))
    c3 = ms_alpha_ld / (1.0d0 + 0.5d0 * (ms_gamma_ld * dt))

    do iimacro = 1, ms%nmacro
        iix = ms%ixyz_tbl(1, iimacro)
        iiy = ms%ixyz_tbl(2, iimacro)
        iiz = ms%ixyz_tbl(3, iimacro)

        fw%vec_j_em_new%v(1:3, iix, iiy, iiz) &
            & = 2.0d0 * c1 * fw%vec_j_em%v(1:3, iix, iiy, iiz) &
            & - c2 * fw%vec_j_em_old%v(1:3, iix, iiy, iiz) &
            & - c3 * (fw%vec_Ac_new%v(1:3, iix, iiy, iiz) &
                & - 2.0d0 * fw%vec_Ac%v(1:3, iix, iiy, iiz) &
                & + fw%vec_Ac_old%v(1:3, iix, iiy, iiz) &
                & )
    end do
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (comm_is_root(ms%id_ms_world)) then
        if (mod(itt, out_ms_step) == 0) write(*, *) "# Lorentz-Drude test", itt
    end if

    end if



    ! Experimental implementation
    if (comm_is_root(ms%id_ms_world)) then
        if (mod(itt, 2) == 0) then
            if (trim(fdtddim) == '1d' .or. trim(fdtddim) == '1D') then
                call write_wave_data_file()
            end if
        end if
    end if

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

        if (comm_is_root(ms%id_ms_world)) then
            call create_directory(trim(idir))
            write(*, *) "Checkpointing macropoint:", i
        end if

        call comm_sync_all()

        if (yn_ms_ld .ne. 'y') then
        if (macropoint_in_mygroup(i)) then
            if (mod(itt,2)==1) then
                call checkpoint_rt(lg,mg,system,info,spsi_out,itt,rt,Vh_stock1,Vh_stock2,singlescale,idir)
            else
                call checkpoint_rt(lg,mg,system,info,spsi_in, itt,rt,Vh_stock1,Vh_stock2,singlescale,idir)
            endif
        end if
        end if
    end do

    if (comm_is_root(ms%id_ms_world)) then
        fh_bin = get_filehandle()

        open(fh_bin,file=trim(idir) // '../vec_Ac_new.bin',form='unformatted')
        write(fh_bin) fw%vec_Ac_new%v
        close(fh_bin)

        open(fh_bin,file=trim(idir) // '../vec_Ac.bin',form='unformatted')
        write(fh_bin) fw%vec_Ac%v
        close(fh_bin)

        open(fh_bin,file=trim(idir) // '../vec_Ac_old.bin',form='unformatted')
        write(fh_bin) fw%vec_Ac_old%v
        close(fh_bin)

        open(fh_bin,file=trim(idir) // '../vec_j_em_old.bin',form='unformatted')
        write(fh_bin) fw%vec_j_em_old%v
        close(fh_bin)

        open(fh_bin,file=trim(idir) // '../vec_j_em.bin',form='unformatted')
        write(fh_bin) fw%vec_j_em%v
        close(fh_bin)

        open(fh_bin,file=trim(idir) // '../vec_j_em_new.bin',form='unformatted')
        write(fh_bin) fw%vec_j_em_new%v
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
    implicit none
    integer ::  iix, iiy, iiz
    integer :: fh_ac_data   
    character(256) :: file_ac_data
    
    if (comm_is_root(ms%id_ms_world)) then
        fh_ac_data = get_filehandle()
        write(file_ac_data, '(a,a,"_Ac_",i6.6,".data")') trim(ms%base_directory_RT_Ac), trim(sysname), itt
        open(fh_ac_data, file=trim(file_ac_data))
        write(fh_ac_data, '(a)') "# Multiscale TDDFT calculation"
        write(fh_ac_data, '(a)') "# IX, IY, IZ: FDTD Grid index"
        write(fh_ac_data, '(a)') "# x, y, z: Coordinates"
        write(fh_ac_data, '(a)') "# Ac: Vector potential field"
        write(fh_ac_data, '(a)') "# E: Electric field"
        write(fh_ac_data, '(a)') "# J_em: Electromagnetic current density"
        write(fh_ac_data, '("#",99(1X,I0,":",A,"[",A,"]"))') &
            & 1, "IX", "none", &
            & 2, "IY", "none", &
            & 3, "IZ", "none", &
            & 4, "Ac_x", trim(t_unit_ac%name), &
            & 5, "Ac_y", trim(t_unit_ac%name), &
            & 6, "Ac_z", trim(t_unit_ac%name), &
            & 7, "E_x", trim(t_unit_elec%name), &
            & 8, "E_y", trim(t_unit_elec%name), &
            & 9, "E_z", trim(t_unit_elec%name), &
            & 10, "B_x", "a.u.", &
            & 11, "B_y", "a.u.", &
            & 12, "B_z", "a.u.", &
            & 13, "Jem_x", trim(t_unit_current%name), &
            & 14, "Jem_y", trim(t_unit_current%name), &
            & 15, "Jem_z", trim(t_unit_current%name), &
            & 16, "E_em", trim(t_unit_energy%name) // "/vol", &
            & 17, "E_abs", trim(t_unit_energy%name) //  "/vol"

        do iiz = fs%mg%is(3), fs%mg%ie(3)
        do iiy = fs%mg%is(2), fs%mg%ie(2)
        do iix = fs%mg%is(1), fs%mg%ie(1)
            write(fh_ac_data, '(3(i6, 1x), 14(e23.15e3, 1x))')  &
                & iix, & 
                & iiy, & 
                & iiz, & 
                & fw%vec_Ac_new%v(1, iix, iiy, iiz) * t_unit_ac%conv, &
                & fw%vec_Ac_new%v(2, iix, iiy, iiz) * t_unit_ac%conv, &
                & fw%vec_Ac_new%v(3, iix, iiy, iiz) * t_unit_ac%conv, &
                & fw%vec_E%v(1, iix, iiy, iiz) * t_unit_elec%conv, &
                & fw%vec_E%v(2, iix, iiy, iiz) * t_unit_elec%conv, &
                & fw%vec_E%v(3, iix, iiy, iiz) * t_unit_elec%conv, &
                & fw%vec_H%v(1, iix, iiy, iiz), &
                & fw%vec_H%v(2, iix, iiy, iiz), &
                & fw%vec_H%v(3, iix, iiy, iiz), &
                & fw%vec_j_em_new%v(1, iix, iiy, iiz) * t_unit_current%conv, &
                & fw%vec_j_em_new%v(2, iix, iiy, iiz) * t_unit_current%conv, &
                & fw%vec_j_em_new%v(3, iix, iiy, iiz) * t_unit_current%conv, &
                & fw%edensity_emfield%f(iix, iiy, iiz) * t_unit_energy%conv / t_unit_length%conv ** 3, &
                & fw%edensity_absorb%f(iix, iiy, iiz) * t_unit_energy%conv / t_unit_length%conv ** 3
        end do
        end do
        end do
        close(fh_ac_data)
    end if

end subroutine write_RT_Ac_file


subroutine incident()
    use em_field, only: calc_Ac_ext_t
    implicit none
    integer :: iix, iiy, iiz
    real(8), allocatable :: Ac(:, :), Ac_new(:, :)
    real(8), allocatable :: Ac_old(:, :)
    integer :: fh_bin

    fw%vec_Ac_new%v = 0d0
    fw%vec_Ac%v = 0d0
    fw%vec_Ac_old%v = 0d0


    call calc_Ac_ext_t(-(fs%mg%is(1)-0.5d0)*fs%hgs(1) / cspeed_au, fw%dt, &
        & -1, nt+2, Ac_inc)

    if (yn_restart == 'y') then

        itt = Mit

        if (comm_is_root(ms%id_ms_world)) then
            fh_bin = get_filehandle()

            open(fh_bin,file=trim(ms%directory_read_data) // 'vec_Ac_new.bin',form='unformatted')
            read(fh_bin) fw%vec_Ac_new%v
            close(fh_bin)

            open(fh_bin,file=trim(ms%directory_read_data) // 'vec_Ac.bin',form='unformatted')
            read(fh_bin) fw%vec_Ac%v
            close(fh_bin)

            open(fh_bin,file=trim(ms%directory_read_data) // 'vec_Ac_old.bin',form='unformatted')
            read(fh_bin) fw%vec_Ac_old%v
            close(fh_bin)

            open(fh_bin,file=trim(ms%directory_read_data) // 'vec_j_em_new.bin',form='unformatted')
            read(fh_bin) fw%vec_j_em_new%v
            close(fh_bin)
            
            open(fh_bin,file=trim(ms%directory_read_data) // 'vec_j_em.bin',form='unformatted')
            read(fh_bin) fw%vec_j_em%v
            close(fh_bin)

            open(fh_bin,file=trim(ms%directory_read_data) // 'vec_j_em_old.bin',form='unformatted')
            read(fh_bin) fw%vec_j_em_old%v
            close(fh_bin)
        end if

        call comm_bcast(fw%vec_Ac_new%v, ms%icomm_ms_world)
        call comm_bcast(fw%vec_Ac_old%v, ms%icomm_ms_world)
        call comm_bcast(fw%vec_Ac%v, ms%icomm_ms_world)
        call comm_bcast(fw%vec_j_em_new%v, ms%icomm_ms_world)
        call comm_bcast(fw%vec_j_em_old%v, ms%icomm_ms_world)
        call comm_bcast(fw%vec_j_em%v, ms%icomm_ms_world)

    else
        itt = -1  ! Prepare itt=-1 Ac_vec field profile
        allocate(Ac_new(1:3, fs%mg%is_overlap(1):0))
        allocate(Ac(1:3, fs%mg%is_overlap(1):0))
        allocate(Ac_old(1:3, fs%mg%is_overlap(1):0)) 
        call calc_Ac_ext_t((itt)*fw%dt, -fs%hgs(1) / cspeed_au, &
            & fs%mg%is_overlap(1), 0, Ac_new)   ! Ac_new = Ac(itt) with itt=-1
        call calc_Ac_ext_t((itt-1)*fw%dt, -fs%hgs(1) / cspeed_au, &
            & fs%mg%is_overlap(1), 0, Ac)       ! Ac = Ac(itt-1) with itt=-1
        call calc_Ac_ext_t((itt-2)*fw%dt, -fs%hgs(1) / cspeed_au, &
            & fs%mg%is_overlap(1), 0, Ac_old)   ! Ac = Ac(itt-2) with itt=-1
        do iiz = fs%mg%is_overlap(3), fs%mg%ie_overlap(3)
            do iiy = fs%mg%is_overlap(2), fs%mg%ie_overlap(2)
                fw%vec_Ac_new%v(1:3, fs%mg%is_overlap(1):0, iiy, iiz) = Ac_new(1:3, fs%mg%is_overlap(1):0)
                fw%vec_Ac%v(1:3, fs%mg%is_overlap(1):0, iiy, iiz) = Ac(1:3, fs%mg%is_overlap(1):0)
                fw%vec_Ac_old%v(1:3, fs%mg%is_overlap(1):0, iiy, iiz) = Ac_old(1:3, fs%mg%is_overlap(1):0)
            end do
        end do
        deallocate(Ac_new, Ac, Ac_old)

        itt = 0 ! Update to itt=0 Ac_vec field prodile
        fw%Ac_inc(:) = Ac_inc(:, itt-1)
        fw%Ac_inc_new(:) = Ac_inc(:, itt)
        call weyl_calc(fs, fw)

    end if

    ! Write out initial field profile:
    call write_RT_Ac_file()

    return
 end subroutine incident


 ! Experimetal Implementation of Incident/Reflection/Transmit field output
 subroutine open_wave_data_file()
    use filesystem, only: open_filehandle
    implicit none
    fh_wave = open_filehandle(trim(ms%base_directory) // trim(sysname) // "_wave.data")
    write(fh_wave, '(a)') "# 1D multiscale calculation:"
    write(fh_wave, '(a)') "# E_inc: E-field amplitude of incident wave"
    write(fh_wave, '(a)') "# E_ref: E-field amplitude of reflected wave"
    write(fh_wave, '(a)') "# E_tra: E-field amplitude of transmitted wave"
    write(fh_wave, '("#",99(1X,I0,":",A,"[",A,"]"))') &
        & 1, "Time", trim(t_unit_time%name), &
        & 2, "E_inc_x", trim(t_unit_elec%name), &
        & 3, "E_inc_y", trim(t_unit_elec%name), &
        & 4, "E_inc_z", trim(t_unit_elec%name), &
        & 5, "E_ref_x", trim(t_unit_elec%name), &
        & 6, "E_ref_y", trim(t_unit_elec%name), &
        & 7, "E_ref_z", trim(t_unit_elec%name), &
        & 8, "E_tra_x", trim(t_unit_elec%name), &
        & 9, "E_tra_y", trim(t_unit_elec%name), &
        & 10, "E_tra_z", trim(t_unit_elec%name)
end subroutine open_wave_data_file

 ! Experimetal Implementation of Incident/Reflection/Transmit field output
subroutine write_wave_data_file()
    implicit none
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

    write(fh_wave, '(99(e23.15e3, 1x))')  &
        & itt * dt * t_unit_time%conv, &
        & e_inc(1) * t_unit_elec%conv, &
        & e_inc(2) * t_unit_elec%conv, &
        & e_inc(3) * t_unit_elec%conv, &
        & e_ref(1) * t_unit_elec%conv, &
        & e_ref(2) * t_unit_elec%conv, &
        & e_ref(3) * t_unit_elec%conv, &
        & e_tra(1) * t_unit_elec%conv, &
        & e_tra(2) * t_unit_elec%conv, &
        & e_tra(3) * t_unit_elec%conv
    return
end subroutine write_wave_data_file
    
    
 ! Experimetal Implementation of Incident/Reflection/Transmit field output
subroutine close_wave_data_file()
    implicit none
    close(fh_wave)
    return
end subroutine close_wave_data_file
    



subroutine calc_weight(iz_sta, iz_end, w, nz, ns)
    implicit none
    integer, intent(in) :: iz_sta
    integer, intent(in) :: iz_end
    real(8), intent(inout) :: w(iz_sta:iz_end)
    integer, intent(in) :: nz
    integer, intent(in) :: ns
    integer :: iz
    w = 0.0d0
    do iz = 1, nz
        if (iz < ns) then
            w(iz) = weight((dble(iz) - dble(0.5)) / dble(ns))
        else if (nz - ns - 1 < iz) then
            w(iz) = weight((dble(nz) - dble(iz) - dble(0.5)) / dble(ns))
        else
            w(iz) = 1.0d0
        end if
    end do
contains
    
    real(8) function weight(tt)
        implicit none
        real(8) :: tt
        weight = -2.0d0 * tt ** 3 + 3.0d0 * tt ** 2
        return
    end function weight
end subroutine


end subroutine main_ms

