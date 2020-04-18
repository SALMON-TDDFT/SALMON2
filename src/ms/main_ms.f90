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
use inputoutput, only: nx_m, ny_m, nz_m
use communication, only: comm_is_root, comm_sync_all, comm_create_group_byid, comm_get_groupinfo, comm_sync_all
use salmon_xc, only: finalize_xc
use timer
use write_sub, only: write_response_0d,write_response_3d,write_pulse_0d,write_pulse_3d
use initialization_rt_sub
use checkpoint_restart_sub
use fdtd_weyl, only: ls_fdtd_weyl, weyl_init, weyl_calc, weyl_finalize
use parallelization, only: nproc_id_global, nproc_size_global, nproc_group_global

use filesystem, only: create_directory
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

integer :: i

integer, allocatable :: iranklists(:)

character(256) :: file_debug_log
integer :: nmacro_mygroup, isize_mygroup

call timer_begin(LOG_TOTAL)

! Open logfile for debugging
write(file_debug_log, "('ms_debug', i3.3, '.log')") nproc_id_global
open(unit=9999, file=file_debug_log)
write(9999, *) 'logging start'; flush(9999)

! Initialization
call initialization_ms()

call comm_sync_all
call timer_enable_sub
call timer_begin(LOG_RT_ITERATION)
TE : do itt=Mit+1,itotNtime

    write(9999, *) 'Step', itt; flush(9999)

    call time_evolution_step_ms_macro()

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
    write(r, '(a, a, a, i6.6, a)') trim(ms%base_directory), trim(sysname), '_m/', imacro, '/'
    ! write(9999, *) trim(tmp)
    ! flush(9999)
    return
end function base_directory_macro


subroutine initialization_ms()
    implicit none

    ms%base_directory = trim(base_directory)
    ms%nmacro = nx_m * ny_m * nz_m
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
    if (0 < hx_m) fs%hgs(1) = hx_m
    if (0 < hy_m) fs%hgs(2) = hy_m
    if (0 < hz_m) fs%hgs(3) = hz_m
    
    fs%mg%is(1) = -nxvacl_m
    
    write(9999, *) 'fs%mg%ndir:', fs%mg%ndir
    write(9999, *) 'fs%mg%nd:', fs%mg%nd
    write(9999, *) 'fs%hgs(1:3):', fs%hgs(1:3)
    flush(9999)
    
    fw%fdtddim = '1d'
    fs%mg%is(1) = -nxvacl_m
    fs%mg%ie(1) = nx_m+nxvacr_m
    fs%mg%is(2) = 1
    fs%mg%ie(2) = ny_m
    fs%mg%is(3) = 1
    fs%mg%ie(3) = nz_m
    fs%mg%is_overlap(1) = 0
    fs%mg%ie_overlap(1) = nx_m + 1
    fs%mg%is_overlap(2) = 0
    fs%mg%ie_overlap(2) = ny_m + 1
    fs%mg%is_overlap(3) = 0
    fs%mg%ie_overlap(3) = nz_m + 1
    fs%mg%is_array(1) = 0
    fs%mg%ie_array(1) = nx_m + 1
    fs%mg%is_array(2) = 0
    fs%mg%ie_array(2) = ny_m + 1
    fs%mg%is_array(3) = 0
    fs%mg%ie_array(3) = nz_m + 1
    allocate(fs%mg%coordinate(minval(fs%mg%is_overlap):maxval(fs%mg%ie_overlap),1:3))
    fs%rlsize(1) = 1d0
    fs%rlsize(2) = 1d0
    fs%rlsize(3) = 1d0
    fs%origin(1) = 0d0
    fs%origin(2) = 0d0
    fs%origin(3) = 0d0
    fs%a_bc(1,1:2) = 'pbc'
    fs%a_bc(2,1:2) = 'pbc'
    fs%a_bc(3,1:2) = 'pbc'
    fs%imedia(:,:,:) = 0

    write(9999, *) 'Initialization FDTD Weyl start';flush(9999)
    call Weyl_init(fs, fw)
    write(9999, *) 'Initialization FDTD Weyl end';flush(9999)
    
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

        
    write(9999, *) 'Initialization Complete'
    flush(9999)
    return    
end subroutine initialization_ms


subroutine time_evolution_step_ms_macro
    implicit none
    ! Override Global Variables
    nproc_group_global = ms%icomm_macropoint
    nproc_id_global = ms%id_macropoint
    nproc_size_global = ms%isize_macropoint

    call weyl_calc(fs, fw)

    if(mod(itt,2)==1)then
        call time_evolution_step(Mit,itotNtime,itt,lg,mg,system,rt,info,stencil,xc_func &
         & ,srg,srg_scalar,pp,ppg,ppn,spsi_in,spsi_out,tpsi,srho,srho_s,V_local,Vbox,sVh,sVh_stock1,sVh_stock2,sVxc &
         & ,sVpsl,dmat,fg,energy,ewald,md,ofl,poisson,singlescale)
      else
        call time_evolution_step(Mit,itotNtime,itt,lg,mg,system,rt,info,stencil,xc_func &
         & ,srg,srg_scalar,pp,ppg,ppn,spsi_out,spsi_in,tpsi,srho,srho_s,V_local,Vbox,sVh,sVh_stock1,sVh_stock2,sVxc &
         & ,sVpsl,dmat,fg,energy,ewald,md,ofl,poisson,singlescale)
    end if

    ! Override Global Variables (Repair)
    nproc_group_global = ms%icomm_ms_world
    nproc_id_global = ms%id_ms_world
    nproc_size_global = ms%isize_ms_world

    return
end subroutine time_evolution_step_ms_macro


end subroutine main_ms

