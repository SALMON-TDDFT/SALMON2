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
use communication, only: comm_is_root, comm_sync_all, comm_create_group_byid, comm_get_groupinfo
use salmon_xc, only: finalize_xc
use timer
use write_sub, only: write_response_0d,write_response_3d,write_pulse_0d,write_pulse_3d
use initialization_rt_sub
use checkpoint_restart_sub
use fdtd_weyl, only: ls_fdtd_weyl, weyl_init, weyl_calc, weyl_finalize
use parallelization, only: nproc_id_global, nproc_size_global, nproc_group_global
implicit none

type(s_rgrid) :: lg
type(s_rgrid) :: mg
type(s_dft_system)  :: system
type(s_rt) :: rt
type(s_parallel_info) :: info
type(s_process_info) :: pinfo
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
type(s_sendrecv_grid) :: srg, srg_scalar
type(s_pp_info) :: pp
type(s_pp_grid) :: ppg
type(s_pp_nlcc) :: ppn
type(s_ofile) :: ofile

type(s_fdtd_system) :: fs
type(ls_fdtd_weyl) :: fw
type(s_multiscale) :: ms

integer :: nmacro
integer :: nproc_group_macropoint
integer :: nproc_size_macropoint
integer :: nproc_id_macropoint
integer :: nmacro_mygrp
integer :: imacro_mygrp_s
integer :: imacro_mygrp_e


integer :: Mit, itt, itotNtime
integer :: nntime

integer :: i

integer, allocatable :: iranklists(:)

character(256) :: file_debug_log
integer :: nmacro_mygroup, isize_mygroup

call timer_begin(LOG_TOTAL)

! Open logfile for debugging
write(file_debug_log, "('debug_log_', i3.3, '.txt')") nproc_id_global
open(unit=9999, file=file_debug_log)
write(9999, *) 'logging start'; flush(9999)

ms%nmacro = nx_m * ny_m * nz_m
ms%icomm_multiscale = nproc_group_global
ms%isize_multiscale = nproc_size_global
ms%id_multiscale = nproc_id_global

if (ms%nmacro <= ms%isize_multiscale) then
    if (mod(ms%isize_multiscale, ms%nmacro) == 0) then
        nmacro_mygroup = 1
        isize_mygroup = ms%isize_multiscale / ms%nmacro
        ms%imacro_mygroup_s = ms%id_multiscale / isize_mygroup + 1
        ms%imacro_mygroup_e = ms%imacro_mygroup_s
        ms%id_mygroup_s = (ms%imacro_mygroup_s - 1) * isize_mygroup
        ms%id_mygroup_e = ms%id_mygroup_s + isize_mygroup - 1
    else
        stop "mod(number of procs, number of points) must be 0!"
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
flush(9999)

if (ms%nmacro < 1) stop "Invalid macropoint number"

allocate(iranklists(isize_mygroup))
do i = 1, isize_mygroup
    iranklists(i) = ms%id_mygroup_s + (i - 1)
end do

write(9999, *) 'iranklists:', iranklists(:)
flush(9999)

ms%icomm_macropoint = comm_create_group_byid(ms%icomm_multiscale, iranklists)
call comm_get_groupinfo(ms%icomm_macropoint, ms%id_macropoint, ms%isize_macropoint)

write(9999, *) 'ms%icomm_multiscale:', ms%icomm_multiscale
write(9999, *) 'ms%id_multiscale:', ms%id_multiscale
write(9999, *) 'ms%isize_multiscale:', ms%isize_multiscale
write(9999, *) 'ms%icomm_macropoint:', ms%icomm_macropoint
write(9999, *) 'ms%id_macropoint:', ms%id_macropoint
write(9999, *) 'ms%isize_macropoint:', ms%isize_macropoint
flush(9999)





fs%lg%ndir = 3
fs%lg%nd = 1
fs%lg%is(1) = 1
fs%lg%ie(1) = nx_m
fs%lg%is(2) = 1
fs%lg%ie(2) = ny_m
fs%lg%is(3) = 1
fs%lg%ie(3) = nz_m
fs%lg%is_overlap(1) = 0
fs%lg%ie_overlap(1) = nx_m + 1
fs%lg%is_overlap(2) = 0
fs%lg%ie_overlap(2) = ny_m + 1
fs%lg%is_overlap(3) = 0
fs%lg%ie_overlap(3) = nz_m + 1
fs%lg%is_array(1) = 0
fs%lg%ie_array(1) = nx_m + 1
fs%lg%is_array(2) = 0
fs%lg%ie_array(2) = ny_m + 1
fs%lg%is_array(3) = 0
fs%lg%ie_array(3) = nz_m + 1
! allocate(fs%lg%idx(fs%lg%is_array(1):fs%lg%ie_array(1)))
! allocate(fs%lg%idy(fs%lg%is_array(2):fs%lg%ie_array(2)))
! allocate(fs%lg%idy(fs%lg%is_array(3):fs%lg%ie_array(3)))
allocate(fs%lg%coordinate(minval(fs%lg%is_overlap):maxval(fs%lg%ie_overlap),1:3))
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


! call initialization_ms( Mit, system, energy, ewald, rt, md, singlescale,  &
!                         stencil, fg, poisson,  &
!                         lg, mg, ng,  &
!                         info, info_field,  &
!                         xc_func, dmat, ofl, j_e,  &
!                         srg, srg_scalar,  &
!                         spsi_in, spsi_out, tpsi, srho, srho_s,  &
!                         V_local, Vbox, sVh, sVh_stock1, sVh_stock2, sVxc, sVpsl,&
!                         pp, ppg, ppn )

! call initialization_rt( Mit, itotNtime, system, energy, ewald, rt, md, &
!                         singlescale,  &
!                         stencil, fg, poisson,  &
!                         lg, mg,   &
!                         info,pinfo,  &
!                         xc_func, dmat, ofl,  &
!                         srg, srg_scalar,  &
!                         spsi_in, spsi_out, tpsi, srho, srho_s,  &
!                         V_local, Vbox, sVh, sVh_stock1, sVh_stock2, sVxc, sVpsl,&
!                         pp, ppg, ppn )



! call comm_sync_all
! call timer_enable_sub
! call timer_begin(LOG_RT_ITERATION)
! TE : do itt=Mit+1,itotNtime

!   if(mod(itt,2)==1)then
!     call time_evolution_step(Mit,itotNtime,itt,lg,mg,system,rt,info,pinfo,stencil,xc_func &
!      & ,srg,srg_scalar,pp,ppg,ppn,spsi_in,spsi_out,tpsi,srho,srho_s,V_local,Vbox,sVh,sVh_stock1,sVh_stock2,sVxc &
!      & ,sVpsl,dmat,fg,energy,ewald,md,ofl,poisson,singlescale)
!   else
!     call time_evolution_step(Mit,itotNtime,itt,lg,mg,system,rt,info,pinfo,stencil,xc_func &
!      & ,srg,srg_scalar,pp,ppg,ppn,spsi_out,spsi_in,tpsi,srho,srho_s,V_local,Vbox,sVh,sVh_stock1,sVh_stock2,sVxc &
!      & ,sVpsl,dmat,fg,energy,ewald,md,ofl,poisson,singlescale)
!   end if

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

! end do TE
! call timer_end(LOG_RT_ITERATION)
! call timer_disable_sub

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
! call timer_end(LOG_TOTAL)

! if(write_rt_wfn_k=='y')then
!   call checkpoint_rt(lg,mg,system,info,spsi_out,Mit,sVh_stock1,sVh_stock2,singlescale,ofl%dir_out_restart)
! end if

! call finalize_xc(xc_func)

close(9999)

end subroutine main_ms

