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
MODULE global_variables_ms

use inputoutput
use allocate_mat_sub
use deallocate_mat_sub
implicit none

END MODULE global_variables_ms

!=======================================================================

subroutine main_ms
use math_constants, only: pi
use structures
use inputoutput, only: nx_m, ny_m, nz_m
use salmon_parallel, only: nproc_group_global, nproc_size_global, nproc_id_global
use salmon_communication, only: comm_is_root, comm_sync_all, comm_create_group, comm_get_groupinfo
use salmon_xc, only: finalize_xc
use timer
use global_variables_ms
use write_sub, only: write_response_0d,write_response_3d,write_pulse_0d,write_pulse_3d
use initialization_rt_sub
use fdtd_coulomb_gauge, only: ls_singlescale
use checkpoint_restart_sub
implicit none

type(s_rgrid) :: lg
type(s_rgrid) :: mg
type(s_rgrid) :: ng
type(s_dft_system)  :: system
type(s_rt) :: rt
type(s_orbital_parallel) :: info
type(s_field_parallel) :: info_field
type(s_poisson) :: poisson
type(s_stencil) :: stencil
type(s_xc_functional) :: xc_func
type(s_reciprocal_grid) :: fg
type(s_dft_energy) :: energy
type(s_md) :: md
type(s_ofile) :: ofl
type(s_scalar) :: sVpsl
type(s_scalar) :: srho,sVh,sVh_stock1,sVh_stock2,Vbox
type(s_scalar),allocatable :: srho_s(:),V_local(:),sVxc(:)
type(s_dmatrix) :: dmat
type(s_orbital) :: spsi_in,spsi_out
type(s_orbital) :: tpsi ! temporary wavefunctions
type(s_sendrecv_grid) :: srg,srg_ng
type(s_pp_info) :: pp
type(s_pp_grid) :: ppg
type(s_pp_nlcc) :: ppn
type(s_vector)  :: j_e ! microscopic electron number current density
type(ls_singlescale) :: singlescale
type(s_ofile) :: ofile

integer :: nmacro
integer :: nproc_group_macropoint
integer :: nproc_size_macropoint
integer :: nproc_id_macropoint
integer :: nmacro_mygrp
integer :: imacro_mygrp_s
integer :: imacro_mygrp_e


integer :: Mit
integer :: nntime

call timer_begin(LOG_TOTAL)

nmacro = nx_m * ny_m * nz_m




! call initialization_ms( Mit, system, energy, rt, md, singlescale,  &
!                         stencil, fg, poisson,  &
!                         lg, mg, ng,  &
!                         info, info_field,  &
!                         xc_func, dmat, ofl, j_e,  &
!                         srg, srg_ng,  &
!                         spsi_in, spsi_out, tpsi, srho, srho_s,  &
!                         V_local, Vbox, sVh, sVh_stock1, sVh_stock2, sVxc, sVpsl,&
!                         pp, ppg, ppn )

! call initialization_rt( Mit, system, energy, rt, md, singlescale,  &
!                         stencil, fg, poisson,  &
!                         lg, mg, ng,  &
!                         info, info_field,  &
!                         xc_func, dmat, ofl, j_e,  &
!                         srg, srg_ng,  &
!                         spsi_in, spsi_out, tpsi, srho, srho_s,  &
!                         V_local, Vbox, sVh, sVh_stock1, sVh_stock2, sVxc, sVpsl,&
!                         pp, ppg, ppn )


! call timer_begin(LOG_RT_ITERATION)
! TE : do itt=Mit+1,itotNtime

!   if(mod(itt,2)==1)then
!     call time_evolution_step(Mit,itotNtime,itt,lg,mg,ng,system,rt,info,info_field,stencil,xc_func &
!      & ,srg,srg_ng,pp,ppg,ppn,spsi_in,spsi_out,tpsi,srho,srho_s,V_local,Vbox,sVh,sVh_stock1,sVh_stock2,sVxc &
!      & ,sVpsl,dmat,fg,energy,md,ofl,poisson,j_e,singlescale)
!   else
!     call time_evolution_step(Mit,itotNtime,itt,lg,mg,ng,system,rt,info,info_field,stencil,xc_func &
!      & ,srg,srg_ng,pp,ppg,ppn,spsi_out,spsi_in,tpsi,srho,srho_s,V_local,Vbox,sVh,sVh_stock1,sVh_stock2,sVxc &
!      & ,sVpsl,dmat,fg,energy,md,ofl,poisson,j_e,singlescale)
!   end if

!   if((checkpoint_interval >= 1) .and. (mod(itt,checkpoint_interval) == 0)) then
!     call timer_begin(LOG_CHECKPOINT_SYNC)
!     call timer_begin(LOG_CHECKPOINT_SELF)
!     if (mod(itt,2)==1) then
!       call checkpoint_rt(lg,mg,ng,system,info,spsi_out,itt,sVh_stock1=sVh_stock1,sVh_stock2=sVh_stock2)
!     else
!       call checkpoint_rt(lg,mg,ng,system,info,spsi_in, itt,sVh_stock1=sVh_stock1,sVh_stock2=sVh_stock2)
!     endif
!     call timer_end(LOG_CHECKPOINT_SELF)
!     call comm_sync_all
!     call timer_end(LOG_CHECKPOINT_SYNC)
!   endif

! end do TE
! call timer_end(LOG_RT_ITERATION)

! close(030) ! laser


! !--------------------------------- end of time-evolution


! !------------ Writing part -----------

! ! write RT: 
! call timer_begin(LOG_WRITE_RT_RESULTS)

! !
! select case(iperiodic)
! case(0)
!   if(theory=="TDDFT_response")then
!     call write_response_0d(ofl,rt)
!   else
!     call write_pulse_0d(ofl,rt)
!   end if
! case(3)
!   if(theory=="TDDFT_response")then
!     call write_response_3d(ofl,rt)
!   else
!     call write_pulse_3d(ofl,rt)
!   end if
! end select

! call timer_end(LOG_WRITE_RT_RESULTS)
! call timer_end(LOG_TOTAL)

! if(write_rt_wfn_k=='y')then
!   call write_bin(ofile%dir_out_restart,lg,mg,ng,system,info,spsi_out,Mit,sVh_stock1=sVh_stock1,sVh_stock2=sVh_stock2)
! end if

! call deallocate_mat

! call finalize_xc(xc_func)

end subroutine main_ms

