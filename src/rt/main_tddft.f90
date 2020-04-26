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

subroutine main_tddft
use math_constants, only: pi
use salmon_global
use structures
use communication, only: comm_is_root, comm_sync_all
use salmon_xc, only: finalize_xc
use timer
use write_sub, only: write_response_0d,write_response_3d,write_pulse_0d,write_pulse_3d
use initialization_rt_sub
use checkpoint_restart_sub
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

integer :: Mit, itt,itotNtime

call timer_begin(LOG_TOTAL)

call print_header()

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


#ifdef __FUJITSU
call fapp_start('time_evol',1,0) ! performance profiling
#endif

call comm_sync_all
call timer_enable_sub
call timer_begin(LOG_RT_ITERATION)
TE : do itt=Mit+1,itotNtime

  if(mod(itt,2)==1)then
    call time_evolution_step(Mit,itotNtime,itt,lg,mg,system,rt,info,stencil,xc_func &
     & ,srg,srg_scalar,pp,ppg,ppn,spsi_in,spsi_out,tpsi,srho,srho_s,V_local,Vbox,sVh,sVh_stock1,sVh_stock2,sVxc &
     & ,sVpsl,dmat,fg,energy,ewald,md,ofl,poisson,singlescale)
  else
    call time_evolution_step(Mit,itotNtime,itt,lg,mg,system,rt,info,stencil,xc_func &
     & ,srg,srg_scalar,pp,ppg,ppn,spsi_out,spsi_in,tpsi,srho,srho_s,V_local,Vbox,sVh,sVh_stock1,sVh_stock2,sVxc &
     & ,sVpsl,dmat,fg,energy,ewald,md,ofl,poisson,singlescale)
  end if

  if((checkpoint_interval >= 1) .and. (mod(itt,checkpoint_interval) == 0)) then
    call timer_begin(LOG_CHECKPOINT_SYNC)
    call timer_begin(LOG_CHECKPOINT_SELF)
    if (mod(itt,2)==1) then
      call checkpoint_rt(lg,mg,system,info,spsi_out,itt,sVh_stock1,sVh_stock2,singlescale)
    else
      call checkpoint_rt(lg,mg,system,info,spsi_in, itt,sVh_stock1,sVh_stock2,singlescale)
    endif
    call timer_end(LOG_CHECKPOINT_SELF)
    call comm_sync_all
    call timer_end(LOG_CHECKPOINT_SYNC)
  endif

end do TE
call timer_end(LOG_RT_ITERATION)
call timer_disable_sub

#ifdef __FUJITSU
call fapp_stop('time_evol',1,0) ! performance profiling
#endif

close(030) ! laser


!--------------------------------- end of time-evolution


!------------ Writing part -----------

! write RT: 
call timer_begin(LOG_WRITE_RT_RESULTS)

!
select case(iperiodic)
case(0)
  if(theory=="tddft_response")then
    call write_response_0d(ofl,rt)
  else
    call write_pulse_0d(ofl,rt)
  end if
case(3)
  if(theory=="tddft_response")then
    call write_response_3d(ofl,rt)
  else
    call write_pulse_3d(ofl,rt)
  end if
end select

call timer_end(LOG_WRITE_RT_RESULTS)
call timer_end(LOG_TOTAL)

if(write_rt_wfn_k=='y')then
  call checkpoint_rt(lg,mg,system,info,spsi_out,Mit,sVh_stock1,sVh_stock2,singlescale,ofl%dir_out_restart)
end if

call finalize_xc(xc_func)

contains

subroutine print_header()
  use parallelization, only: nproc_id_global
  implicit none
  !(header of standard output)
  if(comm_is_root(nproc_id_global))then
    write(*,*)
    select case(iperiodic)
    case(0)
      write(*,'(1x,a10,a11,a48,a15,a18,a10)') &
                  "time-step ", "time[fs]",   &
                  "Dipole moment(xyz)[A]"     &
                ,"electrons", "Total energy[eV]", "iterVh"
    case(3)
      write(*,'(1x,a10,a11,a48,a15,a18)')   &
                  "time-step", "time[fs] ", &
                  "Current(xyz)[a.u.]",     &
                  "electrons", "Total energy[eV] "
    end select
    write(*,'("#",7("----------"))')
  endif
end subroutine print_header

end subroutine main_tddft

