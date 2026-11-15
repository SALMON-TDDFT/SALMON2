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

! SALMONpy one-point embedded TDDFT. Native physics uses existing RT kernels.
#include "config.h"
module salmonpy_tddft_context
use tddft_worker
use density_matrix, only: calc_current
use hamiltonian, only: update_kvector_nonlocalpt
use math_constants, only: pi
use salmon_global
use structures
use parallelization, only: adjust_elapse_time, nproc_group_global
use communication, only: comm_is_root, comm_sync_all, comm_bcast
use salmon_xc, only: finalize_xc
use timer
use write_sub, only: write_response_0d,write_response_3d,write_pulse_0d,write_pulse_3d
use initialization_rt_sub
use checkpoint_restart_sub
use jellium, only: check_condition_jm
use nvtx_wrapper
use parallelization, only: nproc_id_global
use iso_c_binding
use sendrecv_grid, only: dealloc_cache
implicit none
private
public :: sp_tddft_init,sp_tddft_step,sp_tddft_close
type :: rt_context
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
type(s_unfold) :: unfold

real(8) :: worker_A(3),worker_current(3,2)
integer :: Mit, itt
logical :: is_checkpoint_iter, is_shutdown_time, is_checkpoint

end type
type(rt_context),allocatable,save :: state
logical,save :: used=.false.
real(8),save :: accepted_A(3)=0d0
contains
integer(c_int) function sp_tddft_init(path, timestep, capacity, A, decomposition, J) bind(C)
  use inputoutput, only: read_input
  use parallelization, only: setup_parallel
  use math_constants, only: set_math_constants
  use sym_sub, only: read_sw_symmetry
  character(c_char),intent(in) :: path(*)
  real(c_double),value :: timestep
  integer(c_int),value :: capacity
  real(c_double),intent(in) :: A(3)
  integer(c_int),intent(in) :: decomposition(5)
  real(c_double),intent(out) :: J(3)
  character(4096) :: filename
  integer :: i
  sp_tddft_init=1
  if(used) return
  used=.true.
  filename=''
  do i=1,len(filename)
    if(path(i)==c_null_char) exit
    filename(i:i)=path(i)
  end do
  if(i>len(filename)) return
  worker_embedded=.true.
  worker_initial_A=A; embedded_dt=timestep; embedded_steps=capacity
  embedded_decomposition=decomposition
  call set_math_constants
  call setup_parallel
  call read_input(trim(filename))
  calc_mode='RT'
  call read_sw_symmetry(yn_symmetry)
  call timer_initialize
  quiet=.false.
  allocate(state)
call timer_begin(LOG_TOTAL)

call initialization_rt( state%Mit, state%system, state%energy, state%ewald, state%rt, state%md, &
                        state%singlescale,  &
                        state%stencil, state%fg, state%poisson,  &
                        state%lg, state%mg,   &
                        state%info,  &
                        state%xc_func, state%ofl,  &
                        state%srg, state%srg_scalar,  &
                        state%spsi_in, state%spsi_out, state%tpsi, state%rho, state%rho_jm, state%rho_s,  &
                        state%V_local, state%Vbox, state%Vh, state%Vh_stock1, state%Vh_stock2, state%Vxc, state%Vpsl,&
                        state%pp, state%ppg, state%ppn, state%unfold )

if(worker_requested()) then
  state%system%vec_Ac=worker_initial_A
  call update_kvector_nonlocalpt(state%info%ik_s,state%info%ik_e,state%system,state%ppg)
  call calc_current(state%system,state%mg,state%stencil,state%info,state%srg,state%spsi_in,state%ppg, &
     & state%worker_current(:,1:state%system%nspin))
  state%spsi_in%update_zwf_overlap=.true.
  J=sum(state%worker_current(:,1:state%system%nspin),dim=2)
end if

  state%itt=state%Mit
  call comm_sync_all
  call timer_enable_sub
  call timer_begin(LOG_RT_ITERATION)

  accepted_A=A
  sp_tddft_init=0
end function
integer(c_int) function sp_tddft_step(Aold,Anew,J) bind(C)
  real(c_double),intent(in) :: Aold(3),Anew(3)
  real(c_double),intent(out) :: J(3)
  sp_tddft_step=1
  if(.not.allocated(state)) return
  if(any(Aold/=accepted_A).or.state%itt>=nt) return
  state%itt=state%itt+1
  state%rt%Ac_ext(:,state%itt)=Anew
  state%rt%Ac_tot(:,state%itt)=Anew
  if(mod(state%itt,2)==1)then
    call time_evolution_step(state%Mit,nt,state%itt,state%lg,state%mg,state%system,state%rt,state%info,state%stencil,state%xc_func &
     & ,state%srg,state%srg_scalar,state%pp,state%ppg,state%ppn,state%spsi_in,state%spsi_out,state%tpsi &
     & ,state%rho,state%rho_jm,state%rho_s,state%V_local,state%Vbox,state%Vh,state%Vh_stock1,state%Vh_stock2,state%Vxc &
     & ,state%Vpsl,state%fg,state%energy,state%ewald,state%md,state%ofl,state%poisson,state%singlescale,state%unfold)
  else
    call time_evolution_step(state%Mit,nt,state%itt,state%lg,state%mg,state%system,state%rt,state%info,state%stencil,state%xc_func &
     & ,state%srg,state%srg_scalar,state%pp,state%ppg,state%ppn,state%spsi_out,state%spsi_in,state%tpsi &
     & ,state%rho,state%rho_jm,state%rho_s,state%V_local,state%Vbox,state%Vh,state%Vh_stock1,state%Vh_stock2,state%Vxc &
     & ,state%Vpsl,state%fg,state%energy,state%ewald,state%md,state%ofl,state%poisson,state%singlescale,state%unfold)
  end if

  J=state%rt%curr(:,state%itt)

  accepted_A=Anew
  sp_tddft_step=0
end function
integer(c_int) function sp_tddft_close() bind(C)
  sp_tddft_close=1
  if(.not.allocated(state)) return
  call timer_end(LOG_RT_ITERATION)
  call timer_disable_sub
  close(030)
  if(comm_is_root(nproc_id_global)) close(state%ofl%fh_rt)
  call finalize_xc(state%xc_func)
  call dealloc_cache(state%srg)
  call dealloc_cache(state%srg_scalar)
  deallocate(state)
  ! Other SALMON module configuration remains process-owned. Do not reopen.
  sp_tddft_close=0
end function
end module
