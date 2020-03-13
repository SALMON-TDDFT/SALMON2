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
module scf_iteration_sub
  implicit none

contains

subroutine scf_iteration_step(lg,mg,ng,system,info,pinfo,stencil, &
               srg,srg_ng,spsi,shpsi,srho,srho_s, &
               cg,ppg,vlocal,  &
               miter,   &
               iditer_nosubspace_diag,mixing,iter, &
               poisson,fg,sVh,xc_func,ppn,sVxc,energy)
  use salmon_global, only: calc_mode,method_mixing,mixrate &
                        ,yn_subspace_diagonalization,ncg,ncg_init
  use structures
  use timer
  use gram_schmidt_orth, only: gram_schmidt
  use Conjugate_Gradient, only: gscg_zwf,gscg_rwf
  use subspace_diagonalization, only: ssdg
  use density_matrix, only: calc_density
  use mixing_sub
  use hartree_sub, only: hartree
  use salmon_xc
  implicit none

  type(s_rgrid),         intent(in)    :: lg
  type(s_rgrid),         intent(in)    :: mg
  type(s_rgrid),         intent(in)    :: ng
  type(s_dft_system),    intent(in)    :: system
  type(s_parallel_info), intent(in)    :: info
  type(s_process_info),  intent(in)    :: pinfo
  type(s_orbital),       intent(inout) :: spsi,shpsi
  type(s_scalar),        intent(inout) :: srho
  type(s_scalar),        intent(inout) :: srho_s(system%nspin)
  type(s_stencil),       intent(in)    :: stencil
  type(s_sendrecv_grid), intent(inout) :: srg
  type(s_sendrecv_grid), intent(inout) :: srg_ng
  type(s_pp_grid),       intent(in)    :: ppg
  type(s_cg),            intent(inout) :: cg
  type(s_scalar),        intent(in)    :: vlocal(system%nspin)
  integer,               intent(in)    :: miter
  integer,               intent(in)    :: iditer_nosubspace_diag
  type(s_mixing),        intent(inout) :: mixing
  integer,               intent(in)    :: iter
  type(s_poisson),       intent(inout) :: poisson
  type(s_reciprocal_grid),intent(inout):: fg
  type(s_scalar),        intent(inout) :: sVh
  type(s_xc_functional), intent(in)    :: xc_func
  type(s_pp_nlcc),       intent(in)    :: ppn
  type(s_scalar),        intent(inout) :: sVxc(system%nspin)
  type(s_dft_energy),    intent(inout) :: energy
  !
  integer :: j,nncg

  if(miter==1) then
    nncg = ncg_init
  else
    nncg = ncg
  end if

! solve Kohn-Sham equation by minimization techniques
  call timer_begin(LOG_CALC_MINIMIZATION)
  if(system%if_real_orbital) then
    call gscg_rwf(nncg,mg,system,info,stencil,ppg,vlocal,srg,spsi,cg)
  else
    call gscg_zwf(nncg,mg,system,info,stencil,ppg,vlocal,srg,spsi,cg)
  end if
  call timer_end(LOG_CALC_MINIMIZATION)

! Gram Schmidt orghonormalization
  call gram_schmidt(system, mg, info, spsi, pinfo)

! subspace diagonalization
  call timer_begin(LOG_CALC_SUBSPACE_DIAG)
  if(yn_subspace_diagonalization == 'y')then
    if(miter>iditer_nosubspace_diag)then
      call ssdg(mg,system,info,pinfo,stencil,spsi,shpsi,ppg,vlocal,srg)
    end if
  end if
  call timer_end(LOG_CALC_SUBSPACE_DIAG)

  if(mixing%flag_mix_zero) return

! density
  call timer_begin(LOG_CALC_RHO)

  call calc_density(system,srho_s,spsi,info,mg)

  select case(method_mixing)
    case ('simple') ; call simple_mixing(ng,system,1.d0-mixrate,mixrate,srho_s,mixing)
    case ('broyden'); call wrapper_broyden(info%icomm_r,ng,system,srho_s,iter,mixing)
    case ('pulay')  ; call pulay(mg,info,system,srho_s,iter,mixing)
  end select
  call timer_end(LOG_CALC_RHO)

  srho%f = 0d0
  do j=1,system%nspin
    srho%f = srho%f + srho_s(j)%f
  end do

  if(calc_mode/='DFT_BAND')then

    call timer_begin(LOG_CALC_HARTREE)
    call hartree(lg,mg,info,system,fg,poisson,srg_ng,stencil,srho,sVh)
    call timer_end(LOG_CALC_HARTREE)

    call timer_begin(LOG_CALC_EXC_COR)
    call exchange_correlation(system,xc_func,ng,mg,srg_ng,srg,srho_s,ppn,info,spsi,stencil,sVxc,energy%E_xc)
    call timer_end(LOG_CALC_EXC_COR)

  end if

end subroutine scf_iteration_step

end module scf_iteration_sub
