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

subroutine scf_iteration_step(lg,mg,ng,system,info,info_field,pinfo,stencil, &
               srg,srg_ng,spsi,shpsi,srho,srho_s, &
               cg,ppg,vlocal,  &
               miter,   &
               iditer_nosubspace_diag,mixing,iter, &
               poisson,fg,sVh,xc_func,ppn,sVxc,energy)
  use inputoutput, only: calc_mode,iperiodic,method_mixing,mixrate &
                        ,yn_subspace_diagonalization
  use structures
  use timer
  use gram_schmidt_orth, only: gram_schmidt
  use Conjugate_Gradient, only: gscg_isolated,gscg_periodic
  use subspace_diagonalization, only: ssdg_isolated,ssdg_periodic
  use density_matrix, only: calc_density
  use mixing_sub
  use hartree_sub, only: hartree
  use salmon_xc
  implicit none

  type(s_rgrid),         intent(in)    :: lg
  type(s_rgrid),         intent(in)    :: mg
  type(s_rgrid),         intent(in)    :: ng
  type(s_dft_system),    intent(in)    :: system
  type(s_orbital_parallel),intent(in)  :: info
  type(s_field_parallel),intent(in)    :: info_field
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
  integer                              :: j

! solve Kohn-Sham equation by minimization techniques
  call timer_begin(LOG_CALC_MINIMIZATION)

  select case(iperiodic)
  case(0)
    call gscg_isolated(mg,system,info,stencil,ppg,vlocal,srg,spsi,cg)
  case(3)
    call gscg_periodic(mg,system,info,stencil,ppg,vlocal,srg,spsi,cg)
  end select

  call timer_end(LOG_CALC_MINIMIZATION)

! Gram Schmidt orghonormalization
  call gram_schmidt(system, mg, info, spsi)

! subspace diagonalization
  if(yn_subspace_diagonalization == 'y')then
    if(miter>iditer_nosubspace_diag)then
      select case(iperiodic)
      case(0)      
        call ssdg_isolated(mg,system,info,pinfo,stencil,spsi,shpsi,ppg,vlocal,srg)

      case(3)
        call ssdg_periodic(mg,system,info,stencil,spsi,shpsi,ppg,vlocal,srg)
      end select
    end if
  end if
  
! density
  call timer_begin(LOG_CALC_RHO)

  call calc_density(system,srho_s,spsi,info,mg)

  select case(method_mixing)
    case ('simple') ; call simple_mixing(ng,system,1.d0-mixrate,mixrate,srho_s,mixing)
    case ('broyden'); call wrapper_broyden(ng,system,srho_s,iter,mixing)
  end select
  call timer_end(LOG_CALC_RHO)

  srho%f = 0d0
  do j=1,system%nspin
    srho%f = srho%f + srho_s(j)%f
  end do


  if(calc_mode/='DFT_BAND')then

    call timer_begin(LOG_CALC_HARTREE)
    call hartree(lg,mg,ng,info_field,system,poisson,srg_ng,stencil,srho,sVh,fg)
    call timer_end(LOG_CALC_HARTREE)
  
    call timer_begin(LOG_CALC_EXC_COR)
    call exchange_correlation(system,xc_func,ng,srg_ng,srho_s,ppn,info_field%icomm_all,sVxc,energy%E_xc)
    call timer_end(LOG_CALC_EXC_COR)

  end if

end subroutine scf_iteration_step

end module scf_iteration_sub
