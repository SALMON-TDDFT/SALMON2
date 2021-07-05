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
module scf_iteration_sub
  implicit none

contains

subroutine scf_iteration_step(lg,mg,system,info,stencil, &
               srg,srg_scalar,spsi,shpsi,rho,rho_jm,rho_s, &
               cg,ppg,vlocal,  &
               miter,   &
               nscf_init_no_diagonal, mixing, iter, &
               poisson,fg,Vh,xc_func,ppn,Vxc,energy )
  use salmon_global, only: calc_mode,method_mixing  &
                        ,yn_subspace_diagonalization,ncg,ncg_init,yn_jm,yn_spinorbit
  use structures
  use timer
  use gram_schmidt_orth, only: gram_schmidt
  use Conjugate_Gradient, only: gscg_zwf,gscg_rwf
  use subspace_diagonalization, only: ssdg
  use density_matrix, only: calc_density
  use mixing_sub
  use hartree_sub, only: hartree
  use salmon_xc
  use noncollinear_module, only: simple_mixing_so
  implicit none

  type(s_rgrid),          intent(in)    :: lg
  type(s_rgrid),          intent(in)    :: mg
  type(s_dft_system),     intent(in)    :: system
  type(s_parallel_info),  intent(in)    :: info
  type(s_orbital),        intent(inout) :: spsi,shpsi
  type(s_scalar),         intent(inout) :: rho
  type(s_scalar),         intent(in)    :: rho_jm
  type(s_scalar),         intent(inout) :: rho_s(system%nspin)
  type(s_stencil),        intent(in)    :: stencil
  type(s_sendrecv_grid),  intent(inout) :: srg
  type(s_sendrecv_grid),  intent(inout) :: srg_scalar
  type(s_pp_grid),        intent(in)    :: ppg
  type(s_cg),             intent(inout) :: cg
  type(s_scalar),         intent(in)    :: vlocal(system%nspin)
  integer,                intent(in)    :: miter
  integer,                intent(in)    :: nscf_init_no_diagonal
  type(s_mixing),         intent(inout) :: mixing
  integer,                intent(in)    :: iter
  type(s_poisson),        intent(inout) :: poisson
  type(s_reciprocal_grid),intent(inout) :: fg
  type(s_scalar),         intent(inout) :: Vh
  type(s_xc_functional),  intent(in)    :: xc_func
  type(s_pp_nlcc),        intent(in)    :: ppn
  type(s_scalar),         intent(inout) :: Vxc(system%nspin)
  type(s_dft_energy),     intent(inout) :: energy
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
  call gram_schmidt(system, mg, info, spsi)

! subspace diagonalization
  call timer_begin(LOG_CALC_SUBSPACE_DIAG)
  if(yn_subspace_diagonalization == 'y')then
    if(miter > nscf_init_no_diagonal)then
      call ssdg(mg,system,info,stencil,spsi,shpsi,ppg,vlocal,srg)
    end if
  end if
  call timer_end(LOG_CALC_SUBSPACE_DIAG)

  if(mixing%flag_mix_zero) return

  if(calc_mode/='DFT_BAND')then

! density
    call timer_begin(LOG_CALC_RHO)

    call calc_density(system,rho_s,spsi,info,mg)

    select case(method_mixing)
    case ('simple')
      call simple_mixing(mg,system,1.d0-mixing%mixrate,mixing%mixrate,rho_s,mixing)
    case ('simple_dm')
      if(yn_spinorbit=='n') stop 'yn_spinorbit must be y when method_mixing=simple_dm'
      call simple_mixing_so(mg,system,1.d0-mixing%mixrate,mixing%mixrate,rho_s,mixing)
    case ('broyden')
      call wrapper_broyden(info%icomm_r,mg,system,rho_s,iter,mixing)
    case ('pulay')
      call pulay(mg,info,system,rho_s,iter,mixing)
    case default
      stop 'Invalid method_mixing. Specify any one of "simple" or "broyden" or "pulay" for method_mixing.'
    end select
    call timer_end(LOG_CALC_RHO)

    rho%f = 0d0
    do j=1,system%nspin
      rho%f = rho%f + rho_s(j)%f
    end do
  
    if(yn_jm=='y') rho%f = rho%f + rho_jm%f

    call timer_begin(LOG_CALC_HARTREE)
    call hartree(lg,mg,info,system,fg,poisson,srg_scalar,stencil,rho,Vh)
    call timer_end(LOG_CALC_HARTREE)

    call timer_begin(LOG_CALC_EXC_COR)
    call exchange_correlation(system,xc_func,mg,srg_scalar,srg,rho_s,ppn,info,spsi,stencil,Vxc,energy%E_xc)
    call timer_end(LOG_CALC_EXC_COR)

  end if

end subroutine scf_iteration_step

end module scf_iteration_sub
