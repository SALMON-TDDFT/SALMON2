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

subroutine solve_orbitals(mg,system,info,stencil,spsi,shpsi,srg,cg,ppg,vlocal,  &
            &   miter,nscf_init_no_diagonal)
  use salmon_global, only: yn_subspace_diagonalization,ncg,ncg_init
  use structures
  use timer
  use gram_schmidt_orth, only: gram_schmidt
  use Conjugate_Gradient, only: gscg_zwf,gscg_rwf
  use subspace_diagonalization, only: ssdg
  implicit none
  type(s_rgrid),          intent(in)    :: mg
  type(s_dft_system),     intent(in)    :: system
  type(s_parallel_info),  intent(in)    :: info
  type(s_stencil),        intent(in)    :: stencil
  type(s_orbital),        intent(inout) :: spsi,shpsi
  type(s_sendrecv_grid),  intent(inout) :: srg
  type(s_pp_grid),        intent(in)    :: ppg
  type(s_cg),             intent(inout) :: cg
  type(s_scalar),         intent(in)    :: vlocal(system%nspin)
  integer,                intent(in)    :: miter
  integer,                intent(in)    :: nscf_init_no_diagonal
  !
  integer :: nncg

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

end subroutine solve_orbitals

subroutine update_density_and_potential(lg,mg,system,info,stencil,xc_func,pp,ppn,iter, &
               spsi,srg,srg_scalar,poisson,fg,rho,rho_s,rho_jm,Vpsl,Vh,Vxc,vlocal,mixing,energy )
  use structures
  use salmon_global, only: method_mixing,yn_jm,yn_spinorbit
  use timer
  use mixing_sub
  use hartree_sub, only: hartree
  use salmon_xc, only: exchange_correlation
  use noncollinear_module, only: simple_mixing_so
  use hamiltonian, only: update_vlocal
  implicit none
  type(s_rgrid),          intent(in)    :: lg,mg
  type(s_dft_system),     intent(in)    :: system
  type(s_parallel_info),  intent(in)    :: info
  type(s_stencil),        intent(in)    :: stencil
  type(s_xc_functional),  intent(in)    :: xc_func
  type(s_pp_info),        intent(in)    :: pp
  type(s_pp_nlcc),        intent(in)    :: ppn
  integer,                intent(in)    :: iter
  type(s_orbital),        intent(inout) :: spsi
  type(s_sendrecv_grid),  intent(inout) :: srg,srg_scalar
  type(s_poisson),        intent(inout) :: poisson
  type(s_reciprocal_grid),intent(inout) :: fg
  type(s_scalar),         intent(inout) :: rho,rho_s(system%nspin)
  type(s_scalar),         intent(in)    :: rho_jm,Vpsl
  type(s_scalar),         intent(inout) :: Vh,Vxc(system%nspin),vlocal(system%nspin)
  type(s_mixing),         intent(inout) :: mixing
  type(s_dft_energy),     intent(inout) :: energy
  !
  integer :: j

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
  case ('simple_potential')
    ! Nothing is done here since Hartree and XC potentials are mixed instead of density
  case default
    stop 'Invalid method_mixing. Specify any one of "simple" or "broyden" or "pulay" for method_mixing.'
  end select

  rho%f = 0d0
  do j=1,system%nspin
    rho%f = rho%f + rho_s(j)%f
  end do

  if(yn_jm=='y') rho%f = rho%f + rho_jm%f

  call timer_begin(LOG_CALC_HARTREE)
  call hartree(lg,mg,info,system,fg,poisson,srg_scalar,stencil,rho,Vh)
  call timer_end(LOG_CALC_HARTREE)

    call timer_begin(LOG_CALC_EXC_COR)
    call exchange_correlation(system,xc_func,mg,srg_scalar,srg,rho_s,pp,ppn,info,spsi,stencil,Vxc,energy%E_xc)
    call timer_end(LOG_CALC_EXC_COR)


  if(method_mixing=='simple_potential')then
    call simple_mixing_potential(mg,system,1.d0-mixing%mixrate,mixing%mixrate,Vh,Vxc,mixing)
  end if
  
  call update_vlocal(mg,system%nspin,Vh,Vpsl,Vxc,Vlocal)
  
end subroutine update_density_and_potential

end module scf_iteration_sub
