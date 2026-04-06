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

subroutine solve_orbitals(mg,system,info,stencil,spsi,shpsi,sttpsi,srg,cg,ppg,vlocal,  &
            &   miter,nscf_init_no_diagonal)
  use salmon_global, only: yn_subspace_diagonalization,ncg,ncg_init
  use structures
  use timer
  use gram_schmidt_orth, only: gram_schmidt
  use Conjugate_Gradient, only: gscg_zwf,gscg_rwf
  use subspace_diagonalization, only: ssdg
  implicit none
  type(s_rgrid),          intent(in)    :: mg
  type(s_dft_system),     intent(inout) :: system
  type(s_parallel_info),  intent(in)    :: info
  type(s_stencil),        intent(in)    :: stencil
  type(s_orbital),        intent(inout) :: spsi,shpsi,sttpsi
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
  
! subspace diagonalization
  call timer_begin(LOG_CALC_SUBSPACE_DIAG)
  if(yn_subspace_diagonalization == 'y')then
    if(miter > nscf_init_no_diagonal)then
      call ssdg(mg,system,info,stencil,spsi,shpsi,ppg,vlocal,srg)
    end if
  end if
  call timer_end(LOG_CALC_SUBSPACE_DIAG)

! conjugate gradient method
  call timer_begin(LOG_CALC_MINIMIZATION)
  if(system%if_real_orbital) then
    call gscg_rwf(nncg,mg,system,info,stencil,ppg,vlocal,srg,spsi,shpsi,sttpsi,cg)
  else
    call gscg_zwf(nncg,mg,system,info,stencil,ppg,vlocal,srg,spsi,shpsi,sttpsi,cg)
  end if
  call timer_end(LOG_CALC_MINIMIZATION)

! Gram Schmidt orghonormalization
  call gram_schmidt(system, mg, info, spsi)

end subroutine solve_orbitals

subroutine update_density_and_potential(lg,mg,system,info,stencil,xc_func,pp,ppn,iter, &
               spsi,srg,srg_scalar,poisson,fg,rho,rho_s,rho_jm,Vpsl,Vh,Vxc,vlocal,mixing,energy )
  use structures
  use salmon_global, only: method_mixing,yn_jm,yn_spinorbit,yn_dc
  use timer
  use mixing_sub
  use hartree_sub, only: hartree
  use salmon_xc, only: exchange_correlation, calc_tau_from_orbitals, salmon_xctype_tbmbj
  use noncollinear_module, only: simple_mixing_so
  use hamiltonian, only: update_vlocal
  implicit none
  type(s_rgrid),          intent(in)    :: lg,mg
  type(s_dft_system),     intent(inout) :: system
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
  logical :: do_tau_mixing
  logical :: do_j_mixing
  logical :: aux_work_ready

  do_tau_mixing = mixing%use_aux_tau .and. xc_func%use_kinetic_energy
  do_j_mixing = mixing%use_aux_j .and. xc_func%xctype(1) == salmon_xctype_tbmbj
  aux_work_ready = .false.

  if (mixing%use_aux_mixing .and. (mixing%tau_mixrate > 0d0 .or. mixing%j_mixrate > 0d0) .and. .not. xc_func%use_kinetic_energy) then
    stop "aux mixing requested but active XC does not use auxiliary fields"
  end if

  if (do_tau_mixing .or. do_j_mixing) then
    if (system%nspin /= 1) stop "aux mixing currently supports only nspin=1"
    call ensure_aux_mixing_storage(mg, mixing)
    if (method_mixing == 'broyden' .or. method_mixing == 'pulay') then
      if (do_tau_mixing .and. do_j_mixing) then
        call calc_tau_from_orbitals(system,mg,info,srg,stencil,spsi,tau=mixing%tau_work%f,rj=mixing%aux_work%j%v)
      else if (do_j_mixing) then
        call calc_tau_from_orbitals(system,mg,info,srg,stencil,spsi,rj=mixing%aux_work%j%v)
      else
        call calc_tau_from_orbitals(system,mg,info,srg,stencil,spsi,tau=mixing%tau_work%f)
      end if
      mixing%aux_work%use_tau = do_tau_mixing
      if (do_tau_mixing) mixing%aux_work%tau%f = mixing%tau_work%f
      mixing%aux_work%use_j = do_j_mixing
      if (do_j_mixing) then
        do j=1,3
          mixing%j_work(j)%f = mixing%aux_work%j%v(j,:,:,:)
        end do
      end if
      aux_work_ready = .true.
    end if
  end if

  select case(method_mixing)
  case ('simple')
    call simple_mixing(mg,system,1.d0-mixing%mixrate,mixing%mixrate,rho_s,mixing)
  case ('simple_dm')
    if(yn_spinorbit=='n') stop 'yn_spinorbit must be y when method_mixing=simple_dm'
    call simple_mixing_so(mg,system,1.d0-mixing%mixrate,mixing%mixrate,rho_s,mixing)
  case ('broyden')
    if (do_tau_mixing .and. do_j_mixing) then
      call wrapper_broyden(info%icomm_r,mg,system,rho_s,tau=mixing%tau_work,j=mixing%j_work,iter=iter,mixing=mixing)
    else if (do_tau_mixing) then
      call wrapper_broyden(info%icomm_r,mg,system,rho_s,tau=mixing%tau_work,iter=iter,mixing=mixing)
    else if (do_j_mixing) then
      call wrapper_broyden(info%icomm_r,mg,system,rho_s,j=mixing%j_work,iter=iter,mixing=mixing)
    else
      call wrapper_broyden(info%icomm_r,mg,system,rho_s,iter=iter,mixing=mixing)
    end if
  case ('pulay')
    if (do_tau_mixing .and. do_j_mixing) then
      call pulay(mg,info,system,rho_s,tau=mixing%tau_work,j=mixing%j_work,iter=iter,mixing=mixing)
    else if (do_tau_mixing) then
      call pulay(mg,info,system,rho_s,tau=mixing%tau_work,iter=iter,mixing=mixing)
    else if (do_j_mixing) then
      call pulay(mg,info,system,rho_s,j=mixing%j_work,iter=iter,mixing=mixing)
    else
      call pulay(mg,info,system,rho_s,iter=iter,mixing=mixing)
    end if
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

  if ((do_tau_mixing .or. do_j_mixing) .and. .not. aux_work_ready) then
    if (do_tau_mixing .and. do_j_mixing) then
      call calc_tau_from_orbitals(system,mg,info,srg,stencil,spsi,tau=mixing%tau_work%f,rj=mixing%aux_work%j%v)
    else if (do_j_mixing) then
      call calc_tau_from_orbitals(system,mg,info,srg,stencil,spsi,rj=mixing%aux_work%j%v)
    else
      call calc_tau_from_orbitals(system,mg,info,srg,stencil,spsi,tau=mixing%tau_work%f)
    end if
    if (do_tau_mixing) call simple_mixing_tau(mg,1.d0-mixing%tau_mixrate,mixing%tau_mixrate,mixing%tau_work,mixing)
    if (do_j_mixing) call simple_mixing_j(mg,1.d0-mixing%j_mixrate,mixing%j_mixrate,mixing%j_work,mixing)
    mixing%aux_work%use_tau = do_tau_mixing
    if (do_tau_mixing) mixing%aux_work%tau%f = mixing%tau_work%f
    mixing%aux_work%use_j = do_j_mixing
    if (do_j_mixing) then
      do j=1,3
        mixing%j_work(j)%f = mixing%aux_work%j%v(j,:,:,:)
      end do
    end if
    aux_work_ready = .true.
  end if
  
  if(yn_dc=='n') then

    call timer_begin(LOG_CALC_EXC_COR)
    if (do_tau_mixing .or. do_j_mixing) then
      mixing%aux_work%use_tau = do_tau_mixing
      if (do_tau_mixing) mixing%aux_work%tau%f = mixing%tau_work%f
      mixing%aux_work%use_j = do_j_mixing
      if (do_j_mixing) then
        do j=1,3
          mixing%aux_work%j%v(j,:,:,:) = mixing%j_work(j)%f
        end do
      end if
      call exchange_correlation(system,xc_func,mg,srg_scalar,srg,rho_s,pp,ppn,info,spsi,stencil,Vxc,energy%E_xc, &
                                aux_override=mixing%aux_work)
    else
      call exchange_correlation(system,xc_func,mg,srg_scalar,srg,rho_s,pp,ppn,info,spsi,stencil,Vxc,energy%E_xc)
    end if
    call timer_end(LOG_CALC_EXC_COR)

    if(method_mixing=='simple_potential')then
      call simple_mixing_potential(mg,system,1.d0-mixing%mixrate,mixing%mixrate,Vh,Vxc,mixing)
    end if
    
    call update_vlocal(mg,system%nspin,Vh,Vpsl,Vxc,Vlocal)
    
  end if

end subroutine update_density_and_potential

end module scf_iteration_sub
