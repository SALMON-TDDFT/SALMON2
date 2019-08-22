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

subroutine scf_iteration(mg,ng,system,info,stencil,srg,spsi,shpsi,srho,srho_s,itotmst,mst, &
               cg,ppg,vlocal,  &
               iflag_diisjump,energy, &
               norm_diff_psi_stock,  &
               miter,iditerybcg,   &
               iflag_subspace_diag,iditer_nosubspace_diag,ifmst,mixing,iter)
  use inputoutput, only: iperiodic,method_min,method_mixing,mixrate
  use structures
  use timer
  use rmmdiis_sub
  use gram_schmidt_orth, only: gram_schmidt
  use Conjugate_Gradient, only: gscg_isolated,gscg_periodic
  use subspace_diagonalization, only: ssdg_isolated,ssdg_periodic
  use density_matrix, only: calc_density
  use mixing_sub
  implicit none

  type(s_rgrid),         intent(in)    :: mg
  type(s_rgrid),         intent(in)    :: ng
  type(s_dft_system),    intent(in)    :: system
  type(s_orbital_parallel),intent(in)  :: info
  type(s_orbital),       intent(inout) :: spsi,shpsi
  type(s_scalar),        intent(inout) :: srho
  type(s_scalar),        intent(inout) :: srho_s(system%nspin)
  type(s_stencil),       intent(in)    :: stencil
  type(s_sendrecv_grid), intent(inout) :: srg
  type(s_pp_grid),       intent(in)    :: ppg
  integer,               intent(in)    :: itotmst
  integer,               intent(in)    :: mst(2)
  type(s_cg),            intent(inout) :: cg
  type(s_scalar),        intent(in)    :: vlocal(system%nspin)
  integer,               intent(inout) :: iflag_diisjump
  type(s_dft_energy),    intent(inout) :: energy
  real(8),               intent(out)   :: norm_diff_psi_stock(itotmst,1)
  integer,               intent(in)    :: miter
  integer,               intent(in)    :: iditerybcg
  integer,               intent(in)    :: iflag_subspace_diag
  integer,               intent(in)    :: iditer_nosubspace_diag
  integer,               intent(in)    :: ifmst(2)
  type(s_mixing),        intent(inout) :: mixing
  integer,               intent(in)    :: iter
  integer                              :: j

! solve Kohn-Sham equation by minimization techniques
  call timer_begin(LOG_CALC_MINIMIZATION)

  if( method_min == 'cg' .or.       &
    ( method_min == 'cg-diis' .and. Miter <= iDiterYBCG) ) then
    select case(iperiodic)
    case(0)
      call gscg_isolated(mg,system,info,stencil,ppg,vlocal,srg,spsi,cg)
    case(3)
      call gscg_periodic(mg,system,info,stencil,ppg,vlocal,srg,spsi,cg)
    end select
  else if( method_min  == 'diis' .or. method_min == 'cg-diis' ) then
    select case(iperiodic)
    case(0)
      stop "rmmdiis method is not implemented."
!      call rmmdiis(mg,system,info,stencil,srg_ob_1,spsi,energy,itotmst,mst,   &
!                   iflag_diisjump, &
!                   norm_diff_psi_stock,info_ob,ppg,vlocal)
    case(3)
      stop "rmmdiis method is not implemented for periodic systems."
    end select
  end if

  call timer_end(LOG_CALC_MINIMIZATION)

! Gram Schmidt orghonormalization
  call gram_schmidt(system, mg, info, spsi)

! subspace diagonalization
  if(iflag_subspace_diag==1)then
    if(miter>iditer_nosubspace_diag)then
      select case(iperiodic)
      case(0)      
        call ssdg_isolated(mg,system,info,stencil,spsi,shpsi,ppg,vlocal,srg)

      case(3)
        call ssdg_periodic(mg,system,info,stencil,spsi,shpsi,ppg,vlocal,srg)
      end select
    end if
  end if
  
! density
  call timer_begin(LOG_CALC_RHO)

  call calc_density(srho_s,spsi,info,mg,system%nspin)

  select case(method_mixing)
    case ('simple') ; call simple_mixing(ng,system,1.d0-mixrate,mixrate,srho_s,mixing)
    case ('broyden'); call wrapper_broyden(ng,system,srho_s,mst,ifmst,iter,mixing)
  end select
  call timer_end(LOG_CALC_RHO)

  srho%f = 0d0
  do j=1,system%nspin
    srho%f = srho%f + srho_s(j)%f
  end do

end subroutine scf_iteration

end module scf_iteration_sub
