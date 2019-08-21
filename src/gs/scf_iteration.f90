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

subroutine scf_iteration(mg,ng,system,info,stencil,srg,srg_ob_1,spsi,shpsi,srho,srho_s,iflag,itotmst,mst,ilsda,nproc_ob, &
               cg,   &
               info_ob,ppg,vlocal,  &
               iflag_diisjump,energy, &
               norm_diff_psi_stock,  &
               miter,iditerybcg,   &
               iflag_subspace_diag,iditer_nosubspace_diag,iobnum,ifmst,mixing,iter)
  use inputoutput, only: iperiodic,method_min,gscg,method_mixing,mixrate
  use structures
  use timer
  use dtcg_sub
  use gscg_sub
  use dtcg_periodic_sub
  use gscg_periodic_sub
  use rmmdiis_sub
  use gram_schmidt_orth, only: gram_schmidt 
  use subspace_diag_sub
  use subspace_diag_periodic_sub
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
  type(s_sendrecv_grid), intent(inout) :: srg,srg_ob_1
  type(s_pp_grid),       intent(in)    :: ppg
  integer,               intent(inout) :: iflag
  integer,               intent(in)    :: itotmst
  integer,               intent(in)    :: mst(2)
  integer,               intent(in)    :: ilsda
  integer,               intent(in)    :: nproc_ob
  type(s_cg),            intent(inout) :: cg
  type(s_orbital_parallel),intent(in)  :: info_ob
  type(s_scalar),        intent(in)    :: vlocal(system%nspin)
  integer,               intent(inout) :: iflag_diisjump
  type(s_dft_energy),    intent(inout) :: energy
  real(8),               intent(out)   :: norm_diff_psi_stock(itotmst,1)
  integer,               intent(in)    :: miter
  integer,               intent(in)    :: iditerybcg
  integer,               intent(in)    :: iflag_subspace_diag
  integer,               intent(in)    :: iditer_nosubspace_diag
  integer,               intent(in)    :: iobnum
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
      select case(gscg)
      case('y')
        call sgscg(mg,system,info,stencil,srg_ob_1,spsi,iflag,itotmst,mst,ilsda,nproc_ob,cg, &
                   info_ob,ppg,vlocal)
      case('n')
        call dtcg(mg,system,info,stencil,srg_ob_1,spsi,iflag,itotmst,mst,ilsda,nproc_ob,   &
                  info_ob,ppg,vlocal)
      end select
    case(3)
      select case(gscg)
      case('y')
        call gscg_periodic(mg,system,info,stencil,ppg,vlocal,srg,spsi,iflag,cg)
      case('n')
        call dtcg_periodic(mg,system,info,stencil,srg_ob_1,spsi,iflag,itotmst,mst,ilsda,nproc_ob,   &
                           info_ob,ppg,vlocal)
      end select
    end select
  else if( method_min  == 'diis' .or. method_min == 'cg-diis' ) then
    select case(iperiodic)
    case(0)
      call rmmdiis(mg,system,info,stencil,srg_ob_1,spsi,energy,itotmst,mst,   &
                   iflag_diisjump, &
                   norm_diff_psi_stock,info_ob,ppg,vlocal)
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
        call subspace_diag(mg,system,info,stencil,srg_ob_1,spsi,ilsda,nproc_ob,iobnum,itotmst,   &
                           mst,ifmst,info_ob,ppg,vlocal)

      case(3)
        call subspace_diag_periodic(mg,system,info,stencil,spsi,shpsi,ppg,vlocal,srg)
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
