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

subroutine scf_iteration(mg,system,info,stencil,srg_ob_1,spsi,iflag,itotmst,mst,ilsda,nproc_ob,iparaway_ob, &
               num_kpoints_rd,k_rd,   &
               rxk_ob,rhxk_ob,rgk_ob,rpk_ob,   &
               zxk_ob,zhxk_ob,zgk_ob,zpk_ob,zpko_ob,zhtpsi_ob,   &
               info_ob,ppg,vlocal,  &
               iflag_diisjump,energy, &
               norm_diff_psi_stock,  &
               miter,iditerybcg,   &
               iflag_subspace_diag,iditer_nosubspace_diag,cnmat,bnmat,iobnum,ifmst,k_sta,k_end)
  use inputoutput, only: iperiodic,ispin,amin_routine,gscg
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
  implicit none

  type(s_rgrid),         intent(in)    :: mg
  type(s_system),        intent(in)    :: system
  type(s_wf_info),       intent(in)    :: info
  type(s_wavefunction),  intent(inout) :: spsi
  type(s_stencil),       intent(in)    :: stencil
  type(s_sendrecv_grid), intent(inout) :: srg_ob_1
  type(s_pp_grid),       intent(in)    :: ppg
  integer,               intent(inout) :: iflag
  integer,               intent(in)    :: itotmst
  integer,               intent(in)    :: mst(2)
  integer,               intent(in)    :: ilsda
  integer,               intent(in)    :: nproc_ob
  integer,               intent(in)    :: iparaway_ob
  integer,               intent(in)    :: num_kpoints_rd
  real(8),               intent(in)    :: k_rd(3,num_kpoints_rd)
  real(8),               intent(inout) :: rxk_ob(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),1:system%nspin*info%numo)
  real(8),               intent(inout) :: rhxk_ob(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),1:system%nspin*info%numo)
  real(8),               intent(inout) :: rgk_ob(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),1:system%nspin*info%numo)
  real(8),               intent(inout) :: rpk_ob(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),1:system%nspin*info%numo)
  complex(8),            intent(inout)   :: zxk_ob(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),1:system%nspin*info%numo)
  complex(8),            intent(inout)   :: zhxk_ob(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),1:system%nspin*info%numo)
  complex(8),            intent(inout)   :: zgk_ob(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),1:system%nspin*info%numo)
  complex(8),            intent(inout)   :: zpk_ob(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),1:system%nspin*info%numo)
  complex(8),            intent(inout)   :: zpko_ob(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),1:system%nspin*info%numo)
  complex(8),            intent(inout)   :: zhtpsi_ob(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),1:system%nspin*info%numo)
  type(s_wf_info),       intent(in)    :: info_ob
  real(8),               intent(in)    :: vlocal(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),ispin+1)
  integer,               intent(inout) :: iflag_diisjump
  type(s_energy),        intent(inout) :: energy
  real(8),               intent(out)   :: norm_diff_psi_stock(itotmst,1)
  integer,               intent(in)    :: miter
  integer,               intent(in)    :: iditerybcg
  integer,               intent(in)    :: iflag_subspace_diag
  integer,               intent(in)    :: iditer_nosubspace_diag
  real(8),               intent(in)    :: cnmat(0:12,0:12),bnmat(0:12,0:12)
  integer,               intent(in)    :: iobnum
  integer,               intent(in)    :: ifmst(2)
  integer,               intent(in)    :: k_sta,k_end

! solve Kohn-Sham equation by minimization techniques
  call timer_begin(LOG_CALC_MINIMIZATION)

  if( amin_routine == 'cg' .or.       &
    ( amin_routine == 'cg-diis' .and. Miter <= iDiterYBCG) ) then
    select case(iperiodic)
    case(0)
      select case(gscg)
      case('y')
        call sgscg(mg,system,info,stencil,srg_ob_1,spsi,iflag,itotmst,mst,ilsda,nproc_ob,iparaway_ob, &
                   rxk_ob,rhxk_ob,rgk_ob,rpk_ob,   &
                   info_ob,ppg,vlocal)
      case('n')
        call dtcg(mg,system,info,stencil,srg_ob_1,spsi,iflag,itotmst,mst,ilsda,nproc_ob,iparaway_ob,   &
                  info_ob,ppg,vlocal)
      end select
    case(3)
      select case(gscg)
      case('y')
        call gscg_periodic(mg,system,info,stencil,srg_ob_1,spsi,iflag,itotmst,mst,ilsda,nproc_ob,iparaway_ob,   &
                           zxk_ob,zhxk_ob,zgk_ob,zpk_ob,zpko_ob,zhtpsi_ob,   &
                           info_ob,ppg,vlocal,num_kpoints_rd,k_rd)
      case('n')
        call dtcg_periodic(mg,system,info,stencil,srg_ob_1,spsi,iflag,itotmst,mst,ilsda,nproc_ob,iparaway_ob,   &
                           info_ob,ppg,vlocal,num_kpoints_rd,k_rd)
      end select
    end select
  else if( amin_routine  == 'diis' .or. amin_routine == 'cg-diis' ) then
    select case(iperiodic)
    case(0)
      call rmmdiis(mg,system,info,stencil,srg_ob_1,spsi,energy,itotmst,mst,   &
                   iflag_diisjump, &
                   norm_diff_psi_stock,info_ob,ppg,vlocal,iparaway_ob)
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
        call subspace_diag(mg,info,stencil,srg_ob_1,spsi,ilsda,nproc_ob,iparaway_ob,iobnum,itotmst,k_sta,k_end,   &
                           mst,ifmst,system%hvol,info_ob,bnmat,cnmat,system%hgs,ppg,vlocal)

      case(3)
        call subspace_diag_periodic(mg,info,stencil,srg_ob_1,spsi,ilsda,nproc_ob,iparaway_ob,  &
                                    iobnum,itotmst,k_sta,k_end,mst,ifmst,system%hvol,   &
                                    info_ob,bnmat,cnmat,system%hgs,ppg,vlocal,num_kpoints_rd,k_rd)
      end select
    end if
  end if
  
end subroutine scf_iteration

end module scf_iteration_sub
