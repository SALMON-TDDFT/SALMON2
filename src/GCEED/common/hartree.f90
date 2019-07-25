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
!============================ Hartree potential (Solve Poisson equation)
SUBROUTINE Hartree_ns(lg,mg,ng,Brl,srg_ng,stencil,srho,sVh)
  use structures, only: s_rgrid,s_sendrecv_grid,s_stencil,s_scalar
  use hartree_cg_sub
  use hartree_periodic_sub
  use hartree_ffte_sub
  use scf_data
  use new_world_sub
  use allocate_mat_sub
  implicit none
  type(s_rgrid),intent(in) :: lg
  type(s_rgrid),intent(in) :: mg
  type(s_rgrid),intent(in) :: ng
  real(8)      ,intent(in) :: Brl(3,3)
  type(s_sendrecv_grid),intent(inout) :: srg_ng
  type(s_stencil),intent(in) :: stencil
  type(s_scalar) ,intent(in) :: srho
  type(s_scalar)             :: sVh

  select case(iperiodic)
  case(0)
    call Hartree_cg(lg,mg,ng,srho%f,sVh%f,srg_ng,stencil,hconv,itervh,wkbound_h,wk2bound_h,   &
                    meo,lmax_meo,igc_is,igc_ie,gridcoo,hvol,iflag_ps,num_pole,inum_mxin_s,   &
                    iamax,maxval_pole,num_pole_myrank,icorr_polenum,icount_pole,icorr_xyz_pole, &
                    ibox_icoobox_bound,icoobox_bound)
  case(3)
    select case(iflag_hartree)
    case(2)
      call Hartree_periodic(lg,mg,ng,srho%f,sVh%f,hgs,   &
                 ff1,ff1x,ff1y,ff1z,ff2,ff2x,ff2y,ff2z,rhoe_g_tmp,rhoe_g,trho2z,trho3z, &
                 egx,egxc,egy,egyc,egz,egzc,Brl)
    case(4)
      call Hartree_FFTE(lg,mg,ng,srho%f,sVh%f,hgs,npuw,npuy,npuz,   &
                        a_ffte,b_ffte,rhoe_g,coef_poisson)
    end select
  end select

return

END SUBROUTINE Hartree_ns
