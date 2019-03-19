!
!  Copyright 2017 SALMON developers
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
SUBROUTINE Hartree_ns(lg,mg,ng,Brl)
use structures, only: s_rgrid
use hartree_periodic_sub
use hartree_ffte_sub
use scf_data
use new_world_sub
use allocate_mat_sub
implicit none
type(s_rgrid),intent(in) :: lg
type(s_rgrid),intent(in) :: mg
type(s_rgrid),intent(in) :: ng
real(8),intent(in),optional :: Brl(3,3)

if(iSCFRT==1)then
  select case(iperiodic)
  case(0)
    call Hartree_cg(lg,mg,ng,rho,Vh)
  case(3)
    select case(iflag_hartree)
    case(2)
      call Hartree_periodic(lg,mg,ng,rho,Vh,hgs,iscfrt,itcalc_ene,itt,   &
                 ff1,ff1x,ff1y,ff1z,ff2,ff2x,ff2y,ff2z,rhoe_g_tmp,rhoe_g,trho2z,trho3z, &
                 egx,egxc,egy,egyc,egz,egzc)
    case(4)
      call Hartree_FFTE(lg,mg,ng,rho,Vh,icheck_ascorder,hgs,npuw,npuy,npuz,   &
                        a_ffte,b_ffte,rhoe_g,coef_poisson,matbox_l,matbox_l2)
    end select
  end select
else if(iSCFRT==2)then
  select case(iperiodic)
  case(0)
    if(mod(itt,2)==1)then
      call Hartree_cg(lg,mg,ng,rho,Vh_stock2)
    else
      call Hartree_cg(lg,mg,ng,rho,Vh_stock1)
    end if
  case(3)
    if(mod(itt,2)==1)then
      select case(iflag_hartree)
      case(2)
        call Hartree_periodic(lg,mg,ng,rho,Vh_stock2,hgs,iscfrt,itcalc_ene,itt,   &
                 ff1,ff1x,ff1y,ff1z,ff2,ff2x,ff2y,ff2z,rhoe_g_tmp,rhoe_g,trho2z,trho3z, &
                 egx,egxc,egy,egyc,egz,egzc,Brl)
      case(4)
        call Hartree_FFTE(lg,mg,ng,rho,Vh_stock2,icheck_ascorder,hgs,npuw,npuy,npuz,   &
                          a_ffte,b_ffte,rhoe_g,coef_poisson,matbox_l,matbox_l2)
      end select
    else
      select case(iflag_hartree)
      case(2)
        call Hartree_periodic(lg,mg,ng,rho,Vh_stock1,hgs,iscfrt,itcalc_ene,itt,   &
                 ff1,ff1x,ff1y,ff1z,ff2,ff2x,ff2y,ff2z,rhoe_g_tmp,rhoe_g,trho2z,trho3z, &
                 egx,egxc,egy,egyc,egz,egzc,Brl)
      case(4)
        call Hartree_FFTE(lg,mg,ng,rho,Vh_stock1,icheck_ascorder,hgs,npuw,npuy,npuz,   &
                          a_ffte,b_ffte,rhoe_g,coef_poisson,matbox_l,matbox_l2)
      end select
    end if
  end select
end if

return

END SUBROUTINE Hartree_ns
