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
SUBROUTINE Hartree_ns(lg,mg,ng,info_field,system,poisson_cg,srg_ng,stencil,srho,sVh,fg)
  use structures, only: s_rgrid,s_dft_system,s_field_parallel,s_poisson_cg,  &
                        s_sendrecv_grid,s_stencil,s_scalar,s_reciprocal_grid
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
  type(s_field_parallel),intent(in) :: info_field
  type(s_dft_system),intent(in) :: system
  type(s_poisson_cg),intent(in) :: poisson_cg
  type(s_sendrecv_grid),intent(inout) :: srg_ng
  type(s_stencil),intent(in) :: stencil
  type(s_scalar) ,intent(in) :: srho
  type(s_scalar)             :: sVh
  type(s_reciprocal_grid)    :: fg
  !
  integer :: ix,iy,iz,n,nn

  select case(iperiodic)
  case(0)
    call Hartree_cg(lg,mg,ng,info_field,system,poisson_cg,srho%f,sVh%f,srg_ng,stencil,itervh,   &
                    igc_is,igc_ie,gridcoo,iflag_ps)
  case(3)
    select case(iflag_hartree)
    case(2)
      call Hartree_periodic(lg,mg,ng,system,info_field,srho,sVh,fg)
    case(4)
      call Hartree_FFTE(lg,mg,ng,srho%f,sVh%f,hgs,npuw,npuy,npuz,   &
                        a_ffte,b_ffte,rhoe_g,coef_poisson)
!$omp parallel do collapse(2) default(none) &
!$omp             private(iz,iy,ix,n,nn) &
!$omp             shared(lg,ng,fg,rhoe_G,NPUZ,NPUY)
      do iz=1,lg%num(3)/NPUZ
      do iy=1,lg%num(2)/NPUY
      do ix=ng%is(1)-lg%is(1)+1,ng%ie(1)-lg%is(1)+1
        n=(iz-1)*lg%num(2)/NPUY*lg%num(1)+(iy-1)*lg%num(1)+ix
        nn=ix-(ng%is(1)-lg%is(1)+1)+1+(iy-1)*ng%num(1)+(iz-1)*lg%num(2)/NPUY*ng%num(1)+fg%ig_s-1
        fg%zrhoG_ele(nn) = rhoe_G(n)
      enddo
      enddo
      enddo
!$omp end parallel do
    end select
  end select

return

END SUBROUTINE Hartree_ns
