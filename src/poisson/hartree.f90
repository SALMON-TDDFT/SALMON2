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
subroutine hartree_ns(lg,mg,ng,info_field,system,poisson,srg_ng,stencil,srho,sVh,fg)
  use inputoutput, only: iperiodic,yn_ffte
  use structures, only: s_rgrid,s_dft_system,s_field_parallel,s_poisson,  &
                        s_sendrecv_grid,s_stencil,s_scalar,s_reciprocal_grid
  use poisson_cg_sub
  use poisson_periodic_sub
  use poisson_ffte_sub
  implicit none
  type(s_rgrid),intent(in) :: lg
  type(s_rgrid),intent(in) :: mg
  type(s_rgrid),intent(in) :: ng
  type(s_field_parallel),intent(in) :: info_field
  type(s_dft_system),intent(in) :: system
  type(s_poisson),intent(inout) :: poisson
  type(s_sendrecv_grid),intent(inout) :: srg_ng
  type(s_stencil),intent(in) :: stencil
  type(s_scalar) ,intent(in) :: srho
  type(s_scalar)             :: sVh
  type(s_reciprocal_grid)    :: fg

  select case(iperiodic)
  case(0)
    call poisson_cg(lg,mg,ng,info_field,system,poisson,srho%f,sVh%f,srg_ng,stencil)
  case(3)
    select case(yn_ffte)
    case('n')
      call poisson_periodic(lg,mg,ng,system,info_field,srho,sVh,fg)
    case('y')
      call poisson_FFTE(lg,mg,ng,info_field,srho%f,sVh%f,system%hgs,fg,poisson)
    end select
  end select

return

end subroutine hartree_ns
