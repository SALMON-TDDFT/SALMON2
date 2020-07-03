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
!============================ Hartree potential (Solve Poisson equation)
module hartree_sub
  implicit none

contains

!===================================================================================================================================
subroutine hartree(lg,mg,info,system,fg,poisson,srg_scalar,stencil,rho,Vh)
  use inputoutput, only: iperiodic,yn_ffte
  use structures, only: s_rgrid,s_dft_system,s_parallel_info,s_poisson,  &
                        s_sendrecv_grid,s_stencil,s_scalar,s_reciprocal_grid
  use poisson_isolated
  use poisson_periodic
  implicit none
  type(s_rgrid)          ,intent(in)    :: lg
  type(s_rgrid)          ,intent(in)    :: mg
  type(s_parallel_info)  ,intent(in)    :: info
  type(s_dft_system)     ,intent(in)    :: system
  type(s_reciprocal_grid),intent(in)    :: fg
  type(s_poisson)        ,intent(inout) :: poisson
  type(s_sendrecv_grid)  ,intent(inout) :: srg_scalar
  type(s_stencil)        ,intent(in)    :: stencil
  type(s_scalar)         ,intent(in)    :: rho
  type(s_scalar)         ,intent(inout) :: Vh

  select case(iperiodic)
  case(0)
    call poisson_cg(lg,mg,info,system,poisson,rho%f,Vh%f,srg_scalar,stencil)
  case(3)
    select case(yn_ffte)
    case('n')
      call poisson_ft(lg,mg,info,fg,rho,Vh,poisson)
    case('y')
      call poisson_ffte(lg,mg,info,fg,rho,Vh,poisson)
    end select
  end select

return

end subroutine hartree

end module hartree_sub
