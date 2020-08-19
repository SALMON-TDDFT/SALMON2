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
module nondiagonal_so_sub

  use noncollinear_module, only: op_xc_noncollinear

  implicit none

  private
  public :: nondiagonal_so

contains

  subroutine nondiagonal_so( tpsi, hpsi, info, mg )
    use structures, only: s_parallel_info, s_rgrid, s_orbital
    implicit none
    type(s_parallel_info),intent(in) :: info
    type(s_rgrid),intent(in) :: mg
    type(s_orbital),intent(in) :: tpsi
    type(s_orbital),intent(inout) :: hpsi
    !write(*,*) "-------------- nondiagonal_so"
    call op_xc_noncollinear( tpsi, hpsi, info, mg )
    return
  end subroutine nondiagonal_so

end module nondiagonal_so_sub
