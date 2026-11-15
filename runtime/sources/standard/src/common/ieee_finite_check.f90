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
module ieee_finite_check
  implicit none

contains

  ! nvfortran's -acc front end can misresolve the ieee_is_finite intrinsic to
  ! a device-only symbol from ordinary host code, depending on the calling
  ! subroutine's shape. Callers import this under the name ieee_is_finite
  ! (use ieee_finite_check, only: ieee_is_finite => is_finite) so every call
  ! site -- scalar or, by elemental broadcast, any array shape -- is unchanged.
  elemental logical function is_finite(x) result(r)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    real(8), intent(in) :: x
    r = ieee_is_finite(x)
  end function is_finite

end module ieee_finite_check
