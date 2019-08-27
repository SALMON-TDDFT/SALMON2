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
subroutine set_gridcoo(lg)
use structures, only: s_rgrid
use scf_data
!$use omp_lib
implicit none
type(s_rgrid),intent(inout) :: lg
integer :: ix,iy,iz

allocate(lg%coordinate(minval(lg%is_overlap(1:3)):maxval(lg%ie_overlap(1:3)),3))

select case(iperiodic)
case(0)
  select case(imesh_oddeven(1))
    case(1)
!$OMP parallel do
      do ix=lg%is_overlap(1),lg%ie_overlap(1)
        lg%coordinate(ix,1)=dble(ix)*Hgs(1)
      end do
    case(2)
!$OMP parallel do
      do ix=lg%is_overlap(1),lg%ie_overlap(1)
        lg%coordinate(ix,1)=(dble(ix)-0.5d0)*Hgs(1)
      end do
  end select

  select case(imesh_oddeven(2))
    case(1)
!$OMP parallel do
      do iy=lg%is_overlap(2),lg%ie_overlap(2)
        lg%coordinate(iy,2)=dble(iy)*Hgs(2)
      end do
    case(2)
!$OMP parallel do
    do iy=lg%is_overlap(2),lg%ie_overlap(2)
      lg%coordinate(iy,2)=(dble(iy)-0.5d0)*Hgs(2)
    end do
  end select
  
  select case(imesh_oddeven(3))
    case(1)
!$OMP parallel do
      do iz=lg%is_overlap(3),lg%ie_overlap(3)
        lg%coordinate(iz,3)=dble(iz)*Hgs(3)
      end do
    case(2)
!$OMP parallel do
      do iz=lg%is_overlap(3),lg%ie_overlap(3)
        lg%coordinate(iz,3)=(dble(iz)-0.5d0)*Hgs(3)
      end do
  end select
case(3)
!$OMP parallel do
  do ix=lg%is_overlap(1),lg%ie_overlap(1)
    lg%coordinate(ix,1)=dble(ix-1)*Hgs(1)
  end do
!$OMP parallel do
  do iy=lg%is_overlap(2),lg%ie_overlap(2)
    lg%coordinate(iy,2)=dble(iy-1)*Hgs(2)
  end do
!$OMP parallel do
  do iz=lg%is_overlap(3),lg%ie_overlap(3)
    lg%coordinate(iz,3)=dble(iz-1)*Hgs(3)
  end do
end select

end subroutine set_gridcoo
