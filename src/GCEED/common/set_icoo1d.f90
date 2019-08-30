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
subroutine set_icoo1d(lg)
  use structures, only: s_rgrid
  use scf_data
!$ use omp_lib
  implicit none
  type(s_rgrid), intent(in) :: lg
  integer :: ix,iy,iz
  integer :: icount

!$OMP parallel do private(iz,iy,ix,icount) 
  do iz=lg%is(3),lg%ie(3)
  do iy=lg%is(2),lg%ie(2)
  do ix=lg%is(1),lg%ie(1)
    icount=(iz-lg%is(3))*lg%num(2)*lg%num(1)+(iy-lg%is(2))*lg%num(1)+ix-lg%is(1)+1
    icoo1d(1,icount)=ix
    icoo1d(2,icount)=iy
    icoo1d(3,icount)=iz
  end do
  end do
  end do

end subroutine set_icoo1d 
