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
subroutine inner_product4(mg,info,zmatbox1,zmatbox2,zbox2,hvol)
  use structures, only: s_rgrid, s_orbital_parallel
  use salmon_communication, only: comm_summation
  !$ use omp_lib
  implicit none
  type(s_rgrid),intent(in) :: mg
  type(s_orbital_parallel),intent(in) :: info
  complex(8),intent(in)  :: zmatbox1(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
  complex(8),intent(in)  :: zmatbox2(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
  complex(8),intent(out) :: zbox2
  real(8),intent(in) :: hvol
  integer :: ix,iy,iz
  complex(8) :: zbox
  
  zbox=0.d0
  !$omp parallel do private(iz,iy,ix) reduction(+ : zbox)
  do iz=mg%is(3),mg%ie(3)
  do iy=mg%is(2),mg%ie(2)
  do ix=mg%is(1),mg%ie(1)
    zbox=zbox+conjg(zmatbox1(ix,iy,iz))*zmatbox2(ix,iy,iz)
  end do
  end do
  end do
  zbox=zbox*hvol
  call comm_summation(zbox,zbox2,info%icomm_r)
  
end subroutine inner_product4
