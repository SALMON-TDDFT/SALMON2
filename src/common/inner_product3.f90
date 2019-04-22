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
subroutine inner_product3(mg,rmatbox1,rmatbox2,rbox2)
  use structures, only: s_rgrid
  use salmon_parallel, only: nproc_group_korbital
  use salmon_communication, only: comm_summation
  use timer
  !$ use omp_lib
  implicit none
  type(s_rgrid),intent(in) :: mg
  real(8),intent(in) :: rmatbox1(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
  real(8),intent(in) :: rmatbox2(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
  real(8),intent(out) :: rbox2
  integer :: ix,iy,iz
  real(8) :: rbox
  
  rbox=0.d0
  !$omp parallel do reduction(+ : rbox) private(iz,iy,ix)
  do iz=mg%is(3),mg%ie(3)
  do iy=mg%is(2),mg%ie(2)
  do ix=mg%is(1),mg%ie(1)
    rbox=rbox+rmatbox1(ix,iy,iz)*rmatbox2(ix,iy,iz)
  end do
  end do
  end do
  
  call timer_begin(LOG_ALLREDUCE_INNER_PRODUCT3)
  call comm_summation(rbox,rbox2,nproc_group_korbital)
  call timer_end(LOG_ALLREDUCE_INNER_PRODUCT3)
  
end subroutine inner_product3
