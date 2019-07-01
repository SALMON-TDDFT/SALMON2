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
subroutine sendrecv_copy(mg,tpsi)
!$ use omp_lib
use scf_data
use structures, only: s_rgrid
implicit none
integer :: ix,iy,iz,iob,iik
type(s_rgrid),intent(in) :: mg
real(8) :: tpsi(mg%is_overlap(1):mg%ie_overlap(1) &
&              ,mg%is_overlap(2):mg%ie_overlap(2) &
&              ,mg%is_overlap(3):mg%ie_overlap(3), 1:iobnum, k_sta:k_end)

!$omp parallel default(none) &
!$omp          shared(k_sta,k_end,iobnum,mg,tpsi,psi) &
!$omp          private(iik,iob,iz,iy,ix)
do iik=k_sta,k_end
do iob=1,iobnum
!$omp do collapse(2)
  do iz=mg%is_overlap(3),mg%ie_overlap(3)
  do iy=mg%is_overlap(2),mg%ie_overlap(2)
  do ix=mg%is_overlap(1),mg%ie_overlap(1)
    tpsi(ix,iy,iz,iob,iik)=0.d0
  end do
  end do
  end do
!$omp end do
!$omp do collapse(2)
  do iz=mg%is(3),mg%ie(3)
  do iy=mg%is(2),mg%ie(2)
  do ix=mg%is(1),mg%ie(1)
    tpsi(ix,iy,iz,iob,iik)=psi(ix,iy,iz,iob,iik)
  end do
  end do
  end do
!$omp end do
end do
end do
!$omp end parallel

end subroutine sendrecv_copy
