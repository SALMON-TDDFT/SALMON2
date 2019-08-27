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
MODULE gradient_sub

use scf_data
use structures, only: s_rgrid, s_sendrecv_grid
use sendrecv_grid, only: update_overlap_complex8, update_overlap_real8
use pack_unpack, only: copy_data


implicit none 
INTERFACE calc_gradient

  MODULE PROCEDURE R_calc_gradient,C_calc_gradient

END INTERFACE

CONTAINS

!=======================================================================
!=======================================================================

SUBROUTINE R_calc_gradient(mg, srg, is_array_wk, ie_array_wk, wk, grad_wk)
!$ use omp_lib

implicit none
type(s_rgrid), intent(in) :: mg
type(s_sendrecv_grid), intent(inout) :: srg
integer, intent(in) :: is_array_wk(1:3)
integer, intent(in) :: ie_array_wk(1:3)
real(8), intent(in) :: wk( &
  & is_array_wk(1):ie_array_wk(1), &
  & is_array_wk(2):ie_array_wk(2), &
  & is_array_wk(3):ie_array_wk(3))
real(8), intent(out) :: grad_wk(3,  &
  & is_array_wk(1):ie_array_wk(1), &
  & is_array_wk(2):ie_array_wk(2), &
  & is_array_wk(3):ie_array_wk(3))

real(8) :: tmp( &
  & mg%is_array(1):mg%ie_array(1), &
  & mg%is_array(2):mg%ie_array(2), &
  & mg%is_array(3):mg%ie_array(3))
real(8) :: grad_tmp(3, &
  & mg%is_array(1):mg%ie_array(1), &
  & mg%is_array(2):mg%ie_array(2), &
  & mg%is_array(3):mg%ie_array(3))

integer :: ix,iy,iz

call copy_data( &
  & wk( &
    & is_array_wk(1):ie_array_wk(1), &
    & is_array_wk(2):ie_array_wk(2), &
    & is_array_wk(3):ie_array_wk(3)), &
  & tmp( &
    & is_array_wk(1):ie_array_wk(1), &
    & is_array_wk(2):ie_array_wk(2), &
    & is_array_wk(3):ie_array_wk(3)))
  
call update_overlap_real8(srg, mg, tmp)

if(Nd<=3)then
!$OMP parallel do private(iz,iy,ix) 
  do iz=mg%is(3),mg%ie(3)
  do iy=mg%is(2),mg%ie(2)
  do ix=mg%is(1),mg%ie(1)
    grad_tmp(1,ix,iy,iz)=      &
      (tmp(ix+1,iy,iz)-tmp(ix-1,iy,iz))/2.d0/Hgs(1)
    grad_tmp(2,ix,iy,iz)=      &
      (tmp(ix,iy+1,iz)-tmp(ix,iy-1,iz))/2.d0/Hgs(2)
    grad_tmp(3,ix,iy,iz)=      &
      (tmp(ix,iy,iz+1)-tmp(ix,iy,iz-1))/2.d0/Hgs(3)
  end do
  end do
  end do
else if(Nd==4)then
!$OMP parallel do private(iz,iy,ix)
  do iz=mg%is(3),mg%ie(3)
  do iy=mg%is(2),mg%ie(2)
  do ix=mg%is(1),mg%ie(1)
    grad_tmp(1,ix,iy,iz)=      &
      (bN1*(tmp(ix+1,iy,iz)-tmp(ix-1,iy,iz))      &
      +bN2*(tmp(ix+2,iy,iz)-tmp(ix-2,iy,iz))      &
      +bN3*(tmp(ix+3,iy,iz)-tmp(ix-3,iy,iz))      &
      +bN4*(tmp(ix+4,iy,iz)-tmp(ix-4,iy,iz)))/Hgs(1) 
    grad_tmp(2,ix,iy,iz)=      &
      (bN1*(tmp(ix,iy+1,iz)-tmp(ix,iy-1,iz))      &
      +bN2*(tmp(ix,iy+2,iz)-tmp(ix,iy-2,iz))      &
      +bN3*(tmp(ix,iy+3,iz)-tmp(ix,iy-3,iz))      &
      +bN4*(tmp(ix,iy+4,iz)-tmp(ix,iy-4,iz)))/Hgs(2) 
    grad_tmp(3,ix,iy,iz)=      &
      (bN1*(tmp(ix,iy,iz+1)-tmp(ix,iy,iz-1))      &
      +bN2*(tmp(ix,iy,iz+2)-tmp(ix,iy,iz-2))      &
      +bN3*(tmp(ix,iy,iz+3)-tmp(ix,iy,iz-3))      &
      +bN4*(tmp(ix,iy,iz+4)-tmp(ix,iy,iz-4)))/Hgs(3) 
  end do
  end do
  end do
end if

call copy_data( &
  & grad_tmp(1:3, &
    & is_array_wk(1):ie_array_wk(1), &
    & is_array_wk(2):ie_array_wk(2), &
    & is_array_wk(3):ie_array_wk(3)), &
  & grad_wk(1:3, &
    & is_array_wk(1):ie_array_wk(1), &
    & is_array_wk(2):ie_array_wk(2), &
    & is_array_wk(3):ie_array_wk(3)))

return

END SUBROUTINE R_calc_gradient

!=======================================================================

SUBROUTINE C_calc_gradient(mg, srg, is_array_wk, ie_array_wk, wk, grad_wk)
!$ use omp_lib

implicit none
type(s_rgrid), intent(in) :: mg
type(s_sendrecv_grid), intent(inout) :: srg
integer, intent(in) :: is_array_wk(1:3)
integer, intent(in) :: ie_array_wk(1:3)
complex(8), intent(in) :: wk( &
  & is_array_wk(1):ie_array_wk(1), &
  & is_array_wk(2):ie_array_wk(2), &
  & is_array_wk(3):ie_array_wk(3))
complex(8), intent(out) :: grad_wk(3,  &
  & is_array_wk(1):ie_array_wk(1), &
  & is_array_wk(2):ie_array_wk(2), &
  & is_array_wk(3):ie_array_wk(3))
complex(8) :: tmp( &
  & mg%is_array(1):mg%ie_array(1), &
  & mg%is_array(2):mg%ie_array(2), &
  & mg%is_array(3):mg%ie_array(3))
complex(8) :: grad_tmp(1:3, &
  & mg%is_array(1):mg%ie_array(1), &
  & mg%is_array(2):mg%ie_array(2), &
  & mg%is_array(3):mg%ie_array(3))
integer :: ix, iy, iz

  
  call copy_data( &
  & wk( &
    & is_array_wk(1):ie_array_wk(1), &
    & is_array_wk(2):ie_array_wk(2), &
    & is_array_wk(3):ie_array_wk(3)), &
  & tmp( &
    & is_array_wk(1):ie_array_wk(1), &
    & is_array_wk(2):ie_array_wk(2), &
    & is_array_wk(3):ie_array_wk(3)))
  
call update_overlap_complex8(srg, mg, tmp)

if(Nd<=3)then
!$OMP parallel do private(iz,iy,ix) 
  do iz=mg%is(3),mg%ie(3)
  do iy=mg%is(2),mg%ie(2)
  do ix=mg%is(1),mg%ie(1)
    grad_tmp(1,ix,iy,iz)=      &
      (tmp(ix+1,iy,iz)-tmp(ix-1,iy,iz))/2.d0/Hgs(1)  
    grad_tmp(2,ix,iy,iz)=      &
      (tmp(ix,iy+1,iz)-tmp(ix,iy-1,iz))/2.d0/Hgs(2)  
    grad_tmp(3,ix,iy,iz)=      &
      (tmp(ix,iy,iz+1)-tmp(ix,iy,iz-1))/2.d0/Hgs(3)  
  end do
  end do
  end do
else if(Nd==4)then
!$OMP parallel do private(iz,iy,ix) 
  do iz=mg%is(3),mg%ie(3)
  do iy=mg%is(2),mg%ie(2)
  do ix=mg%is(1),mg%ie(1)
    grad_tmp(1,ix,iy,iz)=      &
      (bN1*(tmp(ix+1,iy,iz)-tmp(ix-1,iy,iz))      &
      +bN2*(tmp(ix+2,iy,iz)-tmp(ix-2,iy,iz))      &
      +bN3*(tmp(ix+3,iy,iz)-tmp(ix-3,iy,iz))      &
      +bN4*(tmp(ix+4,iy,iz)-tmp(ix-4,iy,iz)))/Hgs(1) 
    grad_tmp(2,ix,iy,iz)=      &
      (bN1*(tmp(ix,iy+1,iz)-tmp(ix,iy-1,iz))      &
      +bN2*(tmp(ix,iy+2,iz)-tmp(ix,iy-2,iz))      &
      +bN3*(tmp(ix,iy+3,iz)-tmp(ix,iy-3,iz))      &
      +bN4*(tmp(ix,iy+4,iz)-tmp(ix,iy-4,iz)))/Hgs(2) 
    grad_tmp(3,ix,iy,iz)=      &
      (bN1*(tmp(ix,iy,iz+1)-tmp(ix,iy,iz-1))      &
      +bN2*(tmp(ix,iy,iz+2)-tmp(ix,iy,iz-2))      &
      +bN3*(tmp(ix,iy,iz+3)-tmp(ix,iy,iz-3))      &
      +bN4*(tmp(ix,iy,iz+4)-tmp(ix,iy,iz-4)))/Hgs(3) 
  end do
  end do
  end do
end if

call copy_data( &
  & grad_tmp(1:3, &
    & is_array_wk(1):ie_array_wk(1), &
    & is_array_wk(2):ie_array_wk(2), &
    & is_array_wk(3):ie_array_wk(3)), &
  & grad_wk(1:3, &
    & is_array_wk(1):ie_array_wk(1), &
    & is_array_wk(2):ie_array_wk(2), &
    & is_array_wk(3):ie_array_wk(3)))

return

END SUBROUTINE C_calc_gradient

END MODULE gradient_sub


