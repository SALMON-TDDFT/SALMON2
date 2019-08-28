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
use sendrecv_grid, only: update_overlap_real8, update_overlap_complex8, dealloc_cache
use pack_unpack, only: copy_data
implicit none 
INTERFACE calc_gradient

  MODULE PROCEDURE R_calc_gradient,C_calc_gradient

END INTERFACE

CONTAINS

!=======================================================================
!=======================================================================

SUBROUTINE R_calc_gradient(mg, srg_ob_1, is, ie, wk, grad_wk)
  !$ use omp_lib
  implicit none
  type(s_rgrid), intent(in) :: mg
  type(s_sendrecv_grid), intent(inout) :: srg_ob_1
  integer, intent(in) :: is(3), ie(3)
  real(8), intent(in) :: wk( &
    & is(1):ie(1), &
    & is(2):ie(2), &
    & is(3):ie(3))
  real(8), intent(out) :: grad_wk( &
    & 1:3, &
    & is(1):ie(1), &
    & is(2):ie(2), &
    & is(3):ie(3))

  real(8) :: tmp(&
    & mg%is_array(1):mg%ie_array(1), &
    & mg%is_array(2):mg%ie_array(2), &
    & mg%is_array(3):mg%ie_array(3))
  real(8) :: grad_tmp( &
    & 1:3, &
    & mg%is_array(1):mg%ie_array(1), &
    & mg%is_array(2):mg%ie_array(2), &
    & mg%is_array(3):mg%ie_array(3))
  integer :: ix,iy,iz

  tmp = 0d0

  call copy_data( &
    & wk( &
      & is(1):ie(1), &
      & is(2):ie(2), &
      & is(3):ie(3)), &
    & tmp(&
      & is(1):ie(1), &
      & is(2):ie(2), &
      & is(3):ie(3)))

  call update_overlap_real8(srg_ob_1, mg, tmp)

  if(Nd<=3)then
    !$OMP parallel do private(iz,iy,ix) 
      do iz=is(3),ie(3)
      do iy=is(2),ie(2)
      do ix=is(1),ie(1)
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
      do iz=is(3),ie(3)
      do iy=is(2),ie(2)
      do ix=is(1),ie(1)
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
    & grad_tmp( &
      & 1:3, &
      & is(1):ie(1), &
      & is(2):ie(2), &
      & is(3):ie(3)), &
    & grad_wk( &
      & 1:3, &
      & is(1):ie(1), &
      & is(2):ie(2), &
      & is(3):ie(3)))

  return

END SUBROUTINE R_calc_gradient

!=======================================================================

SUBROUTINE C_calc_gradient(mg, srg_ob_1, is, ie, wk, grad_wk)
!$ use omp_lib
  implicit none
  type(s_rgrid), intent(in) :: mg
  type(s_sendrecv_grid), intent(inout) :: srg_ob_1
  integer, intent(in) :: is(3), ie(3)
  complex(8), intent(in) :: wk( &
    & is(1):ie(1), &
    & is(2):ie(2), &
    & is(3):ie(3))
  complex(8), intent(out) :: grad_wk( &
    & 1:3, &
    & is(1):ie(1), &
    & is(2):ie(2), &
    & is(3):ie(3))
  complex(8) :: tmp(&
    & mg%is_array(1):mg%ie_array(1), &
    & mg%is_array(2):mg%ie_array(2), &
    & mg%is_array(3):mg%ie_array(3))
  complex(8) :: grad_tmp( &
    & 1:3, &
    & mg%is_array(1):mg%ie_array(1), &
    & mg%is_array(2):mg%ie_array(2), &
    & mg%is_array(3):mg%ie_array(3))
  integer :: ix,iy,iz

  tmp = 0d0

  call copy_data( &
    & wk( &
      & is(1):ie(1), &
      & is(2):ie(2), &
      & is(3):ie(3)), &
    & tmp(&
      & is(1):ie(1), &
      & is(2):ie(2), &
      & is(3):ie(3)))
    
  call update_overlap_complex8(srg_ob_1, mg, tmp)

  if(Nd<=3)then
    !$OMP parallel do private(iz,iy,ix) 
      do iz=is(3),ie(3)
      do iy=is(2),ie(2)
      do ix=is(1),ie(1)
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
      do iz=is(3),ie(3)
      do iy=is(2),ie(2)
      do ix=is(1),ie(1)
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
    & grad_tmp( &
      & 1:3, &
      & is(1):ie(1), &
      & is(2):ie(2), &
      & is(3):ie(3)), &
    & grad_wk( &
      & 1:3, &
      & is(1):ie(1), &
      & is(2):ie(2), &
      & is(3):ie(3)))

  return

END SUBROUTINE C_calc_gradient

END MODULE gradient_sub


