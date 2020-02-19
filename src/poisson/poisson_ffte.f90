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
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------
module poisson_ffte_sub
  implicit none

contains

subroutine poisson_ffte(lg,mg,ng,info_field,trho,tvh,trhoG_ele,hgs,poisson)
  use structures, only: s_rgrid,s_field_parallel,s_reciprocal_grid,s_poisson
  use communication, only: comm_summation
  use communication, only: comm_is_root
  use math_constants, only : pi
!$  use omp_lib
  implicit none
  type(s_rgrid),intent(in) :: lg
  type(s_rgrid),intent(in) :: mg
  type(s_rgrid),intent(in) :: ng
  type(s_field_parallel),intent(in) :: info_field
  real(8),intent(in)       :: hgs(3)
  type(s_poisson),intent(inout)         :: poisson
  integer :: ix,iy,iz
  integer :: iix,iiy,iiz
  real(8) :: trho(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
  real(8) :: tvh(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
  complex(8) :: trhoG_ele(lg%num(1),lg%num(2),lg%num(3))
  real(8) :: inv_lgnum3
  complex(8),parameter :: zI=(0.d0,1.d0)
  real(8) :: bLx,bLy,bLz

  if(.not.allocated(poisson%coef) .or. &
     .not.allocated(poisson%a_ffte) .or. &
     .not.allocated(poisson%b_ffte) .or. &
     .not.allocated(poisson%a_ffte_tmp))then
    stop 'poisson_ffte: array is not allocated'
  end if

  bLx=2.d0*Pi/(Hgs(1)*dble(lg%num(1)))
  bLy=2.d0*Pi/(Hgs(2)*dble(lg%num(2)))
  bLz=2.d0*Pi/(Hgs(3)*dble(lg%num(3)))

  inv_lgnum3=1.d0/(lg%num(1)*lg%num(2)*lg%num(3))

  if(info_field%isize_ffte(1)==1)then
!$OMP parallel do private(iiz,iiy)
    do iz=1,ng%num(3)
    do iy=1,ng%num(2)
      iiz=iz+ng%is(3)-1
      iiy=iy+ng%is(2)-1
      poisson%a_ffte(1:lg%num(1),iy,iz)=trho(1:lg%num(1),iiy,iiz)
    end do
    end do
  else
    poisson%a_ffte_tmp=0.d0
!$OMP parallel do private(iiz,iiy,ix)
    do iz=1,ng%num(3)
    do iy=1,ng%num(2)
      iiz=iz+ng%is(3)-1
      iiy=iy+ng%is(2)-1
      do ix=1,ng%num(1)
        iix=ix+ng%is(1)-1
        poisson%a_ffte_tmp(ix,iy,iz)=trho(iix,iiy,iiz)
      end do
    end do
    end do
    call comm_summation(poisson%a_ffte_tmp,poisson%a_ffte,  &
                        lg%num(1)*ng%num(2)*ng%num(3),info_field%icomm_ffte(1))
  end if

  CALL PZFFT3DV_MOD(poisson%a_ffte,poisson%b_ffte,lg%num(1),lg%num(2),lg%num(3),   &
                    info_field%isize_ffte(2),info_field%isize_ffte(3),0, &
                    info_field%icomm_ffte(2),info_field%icomm_ffte(3))
  CALL PZFFT3DV_MOD(poisson%a_ffte,poisson%b_ffte,lg%num(1),lg%num(2),lg%num(3),   &
                    info_field%isize_ffte(2),info_field%isize_ffte(3),-1, &
                    info_field%icomm_ffte(2),info_field%icomm_ffte(3))

!$omp parallel do collapse(2) default(none) &
!$omp             private(iz,iy,ix,iiy,iiz) &
!$omp             shared(ng,lg,trhoG_ele,poisson,info_field,inv_lgnum3)
  do iz=1,ng%num(3)
  do iy=1,ng%num(2)
  do ix=1,lg%num(1)
    iiz=iz+ng%is(3)-1
    iiy=iy+ng%is(2)-1
    if (ix == 1 .and. iiy == 1 .and. iiz == 1) then
      trhoG_ele(ix,iiy,iiz)=0d0
    else
      trhoG_ele(ix,iiy,iiz)=poisson%b_ffte(ix,iy,iz)*inv_lgnum3
    end if
    poisson%b_ffte(ix,iy,iz)=poisson%b_ffte(ix,iy,iz)*poisson%coef(ix,iy,iz)
  end do
  end do
  end do
!$omp end parallel do

  CALL PZFFT3DV_MOD(poisson%b_ffte,poisson%a_ffte,lg%num(1),lg%num(2),lg%num(3), &
                    info_field%isize_ffte(2),info_field%isize_ffte(3),1, &
                    info_field%icomm_ffte(2),info_field%icomm_ffte(3))

  if(info_field%isize_ffte(1)==1)then
!$OMP parallel do private(iiz,iiy)
    do iz=1,ng%num(3)
    do iy=1,ng%num(2)
      iiz=iz+ng%is(3)-1
      iiy=iy+ng%is(2)-1
      tvh(1:lg%num(1),iiy,iiz)=poisson%a_ffte(1:lg%num(1),iy,iz)
    end do
    end do
  else
!$OMP parallel do private(iiz,iiy,ix)
    do iz=1,ng%num(3)
    do iy=1,ng%num(2)
      iiz=iz+ng%is(3)-1
      iiy=iy+ng%is(2)-1
      do ix=1,ng%num(1)
        iix=ix+ng%is(1)-1
        tvh(iix,iiy,iiz)=poisson%a_ffte(ix,iy,iz)
      end do
    end do
    end do
  end if

  return
end subroutine poisson_ffte
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------

end module poisson_ffte_sub
