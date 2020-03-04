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
  implicit none
  type(s_rgrid),intent(in) :: lg
  type(s_rgrid),intent(in) :: mg
  type(s_rgrid),intent(in) :: ng
  type(s_field_parallel),intent(in) :: info_field
  real(8),intent(in)       :: hgs(3)
  type(s_poisson),intent(inout)         :: poisson
  integer :: ix,iy,iz
  integer :: iiy,iiz,iix
  real(8) :: trho(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
  real(8) :: tvh(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
  complex(8) :: trhoG_ele(lg%num(1),lg%num(2),lg%num(3))
  real(8) :: inv_lgnum3

  if(.not.allocated(poisson%coef) .or. &
     .not.allocated(poisson%a_ffte) .or. &
     .not.allocated(poisson%b_ffte) .or. &
     .not.allocated(poisson%a_ffte_tmp))then
    stop 'poisson_ffte: array is not allocated'
  end if

  inv_lgnum3=1.d0/(lg%num(1)*lg%num(2)*lg%num(3))

  poisson%a_ffte_tmp=0.d0
!$OMP parallel do private(iiz,iiy,ix) collapse(2)
  do iz=1,ng%num(3)
  do iy=1,ng%num(2)
    iiz=iz+ng%is(3)-1
    iiy=iy+ng%is(2)-1
    poisson%a_ffte_tmp(ng%is(1):ng%ie(1),iy,iz)=trho(ng%is(1):ng%ie(1),iiy,iiz)
  end do
  end do
  call comm_summation(poisson%a_ffte_tmp,poisson%a_ffte,size(poisson%a_ffte),info_field%icomm(1))

  CALL PZFFT3DV_MOD(poisson%a_ffte,poisson%b_ffte,lg%num(1),lg%num(2),lg%num(3),   &
                    info_field%isize(2),info_field%isize(3),-1, &
                    info_field%icomm(2),info_field%icomm(3))

  trhoG_ele=0d0
!$omp parallel do collapse(2) default(none) &
!$omp             private(iz,iy,ix,iiy,iiz,iix) &
!$omp             shared(ng,lg,trhoG_ele,poisson,inv_lgnum3)
  do iz=1,ng%num(3)
  do iy=1,ng%num(2)
    do ix=1,ng%num(1)
      iiz=iz+ng%is(3)-1
      iiy=iy+ng%is(2)-1
      iix=ix+ng%is(1)-1
      trhoG_ele(iix,iiy,iiz)=poisson%b_ffte(iix,iy,iz)*inv_lgnum3
    end do

    do ix=1,lg%num(1)
      poisson%b_ffte(ix,iy,iz)=poisson%b_ffte(ix,iy,iz)*poisson%coef(ix,iy,iz)
    end do
  end do
  end do
!$omp end parallel do

  CALL PZFFT3DV_MOD(poisson%b_ffte,poisson%a_ffte,lg%num(1),lg%num(2),lg%num(3), &
                    info_field%isize(2),info_field%isize(3),1, &
                    info_field%icomm(2),info_field%icomm(3))

!$OMP parallel do private(iiz,iiy) collapse(2)
  do iz=1,ng%num(3)
  do iy=1,ng%num(2)
    iiz=iz+ng%is(3)-1
    iiy=iy+ng%is(2)-1
    tvh(ng%is(1):ng%ie(1),iiy,iiz)=poisson%a_ffte(ng%is(1):ng%ie(1),iy,iz)
  end do
  end do

  return
end subroutine poisson_ffte
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------

end module poisson_ffte_sub
