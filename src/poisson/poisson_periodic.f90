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
!===================================================================================================================================

module poisson_periodic
  implicit none

contains

subroutine poisson_ft(lg,mg,info,fg,rho,Vh,poisson)
  use structures
  use communication, only: comm_summation
  use math_constants, only : pi
  implicit none
  type(s_rgrid)          ,intent(in) :: lg
  type(s_rgrid)          ,intent(in) :: mg
  type(s_parallel_info)  ,intent(in) :: info
  type(s_reciprocal_grid),intent(in) :: fg
  type(s_scalar)         ,intent(in) :: rho
  type(s_scalar)                     :: Vh
  type(s_poisson)                    :: poisson
  !
  integer :: ix,iy,iz,kx,ky,kz

!$omp workshare
  poisson%ff1z = 0d0
!$omp end workshare

!$OMP parallel do private(iz,iy,ix)
  do iz=mg%is(3),mg%ie(3)
  do iy=mg%is(2),mg%ie(2)
  do ix=mg%is(1),mg%ie(1)
    poisson%ff1z(ix,iy,iz) = cmplx(rho%f(ix,iy,iz))
  end do
  end do
  end do

  call comm_summation(poisson%ff1z,poisson%ff2z,mg%num(1)*mg%num(2)*lg%num(3),info%icomm_z)

!$omp workshare
  poisson%ff1y = 0d0
!$omp end workshare

!$OMP parallel do private(kz,iy,ix)
  do kz = mg%is(3),mg%ie(3)
  do iy = mg%is(2),mg%ie(2)
  do ix = mg%is(1),mg%ie(1)
    poisson%ff1y(ix,iy,kz) = sum(fg%egzc(kz,:)*poisson%ff2z(ix,iy,:))
  end do
  end do
  end do
  call comm_summation(poisson%ff1y,poisson%ff2y,mg%num(1)*lg%num(2)*mg%num(3),info%icomm_y)

!$omp workshare
  poisson%ff1x = 0.d0
!$omp end workshare

!$OMP parallel do private(kz,ky,ix)
  do kz = mg%is(3),mg%ie(3)
  do ky = mg%is(2),mg%ie(2)
  do ix = mg%is(1),mg%ie(1)
    poisson%ff1x(ix,ky,kz) = sum(fg%egyc(ky,:)*poisson%ff2y(ix,:,kz))
  end do
  end do
  end do

  call comm_summation(poisson%ff1x,poisson%ff2x,lg%num(1)*mg%num(2)*mg%num(3),info%icomm_x)

!$OMP parallel do private(kz,ky,kx)
  do kz = mg%is(3),mg%ie(3)
  do ky = mg%is(2),mg%ie(2)
  do kx = mg%is(1),mg%ie(1)
    poisson%ff1x(kx,ky,kz) = sum(fg%egxc(kx,:)*poisson%ff2x(:,ky,kz))/dble(lg%num(1)*lg%num(2)*lg%num(3))
  end do
  end do
  end do

  call comm_summation(poisson%ff1x,poisson%ff2x,lg%num(1)*mg%num(2)*mg%num(3),info%icomm_x)

!$omp workshare
  poisson%ff1z = 0.d0
!$omp end workshare

!$OMP parallel do private(kz,ky,kx)
  do kz = mg%is(3),mg%ie(3)
  do ky = mg%is(2),mg%ie(2)
  do kx = mg%is(1),mg%ie(1)
    poisson%zrhoG_ele(kx,ky,kz) = poisson%ff2x(kx,ky,kz)
    poisson%ff1z(kx,ky,kz) = fg%coef(kx,ky,kz)*poisson%ff2x(kx,ky,kz)
  end do
  end do
  end do

  call comm_summation(poisson%ff1z,poisson%ff2z,mg%num(1)*mg%num(2)*lg%num(3),info%icomm_z)

!$omp workshare
  poisson%ff1y = 0.d0
!$omp end workshare

!$OMP parallel do private(iz,ky,kx)
  do iz = mg%is(3),mg%ie(3)
  do ky = mg%is(2),mg%ie(2)
  do kx = mg%is(1),mg%ie(1)
    poisson%ff1y(kx,ky,iz) = sum(fg%egz(:,iz)*poisson%ff2z(kx,ky,:))
  end do
  end do
  end do
  call comm_summation(poisson%ff1y,poisson%ff2y,mg%num(1)*lg%num(2)*mg%num(3),info%icomm_y)

!$omp workshare
  poisson%ff1x = 0.d0
!$omp end workshare

!$OMP parallel do private(iz,iy,kx)
  do iz = mg%is(3),mg%ie(3)
  do iy = mg%is(2),mg%ie(2)
  do kx = mg%is(1),mg%ie(1)
    poisson%ff1x(kx,iy,iz) = sum(fg%egy(:,iy)*poisson%ff2y(kx,:,iz))
  end do
  end do
  end do
  call comm_summation(poisson%ff1x,poisson%ff2x,lg%num(1)*mg%num(2)*mg%num(3),info%icomm_x)

!$OMP parallel do private(iz,iy,ix)
  do iz = mg%is(3),mg%ie(3)
  do iy = mg%is(2),mg%ie(2)
  do ix = mg%is(1),mg%ie(1)
    Vh%f(ix,iy,iz) = sum(fg%egx(:,ix)*poisson%ff2x(:,iy,iz))
  end do
  end do
  end do

  return
end subroutine poisson_ft

!===================================================================================================================================

subroutine poisson_ffte(lg,mg,info,fg,rho,Vh,poisson)
  use structures
  use communication, only: comm_summation
  use salmon_global, only: ffte_parallel
  implicit none
  type(s_rgrid)          ,intent(in) :: lg
  type(s_rgrid)          ,intent(in) :: mg
  type(s_parallel_info)  ,intent(in) :: info
  type(s_reciprocal_grid),intent(in) :: fg
  type(s_scalar)         ,intent(in) :: rho
  type(s_scalar)                     :: Vh
  type(s_poisson)                    :: poisson
  !
  integer :: ix,iy,iz
  integer :: iiy,iiz,iix
  real(8) :: inv_lgnum3

  inv_lgnum3=1.d0/(lg%num(1)*lg%num(2)*lg%num(3))

  select case (ffte_parallel)
  case ('xy')

  poisson%b_ffte=0.d0
!$OMP parallel do private(iix,iiy,ix) collapse(2)
  do ix=1,mg%num(1)
  do iy=1,mg%num(2)
    iix=ix+mg%is(1)-1
    iiy=iy+mg%is(2)-1
    poisson%b_ffte(mg%is(3):mg%ie(3),iy,ix) = cmplx(rho%f(iix,iiy,mg%is(3):mg%ie(3)))
  end do
  end do
  call comm_summation(poisson%b_ffte,poisson%a_ffte,size(poisson%a_ffte),info%icomm_z)

  CALL PZFFT3DV_MOD(poisson%a_ffte,poisson%b_ffte,lg%num(3),lg%num(2),lg%num(1),   &
                    info%isize_y,info%isize_x,-1, &
                    info%icomm_y,info%icomm_x)

  poisson%zrhoG_ele=0d0
!$omp parallel do collapse(2) default(none) &
!$omp             private(iz,iy,ix,iiy,iiz,iix) &
!$omp             shared(mg,lg,poisson,inv_lgnum3,fg)
  do ix=1,mg%num(1)
  do iy=1,mg%num(2)
    iix=ix+mg%is(1)-1
    iiy=iy+mg%is(2)-1
    do iz=1,mg%num(3)
      iiz=iz+mg%is(3)-1
      poisson%zrhoG_ele(iix,iiy,iiz) = poisson%b_ffte(iiz,iy,ix)*inv_lgnum3
    end do
    do iz=1,lg%num(3)
      poisson%b_ffte(iz,iy,ix) = poisson%b_ffte(iz,iy,ix) * fg%coef(iix,iiy,iz)
    end do
  end do
  end do
!$omp end parallel do

  CALL PZFFT3DV_MOD(poisson%b_ffte,poisson%a_ffte,lg%num(3),lg%num(2),lg%num(1), &
                    info%isize_y,info%isize_x,1, &
                    info%icomm_y,info%icomm_x)

!$OMP parallel do private(iix,iiy) collapse(2)
  do ix=1,mg%num(1)
  do iy=1,mg%num(2)
    iix=ix+mg%is(1)-1
    iiy=iy+mg%is(2)-1
    Vh%f(iix,iiy,mg%is(3):mg%ie(3)) = poisson%a_ffte(mg%is(3):mg%ie(3),iy,ix)
  end do
  end do

  case ('yz')

  poisson%b_ffte=0.d0
!$OMP parallel do private(iiz,iiy,ix) collapse(2)
  do iz=1,mg%num(3)
  do iy=1,mg%num(2)
    iiz=iz+mg%is(3)-1
    iiy=iy+mg%is(2)-1
    poisson%b_ffte(mg%is(1):mg%ie(1),iy,iz) = cmplx(rho%f(mg%is(1):mg%ie(1),iiy,iiz))
  end do
  end do
  call comm_summation(poisson%b_ffte,poisson%a_ffte,size(poisson%a_ffte),info%icomm_x)

  CALL PZFFT3DV_MOD(poisson%a_ffte,poisson%b_ffte,lg%num(1),lg%num(2),lg%num(3),   &
                    info%isize_y,info%isize_z,-1, &
                    info%icomm_y,info%icomm_z)

  poisson%zrhoG_ele=0d0
!$omp parallel do collapse(2) default(none) &
!$omp             private(iz,iy,ix,iiy,iiz,iix) &
!$omp             shared(mg,lg,poisson,inv_lgnum3,fg)
  do iz=1,mg%num(3)
  do iy=1,mg%num(2)
    iiz=iz+mg%is(3)-1
    iiy=iy+mg%is(2)-1
    do ix=1,mg%num(1)
      iix=ix+mg%is(1)-1
      poisson%zrhoG_ele(iix,iiy,iiz) = poisson%b_ffte(iix,iy,iz)*inv_lgnum3
    end do
    do ix=1,lg%num(1)
      poisson%b_ffte(ix,iy,iz) = poisson%b_ffte(ix,iy,iz) * fg%coef(ix,iiy,iiz)
    end do
  end do
  end do
!$omp end parallel do

  CALL PZFFT3DV_MOD(poisson%b_ffte,poisson%a_ffte,lg%num(1),lg%num(2),lg%num(3), &
                    info%isize_y,info%isize_z,1, &
                    info%icomm_y,info%icomm_z)

!$OMP parallel do private(iiz,iiy) collapse(2)
  do iz=1,mg%num(3)
  do iy=1,mg%num(2)
    iiz=iz+mg%is(3)-1
    iiy=iy+mg%is(2)-1
    Vh%f(mg%is(1):mg%ie(1),iiy,iiz) = poisson%a_ffte(mg%is(1):mg%ie(1),iy,iz)
  end do
  end do

  end select ! ffte_parallel

  return
end subroutine poisson_ffte

end module poisson_periodic
