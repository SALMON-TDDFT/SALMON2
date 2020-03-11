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
module poisson_periodic_sub
  implicit none

contains

subroutine poisson_periodic(lg,ng,info_field,fg,rho,Vh,poisson)
  use structures
  use communication, only: comm_summation
  use math_constants, only : pi
  implicit none
  type(s_rgrid)          ,intent(in) :: lg
  type(s_rgrid)          ,intent(in) :: ng
  type(s_field_parallel) ,intent(in) :: info_field
  type(s_reciprocal_grid),intent(in) :: fg
  type(s_scalar)         ,intent(in) :: rho
  type(s_scalar)                     :: Vh
  type(s_poisson)                    :: poisson
  !
  integer :: ix,iy,iz,kx,ky,kz

!$omp workshare
  poisson%trho2z = 0d0
!$omp end workshare

!$OMP parallel do private(iz,iy,ix)
  do iz=ng%is(3),ng%ie(3)
  do iy=ng%is(2),ng%ie(2)
  do ix=ng%is(1),ng%ie(1)
    poisson%trho2z(ix,iy,iz) = rho%f(ix,iy,iz)
  end do
  end do
  end do

  call comm_summation(poisson%trho2z,poisson%trho3z,ng%num(1)*ng%num(2)*lg%num(3),info_field%icomm(3))

!$omp workshare
  poisson%ff1y = 0d0
!$omp end workshare

!$OMP parallel do private(kz,iy,ix)
  do kz = ng%is(3),ng%ie(3)
  do iy = ng%is(2),ng%ie(2)
  do ix = ng%is(1),ng%ie(1)
    poisson%ff1y(ix,iy,kz) = sum(fg%egzc(kz,:)*poisson%trho3z(ix,iy,:))
  end do
  end do
  end do
  call comm_summation(poisson%ff1y,poisson%ff2y,ng%num(1)*lg%num(2)*ng%num(3),info_field%icomm(2))

!$omp workshare
  poisson%ff1x = 0.d0
!$omp end workshare

!$OMP parallel do private(kz,ky,ix)
  do kz = ng%is(3),ng%ie(3)
  do ky = ng%is(2),ng%ie(2)
  do ix = ng%is(1),ng%ie(1)
    poisson%ff1x(ix,ky,kz) = sum(fg%egyc(ky,:)*poisson%ff2y(ix,:,kz))
  end do
  end do
  end do

  call comm_summation(poisson%ff1x,poisson%ff2x,lg%num(1)*ng%num(2)*ng%num(3),info_field%icomm(1))

!$OMP parallel do private(kz,ky,kx)
  do kz = ng%is(3),ng%ie(3)
  do ky = ng%is(2),ng%ie(2)
  do kx = ng%is(1),ng%ie(1)
    poisson%ff1x(kx,ky,kz) = sum(fg%egxc(kx,:)*poisson%ff2x(:,ky,kz))/dble(lg%num(1)*lg%num(2)*lg%num(3))
  end do
  end do
  end do

  call comm_summation(poisson%ff1x,poisson%ff2x,lg%num(1)*ng%num(2)*ng%num(3),info_field%icomm(1))

!$omp workshare
  poisson%ff1z = 0.d0
!$omp end workshare

!$OMP parallel do private(kz,ky,kx)
  do kz = ng%is(3),ng%ie(3)
  do ky = ng%is(2),ng%ie(2)
  do kx = ng%is(1),ng%ie(1)
    poisson%zrhoG_ele(kx,ky,kz) = poisson%ff2x(kx,ky,kz)
    poisson%ff1z(kx,ky,kz) = fg%coef(kx,ky,kz)*poisson%ff2x(kx,ky,kz)
  end do
  end do
  end do

  call comm_summation(poisson%ff1z,poisson%ff2z,ng%num(1)*ng%num(2)*lg%num(3),info_field%icomm(3))

!$omp workshare
  poisson%ff1y = 0.d0
!$omp end workshare

!$OMP parallel do private(iz,ky,kx)
  do iz = ng%is(3),ng%ie(3)
  do ky = ng%is(2),ng%ie(2)
  do kx = ng%is(1),ng%ie(1)
    poisson%ff1y(kx,ky,iz) = sum(fg%egz(:,iz)*poisson%ff2z(kx,ky,:))
  end do
  end do
  end do
  call comm_summation(poisson%ff1y,poisson%ff2y,ng%num(1)*lg%num(2)*ng%num(3),info_field%icomm(2))

!$omp workshare
  poisson%ff1x = 0.d0
!$omp end workshare

!$OMP parallel do private(iz,iy,kx)
  do iz = ng%is(3),ng%ie(3)
  do iy = ng%is(2),ng%ie(2)
  do kx = ng%is(1),ng%ie(1)
    poisson%ff1x(kx,iy,iz) = sum(fg%egy(:,iy)*poisson%ff2y(kx,:,iz))
  end do
  end do
  end do
  call comm_summation(poisson%ff1x,poisson%ff2x,lg%num(1)*ng%num(2)*ng%num(3),info_field%icomm(1))

!$OMP parallel do private(iz,iy,ix)
  do iz = ng%is(3),ng%ie(3)
  do iy = ng%is(2),ng%ie(2)
  do ix = ng%is(1),ng%ie(1)
    Vh%f(ix,iy,iz) = sum(fg%egx(:,ix)*poisson%ff2x(:,iy,iz))
  end do
  end do
  end do

  return
end subroutine poisson_periodic
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------

end module poisson_periodic_sub
