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

subroutine poisson_periodic(lg,mg,ng,system,info_field,trho,tVh,trhoG_ele,trhoG_ele_tmp,poisson)
  use structures, only: s_rgrid, s_field_parallel, s_dft_system, &
                        s_scalar, s_reciprocal_grid, s_poisson
  use communication, only: comm_summation
  use math_constants, only : pi
  implicit none
  type(s_rgrid),intent(in) :: lg
  type(s_rgrid),intent(in) :: mg
  type(s_rgrid),intent(in) :: ng
  type(s_field_parallel),intent(in) :: info_field
  type(s_dft_system),intent(in) :: system
  real(8),intent(in)    :: trho(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
  real(8),intent(inout) :: tVh (mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
  complex(8)            :: trhoG_ele_tmp(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3))
  complex(8),intent(out):: trhoG_ele(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3))
  type(s_poisson) :: poisson
  !
  integer :: ix,iy,iz,kx,ky,kz,kkx,kky,kkz
  real(8) :: gx,gy,gz
  real(8) :: g2
  real(8) :: B(3,3)
  integer :: n

!$OMP parallel do private(iz,iy,ix)
  do iz=lg%is(3),lg%ie(3)
  do iy=lg%is(2),lg%ie(2)
  do ix=lg%is(1),lg%ie(1)
    poisson%ff1(ix,iy,iz)=0.d0
    trhoG_ele_tmp(ix,iy,iz)=0.d0
  end do
  end do
  end do

!$OMP parallel do private(iz,iy,ix)
  do iz=lg%is(3),lg%ie(3)
  do iy=ng%is(2),ng%ie(2)
  do ix=ng%is(1),ng%ie(1)
    poisson%trho2z(ix,iy,iz)=0.d0
  end do
  end do
  end do

!$OMP parallel do private(iz,iy,ix)
  do iz=ng%is(3),ng%ie(3)
  do iy=ng%is(2),ng%ie(2)
  do ix=ng%is(1),ng%ie(1)
    poisson%trho2z(ix,iy,iz)=trho(ix,iy,iz)
  end do
  end do
  end do

  call comm_summation(poisson%trho2z,poisson%trho3z,ng%num(1)*ng%num(2)*lg%num(3),info_field%icomm(3))

!$OMP parallel do private(iz,iy,ix)
  do iz=ng%is(3),ng%ie(3)
  do iy=lg%is(2),lg%ie(2)
  do ix=ng%is(1),ng%ie(1)
    poisson%ff1y(ix,iy,iz)=0.d0
  end do
  end do
  end do

!$OMP parallel do private(kz,iy,ix)
  do kz = ng%is(3),ng%ie(3)
  do iy = ng%is(2),ng%ie(2)
  do ix = ng%is(1),ng%ie(1)
    poisson%ff1y(ix,iy,kz)=sum(poisson%egzc(kz,:)*poisson%trho3z(ix,iy,:))
  end do
  end do
  end do
  call comm_summation(poisson%ff1y,poisson%ff2y,ng%num(1)*lg%num(2)*ng%num(3),info_field%icomm(2))

!$OMP parallel do private(iz,iy,ix)
  do iz=ng%is(3),ng%ie(3)
  do iy=ng%is(2),ng%ie(2)
  do ix=lg%is(1),lg%ie(1)
    poisson%ff1x(ix,iy,iz)=0.d0
  end do
  end do
  end do

!$OMP parallel do private(kz,ky,ix)
  do kz = ng%is(3),ng%ie(3)
  do ky = ng%is(2),ng%ie(2)
  do ix = ng%is(1),ng%ie(1)
    poisson%ff1x(ix,ky,kz)=sum(poisson%egyc(ky,:)*poisson%ff2y(ix,:,kz))
  end do
  end do
  end do

  call comm_summation(poisson%ff1x,poisson%ff2x,lg%num(1)*ng%num(2)*ng%num(3),info_field%icomm(1))

!$OMP parallel do private(kz,ky,kx)
  do kz = ng%is(3),ng%ie(3)
  do ky = ng%is(2),ng%ie(2)
  do kx = ng%is(1),ng%ie(1)
    poisson%ff1x(kx,ky,kz)=sum(poisson%egxc(kx,:)*poisson%ff2x(:,ky,kz))/dble(lg%num(1)*lg%num(2)*lg%num(3))
  end do
  end do
  end do

  call comm_summation(poisson%ff1x,poisson%ff2x,lg%num(1)*ng%num(2)*ng%num(3),info_field%icomm(1))

  B = system%primitive_b

!$OMP parallel do private(iz,iy,ix)
  do iz=lg%is(3),lg%ie(3)
  do iy=ng%is(2),ng%ie(2)
  do ix=ng%is(1),ng%ie(1)
    poisson%ff1z(ix,iy,iz)=0.d0
  end do
  end do
  end do

!$OMP parallel do private(kz,ky,kx,n,gx,gy,gz,g2,kkx,kky,kkz)
  do kz = ng%is(3),ng%ie(3)
  do ky = ng%is(2),ng%ie(2)
  do kx = ng%is(1),ng%ie(1)
    n=(kz-lg%is(3))*lg%num(2)*lg%num(1)+(ky-lg%is(2))*lg%num(1)+kx-lg%is(1)+1
    kkx=kx-1-lg%num(1)*(1+sign(1,(kx-1-lg%num(1)/2)))/2
    kky=ky-1-lg%num(2)*(1+sign(1,(ky-1-lg%num(2)/2)))/2
    kkz=kz-1-lg%num(3)*(1+sign(1,(kz-1-lg%num(3)/2)))/2
    gx = kkx*B(1,1) + kky*B(1,2) + kkz*B(1,3)
    gy = kkx*B(2,1) + kky*B(2,2) + kkz*B(2,3)
    gz = kkx*B(3,1) + kky*B(3,2) + kkz*B(3,3)
    g2=gx**2+gy**2+gz**2
    if(kx-1==0.and.ky-1==0.and.kz-1==0)then
      trhoG_ele_tmp(kx,ky,kz)=poisson%ff2x(kx,ky,kz) ! iwata
      poisson%ff1z(kx,ky,kz)=0.d0
    else
      trhoG_ele_tmp(kx,ky,kz)=poisson%ff2x(kx,ky,kz) ! iwata
      poisson%ff1z(kx,ky,kz)=4.d0*Pi/g2*poisson%ff2x(kx,ky,kz)
    end if
  end do
  end do
  end do

  call comm_summation(trhoG_ele_tmp,trhoG_ele,lg%num(1)*lg%num(2)*lg%num(3),info_field%icomm_all)
  call comm_summation(poisson%ff1z,poisson%ff2z,ng%num(1)*ng%num(2)*lg%num(3),info_field%icomm(3))

!$OMP parallel do private(iz,ky,kx)
  do iz = ng%is(3),ng%ie(3)
  do ky = ng%is(2),ng%ie(2)
  do kx = ng%is(1),ng%ie(1)
    poisson%ff1y(kx,ky,iz)=sum(poisson%egz(:,iz)*poisson%ff2z(kx,ky,:))
  end do
  end do
  end do
  call comm_summation(poisson%ff1y,poisson%ff2y,ng%num(1)*lg%num(2)*ng%num(3),info_field%icomm(2))

!$OMP parallel do private(iz,iy,kx)
  do iz = ng%is(3),ng%ie(3)
  do iy = ng%is(2),ng%ie(2)
  do kx = ng%is(1),ng%ie(1)
    poisson%ff1(kx,iy,iz)=sum(poisson%egy(:,iy)*poisson%ff2y(kx,:,iz))
  end do
  end do
  end do
  call comm_summation(poisson%ff1,poisson%ff2,lg%num(1)*lg%num(2)*lg%num(3),info_field%icomm_all)

!$OMP parallel do private(iz,iy,ix)
  do iz = mg%is(3),mg%ie(3)
  do iy = mg%is(2),mg%ie(2)
  do ix = mg%is(1),mg%ie(1)
    tVh(ix,iy,iz)=sum(poisson%egx(:,ix)*poisson%ff2(:,iy,iz))
  end do
  end do
  end do

  return
end subroutine poisson_periodic
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------

end module poisson_periodic_sub
