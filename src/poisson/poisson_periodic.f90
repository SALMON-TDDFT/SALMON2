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

subroutine poisson_periodic(lg,mg,ng,system,info_field,srho,sVh,fg)
  use structures, only: s_rgrid, s_field_parallel, s_dft_system, &
                        s_scalar, s_reciprocal_grid
  use salmon_parallel, only: nproc_group_global
  use salmon_communication, only: comm_summation
  use math_constants, only : pi
  implicit none
  type(s_rgrid),intent(in) :: lg
  type(s_rgrid),intent(in) :: mg
  type(s_rgrid),intent(in) :: ng
  type(s_field_parallel),intent(in) :: info_field
  type(s_dft_system),intent(in) :: system
  type(s_scalar),intent(in)    :: srho
  type(s_scalar),intent(inout) :: sVh
  type(s_reciprocal_grid),intent(inout)  :: fg

  integer :: ix,iy,iz,kx,ky,kz,kkx,kky,kkz
  real(8) :: gx,gy,gz
  complex(8),parameter :: zI=(0.d0,1.d0)
  real(8) :: g2
  real(8) :: B(3,3)
  integer :: n

  complex(8),allocatable :: ff1(:,:,:)
  complex(8),allocatable :: ff1x(:,:,:)
  complex(8),allocatable :: ff1y(:,:,:)
  complex(8),allocatable :: ff1z(:,:,:)
  complex(8),allocatable :: ff2(:,:,:)
  complex(8),allocatable :: ff2x(:,:,:)
  complex(8),allocatable :: ff2y(:,:,:)
  complex(8),allocatable :: ff2z(:,:,:)
  complex(8),allocatable :: rhoe_g_tmp(:)
  real(8),allocatable    :: trho2z(:,:,:)
  real(8),allocatable    :: trho3z(:,:,:)
  complex(8),allocatable :: egx(:,:)
  complex(8),allocatable :: egxc(:,:)
  complex(8),allocatable :: egy(:,:)
  complex(8),allocatable :: egyc(:,:)
  complex(8),allocatable :: egz(:,:)
  complex(8),allocatable :: egzc(:,:)

  if(.not.allocated(ff1))then
    allocate(ff1(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3)))
    allocate(ff1x(lg%is(1):lg%ie(1),ng%is(2):ng%ie(2),ng%is(3):ng%ie(3)))
    allocate(ff1y(ng%is(1):ng%ie(1),lg%is(2):lg%ie(2),ng%is(3):ng%ie(3)))
    allocate(ff1z(ng%is(1):ng%ie(1),ng%is(2):ng%ie(2),lg%is(3):lg%ie(3)))
    allocate(ff2(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3)))
    allocate(ff2x(lg%is(1):lg%ie(1),ng%is(2):ng%ie(2),ng%is(3):ng%ie(3)))
    allocate(ff2y(ng%is(1):ng%ie(1),lg%is(2):lg%ie(2),ng%is(3):ng%ie(3)))
    allocate(ff2z(ng%is(1):ng%ie(1),ng%is(2):ng%ie(2),lg%is(3):lg%ie(3)))
    allocate(egx(lg%is(1):lg%ie(1),lg%is(1):lg%ie(1)))
    allocate(egxc(lg%is(1):lg%ie(1),lg%is(1):lg%ie(1)))
    allocate(egy(lg%is(2):lg%ie(2),lg%is(2):lg%ie(2)))
    allocate(egyc(lg%is(2):lg%ie(2),lg%is(2):lg%ie(2)))
    allocate(egz(lg%is(3):lg%ie(3),lg%is(3):lg%ie(3)))
    allocate(egzc(lg%is(3):lg%ie(3),lg%is(3):lg%ie(3)))
    allocate(rhoe_g_tmp(lg%num(1)*lg%num(2)*lg%num(3)))
    allocate(trho2z(ng%is(1):ng%ie(1),ng%is(2):ng%ie(2),lg%is(3):lg%ie(3)))
    allocate(trho3z(ng%is(1):ng%ie(1),ng%is(2):ng%ie(2),lg%is(3):lg%ie(3)))
  end if

!$OMP parallel do private(iz,iy,ix)
  do iz=lg%is(3),lg%ie(3)
  do iy=lg%is(2),lg%ie(2)
  do ix=lg%is(1),lg%ie(1)
    ff1(ix,iy,iz)=0.d0
  end do
  end do
  end do
!$OMP parallel do
  do n=1,lg%num(1)*lg%num(2)*lg%num(3)
    rhoe_G_tmp(n)=0.d0
  end do

!$OMP parallel do private(iz,iy,ix)
  do iz=lg%is(3),lg%ie(3)
  do iy=ng%is(2),ng%ie(2)
  do ix=ng%is(1),ng%ie(1)
    trho2z(ix,iy,iz)=0.d0
  end do
  end do
  end do

!$OMP parallel do private(iz,iy,ix)
  do iz=ng%is(3),ng%ie(3)
  do iy=ng%is(2),ng%ie(2)
  do ix=ng%is(1),ng%ie(1)
    trho2z(ix,iy,iz)=srho%f(ix,iy,iz)
  end do
  end do
  end do

  call comm_summation(trho2z,trho3z,ng%num(1)*ng%num(2)*lg%num(3),info_field%icomm(3))
  
!$OMP parallel do private(ix,kx)
  do ix=lg%is(1),lg%ie(1)
    do kx=lg%is(1),lg%ie(1)
      egx(kx,ix)=exp(zI*(2.d0*Pi*dble((ix-1)*(kx-1))/dble(lg%num(1))))
      egxc(kx,ix)=conjg(egx(kx,ix))
    end do
  end do
!$OMP parallel do private(iy,ky)
  do iy=lg%is(2),lg%ie(2)
    do ky=lg%is(2),lg%ie(2)
      egy(ky,iy)=exp(zI*(2.d0*Pi*dble((iy-1)*(ky-1))/dble(lg%num(2))))
      egyc(ky,iy)=conjg(egy(ky,iy))
    end do
  end do
!$OMP parallel do private(iz,kz)
  do iz=lg%is(3),lg%ie(3)
    do kz=lg%is(3),lg%ie(3)
      egz(kz,iz)=exp(zI*(2.d0*Pi*dble((iz-1)*(kz-1))/dble(lg%num(3))))
      egzc(kz,iz)=conjg(egz(kz,iz))
    end do
  end do

!$OMP parallel do private(iz,iy,ix)
  do iz=ng%is(3),ng%ie(3)
  do iy=lg%is(2),lg%ie(2)
  do ix=ng%is(1),ng%ie(1)
    ff1y(ix,iy,iz)=0.d0
  end do
  end do
  end do

!$OMP parallel do private(kz,iy,ix)
  do kz = ng%is(3),ng%ie(3)
  do iy = ng%is(2),ng%ie(2)
  do ix = ng%is(1),ng%ie(1)
    ff1y(ix,iy,kz)=sum(egzc(kz,:)*trho3z(ix,iy,:))
  end do
  end do
  end do
  call comm_summation(ff1y,ff2y,ng%num(1)*lg%num(2)*ng%num(3),info_field%icomm(2))

!$OMP parallel do private(iz,iy,ix)
  do iz=ng%is(3),ng%ie(3)
  do iy=ng%is(2),ng%ie(2)
  do ix=lg%is(1),lg%ie(1)
    ff1x(ix,iy,iz)=0.d0
  end do
  end do
  end do

!$OMP parallel do private(kz,ky,ix)
  do kz = ng%is(3),ng%ie(3)
  do ky = ng%is(2),ng%ie(2)
  do ix = ng%is(1),ng%ie(1)
    ff1x(ix,ky,kz)=sum(egyc(ky,:)*ff2y(ix,:,kz))
  end do
  end do
  end do

  call comm_summation(ff1x,ff2x,lg%num(1)*ng%num(2)*ng%num(3),info_field%icomm(1))

!$OMP parallel do private(kz,ky,kx)
  do kz = ng%is(3),ng%ie(3)
  do ky = ng%is(2),ng%ie(2)
  do kx = ng%is(1),ng%ie(1)
    ff1x(kx,ky,kz)=sum(egxc(kx,:)*ff2x(:,ky,kz))/dble(lg%num(1)*lg%num(2)*lg%num(3))
  end do
  end do
  end do

  call comm_summation(ff1x,ff2x,lg%num(1)*ng%num(2)*ng%num(3),info_field%icomm(1))

  B = system%primitive_b

!$OMP parallel do private(iz,iy,ix)
  do iz=lg%is(3),lg%ie(3)
  do iy=ng%is(2),ng%ie(2)
  do ix=ng%is(1),ng%ie(1)
    ff1z(ix,iy,iz)=0.d0
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
      rhoe_G_tmp(n)=ff2x(kx,ky,kz) ! iwata
      ff1z(kx,ky,kz)=0.d0
    else
      rhoe_G_tmp(n)=ff2x(kx,ky,kz)
      ff1z(kx,ky,kz)=4.d0*Pi/g2*ff2x(kx,ky,kz)
    end if
  end do
  end do
  end do

  call comm_summation(rhoe_G_tmp,fg%zrhoG_ele,lg%num(1)*lg%num(2)*lg%num(3),nproc_group_global)
  call comm_summation(ff1z,ff2z,ng%num(1)*ng%num(2)*lg%num(3),info_field%icomm(3))

!$OMP parallel do private(iz,ky,kx)
  do iz = ng%is(3),ng%ie(3)
  do ky = ng%is(2),ng%ie(2)
  do kx = ng%is(1),ng%ie(1)
    ff1y(kx,ky,iz)=sum(egz(:,iz)*ff2z(kx,ky,:))
  end do
  end do
  end do
  call comm_summation(ff1y,ff2y,ng%num(1)*lg%num(2)*ng%num(3),info_field%icomm(2))

!$OMP parallel do private(iz,iy,kx)
  do iz = ng%is(3),ng%ie(3)
  do iy = ng%is(2),ng%ie(2)
  do kx = ng%is(1),ng%ie(1)
    ff1(kx,iy,iz)=sum(egy(:,iy)*ff2y(kx,:,iz))
  end do
  end do
  end do
  call comm_summation(ff1,ff2,lg%num(1)*lg%num(2)*lg%num(3),nproc_group_global)

!$OMP parallel do private(iz,iy,ix)
  do iz = mg%is(3),mg%ie(3)
  do iy = mg%is(2),mg%ie(2)
  do ix = mg%is(1),mg%ie(1)
    sVh%f(ix,iy,iz)=sum(egx(:,ix)*ff2(:,iy,iz))
  end do
  end do
  end do

  return
end subroutine poisson_periodic
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------

end module poisson_periodic_sub
