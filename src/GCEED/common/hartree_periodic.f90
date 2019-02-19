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
!SUBROUTINE Hartree_periodic
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------
subroutine Hartree_periodic(lg,mg,ng,trho,tVh)
  use structures, only: s_rgrid
  use salmon_parallel, only: nproc_group_global, nproc_group_bound
  use salmon_communication, only: comm_summation
  use scf_data
  use new_world_sub
  use allocate_mat_sub
  implicit none
  type(s_rgrid),intent(in) :: lg
  type(s_rgrid),intent(in) :: mg
  type(s_rgrid),intent(in) :: ng
  integer :: ix,iy,iz,kx,ky,kz,kkx,kky,kkz
  real(8) :: trho(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
  real(8) :: tVh(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
  real(8) :: Gx,Gy,Gz
  complex(8),parameter :: zI=(0.d0,1.d0)
  real(8) :: G2
  real(8) :: bLx,bLy,bLz
  integer :: n

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
    trho2z(ix,iy,iz)=trho(ix,iy,iz)
  end do
  end do
  end do

  call comm_summation(trho2z,trho3z,ng%num(1)*ng%num(2)*lg%num(3),nproc_group_bound(3))
  
!$OMP parallel do private(ix,kx)
  do ix=lg%is(1),lg%ie(1)
    do kx=lg%is(1),lg%ie(1)
      eGx(kx,ix)=exp(zI*(2.d0*Pi*dble((ix-1)*(kx-1))/dble(lg%num(1))))
      eGxc(kx,ix)=conjg(eGx(kx,ix))
    end do
  end do
!$OMP parallel do private(iy,ky)
  do iy=lg%is(2),lg%ie(2)
    do ky=lg%is(2),lg%ie(2)
      eGy(ky,iy)=exp(zI*(2.d0*Pi*dble((iy-1)*(ky-1))/dble(lg%num(2))))
      eGyc(ky,iy)=conjg(eGy(ky,iy))
    end do
  end do
!$OMP parallel do private(iz,kz)
  do iz=lg%is(3),lg%ie(3)
    do kz=lg%is(3),lg%ie(3)
      eGz(kz,iz)=exp(zI*(2.d0*Pi*dble((iz-1)*(kz-1))/dble(lg%num(3))))
      eGzc(kz,iz)=conjg(eGz(kz,iz))
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
    ff1y(ix,iy,kz)=sum(eGzc(kz,:)*trho3z(ix,iy,:))
  end do
  end do
  end do
  call comm_summation(ff1y,ff2y,ng%num(1)*lg%num(2)*ng%num(3),nproc_group_bound(2))

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
    ff1x(ix,ky,kz)=sum(eGyc(ky,:)*ff2y(ix,:,kz))
  end do
  end do
  end do

  call comm_summation(ff1x,ff2x,lg%num(1)*ng%num(2)*ng%num(3),nproc_group_bound(1))

!$OMP parallel do private(kz,ky,kx)
  do kz = ng%is(3),ng%ie(3)
  do ky = ng%is(2),ng%ie(2)
  do kx = ng%is(1),ng%ie(1)
    ff1x(kx,ky,kz)=sum(eGxc(kx,:)*ff2x(:,ky,kz))/dble(lg%num(1)*lg%num(2)*lg%num(3))
  end do
  end do
  end do

  call comm_summation(ff1x,ff2x,lg%num(1)*ng%num(2)*ng%num(3),nproc_group_bound(1))

  bLx=2.d0*Pi/(Hgs(1)*dble(lg%num(1)))
  bLy=2.d0*Pi/(Hgs(2)*dble(lg%num(2)))
  bLz=2.d0*Pi/(Hgs(3)*dble(lg%num(3)))

!$OMP parallel do private(iz,iy,ix)
  do iz=lg%is(3),lg%ie(3)
  do iy=ng%is(2),ng%ie(2)
  do ix=ng%is(1),ng%ie(1)
    ff1z(ix,iy,iz)=0.d0
  end do
  end do
  end do

!$OMP parallel do private(kz,ky,kx,n,Gx,Gy,Gz,G2,kkx,kky,kkz)
  do kz = ng%is(3),ng%ie(3)
  do ky = ng%is(2),ng%ie(2)
  do kx = ng%is(1),ng%ie(1)
    n=(kz-lg%is(3))*lg%num(2)*lg%num(1)+(ky-lg%is(2))*lg%num(1)+kx-lg%is(1)+1
!    Gx=2.d0*Pi*kx/lg%num(1)
!    Gy=2.d0*Pi*ky/lg%num(2)
!    Gz=2.d0*Pi*kz/lg%num(3)
    kkx=kx-1-lg%num(1)*(1+sign(1,(kx-1-lg%num(1)/2)))/2
    kky=ky-1-lg%num(2)*(1+sign(1,(ky-1-lg%num(2)/2)))/2
    kkz=kz-1-lg%num(3)*(1+sign(1,(kz-1-lg%num(3)/2)))/2
    Gx=kkx*bLx
    Gy=kky*bLy
    Gz=kkz*bLz
    G2=Gx**2+Gy**2+Gz**2
    if(kx-1==0.and.ky-1==0.and.kz-1==0)then
      rhoe_G_tmp(n)=0.d0
      ff1z(kx,ky,kz)=0.d0
    else
      rhoe_G_tmp(n)=ff2x(kx,ky,kz)
      ff1z(kx,ky,kz)=4.d0*Pi/G2*ff2x(kx,ky,kz)
    end if
  end do
  end do
  end do

  if(iSCFRT==1)then
    call comm_summation(rhoe_G_tmp,rhoe_G,lg%num(1)*lg%num(2)*lg%num(3),nproc_group_global)
  else if(iSCFRT==2)then
    if(itt==1.or.mod(itt,itcalc_ene)==0)then
      call comm_summation(rhoe_G_tmp,rhoe_G,lg%num(1)*lg%num(2)*lg%num(3),nproc_group_global)
    end if
  end if
  call comm_summation(ff1z,ff2z,ng%num(1)*ng%num(2)*lg%num(3),nproc_group_bound(3))

!$OMP parallel do private(iz,ky,kx)
  do iz = ng%is(3),ng%ie(3)
  do ky = ng%is(2),ng%ie(2)
  do kx = ng%is(1),ng%ie(1)
    ff1y(kx,ky,iz)=sum(eGz(:,iz)*ff2z(kx,ky,:))
  end do
  end do
  end do
  call comm_summation(ff1y,ff2y,ng%num(1)*lg%num(2)*ng%num(3),nproc_group_bound(2))

!$OMP parallel do private(iz,iy,kx)
  do iz = ng%is(3),ng%ie(3)
  do iy = ng%is(2),ng%ie(2)
  do kx = ng%is(1),ng%ie(1)
    ff1(kx,iy,iz)=sum(eGy(:,iy)*ff2y(kx,:,iz))
  end do
  end do
  end do
  call comm_summation(ff1,ff2,lg%num(1)*lg%num(2)*lg%num(3),nproc_group_global)

!$OMP parallel do private(iz,iy,ix)
  do iz = mg%is(3),mg%ie(3)
  do iy = mg%is(2),mg%ie(2)
  do ix = mg%is(1),mg%ie(1)
    tVh(ix,iy,iz)=sum(eGx(:,ix)*ff2(:,iy,iz))
  end do
  end do
  end do

  return
end subroutine Hartree_periodic
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------
