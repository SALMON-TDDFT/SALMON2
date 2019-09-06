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
module init_reciprocal_grid_sub
  implicit none

contains

!=====================================================================
subroutine init_reciprocal_grid(lg,ng,fg,system,info_field,poisson)
  use inputoutput,     only : nelem,yn_ffte
  use salmon_parallel, only : nproc_group_global,nproc_id_global,nproc_size_global
  use math_constants,  only : pi
  use structures,      only : s_rgrid,s_reciprocal_grid,s_dft_system,s_field_parallel,s_poisson
  implicit none
  type(s_rgrid),intent(in) :: lg
  type(s_rgrid),intent(in) :: ng
  type(s_reciprocal_grid),intent(inout) :: fg
  type(s_dft_system),intent(in) :: system
  type(s_field_parallel),intent(in) :: info_field
  type(s_poisson),intent(inout) :: poisson
  real(8) :: brl(3,3)
  integer :: n
  integer :: ix,iy,iz,jj
  integer :: iix,iiy,iiz
  integer :: ng_sta_2(3),ng_end_2(3),ng_num_2(3)
  integer :: lg_sta_2(3),lg_end_2(3),lg_num_2(3)
  real(8) :: G2
  integer :: kx,ky,kz
  integer :: kkx,kky,kkz
  integer :: ky2,kz2
  integer :: kx_sta,kx_end,ky_sta,ky_end,kz_sta,kz_end
  real(8) :: bLx,bLy,bLz
  integer :: ky_shift,kz_shift
  
  integer :: npuy,npuz

  brl(:,:)=system%primitive_b(:,:)

  if(allocated(fg%Gx))       deallocate(fg%Gx,fg%Gy,fg%Gz)

  select case(yn_ffte)
  case('n')
    jj = system%ngrid/nproc_size_global
    fg%ig_s = nproc_id_global*jj+1
    fg%ig_e = (nproc_id_global+1)*jj
    if(nproc_id_global==nproc_size_global-1) fg%ig_e = system%ngrid
    fg%icomm_G = nproc_group_global
    fg%ng = system%ngrid
    allocate(fg%Gx(fg%ng),fg%Gy(fg%ng),fg%Gz(fg%ng))

    fg%Gx = 0.d0
    fg%Gy = 0.d0
    fg%Gz = 0.d0
  
    n=0
    do iz=lg%is(3),lg%ie(3)
    do iy=lg%is(2),lg%ie(2)
    do ix=lg%is(1),lg%ie(1)
      n=n+1
      if((ix-1)**2+(iy-1)**2+(iz-1)**2 == 0) fg%iGzero=n
      iix=ix-1-lg%num(1)*(1+sign(1,(ix-1-(lg%num(1)+1)/2)))/2
      iiy=iy-1-lg%num(2)*(1+sign(1,(iy-1-(lg%num(2)+1)/2)))/2
      iiz=iz-1-lg%num(3)*(1+sign(1,(iz-1-(lg%num(3)+1)/2)))/2
      fg%Gx(n) = dble(iix)*brl(1,1) + dble(iiy)*brl(1,2) + dble(iiz)*brl(1,3)
      fg%Gy(n) = dble(iix)*brl(2,1) + dble(iiy)*brl(2,2) + dble(iiz)*brl(2,3)
      fg%Gz(n) = dble(iix)*brl(3,1) + dble(iiy)*brl(3,2) + dble(iiz)*brl(3,3)
    enddo
    enddo
    enddo
  case('y')
    npuy=info_field%isize_ffte(2)
    npuz=info_field%isize_ffte(3)

    jj = system%ngrid/nproc_size_global
    fg%ig_s = nproc_id_global*jj+1
    fg%ig_e = (nproc_id_global+1)*jj
    if(nproc_id_global==nproc_size_global-1) fg%ig_e = system%ngrid
    fg%icomm_G = nproc_group_global
    fg%ng = system%ngrid
    allocate(fg%Gx(fg%ng),fg%Gy(fg%ng),fg%Gz(fg%ng))

    fg%Gx = 0.d0
    fg%Gy = 0.d0
    fg%Gz = 0.d0

    if(.not.allocated(poisson%coef))then
      allocate(poisson%coef(lg%num(1),lg%num(2)/npuy,lg%num(3)/npuz))
    end if
    poisson%coef=0.d0

    lg_sta_2(1:3)=lg%is(1:3)
    lg_end_2(1:3)=lg%ie(1:3)
    lg_num_2(1:3)=lg%num(1:3)
    
    ng_sta_2(1:3)=ng%is(1:3)
    ng_end_2(1:3)=ng%ie(1:3)
    ng_num_2(1:3)=ng%num(1:3)

    bLx=2.d0*pi/(system%hgs(1)*dble(lg_num_2(1)))
    bLy=2.d0*pi/(system%hgs(2)*dble(lg_num_2(2)))
    bLz=2.d0*pi/(system%hgs(3)*dble(lg_num_2(3)))
    
    kx_sta=lg_sta_2(1)
    kx_end=lg_end_2(1)
    ky_sta=1
    ky_end=lg_num_2(2)/npuy
    kz_sta=1
    kz_end=lg_num_2(3)/npuz

    ky_shift=info_field%id_ffte(2)*lg_num_2(2)/npuy
    kz_shift=info_field%id_ffte(3)*lg_num_2(3)/npuz

    fg%iGzero = -1
    do kz = kz_sta,kz_end
    do ky = ky_sta,ky_end
    do kx = kx_sta,kx_end
      ky2=ky+ky_shift
      kz2=kz+kz_shift
      n=(kz2-lg_sta_2(3))*lg_num_2(2)*lg_num_2(1)+(ky2-lg_sta_2(2))*lg_num_2(1)+kx-lg_sta_2(1)+1
      kkx=kx-1-lg_num_2(1)*(1+sign(1,(kx-1-lg_num_2(1)/2)))/2
      kky=ky2-1-lg_num_2(2)*(1+sign(1,(ky2-1-lg_num_2(2)/2)))/2
      kkz=kz2-1-lg_num_2(3)*(1+sign(1,(kz2-1-lg_num_2(3)/2)))/2
      fg%Gx(n)=dble(kkx)*bLx
      fg%Gy(n)=dble(kky)*bLy
      fg%Gz(n)=dble(kkz)*bLz
      G2=fg%Gx(n)**2+fg%Gy(n)**2+fg%Gz(n)**2
      if(kx==1.and.ky2==1.and.kz2==1)then
        fg%iGzero = n
        poisson%coef(kx,ky,kz)=0.d0
      else
        poisson%coef(kx,ky,kz)=4.d0*pi/G2
      end if
    end do
    end do
    end do

  end select 
 
  if(allocated(fg%zrhoG_ion)) deallocate(fg%zrhoG_ion,fg%zrhoG_ele,fg%zrhoG_ele_tmp,fg%zdVG_ion)
  allocate(fg%zrhoG_ion(fg%ng),fg%zrhoG_ele(fg%ng),fg%zrhoG_ele_tmp(fg%ng),fg%zdVG_ion(fg%ng,nelem))

  return
end subroutine init_reciprocal_grid
!=====================================================================

end module init_reciprocal_grid_sub
