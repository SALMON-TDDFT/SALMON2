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
module init_poisson_sub
  implicit none

contains

subroutine init_poisson_fft(lg,ng,system,info_field,poisson)
  use math_constants, only : pi
  use structures,     only: s_rgrid,s_dft_system,s_field_parallel,s_poisson
  implicit none
  type(s_rgrid),intent(in) :: lg
  type(s_rgrid),intent(in) :: ng
  type(s_dft_system),intent(in) :: system
  type(s_field_parallel),intent(in) :: info_field
  type(s_poisson),intent(inout) :: poisson
  integer :: ng_sta_2(3),ng_end_2(3),ng_num_2(3)
  integer :: lg_sta_2(3),lg_end_2(3),lg_num_2(3)
  real(8) :: Gx,Gy,Gz
  real(8) :: G2
  integer :: kx,ky,kz
  integer :: kkx,kky,kkz
  integer :: ky2,kz2
  integer :: n
  integer :: kx_sta,kx_end,ky_sta,ky_end,kz_sta,kz_end
  real(8) :: bLx,bLy,bLz
  integer :: ky_shift,kz_shift
  
  integer :: npuy,npuz

  npuy=info_field%isize_ffte(2)
  npuz=info_field%isize_ffte(3)

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

  do kz = kz_sta,kz_end
  do ky = ky_sta,ky_end
  do kx = kx_sta,kx_end
    ky2=ky+ky_shift
    kz2=kz+kz_shift
    n=(kz2-lg_sta_2(3))*lg_num_2(2)*lg_num_2(1)+(ky2-lg_sta_2(2))*lg_num_2(1)+kx-lg_sta_2(1)+1
    kkx=kx-1-lg_num_2(1)*(1+sign(1,(kx-1-lg_num_2(1)/2)))/2
    kky=ky2-1-lg_num_2(2)*(1+sign(1,(ky2-1-lg_num_2(2)/2)))/2
    kkz=kz2-1-lg_num_2(3)*(1+sign(1,(kz2-1-lg_num_2(3)/2)))/2
    Gx=dble(kkx)*bLx
    Gy=dble(kky)*bLy
    Gz=dble(kkz)*bLz
    G2=Gx**2+Gy**2+Gz**2
    if(kx==1.and.ky2==1.and.kz2==1)then
      poisson%coef(kx,ky,kz)=0.d0
    else
      poisson%coef(kx,ky,kz)=4.d0*pi/G2
    end if
  end do
  end do
  end do

  
  return
end subroutine init_poisson_fft

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------

end module init_poisson_sub
