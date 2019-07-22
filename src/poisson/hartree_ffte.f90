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
module hartree_ffte_sub
  implicit none

contains

subroutine hartree_ffte(lg,mg,ng,trho,tvh,hgs,npuw,npuy,npuz,  &
                        a_ffte,b_ffte,rhoe_g,coef_poisson)
  use structures, only: s_rgrid
  use salmon_parallel, only: nproc_id_icommy
  use salmon_parallel, only: nproc_id_icommz
  use salmon_parallel, only: nproc_group_icommw
  use salmon_communication, only: comm_summation
  use salmon_communication, only: comm_is_root
  use math_constants, only : pi
!$  use omp_lib
  implicit none
  type(s_rgrid),intent(in) :: lg
  type(s_rgrid),intent(in) :: mg
  type(s_rgrid),intent(in) :: ng
  real(8),intent(in)       :: hgs(3)
  integer,intent(in)       :: npuw,npuy,npuz
  complex(8),intent(out)   :: a_ffte(lg%num(1),lg%num(2)/npuy,lg%num(3)/npuz)
  complex(8),intent(out)   :: b_ffte(lg%num(1),lg%num(2)/npuy,lg%num(3)/npuz)
  complex(8),intent(out)   :: rhoe_g(lg%num(1)*lg%num(2)*lg%num(3))
  real(8),intent(in)       :: coef_poisson(lg%num(1),lg%num(2)/npuy,lg%num(3)/npuz)
  integer :: ix,iy,iz
  integer :: iix,iiy,iiz
  integer :: iz_sta,iz_end,iy_sta,iy_end
  real(8) :: trho(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
  real(8) :: tvh(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
  real(8) :: inv_lgnum3
  complex(8),parameter :: zI=(0.d0,1.d0)
  integer :: n
  real(8) :: bLx,bLy,bLz
  complex(8) :: a_ffte_tmp(1:lg%num(1),1:lg%num(2)/npuy,1:lg%num(3)/npuz)

  bLx=2.d0*Pi/(Hgs(1)*dble(lg%num(1)))
  bLy=2.d0*Pi/(Hgs(2)*dble(lg%num(2)))
  bLz=2.d0*Pi/(Hgs(3)*dble(lg%num(3)))

  inv_lgnum3=1.d0/(lg%num(1)*lg%num(2)*lg%num(3))

  iz_sta=1
  iz_end=lg%num(3)/npuz
  iy_sta=1
  iy_end=lg%num(2)/npuy
  
!  rhoe_G_tmp=0.d0

  if(npuw==1)then
!$OMP parallel do private(iiz,iiy)
    do iz=iz_sta,iz_end
      iiz=iz+nproc_id_icommz*lg%num(3)/npuz
      do iy=iy_sta,iy_end
        iiy=iy+nproc_id_icommy*lg%num(2)/npuy
        a_ffte(1:lg%ie(1),iy,iz)=trho(1:lg%ie(1),iiy,iiz)
      end do
    end do
  else
    a_ffte_tmp=0.d0
!$OMP parallel do private(iiz,iiy,ix)
    do iz=iz_sta,iz_end
      iiz=iz+nproc_id_icommz*lg%num(3)/npuz
      do iy=iy_sta,iy_end
        iiy=iy+nproc_id_icommy*lg%num(2)/npuy
        do iix=ng%is(1),ng%ie(1)
          ix=iix-lg%is(1)+1
          a_ffte_tmp(ix,iy,iz)=trho(iix,iiy,iiz)
        end do
      end do
    end do
    call comm_summation(a_ffte_tmp,a_ffte,lg%num(1)*lg%num(2)/npuy*lg%num(3)/npuz,nproc_group_icommw)
  end if

  CALL PZFFT3DV_MOD(a_ffte,b_ffte,lg%num(1),lg%num(2),lg%num(3),npuy,npuz,0) 
  CALL PZFFT3DV_MOD(a_ffte,b_ffte,lg%num(1),lg%num(2),lg%num(3),npuy,npuz,-1) 

!$OMP parallel do private(n)
  do iz=iz_sta,iz_end
    do iy=iy_sta,iy_end
      do ix=1,lg%num(1)
        n=(iz-1)*lg%num(2)/npuy*lg%num(1)+(iy-1)*lg%num(1)+ix
        rhoe_G(n)=b_ffte(ix,iy,iz)*inv_lgnum3
        b_ffte(ix,iy,iz)=b_ffte(ix,iy,iz)*coef_poisson(ix,iy,iz)
      end do
    end do
  end do
  if(nproc_id_icommz==0.and.nproc_id_icommy==0)then
    rhoe_G(1)=0.d0
  end if

  CALL PZFFT3DV_MOD(b_ffte,a_ffte,lg%num(1),lg%num(2),lg%num(3),npuy,npuz,1)

  if(npuw==1)then
!$OMP parallel do private(iiz,iiy)
    do iz=iz_sta,iz_end
      iiz=iz+nproc_id_icommz*lg%num(3)/npuz
      do iy=iy_sta,iy_end
        iiy=iy+nproc_id_icommy*lg%num(2)/npuy
        tvh(1:lg%ie(1),iiy,iiz)=a_ffte(1:lg%ie(1),iy,iz)
      end do
    end do
  else
!$OMP parallel do private(iiz,iiy,ix)
    do iz=iz_sta,iz_end
      iiz=iz+nproc_id_icommz*lg%num(3)/npuz
      do iy=iy_sta,iy_end
        iiy=iy+nproc_id_icommy*lg%num(2)/npuy
        do iix=ng%is(1),ng%ie(1)
          ix=iix-lg%is(1)+1
          tvh(iix,iiy,iiz)=a_ffte(ix,iy,iz)
        end do
      end do
    end do
  end if

  return
end subroutine hartree_ffte
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------

end module hartree_ffte_sub
