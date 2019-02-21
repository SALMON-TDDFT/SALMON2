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
subroutine Hartree_FFTE(lg,mg,ng,trho,tVh,icheck_ascorder)
  use structures, only: s_rgrid
  use salmon_parallel, only: nproc_group_global
  use salmon_parallel, only: nproc_id_icommy, nproc_group_icommy
  use salmon_parallel, only: nproc_id_icommz, nproc_group_icommz
  use salmon_parallel, only: nproc_group_icommw
  use salmon_communication, only: comm_summation
  use salmon_parallel, only: nproc_id_global
  use salmon_communication, only: comm_is_root
  use scf_data
  use allocate_mat_sub
!$  use omp_lib
  implicit none
  type(s_rgrid),intent(in) :: lg
  type(s_rgrid),intent(in) :: mg
  type(s_rgrid),intent(in) :: ng
  integer,intent(in) :: icheck_ascorder
  integer :: ix,iy,iz
  integer :: iix,iiy,iiz
  integer :: iz_sta,iz_end,iy_sta,iy_end
  real(8) :: trho(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
  real(8) :: tVh(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
  real(8) :: inv_lgnum3
  complex(8),parameter :: zI=(0.d0,1.d0)
  integer :: n
  real(8) :: bLx,bLy,bLz
  complex(8) :: A_FFTE_tmp(1:lg%num(1),1:lg%num(2)/NPUY,1:lg%num(3)/NPUZ)

  bLx=2.d0*Pi/(Hgs(1)*dble(lg%num(1)))
  bLy=2.d0*Pi/(Hgs(2)*dble(lg%num(2)))
  bLz=2.d0*Pi/(Hgs(3)*dble(lg%num(3)))

  inv_lgnum3=1.d0/(lg%num(1)*lg%num(2)*lg%num(3))

  iz_sta=1
  iz_end=lg%num(3)/NPUZ
  iy_sta=1
  iy_end=lg%num(2)/NPUY
  
!  rhoe_G_tmp=0.d0

  if(icheck_ascorder==1)then
    if(NPUW==1)then
!$OMP parallel do private(iiz,iiy)
      do iz=iz_sta,iz_end
        iiz=iz+nproc_id_icommz*lg%num(3)/NPUZ
        do iy=iy_sta,iy_end
          iiy=iy+nproc_id_icommy*lg%num(2)/NPUY
          A_FFTE(1:lg%ie(1),iy,iz)=trho(1:lg%ie(1),iiy,iiz)
        end do
      end do
    else
      A_FFTE_tmp=0.d0
!$OMP parallel do private(iiz,iiy,ix)
      do iz=iz_sta,iz_end
        iiz=iz+nproc_id_icommz*lg%num(3)/NPUZ
        do iy=iy_sta,iy_end
          iiy=iy+nproc_id_icommy*lg%num(2)/NPUY
          do iix=ng%is(1),ng%ie(1)
            ix=iix-lg%is(1)+1
            A_FFTE_tmp(ix,iy,iz)=trho(iix,iiy,iiz)
          end do
        end do
      end do
      call comm_summation(A_FFTE_tmp,A_FFTE,lg%num(1)*lg%num(2)/NPUY*lg%num(3)/NPUZ,nproc_group_icommw)
    end if
  else
!$OMP parallel do
    do iz = lg%is(3),lg%ie(3)
    do iy = lg%is(2),lg%ie(2)
    do ix = lg%is(1),lg%ie(1)
      matbox_l(ix,iy,iz)=0.d0
    end do
    end do
    end do
!$OMP parallel do
    do iz = ng%is(3),ng%ie(3)
    do iy = ng%is(2),ng%ie(2)
    do ix = ng%is(1),ng%ie(1)
      matbox_l(ix,iy,iz)=trho(ix,iy,iz)
    end do
    end do
    end do

    call comm_summation(matbox_l,matbox_l2,lg%num(1)*lg%num(2)*lg%num(3),nproc_group_global)

!!$OMP parallel do private(iiz,iiy)
    do iz=iz_sta,iz_end
      iiz=iz+nproc_id_icommz*lg%num(3)/NPUZ
      do iy=iy_sta,iy_end
        iiy=iy+nproc_id_icommy*lg%num(2)/NPUY
        A_FFTE(1:lg%ie(1),iy,iz)=matbox_l2(1:lg%ie(1),iiy,iiz)
      end do
    end do
  end if

  CALL PZFFT3DV_MOD(A_FFTE,B_FFTE,lg%num(1),lg%num(2),lg%num(3),NPUY,NPUZ,0) 
  CALL PZFFT3DV_MOD(A_FFTE,B_FFTE,lg%num(1),lg%num(2),lg%num(3),NPUY,NPUZ,-1) 

!$OMP parallel do private(n)
  do iz=iz_sta,iz_end
    do iy=iy_sta,iy_end
      do ix=1,lg%num(1)
        n=(iz-1)*lg%num(2)/NPUY*lg%num(1)+(iy-1)*lg%num(1)+ix
        rhoe_G(n)=B_FFTE(ix,iy,iz)*inv_lgnum3
        B_FFTE(ix,iy,iz)=B_FFTE(ix,iy,iz)*coef_poisson(ix,iy,iz)
      end do
    end do
  end do
  if(nproc_id_icommz==0.and.nproc_id_icommy==0)then
    rhoe_G(1)=0.d0
  end if

  CALL PZFFT3DV_MOD(B_FFTE,A_FFTE,lg%num(1),lg%num(2),lg%num(3),NPUY,NPUZ,1)

  if(icheck_ascorder==1)then
    if(NPUW==1)then
!$OMP parallel do private(iiz,iiy)
      do iz=iz_sta,iz_end
        iiz=iz+nproc_id_icommz*lg%num(3)/NPUZ
        do iy=iy_sta,iy_end
          iiy=iy+nproc_id_icommy*lg%num(2)/NPUY
          tVh(1:lg%ie(1),iiy,iiz)=A_FFTE(1:lg%ie(1),iy,iz)
        end do
      end do
    else
!$OMP parallel do private(iiz,iiy,ix)
      do iz=iz_sta,iz_end
        iiz=iz+nproc_id_icommz*lg%num(3)/NPUZ
        do iy=iy_sta,iy_end
          iiy=iy+nproc_id_icommy*lg%num(2)/NPUY
          do iix=ng%is(1),ng%ie(1)
            ix=iix-lg%is(1)+1
            tVh(iix,iiy,iiz)=A_FFTE(ix,iy,iz)
          end do
        end do
      end do
    end if
  else
!$OMP parallel do
    do iz = lg%is(3),lg%ie(3)
    do iy = lg%is(2),lg%ie(2)
    do ix = lg%is(1),lg%ie(1)
      matbox_l(ix,iy,iz)=0.d0
    end do
    end do
    end do
!!$OMP parallel do private(iiz,iiy)
    do iz=iz_sta,iz_end
      iiz=iz+nproc_id_icommz*lg%num(3)/NPUZ
      do iy=iy_sta,iy_end
        iiy=iy+nproc_id_icommy*lg%num(2)/NPUY
        matbox_l(1:lg%ie(1),iiy,iiz)=A_FFTE(1:lg%ie(1),iy,iz)
      end do
    end do
    call comm_summation(matbox_l,matbox_l2,lg%num(1)*lg%num(2)*lg%num(3),nproc_group_global)
!$OMP parallel do
    do iz = mg%is(3),mg%ie(3)
    do iy = mg%is(2),mg%ie(2)
    do ix = mg%is(1),mg%ie(1)
      tVh(ix,iy,iz)=matbox_l2(ix,iy,iz)
    end do
    end do
    end do
  end if 

  return
end subroutine Hartree_FFTE
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------
