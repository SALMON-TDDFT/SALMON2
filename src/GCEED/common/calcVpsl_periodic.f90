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
subroutine calcVpsl_periodic(lg,matrix_A,brl)
  use structures,      only: s_rgrid
  use salmon_parallel, only: nproc_group_global, nproc_size_global, nproc_id_global
  use salmon_communication, only: comm_bcast, comm_summation, comm_is_root
  use prep_pp_sub, only: calc_vloc,calc_vpsl
  use scf_data
  use new_world_sub
  use allocate_psl_sub
  use allocate_mat_sub
  implicit none
  type(s_rgrid),intent(in) :: lg
  real(8),intent(in) :: matrix_A(3,3),brl(3,3)
  
  integer :: ii,ix,iy,iz,ak
  integer :: iix,iiy,iiz
  integer :: n
  real(8) :: aLxyz
  integer :: NG_s,NG_e
  integer :: NG_l_s_para,NG_l_e_para
  integer :: numtmp
  complex(8),parameter :: zI=(0.d0,1.d0)
  integer :: lx(lg%num(1)*lg%num(2)*lg%num(3))
  integer :: ly(lg%num(1)*lg%num(2)*lg%num(3))
  integer :: lz(lg%num(1)*lg%num(2)*lg%num(3))
  complex(8),allocatable :: dvloc_g_tmp2(:,:)
  complex(8),allocatable :: rhoion_g_tmp2(:)
  real(8),allocatable :: vpsl_ia(:,:)
  real(8),allocatable :: vpsl_tmp2(:)
  integer :: i,ia
  real(8) :: hx,hy,hz
 
  NG_s=1
  NG_e=lg%num(1)*lg%num(2)*lg%num(3)
  
  numtmp=(NG_e-NG_s+1)/nproc_size_global
  
  NG_l_s_para = nproc_id_global*numtmp+1
  NG_l_e_para = (nproc_id_global+1)*numtmp
  if(nproc_id_global==nproc_size_global-1) NG_l_e_para=NG_e

  allocate(dvloc_g_tmp2(ng_l_s_para:ng_l_e_para,MKI))
 
  nGzero=-1
  
  do ak=1,MKI
    do ii=1,Mr(ak)
      vloctbl(ii,ak)=pp%vpp_f(ii,Lref(ak),ak)
    enddo
  end do
  
  n=0
  do iz=lg%is(3),lg%ie(3)
  do iy=lg%is(2),lg%ie(2)
  do ix=lg%is(1),lg%ie(1)
    n=n+1
    if((ix-1)**2+(iy-1)**2+(iz-1)**2 == 0) nGzero=n
    iix=ix-1-lg%num(1)*(1+sign(1,(ix-1-(lg%num(1)+1)/2)))/2
    iiy=iy-1-lg%num(2)*(1+sign(1,(iy-1-(lg%num(2)+1)/2)))/2
    iiz=iz-1-lg%num(3)*(1+sign(1,(iz-1-(lg%num(3)+1)/2)))/2
    Gx(n) = dble(iix)*brl(1,1) + dble(iiy)*brl(1,2) + dble(iiz)*brl(1,3)
    Gy(n) = dble(iix)*brl(2,1) + dble(iiy)*brl(2,2) + dble(iiz)*brl(2,3)
    Gz(n) = dble(iix)*brl(3,1) + dble(iiy)*brl(3,2) + dble(iiz)*brl(3,3)
!    Gx(n)=dble(iix)*bLx
!    Gy(n)=dble(iiy)*bLy
!    Gz(n)=dble(iiz)*bLz
  enddo
  enddo
  enddo

  call calc_vloc(pp,dvloc_g_tmp2,gx,gy,gz,ng_e,ng_l_s_para,ng_l_e_para,nGzero)

  dvloc_g_tmp=0.d0
  do ak=1,MKI
    do n=ng_l_s_para,ng_l_e_para
      dvloc_g_tmp(n,ak)=dvloc_g_tmp2(n,ak)
    end do
  end do

  call comm_summation(dVloc_G_tmp,dVloc_G,(NG_e-NG_s+1)*MKI,nproc_group_global)

  
  hx=Hgs(1) 
  hy=Hgs(2) 
  hz=Hgs(3)
  aLxyz=Hvol*dble(lg%num(1)*lg%num(2)*lg%num(3))

  do iz=1,lg%num(3)
  do iy=1,lg%num(2)
  do ix=1,lg%num(1)
    i=(iz-1)*lg%num(1)*lg%num(2)+(iy-1)*lg%num(1)+ix
    lx(i)=ix-1
    ly(i)=iy-1
    lz(i)=iz-1
  end do
  end do
  end do
 
  allocate(rhoion_g_tmp2(ng_l_s_para:ng_l_e_para))
  allocate(vpsl_ia(lg%num(1)*lg%num(2)*lg%num(3),MI))
  allocate(vpsl_tmp2(ng_s:ng_e))
 
  call calc_vpsl(pp,rhoion_g_tmp2,vpsl_ia,vpsl_tmp2,dvloc_g_tmp2,  &
                     ngzero,gx,gy,gz,ng_e,ng_l_s_para,ng_l_e_para,ng_e,alxyz,lx,ly,lz,hx,hy,hz,matrix_A)

  deallocate(dvloc_g_tmp2)

  rhoion_G_tmp=0.d0
  rhoion_G_tmp(ng_l_s_para:ng_l_e_para)=rhoion_g_tmp2(ng_l_s_para:ng_l_e_para)

  call comm_summation(rhoion_G_tmp,rhoion_G,NG_e-NG_s+1,nproc_group_global)

  allocate(ppg%Vpsl_atom(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),MI))

  do iz=mg_sta(3),mg_end(3)
  do iy=mg_sta(2),mg_end(2)
  do ix=mg_sta(1),mg_end(1)
    i=(iz-1)*lg%num(1)*lg%num(2)+(iy-1)*lg%num(1)+ix
    Vpsl(ix,iy,iz)=vpsl_tmp2(i)
    do ia=1,MI
      ppg%Vpsl_atom(ix,iy,iz,ia) = Vpsl_ia(i,ia)
    end do
  enddo
  enddo
  enddo

  return

end subroutine calcVpsl_periodic
