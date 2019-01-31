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
MODULE hpsi2_sub

use scf_data
use laplacian2_sub
use new_world_sub
use gradient2_sub
use allocate_mat_sub
implicit none

INTERFACE hpsi2
!   MODULE PROCEDURE R_hpsi2,C_hpsi2
  MODULE PROCEDURE hpsi_test2_R,hpsi_test2_C
END INTERFACE

CONTAINS

!-----------------------------------------------------------------------------------------------------------------------------------
subroutine hpsi_test2_R(tpsi0,htpsi0,iob,ik,nn,isub)
  use structures
  use salmon_parallel, only: nproc_group_korbital
  use init_sendrecv_sub
  use hpsi_sub
  implicit none
  real(8) :: tpsi0(iwksta(1):iwkend(1),iwksta(2):iwkend(2),iwksta(3):iwkend(3))
  real(8) :: htpsi0(iwk3sta(1):iwk3end(1),iwk3sta(2):iwk3end(2),iwk3sta(3):iwk3end(3))
  integer :: iob,nn,isub,ik

  integer :: ix,iy,iz,is,i_all,Norb,i,iobmax,Nspin,ind,j
  type(s_wf_info) :: info
  type(s_rgrid)   :: rg
  type(s_stencil) :: stencil
  type(s_wavefunction) :: tpsi, htpsi
  type(s_scalar),allocatable :: V(:)

  stencil%lap0 = 0.5d0*cNmat(0,Nd)*(1.d0/Hgs(1)**2+1.d0/Hgs(2)**2+1.d0/Hgs(3)**2)
  do j=1,3
    do ind=1,4
      stencil%lapt(ind,j) = cNmat(ind,4)/Hgs(j)**2
      stencil%nabt(ind,j) = bNmat(ind,4)/Hgs(j)
    end do
  end do

  Nspin = numspin
  call set_ispin(iob,is)

  rg%is_array(1) = iwk3sta(1)
  rg%ie_array(1) = iwk3end(1)
  rg%is_array(2) = iwk3sta(2)
  rg%ie_array(2) = iwk3end(2)
  rg%is_array(3) = iwk3sta(3)
  rg%ie_array(3) = iwk3end(3)

  rg%is(1) = mg_sta(1)
  rg%ie(1) = mg_end(1)
  rg%is(2) = mg_sta(2)
  rg%ie(2) = mg_end(2)
  rg%is(3) = mg_sta(3)
  rg%ie(3) = mg_end(3)

  rg%is_overlap = rg%is - 4
  rg%ie_overlap = rg%ie + 4

  allocate(rg%idx(rg%is_overlap(1):rg%ie_overlap(1)) &
          ,rg%idy(rg%is_overlap(2):rg%ie_overlap(2)) &
          ,rg%idz(rg%is_overlap(3):rg%ie_overlap(3)))
  do j=rg%is_overlap(1),rg%ie_overlap(1)
    rg%idx(j) = j
  end do
  do j=rg%is_overlap(2),rg%ie_overlap(2)
    rg%idy(j) = j
  end do
  do j=rg%is_overlap(3),rg%ie_overlap(3)
    rg%idz(j) = j
  end do

  info%io_s = 1
  info%io_e = 1
  info%numo = 1
  info%ik_s = 1
  info%ik_e = 1
  info%numk = 1
  info%i1_s = 1
  info%i1_e = 1
  info%num1 = 1
  info%if_divide_rspace = nproc_Mxin_mul.ne.1
  info%irank_overlap(1) = iup_array(1)
  info%irank_overlap(2) = idw_array(1)
  info%irank_overlap(3) = jup_array(1)
  info%irank_overlap(4) = jdw_array(1)
  info%irank_overlap(5) = kup_array(1)
  info%irank_overlap(6) = kdw_array(1)
  info%icomm_overlap = nproc_group_korbital
  info%icomm_pseudo = nproc_group_korbital

  allocate(tpsi%rwf(rg%is_array(1):rg%ie_array(1),rg%is_array(2):rg%ie_array(2),rg%is_array(3):rg%ie_array(3) &
          ,Nspin,1,1,1) &
         ,htpsi%rwf(rg%is_array(1):rg%ie_array(1),rg%is_array(2):rg%ie_array(2),rg%is_array(3):rg%ie_array(3) &
          ,Nspin,1,1,1))
  call set_ispin(iob,is)
  do iz=rg%is(3),rg%ie(3)
  do iy=rg%is(3),rg%ie(3)
  do ix=rg%is(3),rg%ie(3)
    tpsi%zwf(ix,iy,iz,is,1,1,1) = tpsi0(ix,iy,iz)
  end do
  end do
  end do

  allocate(V(Nspin))
  do j=1,Nspin
    allocate(V(j)%f(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)))
    V(j)%f = Vlocal(:,:,:,j)
  end do

  call hpsi(tpsi,htpsi,info,rg,V,Nspin,stencil,ppg)

  htpsi0 = htpsi%rwf(:,:,:,is,1,1,1)

  call deallocate_rgrid(rg)
  call deallocate_stencil(stencil)
  call deallocate_wavefunction(tpsi)
  call deallocate_wavefunction(htpsi)
  do j=1,Nspin
    call deallocate_scalar(V(j))
  end do
  deallocate(V)
end subroutine hpsi_test2_R

subroutine hpsi_test2_C(tpsi0,htpsi0,iob,ik,nn,isub)
  use structures
  use salmon_parallel, only: nproc_group_korbital
  use init_sendrecv_sub
  use hpsi_sub
  implicit none
  complex(8) :: tpsi0(iwksta(1):iwkend(1),iwksta(2):iwkend(2),iwksta(3):iwkend(3))
  complex(8) :: htpsi0(iwk3sta(1):iwk3end(1),iwk3sta(2):iwk3end(2),iwk3sta(3):iwk3end(3))
  integer,intent(in) :: iob,nn,isub,ik

  integer :: ix,iy,iz,is,i_all,Norb,i,iobmax,Nspin,ind,j,iatom,ikoa,jj
  real(8) :: x,y,z
  complex(8),parameter :: zi=(0d0,1d0)
  type(s_wf_info) :: info
  type(s_rgrid)   :: rg
  type(s_stencil) :: stencil
  type(s_wavefunction) :: tpsi, htpsi
  type(s_scalar),allocatable :: V(:)

  complex(8) :: ekr(maxMps,MI)

  allocate(stencil%kAc(ik:ik,3))
  stencil%lap0 = 0.5d0*ksquare(ik) -0.5d0*cNmat(0,Nd)*(1.d0/Hgs(1)**2+1.d0/Hgs(2)**2+1.d0/Hgs(3)**2)
  do j=1,3
    do ind=1,4
      stencil%lapt(ind,j) = cNmat(ind,4)/Hgs(j)**2
      stencil%nabt(ind,j) = bNmat(ind,4)/Hgs(j)
    end do
    stencil%kAc(ik,j) = k_rd(j,ik)
  end do

  Nspin = numspin
  call set_ispin(iob,is)

  rg%is_array(1) = iwk3sta(1)
  rg%ie_array(1) = iwk3end(1)
  rg%is_array(2) = iwk3sta(2)
  rg%ie_array(2) = iwk3end(2)
  rg%is_array(3) = iwk3sta(3)
  rg%ie_array(3) = iwk3end(3)

  rg%is(1) = mg_sta(1)
  rg%ie(1) = mg_end(1)
  rg%is(2) = mg_sta(2)
  rg%ie(2) = mg_end(2)
  rg%is(3) = mg_sta(3)
  rg%ie(3) = mg_end(3)

  rg%is_overlap = rg%is - 4
  rg%ie_overlap = rg%ie + 4

  allocate(rg%idx(rg%is_overlap(1):rg%ie_overlap(1)) &
          ,rg%idy(rg%is_overlap(2):rg%ie_overlap(2)) &
          ,rg%idz(rg%is_overlap(3):rg%ie_overlap(3)))
  do j=rg%is_overlap(1),rg%ie_overlap(1)
    rg%idx(j) = j
  end do
  do j=rg%is_overlap(2),rg%ie_overlap(2)
    rg%idy(j) = j
  end do
  do j=rg%is_overlap(3),rg%ie_overlap(3)
    rg%idz(j) = j
  end do

  info%io_s = 1
  info%io_e = 1
  info%numo = 1
  info%ik_s = ik
  info%ik_e = ik
  info%numk = 1
  info%i1_s = 1
  info%i1_e = 1
  info%num1 = 1
  info%if_divide_rspace = nproc_Mxin_mul.ne.1
  info%irank_overlap(1) = iup_array(1)
  info%irank_overlap(2) = idw_array(1)
  info%irank_overlap(3) = jup_array(1)
  info%irank_overlap(4) = jdw_array(1)
  info%irank_overlap(5) = kup_array(1)
  info%irank_overlap(6) = kdw_array(1)
  info%icomm_overlap = nproc_group_korbital
  info%icomm_pseudo = nproc_group_korbital

  allocate(tpsi%zwf(rg%is_array(1):rg%ie_array(1),rg%is_array(2):rg%ie_array(2),rg%is_array(3):rg%ie_array(3) &
          ,Nspin,1,ik:ik,1) &
         ,htpsi%zwf(rg%is_array(1):rg%ie_array(1),rg%is_array(2):rg%ie_array(2),rg%is_array(3):rg%ie_array(3) &
          ,Nspin,1,ik:ik,1))
  call set_ispin(iob,is)
  do iz=rg%is(3),rg%ie(3)
  do iy=rg%is(3),rg%ie(3)
  do ix=rg%is(3),rg%ie(3)
    tpsi%zwf(ix,iy,iz,is,1,ik,1) = tpsi0(ix,iy,iz)
  end do
  end do
  end do

  allocate(V(Nspin))
  do j=1,Nspin
    allocate(V(j)%f(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)))
    V(j)%f = Vlocal(:,:,:,j)
  end do

  if(iperiodic==3)then
    do iatom=1,MI
      ikoa=Kion(iatom)
      do jj=1,Mps(iatom)
        x=(dble(Jxyz(1,jj,iatom)-1)-dble(Jxxyyzz(1,jj,iatom)*lg_num(1)))*Hgs(1)
        y=(dble(Jxyz(2,jj,iatom)-1)-dble(Jxxyyzz(2,jj,iatom)*lg_num(2)))*Hgs(2)
        z=(dble(Jxyz(3,jj,iatom)-1)-dble(Jxxyyzz(3,jj,iatom)*lg_num(3)))*Hgs(3)
        ekr(jj,iatom)=exp(zi*(k_rd(1,ik)*x+k_rd(2,ik)*y+k_rd(3,ik)*z))
      end do
    end do
  end if
  call convert_pseudo_GCEED(ppg,ik,ik,ekr)

  call hpsi(tpsi,htpsi,info,rg,V,Nspin,stencil,ppg)

  call set_ispin(iob,is)
  htpsi0 = htpsi%zwf(:,:,:,is,1,ik,1)

  call deallocate_rgrid(rg)
  call deallocate_stencil(stencil)
  call deallocate_wavefunction(tpsi)
  call deallocate_wavefunction(htpsi)
  do j=1,Nspin
    call deallocate_scalar(V(j))
  end do
  deallocate(V)

end subroutine hpsi_test2_C

subroutine convert_pseudo_GCEED(ppg,ik_s,ik_e,ekr)
  use structures
  use scf_data, only: MI,Kion,Mlps,uVu,iwk_size,uV_all,Jxyz,uVu,Hvol,Mps,iperiodic ! GCEED
  type(s_pp_grid) :: ppg
  integer :: ik_s,ik_e
  complex(kind=8) :: ekr(maxMps,MI,ik_s:ik_e)
  !
  integer :: jj,iatom,ik,lm,ilma,jm,NI,Nlma,nps

  nps = ppg%nps
  nlma = ppg%nlma

  if(iperiodic==3)then
    if(.not.allocated(ppg%zproj)) allocate(ppg%zproj(Nps,Nlma,ik_s:ik_e))
    do ik=ik_s,ik_e
      do ilma=1,Nlma
        iatom = ppg%ia_tbl(ilma)
        do jj=1,ppg%Mps(iatom)
          ppg%zproj(jj,ilma,ik) = conjg(ekr(jj,iatom,ik)) * ppg%uv(jj,ilma)
        end do
      end do
    end do
  end if

  return
end subroutine convert_pseudo_GCEED
!-----------------------------------------------------------------------------------------------------------------------------------

END MODULE hpsi2_sub

