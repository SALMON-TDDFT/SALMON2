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
subroutine taylor(tzpsi_in,tzpsi_out,htpsi)
use calc_allob_sub
use scf_data
use allocate_mat_sub
use deallocate_mat_sub
implicit none

integer :: nn,ix,iy,iz
complex(8) :: tzpsi_in(mg_sta(1)-Nd:mg_end(1)+Nd+1,    &
                   mg_sta(2)-Nd:mg_end(2)+Nd,    &
                   mg_sta(3)-Nd:mg_end(3)+Nd,1:iobnum,k_sta:k_end)
complex(8) :: htpsi(mg_sta(1)-Nd:mg_end(1)+Nd+1,    &
                    mg_sta(2)-Nd:mg_end(2)+Nd,    &
                    mg_sta(3)-Nd:mg_end(3)+Nd,1:iobnum,k_sta:k_end)
complex(8) :: tzpsi_out(mg_sta(1)-Nd:mg_end(1)+Nd+1,    &
                   mg_sta(2)-Nd:mg_end(2)+Nd,    &
                   mg_sta(3)-Nd:mg_end(3)+Nd,1:iobnum,k_sta:k_end)


iwk_size=2
call make_iwksta_iwkend

if(ilsda==0)then
!$OMP parallel do private(iz,iy,ix)
  do iz=mg_sta(3),mg_end(3)
  do iy=mg_sta(2),mg_end(2)
  do ix=mg_sta(1),mg_end(1)
    rhobox(ix,iy,iz) = 0.d0
  end do
  end do
  end do
else if(ilsda==1)then
!$OMP parallel do private(iz,iy,ix)
  do iz=mg_sta(3),mg_end(3)
  do iy=mg_sta(2),mg_end(2)
  do ix=mg_sta(1),mg_end(1)
    rhobox_s(ix,iy,iz,1) = 0.d0
    rhobox_s(ix,iy,iz,2) = 0.d0
  end do
  end do
  end do
end if

if(ihpsieff==1)then
  if(ilsda==0)then
!$OMP parallel do private(iz,iy,ix)
    do iz=mg_sta(3),mg_end(3)
    do iy=mg_sta(2),mg_end(2)
    do ix=mg_sta(1),mg_end(1)
      Vlocal2(ix,iy,iz,1)=Vlocal(ix,iy,iz,1)+Vbox(ix,iy,iz)
    end do
    end do
    end do
  else if(ilsda==1)then 
!$OMP parallel do private(iz,iy,ix)
    do iz=mg_sta(3),mg_end(3)
    do iy=mg_sta(2),mg_end(2)
    do ix=mg_sta(1),mg_end(1)
      Vlocal2(ix,iy,iz,1)=Vlocal(ix,iy,iz,1)+Vbox(ix,iy,iz)
      Vlocal2(ix,iy,iz,2)=Vlocal(ix,iy,iz,2)+Vbox(ix,iy,iz)
    end do
    end do
    end do
  end if
end if

!do nn=1,N_hamil
!  if(ihpsieff==1)then
!    if(mod(nn,2)==1)then
!      call hpsi_groupob(tzpsi_in,htpsi,tzpsi_out,Vlocal2,nn,1)
!    else
!      call hpsi_groupob(htpsi,tzpsi_in,tzpsi_out,Vlocal2,nn,1)
!    end if
!  else
!    if(mod(nn,2)==1)then
!      call hpsi_groupob(tzpsi_in,htpsi,tzpsi_out,Vlocal,nn,1)
!    else
!      call hpsi_groupob(htpsi,tzpsi_in,tzpsi_out,Vlocal,nn,1)
!    end if
!  end if
!end do

!-----------------------------------------------------------------------------------------------------------------------------------
do nn=1,N_hamil
  if(mod(nn,2)==1) then
    if(ihpsieff==1) then
      call hpsi_test1(tzpsi_in,htpsi,Vlocal2)
    else
      call hpsi_test1(tzpsi_in,htpsi,Vlocal)
    end if
    call mode_add_polynomial(tzpsi_in,htpsi,tzpsi_out,nn)
  else
    if(ihpsieff==1) then
      call hpsi_test1(htpsi,tzpsi_in,Vlocal2)
    else
      call hpsi_test1(htpsi,tzpsi_in,Vlocal)
    end if
    call mode_add_polynomial(htpsi,tzpsi_in,tzpsi_out,nn)
  end if
end do
!-----------------------------------------------------------------------------------------------------------------------------------

contains

!-----------------------------------------------------------------------------------------------------------------------------------
subroutine hpsi_test1(tpsi0,htpsi0,V0)
  use structures
  use salmon_parallel, only: nproc_group_korbital
  use hpsi_sub
  use init_sendrecv_sub
  use hpsi2_sub, only: convert_pseudo_GCEED
  implicit none
  complex(8) :: tpsi0(iwk2sta(1):iwk2end(1)+1,iwk2sta(2):iwk2end(2),iwk2sta(3):iwk2end(3),   &
                   1:iobnum,k_sta:k_end)
  complex(8) :: htpsi0(iwk2sta(1):iwk2end(1)+1,iwk2sta(2):iwk2end(2),iwk2sta(3):iwk2end(3),   &
                   1:iobnum,k_sta:k_end)
  real(8) :: V0(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),numspin)
  !
  integer :: is,i_all,Norb,i,iobmax,Nspin,ik,ind,j,iatom,ikoa,jj
  real(8) :: x,y,z
  complex(8),parameter :: zi=(0d0,1d0)
  complex(8) :: ekr(maxMps,MI,k_sta:k_end)
  type(s_wf_info) :: info
  type(s_rgrid)   :: rg
  type(s_stencil) :: stencil
  type(s_wavefunction) :: tpsi, htpsi
  type(s_scalar),allocatable :: V(:)

  allocate(stencil%kAc(k_sta:k_end,3))
  stencil%lap0 = 0.5d0*ksquare(ik) -0.5d0*cNmat(0,Nd)*(1.d0/Hgs(1)**2+1.d0/Hgs(2)**2+1.d0/Hgs(3)**2)
  do j=1,3
    do ind=1,4
      stencil%lapt(ind,j) = cNmat(ind,4)/Hgs(j)**2
      stencil%nabt(ind,j) = bNmat(ind,4)/Hgs(j)
    end do
    do ik=k_sta,k_end
      stencil%kAc(ik,j) = k_rd(j,ik)
    end do
  end do

  call calc_pmax(iobmax)
  Nspin = numspin

  rg%is_array(1) = iwk2sta(1)
  rg%ie_array(1) = iwk2end(1)+1
  rg%is_array(2) = iwk2sta(2)
  rg%ie_array(2) = iwk2end(2)
  rg%is_array(3) = iwk2sta(3)
  rg%ie_array(3) = iwk2end(3)

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
  info%io_e = iobnum
  info%numo = iobnum
  info%ik_s = k_sta
  info%ik_e = k_end
  info%numk = k_num
  info%im_s = 1
  info%im_e = 1
  info%numm = 1
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
          ,Nspin,iobnum,k_sta:k_end,1) &
         ,htpsi%zwf(rg%is_array(1):rg%ie_array(1),rg%is_array(2):rg%ie_array(2),rg%is_array(3):rg%ie_array(3) &
          ,Nspin,iobnum,k_sta:k_end,1))
  do ik=k_sta,k_end
    do i=1,iobnum
      call calc_allob(i,i_all,iparaway_ob,itotmst,mst,iobnum)
      call set_is(i_all,is)
      tpsi%zwf(:,:,:,is,i,ik,1) = tpsi0(:,:,:,i,ik)
    end do
  end do

  allocate(V(Nspin))
  do j=1,Nspin
    allocate(V(j)%f(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)))
    V(j)%f = V0(:,:,:,j)
  end do

  if(iperiodic==3)then
    do iatom=1,MI
      ikoa=Kion(iatom)
      do jj=1,Mps(iatom)
        x=(dble(Jxyz(1,jj,iatom)-1)-dble(Jxxyyzz(1,jj,iatom)*lg_num(1)))*Hgs(1)
        y=(dble(Jxyz(2,jj,iatom)-1)-dble(Jxxyyzz(2,jj,iatom)*lg_num(2)))*Hgs(2)
        z=(dble(Jxyz(3,jj,iatom)-1)-dble(Jxxyyzz(3,jj,iatom)*lg_num(3)))*Hgs(3)
        do ik=k_sta,k_end
          ekr(jj,iatom,ik)=exp(zi*(k_rd(1,ik)*x+k_rd(2,ik)*y+k_rd(3,ik)*z))
        end do
      end do
    end do
  end if
  call convert_pseudo_GCEED(ppg,k_sta,k_end,ekr)

  call hpsi(tpsi,htpsi,info,rg,V,Nspin,stencil,srg,ppg)

  do ik=k_sta,k_end
    do i=1,iobnum
      call calc_allob(i,i_all,iparaway_ob,itotmst,mst,iobnum)
      call set_is(i_all,is)
      htpsi0(:,:,:,i,ik) = htpsi%zwf(:,:,:,is,i,ik,1)
    end do
  end do

  call deallocate_rgrid(rg)
  call deallocate_stencil(stencil)
  call deallocate_wavefunction(tpsi)
  call deallocate_wavefunction(htpsi)
  do j=1,Nspin
    call deallocate_scalar(V(j))
  end do
  deallocate(V)

  return
end subroutine hpsi_test1

! hpsi_groupob.f90 --> taylor.f90
subroutine mode_add_polynomial(tpsi,htpsi,tpsi_out,nn)
  implicit none
  complex(8) :: tpsi(mg_sta(1)-Nd:mg_end(1)+Nd+1,    &
                     mg_sta(2)-Nd:mg_end(2)+Nd,    &
                     mg_sta(3)-Nd:mg_end(3)+Nd,1:iobnum,k_sta:k_end)
  complex(8) :: htpsi(mg_sta(1)-Nd:mg_end(1)+Nd+1,    &
                      mg_sta(2)-Nd:mg_end(2)+Nd,    &
                      mg_sta(3)-Nd:mg_end(3)+Nd,1:iobnum,k_sta:k_end)
  complex(8) :: tpsi_out(mg_sta(1)-Nd:mg_end(1)+Nd+1,    &
                     mg_sta(2)-Nd:mg_end(2)+Nd,    &
                     mg_sta(3)-Nd:mg_end(3)+Nd,1:iobnum,k_sta:k_end)
  integer :: nn
  !
  integer :: iobmax
  call calc_pmax(iobmax)

  if(N_hamil==1)then
    if(ikind_eext==0)then
      call add_polynomial(tpsi,htpsi,tpsi_out,iobmax,nn,5)
    else
      call add_polynomial(tpsi,htpsi,tpsi_out,iobmax,nn,4)
    end if
  else if(N_hamil==4)then
    if(nn==1)then
      if(ikind_eext==0)then
        call add_polynomial(tpsi,htpsi,tpsi_out,iobmax,nn,1)
      else
        call add_polynomial(tpsi,htpsi,tpsi_out,iobmax,nn,0)
      end if
    else if(nn==3)then
      call add_polynomial(tpsi,htpsi,tpsi_out,iobmax,nn,-1)
    else if(nn==4)then
      call add_polynomial(tpsi,htpsi,tpsi_out,iobmax,nn,6)
    end if
  else
    if(nn==1)then
      if(ikind_eext==0)then
        call add_polynomial(tpsi,htpsi,tpsi_out,iobmax,nn,1)
      else
        call add_polynomial(tpsi,htpsi,tpsi_out,iobmax,nn,0)
      end if
    else if(nn==N_hamil)then
      call add_polynomial(tpsi,htpsi,tpsi_out,iobmax,nn,6)
    else
      call add_polynomial(tpsi,htpsi,tpsi_out,iobmax,nn,2)
    end if
  end if

  return
end subroutine
!-----------------------------------------------------------------------------------------------------------------------------------

end subroutine taylor

