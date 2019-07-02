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
subroutine calc_force(mg)
use salmon_parallel, only: nproc_group_korbital, nproc_group_global
use salmon_communication, only: comm_is_root, comm_summation
use scf_data
use allocate_mat_sub
use read_pslfile_sub
use new_world_sub
use structures, only: s_rgrid
implicit none
type(s_rgrid),intent(in) :: mg
integer :: ix,iy,iz,iob,ikoa,jj,j2,iatom,ia,ib,lm,ikoa2
real(8) :: rbox1,rbox2
real(8),allocatable :: uVpsibox(:,:,:,:),uVpsibox2(:,:,:,:)
real(8) :: rforce1(3,MI),rforce2(3,MI),rforce3(3,MI)
real(8) :: rab
real(8) :: tpsi(mg%is_overlap(1):mg%ie_overlap(1) &
&              ,mg%is_overlap(2):mg%ie_overlap(2) &
&              ,mg%is_overlap(3):mg%ie_overlap(3), 1:iobnum, k_sta:k_end)

do iatom=1,MI
do j2=1,3
  rforce(j2,iatom)=0.d0
  rforce1(j2,iatom)=0.d0
  rforce2(j2,iatom)=0.d0
  rforce3(j2,iatom)=0.d0
end do
end do

! ion-ion
do ia=1,MI
  ikoa=Kion(ia)
  do ib=1,MI
    if(ia/=ib)then
      ikoa2=Kion(ib)
      rab=sqrt((Rion(1,ia)-Rion(1,ib))**2+   &
               (Rion(2,ia)-Rion(2,ib))**2+   &
               (Rion(3,ia)-Rion(3,ib))**2)
      do j2=1,3
        rforce(j2,ia)=rforce(j2,ia)+Zps(ikoa)*Zps(ikoa2)*(Rion(j2,ia)-Rion(j2,ib))/rab**3
        rforce1(j2,ia)=rforce1(j2,ia)+Zps(ikoa)*Zps(ikoa2)*(Rion(j2,ia)-Rion(j2,ib))/rab**3
      end do
    end if
  end do
end do

tpsi=0.d0
do iob=1,iobnum
!$OMP parallel do private(ix,iy,iz)
  do iz=mg%is(3),mg%ie(3)
  do iy=mg%is(2),mg%ie(2)
  do ix=mg%is(1),mg%ie(1)
    tpsi(ix,iy,iz,iob,1)=psi(ix,iy,iz,iob,1)
  end do
  end do
  end do
end do

call calc_gradient_fast(tpsi,rgrad_wk)

! local part of force
do iatom=1,MI
do j2=1,3
  rbox1=0.d0
  do iob=1,iobnum
!$OMP parallel do private(ix,iy,iz) reduction( + : rbox1 )
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      rbox1=rbox1-2.d0*rocc(iob,1)*rgrad_wk(ix,iy,iz,iob,1,j2)*Vpsl_atom(ix,iy,iz,iatom)*psi(ix,iy,iz,iob,1)
    end do
    end do
    end do
  end do
  call comm_summation(rbox1,rbox2,nproc_group_global)
  rforce(j2,iatom)=rforce(j2,iatom)+rbox2*Hvol
  rforce2(j2,iatom)=rbox2*Hvol
end do
end do

! nonlocal part of force

allocate (uVpsibox(1:iobnum,k_sta:k_end,1:maxlm,1:MI))
allocate (uVpsibox2(1:iobnum,k_sta:k_end,1:maxlm,1:MI))

do iatom=1,MI
  do lm=1,maxlm
    do iob=1,iobnum
      uVpsibox(iob,1,lm,iatom)=0.d0
    end do
  end do
end do

do iatom=1,MI
  ikoa=Kion(iatom)
  do iob=1,iobnum
    loop_lm2 : do lm=1,(Mlps(ikoa)+1)**2
      if ( abs(uVu(lm,iatom))<1.d-5 ) cycle loop_lm2
      rbox1=0.d0
!$OMP parallel do reduction( + : rbox1 )
      do jj=1,Mps(iatom)
        rbox1=rbox1+uV(jj,lm,iatom)*  &
                       psi(Jxyz(1,jj,iatom),Jxyz(2,jj,iatom),Jxyz(3,jj,iatom),iob,1)
      end do
      uVpsibox(iob,1,lm,iatom)=rbox1*Hvol/uVu(lm,iatom)
    end do loop_lm2
  end do
end do

call comm_summation(uVpsibox,uVpsibox2,iobnum*maxlm*MI,nproc_group_korbital)

do iatom=1,MI
  ikoa=Kion(iatom)
  do j2=1,3
    rbox1=0.d0
    do iob=1,iobnum
      do jj=1,Mps(iatom)
        do lm=1,(Mlps(ikoa)+1)**2
          rbox1=rbox1-2.d0*rocc(iob,1)*uV(jj,lm,iatom)*   &
                  rgrad_wk(Jxyz(1,jj,iatom),Jxyz(2,jj,iatom),Jxyz(3,jj,iatom),iob,1,j2)* &
                  uVpsibox2(iob,1,lm,iatom)
        end do
      end do
    end do
    call comm_summation(rbox1,rbox2,nproc_group_global)
    rforce(j2,iatom)=rforce(j2,iatom)+rbox2*Hvol
    rforce3(j2,iatom)=rbox2*Hvol
  end do
end do

deallocate(uVpsibox,uVpsibox2)

end subroutine calc_force
