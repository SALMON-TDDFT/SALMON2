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
module taylor_sub
  implicit none

contains

subroutine taylor(mg,itotmst,mst,lg_sta,lg_end,ilsda,info,info_ob,stencil,tspsi_in,tspsi_out,  &
                  ppg,vlocal,vbox,num_kpoints_rd,k_rd,rhobox,rhobox_s,zc,ihpsieff,rocc,wtk,iparaway_ob)
  use inputoutput, only: iperiodic,ispin,natom,n_hamil
  use structures, only: s_rgrid,s_wf_info,s_wavefunction,s_stencil,s_scalar,s_pp_grid
  use hpsi_sub
  implicit none
  integer,parameter     :: nd=4 
  type(s_rgrid),intent(in) :: mg
  integer,intent(in) :: itotmst
  integer,intent(in) :: mst(2)
  integer,intent(in) :: lg_sta(3)
  integer,intent(in) :: lg_end(3)
  integer,intent(in)    :: ilsda
  type(s_wf_info),intent(in) :: info
  type(s_wf_info),intent(inout) :: info_ob
  type(s_stencil),intent(inout) :: stencil
  type(s_wavefunction),intent(inout) :: tspsi_in
  type(s_wavefunction),intent(inout) :: tspsi_out
  type(s_pp_grid),intent(inout) :: ppg
  real(8),intent(in)    :: vlocal(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),ispin+1)
  real(8),intent(in)    :: vbox(lg_sta(1)-nd:lg_end(1)+nd,  &
                                lg_sta(2)-nd:lg_end(2)+nd,  &
                                lg_sta(3)-nd:lg_end(3)+nd)
  integer,intent(in)    :: num_kpoints_rd
  real(8),intent(in)    :: k_rd(3,num_kpoints_rd)
  real(8),intent(out)   :: rhobox(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
  real(8),intent(out)   :: rhobox_s(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),2)
  complex(8),intent(in) :: zc(n_hamil)
  integer,intent(in)    :: ihpsieff
  real(8),intent(in)    :: rocc(itotmst,num_kpoints_rd)
  real(8),intent(in)    :: wtk(num_kpoints_rd)
  integer,intent(in)    :: iparaway_ob
  type(s_wavefunction) :: stpsi_in_ob
  type(s_wavefunction) :: stpsi_out_ob
  type(s_wavefunction) :: shtpsi_ob
  type(s_scalar),allocatable :: v(:)
  real(8)              :: vlocal2(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),ispin+1)
  integer :: nspin
  integer :: nn,ix,iy,iz
  integer :: ik,iob,iob_allob
  complex(8) :: ekr(ppg%nps,natom)
  integer :: a,iatom
  integer :: ilma,j
  real(8) :: x,y,z
  complex(8),parameter :: zi=(0.d0,1.d0)
  
  allocate(stpsi_in_ob%zwf(mg%is_array(1):mg%ie_array(1),  &
                           mg%is_array(2):mg%ie_array(2),  &
                           mg%is_array(3):mg%ie_array(3),1,1,1,1))
  allocate(stpsi_out_ob%zwf(mg%is_array(1):mg%ie_array(1),  &
                            mg%is_array(2):mg%ie_array(2),  &
                            mg%is_array(3):mg%ie_array(3),1,1,1,1))
  allocate(shtpsi_ob%zwf(mg%is_array(1):mg%ie_array(1),  &
                         mg%is_array(2):mg%ie_array(2),  &
                         mg%is_array(3):mg%ie_array(3),1,1,1,1))

  nspin=1
  allocate(v(1))
  allocate(v(1)%f(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))

  if(ilsda==0)then
  !$OMP parallel do private(iz,iy,ix)
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      rhobox(ix,iy,iz) = 0.d0
    end do
    end do
    end do
  else if(ilsda==1)then
  !$OMP parallel do private(iz,iy,ix)
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      rhobox_s(ix,iy,iz,1) = 0.d0
      rhobox_s(ix,iy,iz,2) = 0.d0
    end do
    end do
    end do
  end if
  
  if(iperiodic==0.and.ihpsieff==1)then
    if(ilsda==0)then
  !$OMP parallel do private(iz,iy,ix)
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        Vlocal2(ix,iy,iz,1)=Vlocal(ix,iy,iz,1)+Vbox(ix,iy,iz)
      end do
      end do
      end do
    else if(ilsda==1)then 
  !$OMP parallel do private(iz,iy,ix)
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        Vlocal2(ix,iy,iz,1)=Vlocal(ix,iy,iz,1)+Vbox(ix,iy,iz)
        Vlocal2(ix,iy,iz,2)=Vlocal(ix,iy,iz,2)+Vbox(ix,iy,iz)
      end do
      end do
      end do
    end if
  end if
  

  do ik=info%ik_s,info%ik_e
  if(iperiodic==3)then
    if(.not.allocated(ppg%zproj)) allocate(ppg%zproj(ppg%nps,ppg%nlma,1:1))
    do a=1,natom
      do j=1,ppg%mps(a)
        x=ppg%rxyz(1,j,a)
        y=ppg%rxyz(2,j,a)
        z=ppg%rxyz(3,j,a)
        ekr(j,a)=exp(zi*(k_rd(1,ik)*x+k_rd(2,ik)*y+k_rd(3,ik)*z))
      end do
    end do
    do ilma=1,ppg%nlma
      iatom = ppg%ia_tbl(ilma)
      do j=1,ppg%mps(iatom)
        ppg%zproj(j,ilma,1) = conjg(ekr(j,iatom)) * ppg%uv(j,ilma)
      end do
    end do
    do j=1,3
      stencil%kAc(1,j) = k_rd(j,ik)
    end do
  end if
  do iob=info%io_s,info%io_e
    call calc_allob(iob,iob_allob,iparaway_ob,itotmst,mst,info%numo)
    if(iperiodic==0.and.ihpsieff==1)then
      if(iob_allob<=MST(1))then
!$OMP parallel do private(iz,iy,ix)
        do iz=mg%is(3),mg%ie(3)
        do iy=mg%is(2),mg%ie(2)
        do ix=mg%is(1),mg%ie(1)
          v(1)%f(ix,iy,iz)=vlocal2(ix,iy,iz,1)
        end do
        end do
        end do
      else
!$OMP parallel do private(iz,iy,ix)
        do iz=mg%is(3),mg%ie(3)
        do iy=mg%is(2),mg%ie(2)
        do ix=mg%is(1),mg%ie(1)
          v(1)%f(ix,iy,iz)=vlocal2(ix,iy,iz,2)
        end do
        end do
        end do
      end if
    else
      if(iob_allob<=MST(1))then
!$OMP parallel do private(iz,iy,ix)
        do iz=mg%is(3),mg%ie(3)
        do iy=mg%is(2),mg%ie(2)
        do ix=mg%is(1),mg%ie(1)
          v(1)%f(ix,iy,iz)=vlocal(ix,iy,iz,1)
        end do
        end do
        end do
      else
!$OMP parallel do private(iz,iy,ix)
        do iz=mg%is(3),mg%ie(3)
        do iy=mg%is(2),mg%ie(2)
        do ix=mg%is(1),mg%ie(1)
          v(1)%f(ix,iy,iz)=vlocal(ix,iy,iz,2)
        end do
        end do
        end do
      end if
    end if

!$OMP parallel do private(iz,iy,ix)
    do iz=mg%is_array(3),mg%ie_array(3)
    do iy=mg%is_array(2),mg%ie_array(2)
    do ix=mg%is_array(1),mg%ie_array(1)
      stpsi_in_ob%zwf(ix,iy,iz,1,1,1,1)=tspsi_in%zwf(ix,iy,iz,1,iob,ik,1)
    end do
    end do
    end do
!$OMP parallel do private(iz,iy,ix)
    do iz=mg%is_array(3),mg%ie_array(3)
    do iy=mg%is_array(2),mg%ie_array(2)
    do ix=mg%is_array(1),mg%ie_array(1)
      stpsi_out_ob%zwf(ix,iy,iz,1,1,1,1)=tspsi_in%zwf(ix,iy,iz,1,iob,ik,1)
    end do
    end do
    end do
!$OMP parallel do private(iz,iy,ix)
    do iz=mg%is_array(3),mg%ie_array(3)
    do iy=mg%is_array(2),mg%ie_array(2)
    do ix=mg%is_array(1),mg%ie_array(1)
      shtpsi_ob%zwf(ix,iy,iz,1,1,1,1)=0.d0
    end do
    end do
    end do

    do nn=1,n_hamil
      if(mod(nn,2)==1)then
        call hpsi(stpsi_in_ob,shtpsi_ob,info_ob,mg,v,nspin,stencil,ppg)
!$OMP parallel do private(iz,iy,ix)
        do iz=mg%is_array(3),mg%ie_array(3)
        do iy=mg%is_array(2),mg%ie_array(2)
        do ix=mg%is_array(1),mg%ie_array(1)
          stpsi_out_ob%zwf(ix,iy,iz,1,1,1,1)=stpsi_out_ob%zwf(ix,iy,iz,1,1,1,1)+ &
                                               zc(nn)*shtpsi_ob%zwf(ix,iy,iz,1,1,1,1)
        end do
        end do
        end do
      else
        call hpsi(shtpsi_ob,stpsi_in_ob,info_ob,mg,v,nspin,stencil,ppg)
!$OMP parallel do private(iz,iy,ix)
        do iz=mg%is_array(3),mg%ie_array(3)
        do iy=mg%is_array(2),mg%ie_array(2)
        do ix=mg%is_array(1),mg%ie_array(1)
          stpsi_out_ob%zwf(ix,iy,iz,1,1,1,1)=stpsi_out_ob%zwf(ix,iy,iz,1,1,1,1)+ &
                                               zc(nn)*stpsi_in_ob%zwf(ix,iy,iz,1,1,1,1)
        end do
        end do
        end do
      end if
    end do
!$OMP parallel do private(iz,iy,ix)
    do iz=mg%is_array(3),mg%ie_array(3)
    do iy=mg%is_array(2),mg%ie_array(2)
    do ix=mg%is_array(1),mg%ie_array(1)
      tspsi_out%zwf(ix,iy,iz,1,iob,ik,1)=stpsi_out_ob%zwf(ix,iy,iz,1,1,1,1)
    end do
    end do
    end do

    if(ilsda==0)then
!$OMP parallel do private(iz,iy,ix)
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        rhobox(ix,iy,iz)=rhobox(ix,iy,iz)+stpsi_out_ob%zwf(ix,iy,iz,1,1,1,1)&
                                           *conjg(stpsi_out_ob%zwf(ix,iy,iz,1,1,1,1))*rocc(iob_allob,ik)*wtk(ik)
      end do
      end do
      end do
    else if(ilsda==1)then
      if(iob_allob<=MST(1))then
!$OMP parallel do private(iz,iy,ix)
        do iz=mg%is(3),mg%ie(3)
        do iy=mg%is(2),mg%ie(2)
        do ix=mg%is(1),mg%ie(1)
          rhobox_s(ix,iy,iz,1)=rhobox_s(ix,iy,iz,1)+stpsi_out_ob%zwf(ix,iy,iz,1,1,1,1)&
                                           *conjg(stpsi_out_ob%zwf(ix,iy,iz,1,1,1,1))*rocc(iob_allob,ik)*wtk(ik)
        end do
        end do
        end do
      else
!$OMP parallel do private(iz,iy,ix)
        do iz=mg%is(3),mg%ie(3)
        do iy=mg%is(2),mg%ie(2)
        do ix=mg%is(1),mg%ie(1)
          rhobox_s(ix,iy,iz,2)=rhobox_s(ix,iy,iz,2)+stpsi_out_ob%zwf(ix,iy,iz,1,1,1,1)&
                                           *conjg(stpsi_out_ob%zwf(ix,iy,iz,1,1,1,1))*rocc(iob_allob,ik)*wtk(ik)
        end do
        end do
        end do
      end if
    end if
  end do
  end do

  deallocate(stpsi_in_ob%zwf,stpsi_out_ob%zwf,shtpsi_ob%zwf)
  deallocate(v(1)%f)
  deallocate(v)
  if(allocated(ppg%zproj)) deallocate(ppg%zproj)

end subroutine taylor

end module taylor_sub

