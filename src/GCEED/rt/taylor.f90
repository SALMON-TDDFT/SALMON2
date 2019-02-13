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
subroutine taylor(mg,info,tzpsi_in,tzpsi_out,htpsi)
  use structures, only: s_rgrid,s_wf_info,s_wavefunction,s_stencil,s_scalar
  use salmon_parallel, only: nproc_group_korbital
  use scf_data
  use init_sendrecv_sub, only: iup_array,idw_array,jup_array,jdw_array,kup_array,kdw_array
  use allocate_mat_sub
  use deallocate_mat_sub
  use hpsi_sub
  implicit none
  
  type(s_rgrid),intent(in) :: mg
  type(s_wf_info),intent(in) :: info
  type(s_wavefunction) :: stpsi_in_ob
  type(s_wavefunction) :: stpsi_out_ob
  type(s_wavefunction) :: shtpsi_ob
  type(s_wf_info) :: info_ob
  type(s_stencil) :: stencil
  type(s_scalar),allocatable :: v(:)
  integer :: nspin
  integer :: nn,ix,iy,iz
  integer :: j,ind
  integer :: ik,iob,iob_allob
  complex(8) :: tzpsi_in(mg%is_array(1):mg%ie_array(1),  &
                         mg%is_array(2):mg%ie_array(2),  &
                         mg%is_array(3):mg%ie_array(3),1:info%numo,info%ik_s:info%ik_e)
  complex(8) :: htpsi(mg%is_array(1):mg%ie_array(1),  &
                      mg%is_array(1):mg%ie_array(1),  &
                      mg%is_array(3):mg%ie_array(3),1:info%numo,info%ik_s:info%ik_e)
  complex(8) :: tzpsi_out(mg%is_array(1):mg%ie_array(1),  &
                          mg%is_array(2):mg%ie_array(2),  &
                          mg%is_array(3):mg%ie_array(3),1:info%numo,info%ik_s:info%ik_e)
  complex(8) :: ekr(ppg%nps,natom)
  integer :: a,iatom
  integer :: ilma
  real(8) :: x,y,z
  complex(8),parameter :: zi=(0.d0,1.d0)
  
  info_ob%im_s = 1
  info_ob%im_e = 1
  info_ob%numm = 1
  info_ob%ik_s = 1
  info_ob%ik_e = 1
  info_ob%numk = 1
  info_ob%io_s = 1
  info_ob%io_e = 1
  info_ob%numo = 1
  info_ob%if_divide_rspace = nproc_mxin_mul.ne.1
  info_ob%irank_r(1) = iup_array(1)
  info_ob%irank_r(2) = idw_array(1)
  info_ob%irank_r(3) = jup_array(1)
  info_ob%irank_r(4) = jdw_array(1)
  info_ob%irank_r(5) = kup_array(1)
  info_ob%irank_r(6) = kdw_array(1)
  info_ob%icomm_r = nproc_group_korbital
  
  allocate(stpsi_in_ob%zwf(mg%is_array(1):mg%ie_array(1),  &
                           mg%is_array(2):mg%ie_array(2),  &
                           mg%is_array(3):mg%ie_array(3),1,1,1,1))
  allocate(stpsi_out_ob%zwf(mg%is_array(1):mg%ie_array(1),  &
                            mg%is_array(2):mg%ie_array(2),  &
                            mg%is_array(3):mg%ie_array(3),1,1,1,1))
  allocate(shtpsi_ob%zwf(mg%is_array(1):mg%ie_array(1),  &
                         mg%is_array(2):mg%ie_array(2),  &
                         mg%is_array(3):mg%ie_array(3),1,1,1,1))

  if(iperiodic==3) allocate(stencil%kAc(1:1,3))

  stencil%lap0 = -0.5d0*cNmat(0,nd)*(1.d0/hgs(1)**2+1.d0/hgs(2)**2+1.d0/hgs(3)**2)

  if(iperiodic==0)then
    do j=1,3
      do ind=1,4
        stencil%lapt(ind,j) = cnmat(ind,4)/hgs(j)**2
        stencil%nabt(ind,j) = 0.d0
      end do
    end do
  else if(iperiodic==3)then
    do j=1,3
      do ind=1,4
        stencil%lapt(ind,j) = cnmat(ind,4)/hgs(j)**2
        stencil%nabt(ind,j) = bnmat(ind,4)/hgs(j)
      end do
    end do
  end if

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
    call calc_allob(iob,iob_allob)
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
      stpsi_in_ob%zwf(ix,iy,iz,1,1,1,1)=tzpsi_in(ix,iy,iz,iob,ik)
    end do
    end do
    end do
!$OMP parallel do private(iz,iy,ix)
    do iz=mg%is_array(3),mg%ie_array(3)
    do iy=mg%is_array(2),mg%ie_array(2)
    do ix=mg%is_array(1),mg%ie_array(1)
      stpsi_out_ob%zwf(ix,iy,iz,1,1,1,1)=tzpsi_in(ix,iy,iz,iob,ik)
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

    do nn=1,N_hamil
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
      tzpsi_out(ix,iy,iz,iob,ik)=stpsi_out_ob%zwf(ix,iy,iz,1,1,1,1)
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
  if(iperiodic==3) deallocate(stencil%kAc)
  deallocate(v(1)%f)
  deallocate(v)
  if(allocated(ppg%zproj)) deallocate(ppg%zproj)

end subroutine taylor

