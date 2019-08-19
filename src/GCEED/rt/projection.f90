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
subroutine projection(mg,info,tzpsi)
use salmon_parallel, only: nproc_group_global, nproc_id_global
use salmon_communication, only: comm_is_root, comm_summation, comm_bcast
use calc_allob_sub
use calc_iroot_sub
use calc_myob_sub
use check_corrkob_sub
use scf_data
use new_world_sub
use allocate_mat_sub
use structures, only: s_rgrid, s_orbital_parallel
implicit none
type(s_rgrid),intent(in) :: mg
type(s_orbital_parallel) :: info
complex(8) :: tzpsi(mg%is_array(1):mg%ie_array(1) &
&                  ,mg%is_array(2):mg%ie_array(2) &
&                  ,mg%is_array(3):mg%ie_array(3), 1:iobnum, k_sta:k_end)
integer :: ix,iy,iz,iob,iik
integer :: iob_myob,icorr_p,job,job_allob
complex(8) :: coef_mat(itotMST,itotMST0,num_kpoints_rd,1)
complex(8) :: coef_mat2(itotMST,itotMST0,num_kpoints_rd,1)
real(8) :: coef(itotMST0,num_kpoints_rd,1)
complex(8) :: cbox
integer :: iobmax
integer :: iroot
complex(8),parameter :: zi=(0.d0,1.d0)
character(100) :: projOutFile

call calc_pmax(iobmax)


if(iSCFRT==2)then
  if(iwrite_projnum==1)then
    write(fileNumber, '(i8)') itt
    projOutFile = trim("proj.")//adjustl(fileNumber)
    open(61,file=projOutFile)
  end if
end if

coef_mat=0.d0

do iik=k_sta,k_end
do iob=1,itotMST0
  call calc_myob(iob,iob_myob,ilsda,nproc_ob,itotmst,mst)
  call check_corrkob(iob,iik,icorr_p,ilsda,nproc_ob,k_sta,k_end,mst)
  if(icorr_p==1)then
!$OMP parallel do private(iz,iy,ix)
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      cmatbox_m(ix,iy,iz)=zpsi_t0(ix,iy,iz,iob_myob,iik)
    end do
    end do
    end do
  end if
  call calc_iroot(iob,iroot,ilsda,nproc_ob,itotmst,mst)
  call comm_bcast(cmatbox_m,info%icomm_o,iroot)
  do job=1,iobmax
    cbox=0.d0
!$OMP parallel do reduction(+:cbox) private(iz,iy,ix)
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      cbox=cbox+conjg(tzpsi(ix,iy,iz,job,iik))*cmatbox_m(ix,iy,iz)
    end do
    end do
    end do
    call calc_allob(job,job_allob,itotmst,mst,iobnum)
    coef_mat(job_allob,iob,iik,1)=coef_mat(job_allob,iob,iik,1)+cbox
  end do
end do
end do

call comm_summation(coef_mat,coef_mat2,itotMST*itotMST0*num_kpoints_rd,nproc_group_global)

coef=0.d0
do iik=1,num_kpoints_rd
do iob=1,itotMST0
  do job=1,itotMST
    coef(iob,iik,1)=coef(iob,iik,1)+abs(coef_mat2(job,iob,iik,1)*Hvol)**2
  end do
end do
end do
if(comm_is_root(nproc_id_global))then
  write(41,'(200f14.8)') dble(itt)*dt*2.41888d-2, &
  & (coef(iwrite_projection_ob(iob),iwrite_projection_k(iob),1),iob=1,num_projection),  &
    sum(coef(1:itotMST,:,1)),sum(coef(1:itotMST0,:,1))
end if
if(mod(itt,100)==0)then
  if(comm_is_root(nproc_id_global))then
    do iik=1,num_kpoints_rd
    do iob=1,itotMST0
      write(*,'(a12,3i6,f16.8)') "projection",iob,iik,coef(iob,iik,1)
    end do
    end do
  end if
end if

if(iSCFRT==2)then
  if(iwrite_projnum==1)then
    close(61)
  end if
end if

end subroutine projection
