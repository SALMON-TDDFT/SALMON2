!
!  Copyright 2018 SALMON developers
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
subroutine subspace_diag(mg)

use structures, only: s_rgrid
use salmon_parallel, only: nproc_group_kgrid, nproc_group_global, nproc_group_korbital
use salmon_communication, only: comm_summation, comm_bcast
use misc_routines, only: get_wtime
use scf_data
use inner_product_sub
use copy_psi_mesh_sub
implicit none
integer :: iob,job,ii,jj
integer :: ix,iy,iz,is
real(8),allocatable :: Amat(:,:)
real(8),allocatable :: Amat2(:,:)
real(8),allocatable :: Smat(:,:)
real(8),allocatable :: tpsi(:,:,:),htpsi(:,:,:)
real(8),allocatable :: psi_box(:,:,:,:)
real(8) :: rbox,rbox1
real(8),allocatable :: evec(:,:)
real(8),allocatable :: rmatbox_m(:,:,:)
integer :: ier2
integer :: is_sta,is_end
integer :: job_myob,iroot,icorr_j,iob_allob,job_allob
integer :: iter
integer :: iobsta(2),iobend(2)
type(s_rgrid) :: mg

elp3(301)=get_wtime()

allocate(tpsi(mg%is_array(1):mg%ie_array(1),  &
              mg%is_array(2):mg%ie_array(2),  &
              mg%is_array(3):mg%ie_array(3)))
allocate(htpsi(mg%is(1):mg%ie(1),  &
               mg%is(2):mg%ie(2),  &
               mg%is(3):mg%ie(3)))
allocate(psi_box(mg%is(1):mg%ie(1),  &
                 mg%is(2):mg%ie(2),  &
                 mg%is(3):mg%ie(3),1:iobnum))
allocate(rmatbox_m(mg%is(1):mg%ie(1),  &
                   mg%is(2):mg%ie(2),  &
                   mg%is(3):mg%ie(3)))

iwk_size=2
call make_iwksta_iwkend

call set_isstaend(is_sta,is_end,ilsda,nproc_ob,nproc_ob_spin)

if(ilsda == 0)then
  iobsta(1)=1
  iobend(1)=itotMST
else if(ilsda == 1)then
  iobsta(1)=1
  iobend(1)=MST(1)
  iobsta(2)=MST(1)+1
  iobend(2)=itotMST
end if

!$OMP parallel do private(iz,iy,ix)
do iz=mg%is_array(3),mg%ie_array(3)
do iy=mg%is_array(2),mg%ie_array(2)
do ix=mg%is_array(1),mg%ie_array(1)
  tpsi(ix,iy,iz)=0.d0
end do
end do
end do

do is=is_sta,is_end

  if(ifMST(is)>=1.and.MST(is)>=1)then

    iter=iobend(is)-iobsta(is)+1
    allocate(evec(iter,iter))
    allocate(Amat(iter,iter))
    allocate(Amat2(iter,iter))
    allocate(Smat(iter,iter))
  
    do jj=1,iter
!$OMP parallel do 
      do ii=1,iter
        Amat2(ii,jj)=0.d0
      end do
    end do
  
    do job=iobsta(is),iobend(is)
      call calc_myob(job,job_myob,ilsda,nproc_ob,iparaway_ob,itotmst,nproc_ob_spin,mst)
      call check_corrkob(job,1,icorr_j,ilsda,nproc_ob,iparaway_ob,itotmst,k_sta,k_end,nproc_ob_spin,mst)
      if(icorr_j==1)then
!$OMP parallel do private(iz,iy,ix)
        do iz=mg%is(3),mg%ie(3)
        do iy=mg%is(2),mg%ie(2)
        do ix=mg%is(1),mg%ie(1)
          tpsi(ix,iy,iz)=psi(ix,iy,iz,job_myob,1)
        end do
        end do
        end do
        call r_hpsi2_buf(tpsi,htpsi,job,1,0,0)
      end if
      call calc_iroot(job,iroot,ilsda,nproc_ob,iparaway_ob,itotmst,nproc_ob_spin,mst)
      call comm_bcast(htpsi,nproc_group_kgrid,iroot)
      
      do iob=1,iobnum
        call calc_allob(iob,iob_allob)
        if(iob_allob>=iobsta(is).and.iob_allob<=iobend(is))then
          rbox=0.d0
!$OMP parallel do reduction(+:rbox) private(iz,iy,ix)
          do iz=mg%is(3),mg%ie(3)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
            rbox=rbox+psi(ix,iy,iz,iob,1)*htpsi(ix,iy,iz)
          end do
          end do
          end do
          Amat2(iob_allob-iobsta(is)+1,job-iobsta(is)+1)=rbox*Hvol
        end if
      end do
    end do
    
    call comm_summation(Amat2,Amat,iter*iter,nproc_group_global)
  
    call eigen_subdiag(Amat,evec,iter,ier2)
   
    do job=1,iobnum
      call calc_allob(job,job_allob)
      if(job_allob>=iobsta(is).and.job_allob<=iobend(is))then
!$OMP parallel do private(iz,iy,ix)
        do iz=mg%is(3),mg%ie(3)
        do iy=mg%is(2),mg%ie(2)
        do ix=mg%is(1),mg%ie(1)
          psi_box(ix,iy,iz,job)=psi(ix,iy,iz,job,1)
          psi(ix,iy,iz,job,1)=0.d0
        end do
        end do
        end do
      end if
    end do
     
    do job=iobsta(is),iobend(is)
      call calc_myob(job,job_myob,ilsda,nproc_ob,iparaway_ob,itotmst,nproc_ob_spin,mst)
      call check_corrkob(job,1,icorr_j,ilsda,nproc_ob,iparaway_ob,itotmst,k_sta,k_end,nproc_ob_spin,mst)
      if(icorr_j==1)then
!$OMP parallel do private(iz,iy,ix)
        do iz=mg%is(3),mg%ie(3)
        do iy=mg%is(2),mg%ie(2)
        do ix=mg%is(1),mg%ie(1)
          rmatbox_m(ix,iy,iz)=psi_box(ix,iy,iz,job_myob)
        end do
        end do
        end do
      end if
      call calc_iroot(job,iroot,ilsda,nproc_ob,iparaway_ob,itotmst,nproc_ob_spin,mst)
      call comm_bcast(rmatbox_m,nproc_group_kgrid,iroot)
      do iob=1,iobnum
        call calc_allob(iob,iob_allob)
        if(iob_allob>=iobsta(is).and.iob_allob<=iobend(is))then
!$OMP parallel do private(iz,iy,ix)
          do iz=mg%is(3),mg%ie(3)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
            psi(ix,iy,iz,iob,1)=psi(ix,iy,iz,iob,1)+evec(job-iobsta(is)+1,iob_allob-iobsta(is)+1)*rmatbox_m(ix,iy,iz)
          end do
          end do
          end do
        end if
      end do
    end do
  
    do iob=1,iobnum
      call calc_allob(iob,iob_allob)
      if(iob_allob>=iobsta(is).and.iob_allob<=iobend(is))then
        rbox=0.d0
!$OMP parallel do reduction(+:rbox) private(iz,iy,ix)
        do iz=mg%is(3),mg%ie(3)
        do iy=mg%is(2),mg%ie(2)
        do ix=mg%is(1),mg%ie(1)
          rbox=rbox+abs(psi(ix,iy,iz,iob,1))**2
        end do
        end do
        end do
        call comm_summation(rbox,rbox1,nproc_group_korbital)
!$OMP parallel do private(iz,iy,ix)
        do iz=mg%is(3),mg%ie(3)
        do iy=mg%is(2),mg%ie(2)
        do ix=mg%is(1),mg%ie(1)
          psi(ix,iy,iz,iob,1)=psi(ix,iy,iz,iob,1)/sqrt(rbox1*Hvol)
        end do
        end do
        end do
      end if
    end do
    deallocate(Amat,Amat2,Smat)
    deallocate(evec)

  end if

end do

deallocate(htpsi,psi_box)

end subroutine subspace_diag
