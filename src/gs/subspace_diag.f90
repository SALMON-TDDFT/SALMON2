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
subroutine subspace_diag(mg,spsi,elp3,ilsda,nproc_ob,iparaway_ob,iobnum,itotmst,k_sta,k_end,nproc_ob_spin,mst,ifmst,hvol,  &
                info_ob,bnmat,cnmat,hgs,ppg,vlocal)

  use inputoutput, only: ispin
  use structures, only: s_rgrid,s_wf_info,s_wavefunction,s_stencil,s_scalar,s_pp_grid
  use salmon_parallel, only: nproc_group_kgrid, nproc_group_global, nproc_group_korbital
  use salmon_communication, only: comm_summation, comm_bcast
  use misc_routines, only: get_wtime
  use hpsi_sub
  implicit none
  type(s_rgrid),intent(inout) :: mg
  type(s_wavefunction),intent(inout) :: spsi
  type(s_stencil) :: stencil
  type(s_pp_grid) :: ppg
  real(8),intent(out) :: elp3(3000)
  integer,intent(in)  :: ilsda
  integer,intent(in)  :: nproc_ob
  integer,intent(in)  :: iparaway_ob
  integer,intent(in)  :: iobnum
  integer,intent(in)  :: itotmst
  integer,intent(in)  :: mst(2),ifmst(2)
  integer,intent(in)  :: k_sta,k_end
  integer,intent(in)  :: nproc_ob_spin(2)
  real(8),intent(in)  :: hvol
  type(s_wf_info)       :: info_ob
  real(8),intent(in)    :: cnmat(0:12,0:12),bnmat(0:12,0:12)
  real(8),intent(in)    :: hgs(3)
  real(8),intent(in)    :: vlocal(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),ispin+1)
  integer,parameter :: nd=4
  integer :: j,ind
  integer :: iob,job,ii,jj
  integer :: ix,iy,iz,is
  integer :: nspin
  type(s_wavefunction)  :: stpsi
  type(s_wavefunction)  :: shtpsi
  type(s_scalar),allocatable :: v(:)
  real(8),allocatable :: amat(:,:)
  real(8),allocatable :: amat2(:,:)
  real(8),allocatable :: smat(:,:)
  real(8),allocatable :: htpsi(:,:,:)
  real(8),allocatable :: psi_box(:,:,:,:)
  real(8) :: rbox,rbox1
  real(8),allocatable :: evec(:,:)
  real(8),allocatable :: rmatbox_m(:,:,:)
  integer :: ier2
  integer :: is_sta,is_end
  integer :: job_myob,iroot,icorr_j,iob_allob,job_allob
  integer :: iter
  integer :: iobsta(2),iobend(2)
  
  elp3(301)=get_wtime()

  allocate(stpsi%rwf(mg%is_array(1):mg%ie_array(1),  &
                     mg%is_array(2):mg%ie_array(2),  &
                     mg%is_array(3):mg%ie_array(3),1,1,1,1))
  allocate(shtpsi%rwf(mg%is_array(1):mg%ie_array(1),  &
                      mg%is_array(2):mg%ie_array(2),  &
                      mg%is_array(3):mg%ie_array(3),1,1,1,1))

  stencil%lap0 = 0.5d0*cNmat(0,nd)*(1.d0/hgs(1)**2+1.d0/hgs(2)**2+1.d0/hgs(3)**2)
  do j=1,3
    do ind=1,4
      stencil%lapt(ind,j) = cnmat(ind,4)/hgs(j)**2
      stencil%nabt(ind,j) = bnmat(ind,4)/hgs(j)
    end do
  end do

  mg%is_overlap = mg%is - 4
  mg%ie_overlap = mg%ie + 4

  allocate(mg%idx(mg%is_overlap(1):mg%ie_overlap(1)) &
          ,mg%idy(mg%is_overlap(2):mg%ie_overlap(2)) &
          ,mg%idz(mg%is_overlap(3):mg%ie_overlap(3)))
  do j=mg%is_overlap(1),mg%ie_overlap(1)
    mg%idx(j) = j
  end do
  do j=mg%is_overlap(2),mg%ie_overlap(2)
    mg%idy(j) = j
  end do
  do j=mg%is_overlap(3),mg%ie_overlap(3)
    mg%idz(j) = j
  end do

  nspin=1
  allocate(v(1))
  allocate(v(1)%f(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))
 
  allocate(htpsi(mg%is(1):mg%ie(1),  &
                 mg%is(2):mg%ie(2),  &
                 mg%is(3):mg%ie(3)))
  allocate(psi_box(mg%is(1):mg%ie(1),  &
                   mg%is(2):mg%ie(2),  &
                   mg%is(3):mg%ie(3),1:iobnum))
  allocate(rmatbox_m(mg%is(1):mg%ie(1),  &
                     mg%is(2):mg%ie(2),  &
                     mg%is(3):mg%ie(3)))
  
  call set_isstaend(is_sta,is_end,ilsda,nproc_ob,nproc_ob_spin)
  
  if(ilsda == 0)then
    iobsta(1)=1
    iobend(1)=itotmst
  else if(ilsda == 1)then
    iobsta(1)=1
    iobend(1)=mst(1)
    iobsta(2)=mst(1)+1
    iobend(2)=itotmst
  end if
  
  !$OMP parallel do private(iz,iy,ix)
  do iz=mg%is_array(3),mg%ie_array(3)
  do iy=mg%is_array(2),mg%ie_array(2)
  do ix=mg%is_array(1),mg%ie_array(1)
    stpsi%rwf(ix,iy,iz,1,1,1,1)=0.d0
  end do
  end do
  end do
  
  do is=is_sta,is_end
  
    if(is==1)then
  !$OMP parallel do private(iz,iy,ix) collapse(2)
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        v(1)%f(ix,iy,iz) = vlocal(ix,iy,iz,1)
      end do
      end do
      end do
    else
  !$OMP parallel do private(iz,iy,ix) collapse(2)
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        v(1)%f(ix,iy,iz) = vlocal(ix,iy,iz,2)
      end do
      end do
      end do
    end if

    if(ifmst(is)>=1.and.mst(is)>=1)then
  
      iter=iobend(is)-iobsta(is)+1
      allocate(evec(iter,iter))
      allocate(amat(iter,iter))
      allocate(amat2(iter,iter))
      allocate(smat(iter,iter))
    
      do jj=1,iter
  !$OMP parallel do 
        do ii=1,iter
          amat2(ii,jj)=0.d0
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
            stpsi%rwf(ix,iy,iz,1,1,1,1)=spsi%rwf(ix,iy,iz,1,job_myob,1,1)
          end do
          end do
          end do
          call hpsi(stpsi,shtpsi,info_ob,mg,v,nspin,stencil,ppg)
  !$OMP parallel do private(iz,iy,ix)
          do iz=mg%is(3),mg%ie(3)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
            htpsi(ix,iy,iz)=shtpsi%rwf(ix,iy,iz,1,1,1,1)
          end do
          end do
          end do
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
              rbox=rbox+spsi%rwf(ix,iy,iz,1,iob,1,1)*htpsi(ix,iy,iz)
            end do
            end do
            end do
            amat2(iob_allob-iobsta(is)+1,job-iobsta(is)+1)=rbox*Hvol
          end if
        end do
      end do
      
      call comm_summation(amat2,amat,iter*iter,nproc_group_global)
    
      call eigen_subdiag(amat,evec,iter,ier2)
     
      do job=1,iobnum
        call calc_allob(job,job_allob)
        if(job_allob>=iobsta(is).and.job_allob<=iobend(is))then
  !$OMP parallel do private(iz,iy,ix)
          do iz=mg%is(3),mg%ie(3)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
            psi_box(ix,iy,iz,job)=spsi%rwf(ix,iy,iz,1,job,1,1)
            spsi%rwf(ix,iy,iz,1,job,1,1)=0.d0
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
              spsi%rwf(ix,iy,iz,1,iob,1,1)=spsi%rwf(ix,iy,iz,1,iob,1,1)  &
                                             +evec(job-iobsta(is)+1,iob_allob-iobsta(is)+1)*rmatbox_m(ix,iy,iz)
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
            rbox=rbox+abs(spsi%rwf(ix,iy,iz,1,iob,1,1))**2
          end do
          end do
          end do
          call comm_summation(rbox,rbox1,nproc_group_korbital)
  !$OMP parallel do private(iz,iy,ix)
          do iz=mg%is(3),mg%ie(3)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
            spsi%rwf(ix,iy,iz,1,iob,1,1)=spsi%rwf(ix,iy,iz,1,iob,1,1)/sqrt(rbox1*Hvol)
          end do
          end do
          end do
        end if
      end do
      deallocate(amat,amat2,smat)
      deallocate(evec)
  
    end if
  
  end do
  
  deallocate(htpsi,psi_box)
  deallocate(stpsi%rwf,shtpsi%rwf)
  deallocate(mg%idx,mg%idy,mg%idz)
  deallocate(v(1)%f)
  deallocate(v)

end subroutine subspace_diag
