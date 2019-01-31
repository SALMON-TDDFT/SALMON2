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
subroutine subspace_diag_periodic(mg,spsi,elp3,ilsda,nproc_ob,iparaway_ob,  &
                                     iobnum,itotmst,k_sta,k_end,nproc_ob_spin,mst,ifmst,hvol)

  use structures, only: s_rgrid,s_wavefunction
  use salmon_parallel, only: nproc_group_korbital, nproc_group_k, nproc_group_kgrid
  use salmon_communication, only: comm_bcast, comm_summation
  use misc_routines, only: get_wtime
  implicit none
  type(s_rgrid),intent(in) :: mg
  type(s_wavefunction),intent(inout) :: spsi
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
  integer :: ii,jj,ik
  integer :: ix,iy,iz
  complex(8),allocatable :: amat(:,:)
  complex(8),allocatable :: amat2(:,:)
  complex(8),allocatable :: smat(:,:)
  complex(8):: tpsi(mg%is_array(1):mg%ie_array(1),  &
                    mg%is_array(2):mg%ie_array(2),  &
                    mg%is_array(3):mg%ie_array(3))
  complex(8),allocatable :: htpsi(:,:,:)
  complex(8),allocatable :: htpsi_groupob(:,:,:,:)
  complex(8),allocatable :: ztpsi_groupob(:,:,:,:)
  complex(8) :: cbox
  real(8) :: rbox,rbox1
  complex(8),allocatable :: evec(:,:)
  integer :: iter,ier2
  integer :: iroot
  integer :: j_myob,i_allob,j_allob
  integer :: icorr_j
  integer :: is,is_sta,is_end
  integer :: iobsta(2),iobend(2)
  
  elp3(301)=get_wtime()
  
  allocate(htpsi(mg%is(1):mg%ie(1),  &
                 mg%is(2):mg%ie(2),  &
                 mg%is(3):mg%ie(3)))
  allocate(htpsi_groupob(mg%is(1):mg%ie(1),  &
                         mg%is(2):mg%ie(2),  &
                         mg%is(3):mg%ie(3),1:iobnum))
  allocate(ztpsi_groupob(mg%is(1):mg%ie(1),  &
                         mg%is(2):mg%ie(2),  &
                         mg%is(3):mg%ie(3),1:iobnum))
  
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
    tpsi(ix,iy,iz)=0.d0
  end do
  end do
  end do
  
  elp3(302)=get_wtime()
  elp3(352)=elp3(352)+elp3(302)-elp3(301)
  
  elp3(303)=get_wtime()
  elp3(353)=elp3(353)+elp3(303)-elp3(302)
  
  call set_isstaend(is_sta,is_end,ilsda,nproc_ob,nproc_ob_spin)
  
  do ik=k_sta,k_end
  do is=is_sta,is_end
  
    if(ifmst(is)>=1.and.mst(is)>=1)then
  
      iter=iobend(is)-iobsta(is)+1
    
      allocate(evec(iter,iter))
      allocate(amat(iter,iter))
      allocate(amat2(iter,iter))
      allocate(smat(iter,iter))
    
  !$OMP parallel do private(jj,ii)
      do jj=1,iter
        do ii=1,iter
          amat(ii,jj)=0.d0
        end do
      end do
    
  !do jj=1,itotmst
      do jj=1,iobnum
        call calc_allob(jj,j_allob)
        if(j_allob>=iobsta(is).and.j_allob<=iobend(is))then
  !$OMP parallel do private(iz,iy,ix)
          do iz=mg%is(3),mg%ie(3)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
            tpsi(ix,iy,iz)=spsi%zwf(ix,iy,iz,1,jj,ik,1)
          end do
          end do
          end do
          call c_hpsi2_buf(tpsi,htpsi,j_allob,ik,0,0)
  !$OMP parallel do private(iz,iy,ix)
          do iz=mg%is(3),mg%ie(3)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
            htpsi_groupob(ix,iy,iz,jj)=htpsi(ix,iy,iz)
          end do
          end do
          end do
        end if
      end do
    
      do jj=iobsta(is),iobend(is)
        call calc_myob(jj,j_myob,ilsda,nproc_ob,iparaway_ob,itotmst,nproc_ob_spin,mst)
        call check_corrkob(jj,ik,icorr_j,ilsda,nproc_ob,iparaway_ob,itotmst,k_sta,k_end,nproc_ob_spin,mst)
        if(icorr_j==1)then
  !$OMP parallel do private(iz,iy,ix)
          do iz=mg%is(3),mg%ie(3)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
            htpsi(ix,iy,iz)=htpsi_groupob(ix,iy,iz,j_myob)
          end do
          end do
          end do
        end if
        call calc_iroot(jj,iroot,ilsda,nproc_ob,iparaway_ob,itotmst,nproc_ob_spin,mst)
        call comm_bcast(htpsi,nproc_group_kgrid,iroot)
        do ii=1,iobnum
          call calc_allob(ii,i_allob)
          if(i_allob>=iobsta(is).and.i_allob<=iobend(is))then
            cbox=0.d0
  !$OMP parallel do private(iz,iy,ix) reduction (+ : cbox)
            do iz=mg%is(3),mg%ie(3)
            do iy=mg%is(2),mg%ie(2)
            do ix=mg%is(1),mg%ie(1)
               cbox=cbox+conjg(spsi%zwf(ix,iy,iz,1,ii,ik,1))*htpsi(ix,iy,iz)
            end do
            end do
            end do
            amat(i_allob-iobsta(is)+1,jj-iobsta(is)+1)=cbox*hvol
          end if
        end do
      end do
    
      call comm_summation(amat,amat2,iter*iter,nproc_group_k)
    
      call eigen_subdiag_periodic(amat2,evec,iter,ier2)
    
      do jj=1,iobnum
        call calc_allob(jj,j_allob)
        if(j_allob>=iobsta(is).and.j_allob<=iobend(is))then
  !$OMP parallel do private(iz,iy,ix)
          do iz=mg%is(3),mg%ie(3)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
            ztpsi_groupob(ix,iy,iz,jj)=spsi%zwf(ix,iy,iz,1,jj,ik,1)
            spsi%zwf(ix,iy,iz,1,jj,ik,1)=0.d0
          end do
          end do
          end do
        end if
      end do
      
      do jj=iobsta(is),iobend(is)
        call calc_myob(jj,j_myob,ilsda,nproc_ob,iparaway_ob,itotmst,nproc_ob_spin,mst)
        call check_corrkob(jj,ik,icorr_j,ilsda,nproc_ob,iparaway_ob,itotmst,k_sta,k_end,nproc_ob_spin,mst)
        if(icorr_j==1)then
  !$OMP parallel do private(iz,iy,ix)
          do iz=mg%is(3),mg%ie(3)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
            htpsi(ix,iy,iz)=ztpsi_groupob(ix,iy,iz,j_myob) ! htpsi is making a role of original tpsi
          end do
          end do
          end do
        end if
        call calc_iroot(jj,iroot,ilsda,nproc_ob,iparaway_ob,itotmst,nproc_ob_spin,mst)
        call comm_bcast(htpsi,nproc_group_kgrid,iroot)
        do ii=1,iobnum
          call calc_allob(ii,i_allob)
          if(i_allob>=iobsta(is).and.i_allob<=iobend(is))then
  !$OMP parallel do private(iz,iy,ix)
            do iz=mg%is(3),mg%ie(3)
            do iy=mg%is(2),mg%ie(2)
            do ix=mg%is(1),mg%ie(1)
              spsi%zwf(ix,iy,iz,1,ii,ik,1)=spsi%zwf(ix,iy,iz,1,ii,ik,1)+evec(jj-iobsta(is)+1,i_allob-iobsta(is)+1)*htpsi(ix,iy,iz)
            end do
            end do
            end do
          end if
        end do
      end do
    
      do ii=1,iobnum
        call calc_allob(ii,i_allob)
        if(i_allob>=iobsta(is).and.i_allob<=iobend(is))then
          rbox=0.d0
  !$OMP parallel do private(iz,iy,ix) reduction(+:rbox)
          do iz=mg%is(3),mg%ie(3)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
            rbox=rbox+abs(spsi%zwf(ix,iy,iz,1,ii,ik,1))**2
          end do
          end do
          end do
          call comm_summation(rbox,rbox1,nproc_group_korbital)
  !$OMP parallel do private(iz,iy,ix) 
          do iz=mg%is(3),mg%ie(3)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
            spsi%zwf(ix,iy,iz,1,ii,ik,1)=spsi%zwf(ix,iy,iz,1,ii,ik,1)/sqrt(rbox1*hvol)
          end do
          end do
          end do
        end if
      end do
    
      deallocate(amat,amat2,smat)
      deallocate(evec)
  
    end if
  
  end do
  end do
  
  deallocate(htpsi,htpsi_groupob,ztpsi_groupob)

end subroutine subspace_diag_periodic
