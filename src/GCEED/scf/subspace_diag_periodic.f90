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
subroutine subspace_diag_periodic(mg)

  use structures, only: s_rgrid
  use salmon_parallel, only: nproc_group_korbital, nproc_group_k, nproc_group_kgrid
  use salmon_communication, only: comm_bcast, comm_summation
  use misc_routines, only: get_wtime
  use scf_data
  use new_world_sub
  use hpsi2_sub
  implicit none
  type(s_rgrid),intent(in) :: mg
  integer :: ii,jj,ik
  integer :: ix,iy,iz
  complex(8),allocatable :: Amat(:,:)
  complex(8),allocatable :: Amat2(:,:)
  complex(8),allocatable :: Smat(:,:)
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
  
  iwk_size=2
  call make_iwksta_iwkend
  
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
  
  elp3(302)=get_wtime()
  elp3(352)=elp3(352)+elp3(302)-elp3(301)
  
  elp3(303)=get_wtime()
  elp3(353)=elp3(353)+elp3(303)-elp3(302)
  
  call set_isstaend(is_sta,is_end,ilsda,nproc_ob,nproc_ob_spin)
  
  do ik=k_sta,k_end
  do is=is_sta,is_end
  
    if(ifMST(is)>=1.and.MST(is)>=1)then
  
      iter=iobend(is)-iobsta(is)+1
    
      allocate(evec(iter,iter))
      allocate(Amat(iter,iter))
      allocate(Amat2(iter,iter))
      allocate(Smat(iter,iter))
    
  !$OMP parallel do private(jj,ii)
      do jj=1,iter
        do ii=1,iter
          Amat(ii,jj)=0.d0
        end do
      end do
    
  !do jj=1,itotMST
      do jj=1,iobnum
        call calc_allob(jj,j_allob)
        if(j_allob>=iobsta(is).and.j_allob<=iobend(is))then
  !$OMP parallel do private(iz,iy,ix)
          do iz=mg%is(3),mg%ie(3)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
            tpsi(ix,iy,iz)=zpsi(ix,iy,iz,jj,ik)
          end do
          end do
          end do
          call hpsi2(tpsi,htpsi,j_allob,ik,0,0)
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
               cbox=cbox+conjg(zpsi(ix,iy,iz,ii,ik))*htpsi(ix,iy,iz)
            end do
            end do
            end do
            Amat(i_allob-iobsta(is)+1,jj-iobsta(is)+1)=cbox*Hvol
          end if
        end do
      end do
    
      call comm_summation(Amat,Amat2,iter*iter,nproc_group_k)
    
      call eigen_subdiag_periodic(Amat2,evec,iter,ier2)
    
      do jj=1,iobnum
        call calc_allob(jj,j_allob)
        if(j_allob>=iobsta(is).and.j_allob<=iobend(is))then
  !$OMP parallel do private(iz,iy,ix)
          do iz=mg%is(3),mg%ie(3)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
            ztpsi_groupob(ix,iy,iz,jj)=zpsi(ix,iy,iz,jj,ik)
            zpsi(ix,iy,iz,jj,ik)=0.d0
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
              zpsi(ix,iy,iz,ii,ik)=zpsi(ix,iy,iz,ii,ik)+evec(jj-iobsta(is)+1,i_allob-iobsta(is)+1)*htpsi(ix,iy,iz)
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
            rbox=rbox+abs(zpsi(ix,iy,iz,ii,ik))**2
          end do
          end do
          end do
          call comm_summation(rbox,rbox1,nproc_group_korbital)
  !$OMP parallel do private(iz,iy,ix) 
          do iz=mg%is(3),mg%ie(3)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
            zpsi(ix,iy,iz,ii,ik)=zpsi(ix,iy,iz,ii,ik)/sqrt(rbox1*Hvol)
          end do
          end do
          end do
        end if
      end do
    
      deallocate(Amat,Amat2,Smat)
      deallocate(evec)
  
    end if
  
  end do
  end do
  
  deallocate(htpsi,htpsi_groupob,ztpsi_groupob)

end subroutine subspace_diag_periodic
