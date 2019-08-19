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
module subspace_diag_periodic_sub
  implicit none

contains

subroutine subspace_diag_periodic(mg,system,info,stencil,srg_ob_1,spsi,ilsda,nproc_ob,  &
                                  iobnum,itotmst,mst,ifmst,   &
                                  info_ob,ppg,vlocal)

  use structures, only: s_rgrid,s_dft_system,s_orbital_parallel,s_orbital,s_stencil,s_scalar,s_pp_grid
  use salmon_parallel, only: nproc_group_k
  use salmon_communication, only: comm_bcast, comm_summation
  use timer
  use hpsi_sub
  use eigen_subdiag_periodic_sub
  use calc_allob_sub
  use calc_iroot_sub
  use calc_myob_sub
  use check_corrkob_sub
  use set_isstaend_sub
  use sendrecv_grid, only: s_sendrecv_grid
  implicit none
  type(s_rgrid),intent(in) :: mg
  type(s_dft_system),intent(in) :: system
  type(s_orbital_parallel)     :: info
  type(s_orbital),intent(inout) :: spsi
  type(s_stencil) :: stencil
  type(s_sendrecv_grid),intent(inout) :: srg_ob_1
  type(s_pp_grid) :: ppg
  integer,intent(in)  :: ilsda
  integer,intent(in)  :: nproc_ob
  integer,intent(in)  :: iobnum
  integer,intent(in)  :: itotmst
  integer,intent(in)  :: mst(2),ifmst(2)
  type(s_orbital_parallel)       :: info_ob
  type(s_scalar),intent(in) :: vlocal(system%nspin)
  integer,parameter :: nd=4
  integer :: j
  integer :: ii,jj,ik
  integer :: ix,iy,iz
  complex(8),allocatable :: amat(:,:)
  complex(8),allocatable :: amat2(:,:)
  complex(8),allocatable :: smat(:,:)
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
  integer :: nspin_1
  type(s_orbital)  :: stpsi
  type(s_orbital)  :: shtpsi
  type(s_scalar),allocatable :: v(:)
  complex(8),parameter :: zi=(0.d0,1.d0)
  type(s_dft_system) :: system_spin1 ! temporary
  
  call timer_begin(LOG_DIAG_TOTAL)
  
  call timer_begin(LOG_DIAG_INIT)
  allocate(stpsi%zwf(mg%is_array(1):mg%ie_array(1),  &
                     mg%is_array(2):mg%ie_array(2),  &
                     mg%is_array(3):mg%ie_array(3),1,1,1,1))
  allocate(shtpsi%zwf(mg%is_array(1):mg%ie_array(1),  &
                      mg%is_array(2):mg%ie_array(2),  &
                      mg%is_array(3):mg%ie_array(3),1,1,1,1))

  allocate(stencil%vec_kAc(3,1:1))

  nspin_1=1
  system_spin1%nspin = 1
  allocate(v(nspin_1))
  allocate(v(nspin_1)%f(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))
  
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
    stpsi%zwf(ix,iy,iz,1,1,1,1)=0.d0
  end do
  end do
  end do
  
  call set_isstaend(is_sta,is_end,ilsda)
  call timer_end(LOG_DIAG_INIT)

  
  do ik=info%ik_s,info%ik_e
  do is=is_sta,is_end

    call timer_begin(LOG_DIAG_VLOCAL)
    do j=1,3
      stencil%vec_kAc(j,1) = system%vec_k(j,ik)
    end do
    call update_kvector_nonlocalpt(ppg,stencil%vec_kAc,1,1)

    if(is==1)then
  !$OMP parallel do private(iz,iy,ix) collapse(2)
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        v(1)%f(ix,iy,iz) = vlocal(1)%f(ix,iy,iz)
      end do
      end do
      end do
    else
  !$OMP parallel do private(iz,iy,ix) collapse(2)
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        v(1)%f(ix,iy,iz) = vlocal(2)%f(ix,iy,iz)
      end do
      end do
      end do
    end if
    call timer_end(LOG_DIAG_VLOCAL)


    call timer_begin(LOG_DIAG_AMAT)
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
        call calc_allob(jj,j_allob,itotmst,mst,iobnum)
        if(j_allob>=iobsta(is).and.j_allob<=iobend(is))then
  !$OMP parallel do private(iz,iy,ix)
          do iz=mg%is(3),mg%ie(3)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
            stpsi%zwf(ix,iy,iz,1,1,1,1)=spsi%zwf(ix,iy,iz,is,jj-(is-1)*info%numo,ik,1)
          end do
          end do
          end do
          call hpsi(stpsi,shtpsi,info_ob,mg,v,system_spin1,stencil,srg_ob_1,ppg)
  !$OMP parallel do private(iz,iy,ix)
          do iz=mg%is(3),mg%ie(3)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
            htpsi_groupob(ix,iy,iz,jj)=shtpsi%zwf(ix,iy,iz,1,1,1,1)
          end do
          end do
          end do
        end if
      end do
    
      do jj=iobsta(is),iobend(is)
        call calc_myob(jj,j_myob,ilsda,nproc_ob,itotmst,mst)
        call check_corrkob(jj,ik,icorr_j,ilsda,nproc_ob,info%ik_s,info%ik_e,mst)
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
        call calc_iroot(jj,iroot,ilsda,nproc_ob,itotmst,mst)
        call comm_bcast(htpsi,info%icomm_o,iroot)
        do ii=1,iobnum
          call calc_allob(ii,i_allob,itotmst,mst,iobnum)
          if(i_allob>=iobsta(is).and.i_allob<=iobend(is))then
            cbox=0.d0
  !$OMP parallel do private(iz,iy,ix) reduction (+ : cbox)
            do iz=mg%is(3),mg%ie(3)
            do iy=mg%is(2),mg%ie(2)
            do ix=mg%is(1),mg%ie(1)
               cbox=cbox+conjg(spsi%zwf(ix,iy,iz,is,ii-(is-1)*info%numo,ik,1))*htpsi(ix,iy,iz)
            end do
            end do
            end do
            amat(i_allob-iobsta(is)+1,jj-iobsta(is)+1)=cbox*system%hvol
          end if
        end do
      end do
      call timer_end(LOG_DIAG_AMAT)
    

      call timer_begin(LOG_DIAG_ALLREDUCE)
      call comm_summation(amat,amat2,iter*iter,nproc_group_k)
      call timer_end(LOG_DIAG_ALLREDUCE)
    

      call timer_begin(LOG_DIAG_EIGEN)
      call eigen_subdiag_periodic(amat2,evec,iter,ier2)
      call timer_end(LOG_DIAG_EIGEN)
    

      call timer_begin(LOG_DIAG_SET_ORBITAL)
      do jj=1,iobnum
        call calc_allob(jj,j_allob,itotmst,mst,iobnum)
        if(j_allob>=iobsta(is).and.j_allob<=iobend(is))then
  !$OMP parallel do private(iz,iy,ix)
          do iz=mg%is(3),mg%ie(3)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
            ztpsi_groupob(ix,iy,iz,jj)=spsi%zwf(ix,iy,iz,is,jj-(is-1)*info%numo,ik,1)
            spsi%zwf(ix,iy,iz,is,jj-(is-1)*info%numo,ik,1)=0.d0
          end do
          end do
          end do
        end if
      end do
      call timer_end(LOG_DIAG_SET_ORBITAL)


      call timer_begin(LOG_DIAG_UPDATE)
      do jj=iobsta(is),iobend(is)
        call calc_myob(jj,j_myob,ilsda,nproc_ob,itotmst,mst)
        call check_corrkob(jj,ik,icorr_j,ilsda,nproc_ob,info%ik_s,info%ik_e,mst)
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
        call calc_iroot(jj,iroot,ilsda,nproc_ob,itotmst,mst)
        call comm_bcast(htpsi,info%icomm_o,iroot)
        do ii=1,iobnum
          call calc_allob(ii,i_allob,itotmst,mst,iobnum)
          if(i_allob>=iobsta(is).and.i_allob<=iobend(is))then
  !$OMP parallel do private(iz,iy,ix)
            do iz=mg%is(3),mg%ie(3)
            do iy=mg%is(2),mg%ie(2)
            do ix=mg%is(1),mg%ie(1)
              spsi%zwf(ix,iy,iz,is,ii-(is-1)*info%numo,ik,1)=   &
                spsi%zwf(ix,iy,iz,is,ii-(is-1)*info%numo,ik,1)+evec(jj-iobsta(is)+1,i_allob-iobsta(is)+1)*htpsi(ix,iy,iz)
            end do
            end do
            end do
          end if
        end do
      end do
    
      do ii=1,iobnum
        call calc_allob(ii,i_allob,itotmst,mst,iobnum)
        if(i_allob>=iobsta(is).and.i_allob<=iobend(is))then
          rbox=0.d0
  !$OMP parallel do private(iz,iy,ix) reduction(+:rbox)
          do iz=mg%is(3),mg%ie(3)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
            rbox=rbox+abs(spsi%zwf(ix,iy,iz,is,ii-(is-1)*info%numo,ik,1))**2
          end do
          end do
          end do
          call comm_summation(rbox,rbox1,info%icomm_r)
  !$OMP parallel do private(iz,iy,ix) 
          do iz=mg%is(3),mg%ie(3)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
            spsi%zwf(ix,iy,iz,is,ii-(is-1)*info%numo,ik,1)=   &
              spsi%zwf(ix,iy,iz,is,ii-(is-1)*info%numo,ik,1)/sqrt(rbox1*system%hvol)
          end do
          end do
          end do
        end if
      end do
      call timer_end(LOG_DIAG_UPDATE)
    
      deallocate(amat,amat2,smat)
      deallocate(evec)
  
    end if
  
  end do
  end do
  
  deallocate(htpsi,htpsi_groupob,ztpsi_groupob)

  deallocate(stpsi%zwf,shtpsi%zwf)
  deallocate(stencil%vec_kAc)
  deallocate(v(nspin_1)%f)
  deallocate(v)
  if(allocated(ppg%zekr_uV)) deallocate(ppg%zekr_uV)

  call timer_end(LOG_DIAG_TOTAL)

end subroutine subspace_diag_periodic

end module subspace_diag_periodic_sub
