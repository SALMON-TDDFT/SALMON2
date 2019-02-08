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
                                  iobnum,itotmst,k_sta,k_end,nproc_ob_spin,mst,ifmst,hvol,   &
                                  info_ob,bnmat,cnmat,hgs,ppg,vlocal,num_kpoints_rd,k_rd)

  use inputoutput, only: ispin,natom
  use structures, only: s_rgrid,s_wf_info,s_wavefunction,s_stencil,s_scalar,s_pp_grid
  use salmon_parallel, only: nproc_group_korbital, nproc_group_k, nproc_group_kgrid
  use salmon_communication, only: comm_bcast, comm_summation
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
  integer,intent(in)    :: num_kpoints_rd
  real(8),intent(in)    :: k_rd(3,num_kpoints_rd)
  integer,parameter :: nd=4
  integer :: j,ind
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
  integer :: nspin
  type(s_wavefunction)  :: stpsi
  type(s_wavefunction)  :: shtpsi
  type(s_scalar),allocatable :: v(:)
  complex(8) :: ekr(ppg%nps,natom)
  real(8) :: x,y,z
  integer :: a,iatom,ilma
  complex(8),parameter :: zi=(0.d0,1.d0)
  
  elp3(301)=get_wtime()
  
  allocate(stpsi%zwf(mg%is_array(1):mg%ie_array(1),  &
                     mg%is_array(2):mg%ie_array(2),  &
                     mg%is_array(3):mg%ie_array(3),1,1,1,1))
  allocate(shtpsi%zwf(mg%is_array(1):mg%ie_array(1),  &
                      mg%is_array(2):mg%ie_array(2),  &
                      mg%is_array(3):mg%ie_array(3),1,1,1,1))

  allocate(stencil%kAc(1:1,3))

  stencil%lap0 = -0.5d0*cNmat(0,Nd)*(1.d0/Hgs(1)**2+1.d0/Hgs(2)**2+1.d0/Hgs(3)**2)
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
  
  elp3(302)=get_wtime()
  elp3(352)=elp3(352)+elp3(302)-elp3(301)
  
  elp3(303)=get_wtime()
  elp3(353)=elp3(353)+elp3(303)-elp3(302)
  
  call set_isstaend(is_sta,is_end,ilsda,nproc_ob,nproc_ob_spin)
  
  do ik=k_sta,k_end
  do is=is_sta,is_end
  
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
            stpsi%zwf(ix,iy,iz,1,1,1,1)=spsi%zwf(ix,iy,iz,1,jj,ik,1)
          end do
          end do
          end do
          call hpsi(stpsi,shtpsi,info_ob,mg,v,nspin,stencil,ppg)
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

  deallocate(stpsi%zwf,shtpsi%zwf)
  deallocate(stencil%kAc)
  deallocate(mg%idx,mg%idy,mg%idz)
  deallocate(v(1)%f)
  deallocate(v)
  if(allocated(ppg%zproj)) deallocate(ppg%zproj)


end subroutine subspace_diag_periodic
