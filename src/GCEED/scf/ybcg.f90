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
!=======================================================================
!======================================= Conjugate-Gradient minimization

subroutine dtcg(mg,info,spsi,iflag,itotmst,mst,hvol,ilsda,nproc_ob,nproc_ob_spin,iparaway_ob)
  use inputoutput, only: ncg
  use structures, only: s_rgrid,s_wf_info,s_wavefunction
  use salmon_parallel, only: nproc_group_grid
  use salmon_communication, only: comm_bcast
  use misc_routines, only: get_wtime
  use inner_product_sub
  !$ use omp_lib
  implicit none
  
  type(s_rgrid),intent(in) :: mg
  type(s_wf_info) :: info
  type(s_wavefunction),intent(inout) :: spsi
  integer,intent(inout) :: iflag
  integer,intent(in)    :: itotmst
  integer,intent(in)    :: mst(2)
  real(8),intent(in)    :: hvol
  integer,intent(in)    :: ilsda
  integer,intent(in)    :: nproc_ob
  integer,intent(in)    :: nproc_ob_spin
  integer,intent(in)    :: iparaway_ob
  integer :: iter,iob,job
  integer :: ix,iy,iz
  integer :: is,iobsta(2),iobend(2)
  real(8) :: sum0,xkhxk,xkxk,Rk,gkgk,xkpk,pkpk,pkhxk,pkhpk
  real(8) :: uk,alpha,ak,bk,ck
  real(8) , allocatable :: xk(:,:,:),hxk(:,:,:),gk(:,:,:),pk(:,:,:)
  real(8) , allocatable :: gk2(:,:,:)
  real(8) :: elp2(2000)
  real(8):: tpsi(mg%is_array(1):mg%ie_array(1),  &
                 mg%is_array(2):mg%ie_array(2),  &
                 mg%is_array(3):mg%ie_array(3))
  integer :: iob_myob,job_myob
  integer :: icorr,jcorr               
  integer :: iroot
  integer :: is_sta,is_end
  character(30) :: commname
  
  commname='nproc_group_korbital'
  
  allocate (xk(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))
  allocate (hxk(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))
  allocate (gk(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))
  allocate (gk2(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))
  allocate (pk(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))
  
  call set_isstaend(is_sta,is_end,ilsda,nproc_ob,nproc_ob_spin)
  
  !$OMP parallel do private(iz,iy,ix) 
  do iz=mg%is_array(3),mg%ie_array(3)
  do iy=mg%is_array(2),mg%ie_array(2)
  do ix=mg%is_array(1),mg%ie_array(1)
    tpsi(ix,iy,iz)=0.d0
  end do
  end do
  end do
  
  elp2(:)=0d0
  elp2(1)=get_wtime()
  
  if(ilsda == 0)then
    iobsta(1)=1
    iobend(1)=itotMST
  else if(ilsda == 1)then
    iobsta(1)=1
    iobend(1)=MST(1)
    iobsta(2)=MST(1)+1
    iobend(2)=itotMST
  end if
  
  do is=is_sta,is_end
  
  
  orbital : do iob=iobsta(is),iobend(is)
    call calc_myob(iob,iob_myob,ilsda,nproc_ob,iparaway_ob,itotmst,nproc_ob_spin,mst)
    call check_corrkob(iob,1,icorr,ilsda,nproc_ob,iparaway_ob,itotmst,info%ik_s,info%ik_e,nproc_ob_spin,mst)
    elp2(2)=get_wtime()
  
    if(icorr==1)then
  
  !$OMP parallel do private(iz,iy,ix) 
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        xk(ix,iy,iz)=spsi%rwf(ix,iy,iz,1,iob_myob,1,1)
      end do
      end do
      end do
  
  !$OMP parallel do private(iz,iy,ix) 
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        tpsi(ix,iy,iz)=xk(ix,iy,iz)
      end do
      end do
      end do
  
      call r_hpsi2_buf(tpsi,hxk,iob,1,0,0)
  
      call inner_product(mg,xk,hxk,xkhxk,commname)
  
      xkhxk=xkhxk*hvol ; xkxk=1.d0 ; Rk=xkhxk/xkxk
  
    end if
  
    Iteration : do iter=1,ncg
  
      if(icorr==1)then
  !$OMP parallel do private(iz,iy,ix) 
        do iz=mg%is(3),mg%ie(3)
        do iy=mg%is(2),mg%ie(2)
        do ix=mg%is(1),mg%ie(1)
          gk(ix,iy,iz) = 2*( hxk(ix,iy,iz) - Rk*xk(ix,iy,iz) )
        end do
        end do
        end do
      end if
      call calc_iroot(iob,iroot,ilsda,nproc_ob,iparaway_ob,itotmst,nproc_ob_spin,mst)
      call comm_bcast(gk,nproc_group_grid,iroot)
  
      do job=iobsta(is),iob-1
        sum0=0.d0
        call calc_myob(job,job_myob,ilsda,nproc_ob,iparaway_ob,itotmst,nproc_ob_spin,mst)
        call check_corrkob(job,1,jcorr,ilsda,nproc_ob,iparaway_ob,itotmst,info%ik_s,info%ik_e,nproc_ob_spin,mst)
        if(jcorr==1)then
          call inner_product(mg,spsi%rwf(:,:,:,1,job_myob,1,1),gk(:,:,:),sum0,commname)
          sum0=sum0*hvol
  !$OMP parallel do private(iz,iy,ix)
          do iz=mg%is(3),mg%ie(3)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
            gk(ix,iy,iz)=gk(ix,iy,iz)-sum0*spsi%rwf(ix,iy,iz,1,job_myob,1,1)
          end do
          end do
          end do
        end if
  
        call calc_iroot(job,iroot,ilsda,nproc_ob,iparaway_ob,itotmst,nproc_ob_spin,mst)
        call comm_bcast(gk,nproc_group_grid,iroot)
      end do
  
      if(icorr==1)then
  
        call inner_product(mg,gk,gk,sum0,commname)
        sum0=sum0*hvol
  
        if ( iter==1 ) then
  !$OMP parallel do private(iz,iy,ix)
          do iz=mg%is(3),mg%ie(3)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
            pk(ix,iy,iz) = -gk(ix,iy,iz)
          end do
          end do
          end do
        else
          uk=sum0/gkgk 
  !$OMP parallel do private(iz,iy,ix)
          do iz=mg%is(3),mg%ie(3)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
            pk(ix,iy,iz) = -gk(ix,iy,iz) + uk*pk(ix,iy,iz)
          end do
          end do
          end do
        end if
        gkgk=sum0
  
        xkpk=0.d0 ; pkpk=0.d0 ; pkhxk=0.d0
  
        call inner_product(mg,xk,pk,xkpk,commname)
        xkpk = xkpk*hvol
  
        call inner_product(mg,pk,pk,pkpk,commname)
        pkpk = pkpk*hvol
  
        call inner_product(mg,pk,hxk,pkhxk,commname)
        pkhxk = pkhxk*hvol
  
  !$OMP parallel do private(iz,iy,ix)
        do iz=mg%is(3),mg%ie(3)
        do iy=mg%is(2),mg%ie(2)
        do ix=mg%is(1),mg%ie(1)
          tpsi(ix,iy,iz) = pk(ix,iy,iz)
        end do
        end do
        end do
        call r_hpsi2_buf(tpsi,gk,iob,1,0,0)
  
        call inner_product(mg,pk,gk,pkhpk,commname)
        pkhpk = pkhpk*hvol
  
        ak=pkhpk*xkpk-pkhxk*pkpk
        bk=pkhpk*xkxk-xkhxk*pkpk
        ck=pkhxk*xkxk-xkhxk*xkpk
        alpha=(-bk+sqrt(bk*bk-4*ak*ck))/(2*ak)
  
        xk = xk + alpha*pk
        hxk=hxk + alpha*gk
  
        call inner_product(mg,xk,hxk,xkhxk,commname)
        xkhxk = xkhxk*hvol
  
        call inner_product(mg,xk,xk,xkxk,commname)
        xkxk = xkxk*hvol
      
        Rk=xkhxk/xkxk
  
      end if
  
    end do Iteration
  
    if(icorr==1)then
      call inner_product(mg,xk,xk,sum0,commname)
  !$OMP parallel do private(iz,iy,ix)
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        spsi%rwf(ix,iy,iz,1,iob_myob,1,1)=xk(ix,iy,iz)/sqrt(sum0*hvol)
      end do
      end do
      end do
    end if
  
  end do orbital
  
  end do
  
  if(iflag.eq.1) then
    iflag=0
  end if
  
  deallocate (xk,hxk,gk,pk,gk2)
  return

end subroutine dtcg
