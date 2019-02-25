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
module dtcg_sub
  implicit none

contains

!=======================================================================
!======================================= Conjugate-Gradient minimization

subroutine dtcg(mg,nspin_2,info,stencil,spsi,iflag,itotmst,mst,hvol,ilsda,nproc_ob,iparaway_ob,  &
                info_ob,bnmat,cnmat,hgs,ppg,vlocal)
  use inputoutput, only: ncg,ispin
  use structures, only: s_rgrid,s_wf_info,s_wavefunction,s_stencil,s_scalar,s_pp_grid
  use salmon_parallel, only: nproc_group_grid,nproc_group_korbital
  use salmon_communication, only: comm_bcast,comm_summation
  use misc_routines, only: get_wtime
  use inner_product_sub
  use hpsi_sub
  use calc_iroot_sub
  use calc_myob_sub
  use check_corrkob_sub
  use set_isstaend_sub
  !$ use omp_lib
  implicit none
  
  type(s_rgrid),intent(in) :: mg
  integer,intent(in)    :: nspin_2
  type(s_wf_info) :: info
  type(s_wavefunction),intent(inout) :: spsi
  type(s_stencil) :: stencil
  type(s_pp_grid) :: ppg
  integer,intent(inout) :: iflag
  integer,intent(in)    :: itotmst
  integer,intent(in)    :: mst(2)
  real(8),intent(in)    :: hvol
  integer,intent(in)    :: ilsda
  integer,intent(in)    :: nproc_ob
  integer,intent(in)    :: iparaway_ob
  type(s_wf_info)       :: info_ob
  real(8),intent(in)    :: cnmat(0:12,0:12),bnmat(0:12,0:12)
  real(8),intent(in)    :: hgs(3)
  real(8),intent(in)    :: vlocal(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),ispin+1)
  integer,parameter :: nd=4
  integer :: j,ind
  integer :: iter,iob,job
  integer :: ix,iy,iz
  integer :: is,iobsta(2),iobend(2)
  integer :: nspin
  type(s_wavefunction)  :: stpsi
  type(s_wavefunction)  :: shtpsi
  type(s_scalar),allocatable :: v(:)
  real(8) :: sum0,xkhxk,xkxk,Rk,gkgk,xkpk,pkpk,pkhxk,pkhpk,sum1
  real(8) :: uk,alpha,ak,bk,ck
  real(8) , allocatable :: xk(:,:,:),hxk(:,:,:),gk(:,:,:),pk(:,:,:)
  real(8) , allocatable :: gk2(:,:,:)
  real(8) :: elp2(2000)
  integer :: iob_myob,job_myob
  integer :: icorr,jcorr               
  integer :: iroot
  integer :: is_sta,is_end
  character(30) :: commname

  commname='nproc_group_korbital'

  allocate(stpsi%rwf(mg%is_array(1):mg%ie_array(1),  &
                     mg%is_array(2):mg%ie_array(2),  &
                     mg%is_array(3):mg%ie_array(3),1,1,1,1))
  allocate(shtpsi%rwf(mg%is_array(1):mg%ie_array(1),  &
                      mg%is_array(2):mg%ie_array(2),  &
                      mg%is_array(3):mg%ie_array(3),1,1,1,1))

  nspin=1
  allocate(v(1))
  allocate(v(1)%f(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))

  allocate (xk(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))
  allocate (hxk(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))
  allocate (gk(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))
  allocate (gk2(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))
  allocate (pk(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))
  
  call set_isstaend(is_sta,is_end,ilsda)
  
  !$OMP parallel do private(iz,iy,ix) 
  do iz=mg%is_array(3),mg%ie_array(3)
  do iy=mg%is_array(2),mg%ie_array(2)
  do ix=mg%is_array(1),mg%ie_array(1)
    stpsi%rwf(ix,iy,iz,1,1,1,1)=0.d0
    shtpsi%rwf(ix,iy,iz,1,1,1,1)=0.d0
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
  
  !$OMP parallel do private(iz,iy,ix) 
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      v(1)%f(ix,iy,iz) = vlocal(ix,iy,iz,is)
    end do
    end do
    end do

  orbital : do iob=iobsta(is),iobend(is)
    call calc_myob(iob,iob_myob,ilsda,nproc_ob,iparaway_ob,itotmst,mst,nspin_2*info%numo)
    call check_corrkob(iob,1,icorr,ilsda,nproc_ob,iparaway_ob,info%ik_s,info%ik_e,mst)
    elp2(2)=get_wtime()
  
    if(icorr==1)then
  
  !$OMP parallel do private(iz,iy,ix) 
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        xk(ix,iy,iz)=spsi%rwf(ix,iy,iz,is,iob_myob-(is-1)*info%numo,1,1)
      end do
      end do
      end do
  
  !$OMP parallel do private(iz,iy,ix) 
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        stpsi%rwf(ix,iy,iz,1,1,1,1)=xk(ix,iy,iz)
      end do
      end do
      end do
 
      call hpsi(stpsi,shtpsi,info_ob,mg,v,nspin,stencil,ppg)
 
  !$OMP parallel do private(iz,iy,ix) 
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        hxk(ix,iy,iz)=shtpsi%rwf(ix,iy,iz,1,1,1,1)
      end do
      end do
      end do
  
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
      call calc_iroot(iob,iroot,ilsda,nproc_ob,iparaway_ob,itotmst,mst)
      call comm_bcast(gk,nproc_group_grid,iroot)
  
      do job=iobsta(is),iob-1
        sum0=0.d0
        call calc_myob(job,job_myob,ilsda,nproc_ob,iparaway_ob,itotmst,mst,nspin_2*info%numo)
        call check_corrkob(job,1,jcorr,ilsda,nproc_ob,iparaway_ob,info%ik_s,info%ik_e,mst)
        if(jcorr==1)then
  !$OMP parallel do reduction(+ : sum0) private(iz,iy,ix) 
          do iz=mg%is(3),mg%ie(3)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
            sum0=sum0+spsi%rwf(ix,iy,iz,is,job_myob-(is-1)*info%numo,1,1)*gk(ix,iy,iz)
          end do
          end do
          end do
          call comm_summation(sum0,sum1,nproc_group_korbital)
          sum0=sum1*hvol
  !$OMP parallel do private(iz,iy,ix)
          do iz=mg%is(3),mg%ie(3)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
            gk(ix,iy,iz)=gk(ix,iy,iz)-sum0*spsi%rwf(ix,iy,iz,is,job_myob-(is-1)*info%numo,1,1)
          end do
          end do
          end do
        end if
  
        call calc_iroot(job,iroot,ilsda,nproc_ob,iparaway_ob,itotmst,mst)
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
          stpsi%rwf(ix,iy,iz,1,1,1,1)=pk(ix,iy,iz)
        end do
        end do
        end do

        call hpsi(stpsi,shtpsi,info_ob,mg,v,nspin,stencil,ppg)
  
  !$OMP parallel do private(iz,iy,ix) 
        do iz=mg%is(3),mg%ie(3)
        do iy=mg%is(2),mg%ie(2)
        do ix=mg%is(1),mg%ie(1)
          gk(ix,iy,iz)=shtpsi%rwf(ix,iy,iz,1,1,1,1)
        end do
        end do
        end do

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
        spsi%rwf(ix,iy,iz,is,iob_myob-(is-1)*info%numo,1,1)=xk(ix,iy,iz)/sqrt(sum0*hvol)
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
  
  deallocate(v(1)%f)
  deallocate(v)

  return

end subroutine dtcg

end module dtcg_sub
