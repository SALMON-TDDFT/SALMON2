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
!=======================================================================
module poisson_isolated
  implicit none

contains

!============================ Hartree potential (Solve Poisson equation)
subroutine poisson_cg(lg,mg,info,system,poisson,trho,tVh,srg_scalar,stencil)
  use inputoutput, only: threshold_cg
  use structures, only: s_rgrid,s_parallel_info,s_dft_system,s_poisson,s_sendrecv_grid,s_stencil
  use communication, only: comm_is_root, comm_summation
  use math_constants, only : pi
  use sendrecv_grid, only: update_overlap_real8
  implicit none
  integer,parameter :: ndh=4
  type(s_rgrid),intent(in) :: lg
  type(s_rgrid),intent(in) :: mg
  type(s_parallel_info),intent(in) :: info
  type(s_dft_system),intent(in) :: system
  type(s_poisson),intent(inout) :: poisson
  real(8) :: trho(mg%is(1):mg%ie(1),    &
                  mg%is(2):mg%ie(2),      &
                  mg%is(3):mg%ie(3))
  real(8) :: tVh(mg%is(1):mg%ie(1),    &
                 mg%is(2):mg%ie(2),      &
                 mg%is(3):mg%ie(3))
  type(s_sendrecv_grid),intent(inout) :: srg_scalar
  type(s_stencil),intent(in) :: stencil
  
  integer,parameter :: maxiter=1000
  integer :: ix,iy,iz,iter
  real(8) :: sum1,sum2,ak,ck
  real(8) :: tottmp
  real(8) :: totbox
  real(8) :: rlap_wk(mg%is_array(1):mg%ie_array(1),    &
                     mg%is_array(2):mg%ie_array(2),      &
                     mg%is_array(3):mg%ie_array(3))
  real(8) :: zk(mg%is_array(1):mg%ie_array(1),   &
                mg%is_array(2):mg%ie_array(2),   &
                mg%is_array(3):mg%ie_array(3))
  real(8) :: pk(mg%is_array(1):mg%ie_array(1),   &
                mg%is_array(2):mg%ie_array(2),   &
                mg%is_array(3):mg%ie_array(3))
  
  call poisson_boundary(lg,mg,info,system,poisson,trho,pk)
  
!------------------------- C-G minimization
  
!$omp parallel do private(iz,iy,ix) collapse(2)
  do iz=mg%is_array(3),mg%ie_array(3)
  do iy=mg%is_array(2),mg%ie_array(2)
  do ix=mg%is_array(1),mg%ie_array(1)
    zk(ix,iy,iz)=0.d0
  end do
  end do
  end do
  
!$omp parallel do private(iz,iy,ix) collapse(2)
  do iz=mg%is(3),mg%ie(3)
  do iy=mg%is(2),mg%ie(2)
  do ix=mg%is(1),mg%ie(1)
    pk(ix,iy,iz)=tVh(ix,iy,iz)
    zk(ix,iy,iz)=-4.d0*pi*trho(ix,iy,iz)
  end do
  end do
  end do
  call update_overlap_real8(srg_scalar, mg, pk)
  call laplacian_poisson(mg,pk,rlap_wk,stencil%coef_lap0,stencil%coef_lap)
  
!$omp parallel do private(iz,iy,ix) collapse(2)
  do iz=mg%is_array(3),mg%ie_array(3)
  do iy=mg%is_array(2),mg%ie_array(2)
  do ix=mg%is_array(1),mg%ie_array(1)
    pk(ix,iy,iz)=0.d0
  end do
  end do
  end do
  
!$omp parallel do private(iz,iy,ix) collapse(2)
  do iz=mg%is(3),mg%ie(3)
  do iy=mg%is(2),mg%ie(2)
  do ix=mg%is(1),mg%ie(1)
    zk(ix,iy,iz)=zk(ix,iy,iz)-rlap_wk(ix,iy,iz)
    pk(ix,iy,iz)=zk(ix,iy,iz)
  end do
  end do
  end do
  
  sum1=0.d0
!$omp parallel do reduction(+ : sum1) private(iz,iy,ix) collapse(2)
  do iz=mg%is(3),mg%ie(3)
  do iy=mg%is(2),mg%ie(2)
  do ix=mg%is(1),mg%ie(1)
    sum1=sum1+zk(ix,iy,iz)**2*system%hvol
  end do
  end do
  end do
  
  if(info%isize_r==1)then
  else
    call comm_summation(sum1,sum2,info%icomm_r)
    sum1=sum2
  end if
  
  iteration : do iter=1,maxiter
  
    call update_overlap_real8(srg_scalar, mg, pk)
    call laplacian_poisson(mg,pk,rlap_wk,stencil%coef_lap0,stencil%coef_lap)
  
    totbox=0d0
!$omp parallel do reduction(+ : totbox) private(iz,iy,ix) collapse(2)
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      totbox=totbox+(zk(ix,iy,iz)*rlap_wk(ix,iy,iz))
    end do
    end do
    end do
  
    if(info%isize_r==1)then
      tottmp=totbox
    else
      call comm_summation(totbox,tottmp,info%icomm_r)
    end if
  
    ak=sum1/tottmp/system%hvol
  
!$omp parallel do private(iz,iy,ix) firstprivate(ak) collapse(2)
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
       tVh(ix,iy,iz)=tVh(ix,iy,iz)+ak*pk(ix,iy,iz)
       zk(ix,iy,iz)=zk(ix,iy,iz)-ak*rlap_wk(ix,iy,iz)
    end do
    end do
    end do
  
    totbox=0d0
!$omp parallel do reduction(+ : totbox) private(iz,iy,ix) collapse(2)
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      totbox=totbox+zk(ix,iy,iz)**2
    end do
    end do
    end do
  
    if(info%isize_r==1)then
      tottmp=totbox
    else
      call comm_summation(totbox,tottmp,info%icomm_r)
    end if
  
    sum2=tottmp*system%hvol
  
    if ( abs(sum2) < threshold_cg*dble(lg%num(1)*lg%num(2)*lg%num(3)) ) exit
  
    ck=sum2/sum1 ; sum1=sum2
  
!$omp parallel do private(iz,iy,ix) firstprivate(ck) collapse(2)
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      pk(ix,iy,iz)=zk(ix,iy,iz)+ck*pk(ix,iy,iz)
    end do
    end do
    end do
     
  end do iteration
  
  poisson%iterVh=iter
  if ( poisson%iterVh>maxiter .and. comm_is_root(info%id_r)) then
     write(*,*) "Warning:Vh iteration is not converged"
     write(*,'("||tVh(i)-tVh(i-1)||**2/(# of grids) = ",e15.8)') &
                                sum2/dble(lg%num(1)*lg%num(2)*lg%num(3))
  end if
  
  return

end subroutine poisson_cg

subroutine laplacian_poisson(mg,pk,rlap_wk,lap0,lapt)
  use structures, only: s_rgrid
  implicit none
  type(s_rgrid),intent(in) :: mg
  real(8),intent(in) :: pk(mg%is_array(1):mg%ie_array(1),  &
                           mg%is_array(2):mg%ie_array(2),  &
                           mg%is_array(3):mg%ie_array(3))
  real(8),intent(out) :: rlap_wk(mg%is_array(1):mg%ie_array(1),  &
                                 mg%is_array(2):mg%ie_array(2),  &
                                 mg%is_array(3):mg%ie_array(3))
  real(8),intent(in)  :: lap0,lapt(4,3)
  integer :: ix,iy,iz

!$omp parallel do private(iz,iy,ix) collapse(2)
  do iz=mg%is(3),mg%ie(3)
  do iy=mg%is(2),mg%ie(2)
  do ix=mg%is(1),mg%ie(1)
    rlap_wk(ix,iy,iz)=-2.d0*lap0*pk(ix,iy,iz)+(  &
                      lapt(1,1)*(pk(ix+1,iy,iz) + pk(ix-1,iy,iz)) &
                     +lapt(2,1)*(pk(ix+2,iy,iz) + pk(ix-2,iy,iz)) &
                     +lapt(3,1)*(pk(ix+3,iy,iz) + pk(ix-3,iy,iz)) &
                     +lapt(4,1)*(pk(ix+4,iy,iz) + pk(ix-4,iy,iz)) &
                     +lapt(1,2)*(pk(ix,iy+1,iz) + pk(ix,iy-1,iz)) &
                     +lapt(2,2)*(pk(ix,iy+2,iz) + pk(ix,iy-2,iz)) &
                     +lapt(3,2)*(pk(ix,iy+3,iz) + pk(ix,iy-3,iz)) &
                     +lapt(4,2)*(pk(ix,iy+4,iz) + pk(ix,iy-4,iz)) &
                     +lapt(1,3)*(pk(ix,iy,iz+1) + pk(ix,iy,iz-1)) &
                     +lapt(2,3)*(pk(ix,iy,iz+2) + pk(ix,iy,iz-2)) &
                     +lapt(3,3)*(pk(ix,iy,iz+3) + pk(ix,iy,iz-3)) &
                     +lapt(4,3)*(pk(ix,iy,iz+4) + pk(ix,iy,iz-4)))
  end do
  end do
  end do

  return 

end subroutine laplacian_poisson

subroutine poisson_boundary(lg,mg,info,system,poisson,trho,wk2)
  use inputoutput, only: natom,rion,lmax_lmp,layout_multipole,natom
  use structures, only: s_rgrid,s_parallel_info,s_dft_system,s_poisson
  use salmon_math, only: ylm
  use communication, only: comm_summation
  
  use omp_lib, only: omp_get_num_threads, omp_get_thread_num, omp_get_max_threads
  use misc_routines, only: ceiling_pow2
  
  implicit none
  integer,parameter :: ndh=4
  type(s_rgrid),intent(in) :: lg
  type(s_rgrid),intent(in) :: mg
  type(s_parallel_info),intent(in) :: info
  type(s_dft_system),intent(in) :: system
  type(s_poisson),intent(inout) :: poisson
  real(8) :: trho(mg%is(1):mg%ie(1),    &
                 mg%is(2):mg%ie(2),      &
                 mg%is(3):mg%ie(3))
  real(8) :: wk2(mg%is_array(1):mg%ie_array(1),    &
                 mg%is_array(2):mg%ie_array(2),      &
                 mg%is_array(3):mg%ie_array(3))
  integer,parameter :: maxiter=1000
  integer :: ii,jj,ix,iy,iz,lm,ll,m,icen  !,kk,pl,cl
  integer :: ixbox,iybox,izbox
  integer :: j,k
  integer :: istart(0:info%isize_r-1),iend(0:info%isize_r-1)
  integer :: icount
  integer,allocatable :: itrho(:)
  integer :: num_center
  real(8) :: ylm2(25)
  integer :: l2(25)
  real(8) :: xx,yy,zz,rr,sum1,xxxx,yyyy,zzzz,rrrr,sumbox1,sumbox2,sumbox3
  real(8) :: rholm2box
  real(8),allocatable :: rholm(:,:),rholm2(:,:) !,rholm3(:,:)
 !integer :: tid
  real(8) :: center_trho2(3)
  real(8),allocatable :: center_trho(:,:)
  real(8),allocatable :: center_trho_nume_deno(:,:)
  real(8),allocatable :: center_trho_nume_deno2(:,:)
  real(8) :: xp2,yp2,zp2,xy,yz,xz
  real(8) :: deno(25)
  real(8) :: rinv
  real(8) :: rbox
  real(8),allocatable :: rion2(:,:)
  real(8) :: hvol
  integer,allocatable :: ig_num(:)
  integer,allocatable :: ig(:,:,:)
  real(8),allocatable :: coordinate(:,:)
 !integer :: lmax_lmp_tmp !iwata
  integer :: comm_xyz(3),nproc_xyz(3),myrank_xyz(3)
  
  comm_xyz(1) = info%icomm_x
  comm_xyz(2) = info%icomm_y
  comm_xyz(3) = info%icomm_z
  nproc_xyz(1) = info%isize_x
  nproc_xyz(2) = info%isize_y
  nproc_xyz(3) = info%isize_z
  myrank_xyz(1) = info%id_x
  myrank_xyz(2) = info%id_y
  myrank_xyz(3) = info%id_z
  
  !------------------------- Boundary condition (multipole expansion)

  hvol=system%hvol
  if(allocated(poisson%ig_num)) then
    if(.not.allocated(ig_num)) then
      if(layout_multipole==2)then
        allocate(ig_num(natom))
      else if(layout_multipole==3)then
        allocate(ig_num(poisson%npole_partial))
      end if
    end if
    ig_num=poisson%ig_num
  end if
  if(.not.allocated(ig))then
    allocate(ig(3,maxval(poisson%ig_num(:)),poisson%npole_partial))
  end if
  ig=poisson%ig

  if(.not.allocated(coordinate))then
    allocate(coordinate(minval(lg%is_overlap(1:3)):maxval(lg%ie_overlap(1:3)),3))
  end if
  do j=1,3
    coordinate(lg%is_overlap(j):lg%ie_overlap(j),j)=lg%coordinate(lg%is_overlap(j):lg%ie_overlap(j),j)
  end do
 
  if(.not.allocated(poisson%wkbound))then
    allocate(poisson%wkbound(lg%num(1)*lg%num(2)*lg%num(3)/minval(lg%num(1:3))*6*ndh))
    allocate(poisson%wkbound2(lg%num(1)*lg%num(2)*lg%num(3)/minval(lg%num(1:3))*6*ndh))
  end if
 
  select case( layout_multipole )
  
  case(1)
  
  num_center=1
  allocate (rholm((lmax_lmp+1)**2,1))
  allocate (rholm2((lmax_lmp+1)**2,1))
  allocate(itrho(1))
  allocate(center_trho(3,1))
  do jj=1,3
    center_trho(jj,1)=0.d0
    center_trho2(jj)=0.d0
  end do
  do lm=1,(lmax_lmp+1)**2
    rholm(lm,1)=0.d0
    rholm2(lm,1)=0.d0
  end do
  itrho(1)=1
  
  do ll=0,lmax_lmp
  do m=-ll,ll
    lm=ll*ll+ll+1+m
    rholm2box=0.d0
  !$OMP parallel do reduction ( + : rholm2box)&
  !$OMP private(ix,iy,iz,xx,yy,zz,rr,xxxx,yyyy,zzzz)
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      xx=coordinate(ix,1)-center_trho(1,1)
      yy=coordinate(iy,2)-center_trho(2,1)
      zz=coordinate(iz,3)-center_trho(3,1)
      rr=sqrt(xx*xx+yy*yy+zz*zz)+1.d-50 ; xxxx=xx/rr ; yyyy=yy/rr ; zzzz=zz/rr
      rholm2box=rholm2box+rr**ll*ylm(xxxx,yyyy,zzzz,ll,m)*trho(ix,iy,iz)*hvol
    end do
    end do
    end do
    rholm2(lm,1)=rholm2box
  end do
  end do
  
  case(2)
  
  num_center=natom
  allocate (rholm((lmax_lmp+1)**2,natom))
  allocate (rholm2((lmax_lmp+1)**2,natom))
  allocate(itrho(natom))
  allocate(center_trho(3,natom))
  allocate(rion2(3,natom))
  rion2(:,:)=rion(:,:)
  
  !$OMP parallel do private(icen, lm, jj)
  do icen=1,num_center
    do lm=1,(lmax_lmp+1)**2
      rholm(lm,icen)=0.d0
      rholm2(lm,icen)=0.d0
    end do
    do jj=1,3
      center_trho(jj,icen)=rion2(jj,icen)
    end do
    itrho(icen)=1
  end do
  
  do icen=1,poisson%npole_partial
    do ll=0,lmax_lmp
    do m=-ll,ll
      lm=ll*ll+ll+1+m
      rholm2box=0.d0
  !$OMP parallel do reduction ( + : rholm2box)&
  !$OMP private(jj,ix,iy,iz,xx,yy,zz,rr,xxxx,yyyy,zzzz)
      do jj=1,ig_num(icen)
        ix=ig(1,jj,icen)
        iy=ig(2,jj,icen)
        iz=ig(3,jj,icen)
        xx=coordinate(ix,1)-rion2(1,poisson%ipole_tbl(icen))
        yy=coordinate(iy,2)-rion2(2,poisson%ipole_tbl(icen))
        zz=coordinate(iz,3)-rion2(3,poisson%ipole_tbl(icen))
        rr=sqrt(xx*xx+yy*yy+zz*zz)+1.d-50 ; xxxx=xx/rr ; yyyy=yy/rr ; zzzz=zz/rr
        rholm2box=rholm2box+rr**ll*ylm(xxxx,yyyy,zzzz,ll,m)*trho(ix,iy,iz)*hvol
      end do
      rholm2(lm,poisson%ipole_tbl(icen))=rholm2box
    end do
    end do
  end do
  
  case(3)
  
  num_center=poisson%npole_total
  
  allocate (rholm((lmax_lmp+1)**2,num_center))
  allocate (rholm2((lmax_lmp+1)**2,num_center))
  allocate(itrho(num_center))
  allocate(center_trho(3,num_center))
  allocate(center_trho_nume_deno(4,num_center))
  allocate(center_trho_nume_deno2(4,num_center))
  
  !$OMP parallel do private(icen, jj, lm)
  do icen=1,num_center
    do jj=1,4
      center_trho_nume_deno2(jj,icen)=0.d0
    end do
    do lm=1,(lmax_lmp+1)**2
      rholm(lm,icen)=0.d0
      rholm2(lm,icen)=0.d0
    end do
  end do
  
  do ii=1,poisson%npole_partial
    sum1=0.d0
    sumbox1=0.d0
    sumbox2=0.d0
    sumbox3=0.d0
  !$OMP parallel do reduction (+ : sumbox1, sumbox2, sumbox3, sum1) &
  !$OMP private(jj,ixbox,iybox,izbox,xx,yy,zz)
    do jj=1,ig_num(ii)
      ixbox=ig(1,jj,ii)
      iybox=ig(2,jj,ii)
      izbox=ig(3,jj,ii)
      xx=coordinate(ixbox,1)
      yy=coordinate(iybox,2)
      zz=coordinate(izbox,3)
      sumbox1=sumbox1+trho(ixbox,iybox,izbox)*xx
      sumbox2=sumbox2+trho(ixbox,iybox,izbox)*yy
      sumbox3=sumbox3+trho(ixbox,iybox,izbox)*zz
      sum1=sum1+trho(ixbox,iybox,izbox)
    end do
    center_trho_nume_deno2(1,poisson%ipole_tbl(ii))=sumbox1
    center_trho_nume_deno2(2,poisson%ipole_tbl(ii))=sumbox2
    center_trho_nume_deno2(3,poisson%ipole_tbl(ii))=sumbox3
    center_trho_nume_deno2(4,poisson%ipole_tbl(ii))=sum1
  end do
  
  call comm_summation(center_trho_nume_deno2,center_trho_nume_deno,4*poisson%npole_total,info%icomm_r)
  
  do ii=1,poisson%npole_total
    if(center_trho_nume_deno(4,ii)*hvol>=1.d-12)then
      itrho(ii)=1
      center_trho(1:3,ii)=center_trho_nume_deno(1:3,ii)/center_trho_nume_deno(4,ii)
    else
      itrho(ii)=0
      center_trho(1:3,ii)=0.d0
    end if
  end do
  
!  if(omp_get_max_threads() > 16) then
!  !$omp parallel shared(rholm3,lmax_lmp)
!  !$omp master
!      allocate(rholm3((lmax_lmp+1)**2,0:ceiling_pow2(omp_get_num_threads())-1))
!  !$omp end master
!  !$omp end parallel
!
!    lmax_lmp_tmp=lmax_lmp !iwata
!    rholm2=0.d0
!    rholm3=0.d0
!    do ii=1,poisson%npole_partial
!      pl=poisson%ipole_tbl(ii)
!      cl=ig_num(ii)
!      if(itrho(pl)==1)then
!  !$omp parallel default(none) &
!  !$omp          shared(ig,coordinate,center_trho,trho,rholm3) &
!  !$omp          private(tid,kk,jj,ll,lm,ixbox,iybox,izbox,xx,yy,zz,rr,rinv,xxxx,yyyy,zzzz) &
!  !$omp          firstprivate(ii,pl,cl,lmax_lmp_tmp,hvol) !iwata
!        tid=omp_get_thread_num()
!        rholm3(:,tid)=0.d0
!
!  !$omp do
!        do jj=1,cl
!          ixbox=ig(1,jj,ii)
!          iybox=ig(2,jj,ii)
!          izbox=ig(3,jj,ii)
!          xx=coordinate(ixbox,1)-center_trho(1,pl)
!          yy=coordinate(iybox,2)-center_trho(2,pl)
!          zz=coordinate(izbox,3)-center_trho(3,pl)
!          rr=sqrt(xx*xx+yy*yy+zz*zz)+1.d0-50.d0
!          rinv=1.0d0/rr
!          xxxx=xx*rinv
!          yyyy=yy*rinv
!          zzzz=zz*rinv
!          do ll=0,lmax_lmp_tmp !iwata
!          do m=-ll,ll
!            lm=ll*ll+ll+1+m
!            rholm3(lm,tid)=rholm3(lm,tid)+rr**ll*ylm(xxxx,yyyy,zzzz,ll,m)*trho(ixbox,iybox,izbox)*hvol
!          end do
!          end do
!        end do
!  !$omp end do
!
!        kk = ceiling_pow2(omp_get_num_threads())/2
!        do while(kk > 0)
!          if(tid < kk) then
!            rholm3(:,tid) = rholm3(:,tid) + rholm3(:,tid+kk)
!          end if
!          kk = kk/2
!  !$omp barrier
!        end do
!  !$omp end parallel
!      end if
!      rholm2(:,pl)=rholm3(:,0)
!    end do
!    deallocate(rholm3)
!  else
    rholm2=0.d0
    do ii=1,poisson%npole_partial
      if(itrho(poisson%ipole_tbl(ii))==1)then
        rholm=0.d0
        do ll=0,lmax_lmp
        do m=-ll,ll
          lm=ll*ll+ll+1+m
          rholm2box=0.d0
  !$OMP parallel do reduction ( + : rholm2box)&
  !$OMP private(jj,ixbox,iybox,izbox,xx,yy,zz,rr,xxxx,yyyy,zzzz)
          do jj=1,ig_num(ii)
            ixbox=ig(1,jj,ii)
            iybox=ig(2,jj,ii)
            izbox=ig(3,jj,ii)
            xx=coordinate(ixbox,1)-center_trho(1,poisson%ipole_tbl(ii))
            yy=coordinate(iybox,2)-center_trho(2,poisson%ipole_tbl(ii))
            zz=coordinate(izbox,3)-center_trho(3,poisson%ipole_tbl(ii))
            rr=sqrt(xx*xx+yy*yy+zz*zz)+1.d-50 ; xxxx=xx/rr ; yyyy=yy/rr ; zzzz=zz/rr
            rholm2box=rholm2box+rr**ll*ylm(xxxx,yyyy,zzzz,ll,m)*trho(ixbox,iybox,izbox)*hvol
          end do
          rholm2(lm,poisson%ipole_tbl(ii))=rholm2box
        end do
        end do
      end if
    end do
!  endif
  
  deallocate(center_trho_nume_deno)
  deallocate(center_trho_nume_deno2)
  
  end select
  
  if(info%isize_r==1)then
  !$OMP parallel do
    do icen=1,num_center
      rholm(:,icen)=rholm2(:,icen)
    end do
  else
    call comm_summation(rholm2,rholm,(lmax_lmp+1)**2*num_center,info%icomm_r)
  end if
  
  !$OMP parallel do private(iz,iy,ix) collapse(2)
  do iz=mg%is_array(3),mg%ie_array(3)
  do iy=mg%is_array(2),mg%ie_array(2)
  do ix=mg%is_array(1),mg%ie_array(1)
    wk2(ix,iy,iz)=0.d0
  end do
  end do
  end do
  
  do k=1,3
  
  !$OMP parallel do
    do jj=1,lg%num(1)*lg%num(2)*lg%num(3)/minval(lg%num(1:3))*6*ndh
      poisson%wkbound2(jj)=0.d0
    end do
  
    icount=mg%num(1)*mg%num(2)*mg%num(3)/mg%num(k)*2*ndh

  !$OMP parallel do
    do ii=0,nproc_xyz(k)-1
      istart(ii)=ii*icount/nproc_xyz(k)+1
      iend(ii)=(ii+1)*icount/nproc_xyz(k)
    end do
          
    do ll=0,lmax_lmp
      do lm=ll**2+1,(ll+1)**2
        l2(lm)=ll
      end do
    end do
  
  !$OMP parallel do &
  !$OMP private(xx,yy,zz,rr,xxxx,yyyy,zzzz,lm,ll,sum1,ylm2,rrrr,xp2,yp2,zp2,xy,yz,xz,rinv,rbox,deno,icen)
    do jj=istart(myrank_xyz(k)),iend(myrank_xyz(k))
      do icen=1,num_center
        if(itrho(icen)==1)then
          xx=coordinate(poisson%ig_bound(1,jj,k),1)-center_trho(1,icen)
          yy=coordinate(poisson%ig_bound(2,jj,k),2)-center_trho(2,icen)
          zz=coordinate(poisson%ig_bound(3,jj,k),3)-center_trho(3,icen)
          rr=sqrt(xx**2+yy**2+zz**2)+1.d-50
          rinv=1.d0/rr
  !        xxxx=xx/rr ; yyyy=yy/rr ; zzzz=zz/rr
          xx=xx*rinv ; yy=yy*rinv ; zz=zz*rinv
  
          xp2=xx**2
          yp2=yy**2
          zp2=zz**2
          xy=xx*yy
          yz=yy*zz
          xz=xx*zz
  
          rrrr=(xp2+yp2+zp2)**2
  
          ylm2(1)=1.d0
          ylm2(2)=yy
          ylm2(3)=zz
          ylm2(4)=xx
          ylm2(5)=sqrt(3.d0)*xy                      ! lm=5  (2 -2)
          ylm2(6)=sqrt(3.d0)*yz                      ! lm=6  (2 -1)
          ylm2(7)=(2*zp2-xp2-yp2)/2.d0               ! lm=7  (2 0)
          ylm2(8)=sqrt(3.d0)*xz                      ! lm=8  (2 1)
          ylm2(9)=sqrt(3.d0/4.d0)*(xp2-yp2)          ! lm=9  (2 2)
          ylm2(10)=sqrt(5.d0/8.d0)*yy*(3*xp2-yp2)    ! lm=10 (3 -3)
          ylm2(11)=sqrt(15.d0)*xx*yy*zz               ! lm=11 (3 -2)
          ylm2(12)=sqrt(3.d0/8.d0)*yy*(4*zp2-xp2-yp2)  ! lm=12 (3 -1)
          ylm2(13)=zz*(2*zp2-3*xp2-3*yp2)/2.d0         ! lm=13 (3 0)
          ylm2(14)=sqrt(3.d0/8.d0)*xx*(4*zp2-xp2-yp2)  ! lm=14 (3 1)
          ylm2(15)=sqrt(15.d0/4.d0)*zz*(xp2-yp2)       ! lm=15 (3 2)
          ylm2(16)=sqrt(5.d0/8.d0)*xx*(xp2-3*yp2)      ! lm=16 (3 3)
  
          ylm2(17)=rrrr*sqrt(35.d0)/2.d0*xy*(xp2-yp2)
          ylm2(18)=rrrr*sqrt(35.d0/8.d0)*yz*(3*xp2-yp2)
          ylm2(19)=rrrr*sqrt(5.d0)/2.d0*xy*(7*zp2-1.d0)
          ylm2(20)=rrrr*sqrt(5.d0/8.d0)*yz*(7*zp2-3.d0)
          ylm2(21)=rrrr*(35*zp2**2-30*zp2+3.d0)/8.d0
          ylm2(22)=rrrr*sqrt(5.d0/8.d0)*xz*(7*zp2-3.d0)
          ylm2(23)=rrrr*sqrt(5.d0)/4.d0*(7*zp2-1)*(xp2-yp2)
          ylm2(24)=rrrr*sqrt(35.d0/8.d0)*xz*(xp2-3*yp2)
          ylm2(25)=rrrr*sqrt(35.d0)/8.d0*(xp2**2+yp2**2-6*xp2*yp2)
  
          deno(1)=rinv
          rbox=rinv*rinv
          deno(2:4)=rbox
          rbox=rbox*rinv
          deno(5:9)=rbox
          rbox=rbox*rinv
          deno(10:16)=rbox
          rbox=rbox*rinv
          deno(17:25)=rbox
  
          sum1=0.d0
          do lm=1,(lmax_lmp+1)**2
            sum1=sum1+ylm2(lm)*deno(lm)*rholm(lm,icen)
          end do
          poisson%wkbound2(jj) = poisson%wkbound2(jj) + sum1
        end if
      end do
    end do
  
    if(info%isize_r==1)then
  !$OMP parallel do
      do jj=1,icount
        poisson%wkbound(jj)=poisson%wkbound2(jj)
      end do
    else
      call comm_summation( &
        poisson%wkbound2,              poisson%wkbound,              icount/2, comm_xyz(k), 0                    )
      call comm_summation( &
        poisson%wkbound2(icount/2+1:), poisson%wkbound(icount/2+1:), icount/2, comm_xyz(k), nproc_xyz(k)-1)
    end if
  
    if(myrank_xyz(k)==0) then
  !$OMP parallel do
      do jj=1,icount/2
        wk2(poisson%ig_bound(1,jj,k),poisson%ig_bound(2,jj,k),poisson%ig_bound(3,jj,k))=poisson%wkbound(jj)
      end do
    end if
    if(myrank_xyz(k)==nproc_xyz(k)-1) then
  !$OMP parallel do
      do jj=icount/2+1,icount
        wk2(poisson%ig_bound(1,jj,k),poisson%ig_bound(2,jj,k),poisson%ig_bound(3,jj,k))=poisson%wkbound(jj)
      end do
    end if
  
  end do
  
  deallocate(rholm,rholm2)
  deallocate(center_trho)
  
  return
  
end subroutine poisson_boundary

end module poisson_isolated
