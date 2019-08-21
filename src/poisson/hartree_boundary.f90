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
module hartree_boundary_sub
  implicit none

contains

!============================ Hartree potential (Solve Poisson equation)

subroutine hartree_boundary(lg,mg,ng,info_field,trho,wk2,wkbound_h,wk2bound_h,   &
                            layout_multipole,lmax_lmp,igc_is,igc_ie,gridcoo,hvol,iflag_ps,num_pole,inum_mxin_s,   &
                            iamax,maxval_pole,num_pole_myrank,icorr_polenum,icount_pole,icorr_xyz_pole,   &
                            ibox_icoobox_bound,icoobox_bound)
  use inputoutput, only: natom,rion
  use structures, only: s_rgrid,s_field_parallel
  use salmon_parallel, only: nproc_id_global, nproc_size_global, nproc_group_h
  use salmon_communication, only: comm_summation
  use timer
  
  use omp_lib, only: omp_get_num_threads, omp_get_thread_num, omp_get_max_threads
  use misc_routines, only: ceiling_pow2
  
  implicit none
  integer,parameter :: ndh=4
  type(s_rgrid),intent(in) :: lg
  type(s_rgrid),intent(in) :: mg
  type(s_rgrid),intent(in) :: ng
  type(s_field_parallel),intent(in) :: info_field
  real(8) :: trho(mg%is(1):mg%ie(1),    &
                 mg%is(2):mg%ie(2),      &
                 mg%is(3):mg%ie(3))
  real(8) :: wk2(ng%is(1)-ndh:ng%ie(1)+ndh,    &
                 ng%is(2)-ndh:ng%ie(2)+ndh,      &
                 ng%is(3)-ndh:ng%ie(3)+ndh)
  real(8),intent(out) :: wkbound_h(lg%num(1)*lg%num(2)*lg%num(3)/minval(lg%num(1:3))*6*ndh)
  real(8),intent(out) :: wk2bound_h(lg%num(1)*lg%num(2)*lg%num(3)/minval(lg%num(1:3))*6*ndh)
  integer,intent(in) :: layout_multipole
  integer,intent(in) :: lmax_lmp
  integer,intent(in) :: igc_is
  integer,intent(in) :: igc_ie
  real(8),intent(in) :: gridcoo(igc_is:igc_ie,3)
  real(8),intent(in) :: hvol
  integer,intent(in) :: iflag_ps
  integer,intent(in) :: num_pole
  integer,intent(in) :: inum_mxin_s(3,0:nproc_size_global-1)
  integer,intent(in) :: iamax
  integer,intent(in) :: maxval_pole
  integer,intent(in) :: num_pole_myrank
  integer,intent(in) :: icorr_polenum(iamax)
  integer,intent(in) :: icount_pole(iamax)
  integer,intent(in) :: icorr_xyz_pole(3,maxval_pole,num_pole_myrank)
  integer,intent(in) :: ibox_icoobox_bound
  integer,intent(in) :: icoobox_bound(3,ibox_icoobox_bound,3)
  integer,parameter :: maxiter=1000
  integer :: ii,jj,kk,ix,iy,iz,lm,ll,icen,pl,cl
  integer :: ixbox,iybox,izbox
  integer :: k
  integer :: istart(0:nproc_size_global-1),iend(0:nproc_size_global-1)
  integer :: icount
  integer,allocatable :: itrho(:)
  integer :: num_center
  real(8) :: ylm
  real(8) :: ylm2(25)
  integer :: l2(25)
  real(8) :: xx,yy,zz,rr,sum1,xxxx,yyyy,zzzz,rrrr,sumbox1,sumbox2,sumbox3
  real(8) :: rholm2box
  real(8),allocatable :: rholm(:,:),rholm2(:,:),rholm3(:,:)
  integer :: tid
  real(8) :: center_trho2(3)
  real(8),allocatable :: center_trho(:,:)
  real(8),allocatable :: center_trho_nume_deno(:,:)
  real(8),allocatable :: center_trho_nume_deno2(:,:)
  real(8) :: xp2,yp2,zp2,xy,yz,xz
  real(8) :: deno(25)
  real(8) :: rinv
  real(8) :: rbox
  real(8),allocatable :: rion2(:,:)
  
  !------------------------- Boundary condition (multipole expansion)
  
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
  do lm=ll**2+1,(ll+1)**2
    rholm2box=0.d0
  !$OMP parallel do reduction ( + : rholm2box)&
  !$OMP private(ix,iy,iz,xx,yy,zz,rr,xxxx,yyyy,zzzz,ylm)
    do iz=ng%is(3),ng%ie(3)
    do iy=ng%is(2),ng%ie(2)
    do ix=ng%is(1),ng%ie(1)
      xx=gridcoo(ix,1)-center_trho(1,1)
      yy=gridcoo(iy,2)-center_trho(2,1)
      zz=gridcoo(iz,3)-center_trho(3,1)
      rr=sqrt(xx*xx+yy*yy+zz*zz)+1.d-50 ; xxxx=xx/rr ; yyyy=yy/rr ; zzzz=zz/rr
      call ylm_sub(xxxx,yyyy,zzzz,lm,ylm)
      rholm2box=rholm2box+rr**ll*ylm*trho(ix,iy,iz)*Hvol
    end do
    end do
    end do
    rholm2(lm,1)=rholm2box
  end do
  end do
  
  case(2)
  
  if(iflag_ps==1)then
    num_center=natom
    allocate (rholm((lmax_lmp+1)**2,natom))
    allocate (rholm2((lmax_lmp+1)**2,natom))
    allocate(itrho(natom))
    allocate(center_trho(3,natom))
    allocate(rion2(3,natom))
    rion2(:,:)=rion(:,:)
  end if
  
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
  
  do icen=1,num_pole_myrank
    do ll=0,lmax_lmp
    do lm=ll**2+1,(ll+1)**2
      rholm2box=0.d0
  !$OMP parallel do reduction ( + : rholm2box)&
  !$OMP private(jj,ix,iy,iz,xx,yy,zz,rr,xxxx,yyyy,zzzz,ylm)
      do jj=1,icount_pole(icen)
        ix=icorr_xyz_pole(1,jj,icen)
        iy=icorr_xyz_pole(2,jj,icen)
        iz=icorr_xyz_pole(3,jj,icen)
        xx=gridcoo(ix,1)-rion2(1,icorr_polenum(icen))
        yy=gridcoo(iy,2)-rion2(2,icorr_polenum(icen))
        zz=gridcoo(iz,3)-rion2(3,icorr_polenum(icen))
        rr=sqrt(xx*xx+yy*yy+zz*zz)+1.d-50 ; xxxx=xx/rr ; yyyy=yy/rr ; zzzz=zz/rr
        call ylm_sub(xxxx,yyyy,zzzz,lm,ylm)
        rholm2box=rholm2box+rr**ll*ylm*trho(ix,iy,iz)*Hvol
      end do
      rholm2(lm,icorr_polenum(icen))=rholm2box
    end do
    end do
  end do
  
  case(3)
  
  num_center=num_pole
  
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
  
  do ii=1,num_pole_myrank
    sum1=0.d0
    sumbox1=0.d0
    sumbox2=0.d0
    sumbox3=0.d0
  !$OMP parallel do reduction (+ : sumbox1, sumbox2, sumbox3, sum1) &
  !$OMP private(jj,ixbox,iybox,izbox,xx,yy,zz)
    do jj=1,icount_pole(ii)
      ixbox=icorr_xyz_pole(1,jj,ii)
      iybox=icorr_xyz_pole(2,jj,ii)
      izbox=icorr_xyz_pole(3,jj,ii)
      xx=gridcoo(ixbox,1)
      yy=gridcoo(iybox,2)
      zz=gridcoo(izbox,3)
      sumbox1=sumbox1+trho(ixbox,iybox,izbox)*xx
      sumbox2=sumbox2+trho(ixbox,iybox,izbox)*yy
      sumbox3=sumbox3+trho(ixbox,iybox,izbox)*zz
      sum1=sum1+trho(ixbox,iybox,izbox)
    end do
    center_trho_nume_deno2(1,icorr_polenum(ii))=sumbox1
    center_trho_nume_deno2(2,icorr_polenum(ii))=sumbox2
    center_trho_nume_deno2(3,icorr_polenum(ii))=sumbox3
    center_trho_nume_deno2(4,icorr_polenum(ii))=sum1
  end do
  
  call timer_begin(LOG_ALLREDUCE_HARTREE)
  call comm_summation(center_trho_nume_deno2,center_trho_nume_deno,4*num_pole,nproc_group_h)
  call timer_end(LOG_ALLREDUCE_HARTREE)
  
  do ii=1,num_pole
    if(center_trho_nume_deno(4,ii)*Hvol>=1.d-12)then
      itrho(ii)=1
      center_trho(1:3,ii)=center_trho_nume_deno(1:3,ii)/center_trho_nume_deno(4,ii)
    else
      itrho(ii)=0
      center_trho(1:3,ii)=0.d0
    end if
  end do
  
  if(omp_get_max_threads() > 16) then
  !$omp parallel shared(rholm3,lmax_lmp)
  !$omp master
      allocate(rholm3((lmax_lmp+1)**2,0:ceiling_pow2(omp_get_num_threads())-1))
  !$omp end master
  !$omp end parallel
  
    rholm2=0.d0
    rholm3=0.d0
    do ii=1,num_pole_myrank
      pl=icorr_polenum(ii)
      cl=icount_pole(ii)
      if(itrho(pl)==1)then
  !$omp parallel default(none) &
  !$omp          shared(icorr_xyz_pole,gridcoo,center_trho,trho,rholm3) &
  !$omp          private(tid,kk,jj,ll,lm,ixbox,iybox,izbox,xx,yy,zz,rr,rinv,xxxx,yyyy,zzzz,ylm) &
  !$omp          firstprivate(ii,pl,cl,lmax_lmp,Hvol)
        tid=omp_get_thread_num()
        rholm3(:,tid)=0.d0
  
  !$omp do
        do jj=1,cl
          ixbox=icorr_xyz_pole(1,jj,ii)
          iybox=icorr_xyz_pole(2,jj,ii)
          izbox=icorr_xyz_pole(3,jj,ii)
          xx=gridcoo(ixbox,1)-center_trho(1,pl)
          yy=gridcoo(iybox,2)-center_trho(2,pl)
          zz=gridcoo(izbox,3)-center_trho(3,pl)
          rr=sqrt(xx*xx+yy*yy+zz*zz)+1.d0-50.d0
          rinv=1.0d0/rr
          xxxx=xx*rinv
          yyyy=yy*rinv
          zzzz=zz*rinv
          do ll=0,lmax_lmp
          do lm=ll**2+1,(ll+1)**2
            call ylm_sub(xxxx,yyyy,zzzz,lm,ylm)
            rholm3(lm,tid)=rholm3(lm,tid)+rr**ll*ylm*trho(ixbox,iybox,izbox)*Hvol
          end do
          end do
        end do
  !$omp end do
  
        kk = ceiling_pow2(omp_get_num_threads())/2
        do while(kk > 0)
          if(tid < kk) then
            rholm3(:,tid) = rholm3(:,tid) + rholm3(:,tid+kk)
          end if
          kk = kk/2
  !$omp barrier
        end do
  !$omp end parallel
      end if
      rholm2(:,pl)=rholm3(:,0)
    end do
    deallocate(rholm3)
  else
    rholm2=0.d0
    do ii=1,num_pole_myrank
      if(itrho(icorr_polenum(ii))==1)then
        rholm=0.d0
        do ll=0,lmax_lmp
        do lm=ll**2+1,(ll+1)**2
          rholm2box=0.d0
  !$OMP parallel do reduction ( + : rholm2box)&
  !$OMP private(jj,ixbox,iybox,izbox,xx,yy,zz,rr,xxxx,yyyy,zzzz,ylm)
          do jj=1,icount_pole(ii)
            ixbox=icorr_xyz_pole(1,jj,ii)
            iybox=icorr_xyz_pole(2,jj,ii)
            izbox=icorr_xyz_pole(3,jj,ii)
            xx=gridcoo(ixbox,1)-center_trho(1,icorr_polenum(ii))
            yy=gridcoo(iybox,2)-center_trho(2,icorr_polenum(ii))
            zz=gridcoo(izbox,3)-center_trho(3,icorr_polenum(ii))
            rr=sqrt(xx*xx+yy*yy+zz*zz)+1.d-50 ; xxxx=xx/rr ; yyyy=yy/rr ; zzzz=zz/rr
            call ylm_sub(xxxx,yyyy,zzzz,lm,ylm)
            rholm2box=rholm2box+rr**ll*ylm*trho(ixbox,iybox,izbox)*Hvol
          end do
          rholm2(lm,icorr_polenum(ii))=rholm2box
        end do
        end do
      end if
    end do
  endif
  
  deallocate(center_trho_nume_deno)
  deallocate(center_trho_nume_deno2)
  
  end select
  
  if(nproc_size_global==1)then
  !$OMP parallel do
    do icen=1,num_center
      rholm(:,icen)=rholm2(:,icen)
    end do
  else
    call timer_begin(LOG_ALLREDUCE_HARTREE)
    call comm_summation(rholm2,rholm,(lmax_lmp+1)**2*num_center,nproc_group_h)
    call timer_end(LOG_ALLREDUCE_HARTREE)
  end if
  
  !$OMP parallel do private(iz,iy,ix) collapse(2)
  do iz=ng%is(3)-ndh,ng%ie(3)+ndh
  do iy=ng%is(2)-ndh,ng%ie(2)+ndh
  do ix=ng%is(1)-ndh,ng%ie(1)+ndh
    wk2(ix,iy,iz)=0.d0
  end do
  end do
  end do
  
  do k=1,3
  
  !$OMP parallel do
    do jj=1,lg%num(1)*lg%num(2)*lg%num(3)/minval(lg%num(1:3))*6*ndh
      wk2bound_h(jj)=0.d0
    end do
  
    icount=inum_Mxin_s(1,nproc_id_global)*inum_Mxin_s(2,nproc_id_global)  &
           *inum_Mxin_s(3,nproc_id_global)/inum_Mxin_s(k,nproc_id_global)*2*ndh
  !$OMP parallel do
    do ii=0,info_field%isize(k)-1
      istart(ii)=ii*icount/info_field%isize(k)+1
      iend(ii)=(ii+1)*icount/info_field%isize(k)
    end do
          
    do ll=0,lmax_lmp
      do lm=ll**2+1,(ll+1)**2
        l2(lm)=ll
      end do
    end do
  
  !$OMP parallel do &
  !$OMP private(xx,yy,zz,rr,xxxx,yyyy,zzzz,lm,ll,sum1,ylm2,rrrr,xp2,yp2,zp2,xy,yz,xz,rinv,rbox,deno,icen)
    do jj=istart(info_field%id(k)),iend(info_field%id(k))
      do icen=1,num_center
        if(itrho(icen)==1)then
          xx=gridcoo(icoobox_bound(1,jj,k),1)-center_trho(1,icen)
          yy=gridcoo(icoobox_bound(2,jj,k),2)-center_trho(2,icen)
          zz=gridcoo(icoobox_bound(3,jj,k),3)-center_trho(3,icen)
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
          wk2bound_h(jj) = wk2bound_h(jj) + sum1
        end if
      end do
    end do
  
    if(nproc_size_global==1)then
  !$OMP parallel do
      do jj=1,icount
        wkbound_h(jj)=wk2bound_h(jj)
      end do
    else
      call timer_begin(LOG_ALLREDUCE_HARTREE)
      call comm_summation( &
        wk2bound_h,              wkbound_h,              icount/2, info_field%icomm(k), 0                    )
      call comm_summation( &
        wk2bound_h(icount/2+1:), wkbound_h(icount/2+1:), icount/2, info_field%icomm(k), info_field%isize(k)-1)
      call timer_end(LOG_ALLREDUCE_HARTREE)
    end if
  
    if(info_field%id(k)==0) then
  !$OMP parallel do
      do jj=1,icount/2
        wk2(icoobox_bound(1,jj,k),icoobox_bound(2,jj,k),icoobox_bound(3,jj,k))=wkbound_h(jj)
      end do
    end if
    if(info_field%id(k)==info_field%isize(k)-1) then
  !$OMP parallel do
      do jj=icount/2+1,icount
        wk2(icoobox_bound(1,jj,k),icoobox_bound(2,jj,k),icoobox_bound(3,jj,k))=wkbound_h(jj)
      end do
    end if
  
  end do
  
  deallocate(rholm,rholm2)
  deallocate(center_trho)
  
  return
  
end subroutine hartree_boundary

end module hartree_boundary_sub
