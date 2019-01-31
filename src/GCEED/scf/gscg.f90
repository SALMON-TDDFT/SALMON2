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
!=======================================================================
!======================================= Conjugate-Gradient minimization

subroutine sgscg(mg,spsi,iflag,itotmst,mst,hvol,ilsda,nproc_ob,nproc_ob_spin,iparaway_ob,elp3, &
                 rxk_ob,rhxk_ob,rgk_ob,rpk_ob)
use inputoutput, only: ncg
use structures, only: s_rgrid,s_wavefunction
use salmon_parallel, only: nproc_group_grid, nproc_group_global, nproc_group_korbital
use salmon_communication, only: comm_summation, comm_bcast
use misc_routines, only: get_wtime
!$ use omp_lib
implicit none

type(s_rgrid),intent(in) :: mg
type(s_wavefunction),intent(inout) :: spsi
integer,intent(inout) :: iflag
  integer,intent(in)    :: itotmst
  integer,intent(in)    :: mst(2)
  real(8),intent(in)    :: hvol
  integer,intent(in)    :: ilsda
  integer,intent(in)    :: nproc_ob
  integer,intent(in)    :: nproc_ob_spin
  integer,intent(in)    :: iparaway_ob
  real(8),intent(out)    :: elp3(3000)
  real(8),intent(inout) :: rxk_ob(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),1:spsi%numo)
  real(8),intent(inout) :: rhxk_ob(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),1:spsi%numo)
  real(8),intent(inout) :: rgk_ob(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),1:spsi%numo)
  real(8),intent(inout) :: rpk_ob(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),1:spsi%numo)
integer :: iter,iob,job
integer :: ix,iy,iz
integer :: is,iobsta(2),iobend(2)
real(8) :: sum0,sum1
real(8) :: sum_ob0(itotMST)
real(8) :: sum_obmat0(itotMST,itotMST),sum_obmat1(itotMST,itotMST)
real(8) :: xkHxk_ob(itotMST),xkxk_ob(itotMST),Rk_ob(itotMST)
real(8) :: gkgk_ob(itotMST),pkpk_ob(itotMST),xkpk_ob(itotMST)
real(8) :: pkHxk_ob(itotMST),pkHpk_ob(itotMST)
real(8) :: uk,alpha,Ak,Bk,Ck
real(8) , allocatable :: gk(:,:,:)
real(8) :: elp2(2000)
real(8):: tpsi(mg%is_array(1):mg%ie_array(1),  &
               mg%is_array(2):mg%ie_array(2),  &
               mg%is_array(3):mg%ie_array(3))
real(8):: htpsi(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
real(8):: rmatbox_m(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
integer :: iob_myob,job_myob
integer :: iob_allob
integer :: icorr_iob,icorr_job
integer :: iroot
integer :: is_sta,is_end

allocate (gk(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))

call set_isstaend(is_sta,is_end,ilsda,nproc_ob,nproc_ob_spin)

!$OMP parallel do private(iz,iy,ix) collapse(2)
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

do iob=1,spsi%numo
  call calc_allob(iob,iob_allob)

!$OMP parallel do private(iz,iy,ix) collapse(2)
  do iz=mg%is(3),mg%ie(3)
  do iy=mg%is(2),mg%ie(2)
  do ix=mg%is(1),mg%ie(1)
    rxk_ob(ix,iy,iz,iob)=spsi%rwf(ix,iy,iz,1,1,iob,1)
    tpsi(ix,iy,iz)=rxk_ob(ix,iy,iz,iob)
  end do
  end do
  end do
  call r_hpsi2_buf(tpsi,htpsi,iob_allob,1,0,0)
  
!$OMP parallel do private(iz,iy,ix) collapse(2)
  do iz=mg%is(3),mg%ie(3)
  do iy=mg%is(2),mg%ie(2)
  do ix=mg%is(1),mg%ie(1)
    rhxk_ob(ix,iy,iz,iob)=htpsi(ix,iy,iz)
  end do
  end do
  end do
end do

call inner_product7(mg,itotmst,spsi%numo,rxk_ob,rhxk_ob,xkHxk_ob,elp3,hvol)

xkxk_ob(:)=1.d0 
Rk_ob(:)=xkHxk_ob(:)/xkxk_ob(:)

Iteration : do iter=1,ncg
elp2(2)=get_wtime()
  do iob=1,spsi%numo
    call calc_allob(iob,iob_allob)
!$OMP parallel do private(iz,iy,ix) collapse(2)
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      rgk_ob(ix,iy,iz,iob) = 2*( rhxk_ob(ix,iy,iz,iob) - Rk_ob(iob_allob)*rxk_ob(ix,iy,iz,iob) )
    end do
    end do
    end do
  end do
 if(nproc_ob==1)then
    do is=is_sta,is_end
    do iob=iobsta(is),iobend(is)
      do job=iobsta(is),iob-1
        sum0=0.d0
  !$omp parallel do private(iz,iy,ix) collapse(2) reduction(+ : sum0)
        do iz=mg%is(3),mg%ie(3)
        do iy=mg%is(2),mg%ie(2)
        do ix=mg%is(1),mg%ie(1)
          sum0=sum0+spsi%rwf(ix,iy,iz,1,1,job,1)*rgk_ob(ix,iy,iz,iob)
        end do
        end do
        end do
        sum_obmat0(iob,job)=sum0*Hvol
      end do
    end do 
    end do
    call comm_summation(sum_obmat0,sum_obmat1,itotMST*itotMST,nproc_group_global)
    do is=is_sta,is_end
    do iob=iobsta(is),iobend(is)
      do job=iobsta(is),iob-1
  !$omp parallel do private(iz,iy,ix) collapse(2)
        do iz=mg%is(3),mg%ie(3)
        do iy=mg%is(2),mg%ie(2)
        do ix=mg%is(1),mg%ie(1)
          rgk_ob(ix,iy,iz,iob)=rgk_ob(ix,iy,iz,iob)-sum_obmat1(iob,job)*spsi%rwf(ix,iy,iz,1,1,job,1)
        end do
        end do
        end do
      end do
    end do
    end do
  else
    do iob=iobsta(is),iobend(is)
      call calc_myob(iob,iob_myob,ilsda,nproc_ob,iparaway_ob,itotmst,nproc_ob_spin,mst)
      call check_corrkob(iob,1,icorr_iob,ilsda,nproc_ob,iparaway_ob,itotmst,spsi%ik_s,spsi%ik_e,nproc_ob_spin,mst)
      do job=iobsta(is),iob-1
        call calc_myob(job,job_myob,ilsda,nproc_ob,iparaway_ob,itotmst,nproc_ob_spin,mst)
        call check_corrkob(job,1,icorr_job,ilsda,nproc_ob,iparaway_ob,itotmst,spsi%ik_s,spsi%ik_e,nproc_ob_spin,mst)
        if(icorr_job==1)then
  !$omp parallel do private(iz,iy,ix) collapse(2)
          do iz=mg%is(3),mg%ie(3)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
            rmatbox_m(ix,iy,iz)=spsi%rwf(ix,iy,iz,1,1,job_myob,1)
          end do
          end do
          end do
        end if
        call calc_iroot(job,iroot,ilsda,nproc_ob,iparaway_ob,itotmst,nproc_ob_spin,mst)
        call comm_bcast(rmatbox_m,nproc_group_grid,iroot)
        sum0=0.d0
  !$omp parallel do private(iz,iy,ix) collapse(2) reduction(+ : sum0)
        do iz=mg%is(3),mg%ie(3)
        do iy=mg%is(2),mg%ie(2)
        do ix=mg%is(1),mg%ie(1)
          sum0=sum0+rmatbox_m(ix,iy,iz)*rgk_ob(ix,iy,iz,iob_myob)
        end do
        end do
        end do
        sum0=sum0*Hvol
        call comm_summation(sum0,sum1,nproc_group_korbital)
        do iz=mg%is(3),mg%ie(3)
        do iy=mg%is(2),mg%ie(2)
        do ix=mg%is(1),mg%ie(1)
          rgk_ob(ix,iy,iz,iob_myob)=rgk_ob(ix,iy,iz,iob_myob)-sum1*rmatbox_m(ix,iy,iz)
        end do
        end do
        end do
      end do
    end do
  end if 
 call inner_product7(mg,itotmst,spsi%numo,rgk_ob,rgk_ob,sum_ob0,elp3,hvol)
 if ( iter==1 ) then
    do iob=1,spsi%numo
      call calc_allob(iob,iob_allob)
!$OMP parallel do private(iz,iy,ix) collapse(2)
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        rpk_ob(ix,iy,iz,iob) = -rgk_ob(ix,iy,iz,iob)
      end do
      end do
      end do
    end do
  else
    do iob=1,spsi%numo
      call calc_allob(iob,iob_allob)
      uk=sum_ob0(iob_allob)/gkgk_ob(iob_allob)
!$OMP parallel do private(iz,iy,ix)
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        rpk_ob(ix,iy,iz,iob) = -rgk_ob(ix,iy,iz,iob) + uk*rpk_ob(ix,iy,iz,iob)
      end do
      end do
      end do
    end do
  end if 
  gkgk_ob(:)=sum_ob0(:)
  call inner_product7(mg,itotmst,spsi%numo,rxk_ob,rpk_ob,xkpk_ob,elp3,hvol)
  call inner_product7(mg,itotmst,spsi%numo,rpk_ob,rpk_ob,pkpk_ob,elp3,hvol)
  call inner_product7(mg,itotmst,spsi%numo,rpk_ob,rhxk_ob,pkHxk_ob,elp3,hvol)

  do iob=1,spsi%numo
    call calc_allob(iob,iob_allob)
!$OMP parallel do private(iz,iy,ix) collapse(2)
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      tpsi(ix,iy,iz) = rpk_ob(ix,iy,iz,iob)
    end do
    end do
    end do
    call r_hpsi2_buf(tpsi,htpsi,iob_allob,1,0,0)
    
!$OMP parallel do private(iz,iy,ix) collapse(2)
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
       rgk_ob(ix,iy,iz,iob)=htpsi(ix,iy,iz)
    end do
    end do
    end do
  end do
  call inner_product7(mg,itotmst,spsi%numo,rpk_ob,rgk_ob,pkHpk_ob,elp3,hvol)
 do iob=1,spsi%numo
    call calc_allob(iob,iob_allob)
    Ak=pkHpk_ob(iob_allob)*xkpk_ob(iob_allob)-pkHxk_ob(iob_allob)*pkpk_ob(iob_allob)
    Bk=pkHpk_ob(iob_allob)*xkxk_ob(iob_allob)-xkHxk_ob(iob_allob)*pkpk_ob(iob_allob)
    Ck=pkHxk_ob(iob_allob)*xkxk_ob(iob_allob)-xkHxk_ob(iob_allob)*xkpk_ob(iob_allob)
    alpha=(-Bk+sqrt(Bk*Bk-4*Ak*Ck))/(2*Ak)

!$OMP parallel do private(iz,iy,ix) collapse(2)
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      rxk_ob(ix,iy,iz,iob) = rxk_ob(ix,iy,iz,iob) + alpha*rpk_ob(ix,iy,iz,iob)
      rhxk_ob(ix,iy,iz,iob) = rhxk_ob(ix,iy,iz,iob) + alpha*rgk_ob(ix,iy,iz,iob)
    end do
    end do
    end do
  end do
  call inner_product7(mg,itotmst,spsi%numo,rxk_ob,rhxk_ob,xkHxk_ob,elp3,hvol)
  call inner_product7(mg,itotmst,spsi%numo,rxk_ob,rxk_ob,xkxk_ob,elp3,hvol)
  Rk_ob(:)=xkHxk_ob(:)/xkxk_ob(:)


end do Iteration

call inner_product7(mg,itotmst,spsi%numo,rxk_ob,rxk_ob,sum_ob0,elp3,hvol)
do iob=1,spsi%numo
  call calc_allob(iob,iob_allob)
!$OMP parallel do private(iz,iy,ix) collapse(2)
  do iz=mg%is(3),mg%ie(3)
  do iy=mg%is(2),mg%ie(2)
  do ix=mg%is(1),mg%ie(1)
    spsi%rwf(ix,iy,iz,1,1,iob,1)=rxk_ob(ix,iy,iz,iob)/sqrt(sum_ob0(iob_allob))
  end do
  end do
  end do
end do

if(iflag.eq.1) then
  iflag=0
end if

deallocate (gk)
return

end subroutine sgscg
