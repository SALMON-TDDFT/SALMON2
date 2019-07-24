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
module gscg_sub
  implicit none

contains

!=======================================================================
!======================================= Conjugate-Gradient minimization

subroutine sgscg(mg,system,info,stencil,srg_ob_1,spsi,iflag,itotmst,mst,ilsda,nproc_ob,cg, &
                 info_ob,ppg,vlocal)
  use inputoutput, only: ncg,ispin
  use structures, only: s_rgrid,s_dft_system,s_orbital_parallel,s_orbital,s_stencil,s_scalar,s_pp_grid,s_cg
  use salmon_parallel, only: nproc_group_grid, nproc_group_global, nproc_group_korbital
  use salmon_communication, only: comm_summation, comm_bcast
  use timer
  use hpsi_sub
  use calc_allob_sub
  use calc_iroot_sub
  use calc_myob_sub
  use check_corrkob_sub
  use set_isstaend_sub
  use sendrecv_grid, only: s_sendrecv_grid
  !$ use omp_lib
  implicit none
  
  type(s_rgrid),intent(in) :: mg
  type(s_dft_system),intent(in) :: system
  type(s_orbital_parallel) :: info
  type(s_orbital),intent(inout) :: spsi
  type(s_stencil) :: stencil
  type(s_sendrecv_grid),intent(inout) :: srg_ob_1
  type(s_pp_grid) :: ppg
  integer,intent(inout) :: iflag
  integer,intent(in)    :: itotmst
  integer,intent(in)    :: mst(2)
  integer,intent(in)    :: ilsda
  integer,intent(in)    :: nproc_ob
  type(s_cg),intent(inout)            :: cg
  type(s_orbital_parallel)       :: info_ob
  real(8),intent(in)    :: vlocal(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),ispin+1)
  integer,parameter :: nd=4
  integer :: iter,iob,job
  integer :: ix,iy,iz
  integer :: is,iobsta(2),iobend(2)
  integer :: nspin_1
  type(s_orbital)  :: stpsi
  type(s_orbital)  :: shtpsi
  type(s_scalar),allocatable :: v(:)
  real(8) :: sum0,sum1
  real(8) :: sum_ob0(itotmst)
  real(8) :: sum_obmat0(itotmst,itotmst),sum_obmat1(itotmst,itotmst)
  real(8) :: xkhxk_ob(itotmst),xkxk_ob(itotmst),rk_ob(itotmst)
  real(8) :: gkgk_ob(itotmst),pkpk_ob(itotmst),xkpk_ob(itotmst)
  real(8) :: pkhxk_ob(itotmst),pkHpk_ob(itotmst)
  real(8) :: uk,alpha,ak,bk,ck
  real(8) , allocatable :: gk(:,:,:)
  real(8):: rmatbox_m(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
  integer :: iob_myob,job_myob
  integer :: iob_allob
  integer :: icorr_iob,icorr_job
  integer :: iroot
  integer :: is_sta,is_end

  call timer_begin(LOG_GSCG_TOTAL)

  call timer_begin(LOG_GSCG_INIT)
  allocate(stpsi%rwf(mg%is_array(1):mg%ie_array(1),  &
                     mg%is_array(2):mg%ie_array(2),  &
                     mg%is_array(3):mg%ie_array(3),1,1,1,1))
  allocate(shtpsi%rwf(mg%is_array(1):mg%ie_array(1),  &
                      mg%is_array(2):mg%ie_array(2),  &
                      mg%is_array(3):mg%ie_array(3),1,1,1,1))

  nspin_1=1
  allocate(v(nspin_1))
  allocate(v(nspin_1)%f(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))

  allocate (gk(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))
  
  call set_isstaend(is_sta,is_end,ilsda)
  
  !$OMP parallel do private(iz,iy,ix) collapse(2)
  do iz=mg%is_array(3),mg%ie_array(3)
  do iy=mg%is_array(2),mg%ie_array(2)
  do ix=mg%is_array(1),mg%ie_array(1)
    stpsi%rwf(ix,iy,iz,1,1,1,1)=0.d0
  end do
  end do
  end do
  
  if(ilsda == 0)then
    iobsta(1)=1
    iobend(1)=itotmst
  else if(ilsda == 1)then
    iobsta(1)=1
    iobend(1)=mst(1)
    iobsta(2)=mst(1)+1
    iobend(2)=itotmst
  end if
  call timer_end(LOG_GSCG_INIT)


  call timer_begin(LOG_GSCG_INIT_ITERATION)
  do iob=1,system%nspin*info%numo
    call calc_allob(iob,iob_allob,itotmst,mst,system%nspin*info%numo)
    if(ilsda==0.or.ilsda==1.and.iob<=info%numo)then
      is=1
    else
      is=2
    end if

  !$OMP parallel do private(iz,iy,ix) collapse(2)
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      cg%rxk_ob(ix,iy,iz,iob)=spsi%rwf(ix,iy,iz,is,iob-(is-1)*info%numo,1,1)
      stpsi%rwf(ix,iy,iz,1,1,1,1)=cg%rxk_ob(ix,iy,iz,iob)
    end do
    end do
    end do

    if(iob_allob<=mst(1))then
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

    call hpsi(stpsi,shtpsi,info_ob,mg,v,nspin_1,stencil,srg_ob_1,ppg)
    
  !$OMP parallel do private(iz,iy,ix) collapse(2)
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      cg%rhxk_ob(ix,iy,iz,iob)=shtpsi%rwf(ix,iy,iz,1,1,1,1)
    end do
    end do
    end do
  end do
  
  call inner_product7(mg,itotmst,mst,system%nspin*info%numo,cg%rxk_ob,cg%rhxk_ob,xkhxk_ob,system%hvol)
  
  xkxk_ob(:)=1.d0 
  rk_ob(:)=xkhxk_ob(:)/xkxk_ob(:)
  call timer_end(LOG_GSCG_INIT_ITERATION)

  
  call timer_begin(LOG_GSCG_ITERATION)
  Iteration : do iter=1,ncg
    do iob=1,system%nspin*info%numo
      call calc_allob(iob,iob_allob,itotmst,mst,system%nspin*info%numo)
  !$OMP parallel do private(iz,iy,ix) collapse(2)
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        cg%rgk_ob(ix,iy,iz,iob) = 2*( cg%rhxk_ob(ix,iy,iz,iob) - rk_ob(iob_allob)*cg%rxk_ob(ix,iy,iz,iob) )
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
            sum0=sum0+spsi%rwf(ix,iy,iz,is,job-(is-1)*info%numo,1,1)*cg%rgk_ob(ix,iy,iz,iob)
          end do
          end do
          end do
          sum_obmat0(iob,job)=sum0*system%hvol
        end do
      end do 
      end do
      call comm_summation(sum_obmat0,sum_obmat1,itotmst*itotmst,nproc_group_global)
      do is=is_sta,is_end
      do iob=iobsta(is),iobend(is)
        do job=iobsta(is),iob-1
    !$omp parallel do private(iz,iy,ix) collapse(2)
          do iz=mg%is(3),mg%ie(3)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
            cg%rgk_ob(ix,iy,iz,iob)=cg%rgk_ob(ix,iy,iz,iob)-sum_obmat1(iob,job)*spsi%rwf(ix,iy,iz,is,job-(is-1)*info%numo,1,1)
          end do
          end do
          end do
        end do
      end do
      end do
    else
      do is=is_sta,is_end
      do iob=iobsta(is),iobend(is)
        call calc_myob(iob,iob_myob,ilsda,nproc_ob,itotmst,mst,system%nspin*info%numo)
        call check_corrkob(iob,1,icorr_iob,ilsda,nproc_ob,info%ik_s,info%ik_e,mst)
        do job=iobsta(is),iob-1
          call calc_myob(job,job_myob,ilsda,nproc_ob,itotmst,mst,system%nspin*info%numo)
          call check_corrkob(job,1,icorr_job,ilsda,nproc_ob,info%ik_s,info%ik_e,mst)
          if(icorr_job==1)then
    !$omp parallel do private(iz,iy,ix) collapse(2)
            do iz=mg%is(3),mg%ie(3)
            do iy=mg%is(2),mg%ie(2)
            do ix=mg%is(1),mg%ie(1)
              rmatbox_m(ix,iy,iz)=spsi%rwf(ix,iy,iz,is,job_myob-(is-1)*info%numo,1,1)
            end do
            end do
            end do
          end if
          call calc_iroot(job,iroot,ilsda,nproc_ob,itotmst,mst)
          call comm_bcast(rmatbox_m,nproc_group_grid,iroot)
          sum0=0.d0
    !$omp parallel do private(iz,iy,ix) collapse(2) reduction(+ : sum0)
          do iz=mg%is(3),mg%ie(3)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
            sum0=sum0+rmatbox_m(ix,iy,iz)*cg%rgk_ob(ix,iy,iz,iob_myob)
          end do
          end do
          end do
          sum0=sum0*system%hvol
          call comm_summation(sum0,sum1,nproc_group_korbital)
          do iz=mg%is(3),mg%ie(3)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
            cg%rgk_ob(ix,iy,iz,iob_myob)=cg%rgk_ob(ix,iy,iz,iob_myob)-sum1*rmatbox_m(ix,iy,iz)
          end do
          end do
          end do
        end do
      end do
      end do
    end if 
    call inner_product7(mg,itotmst,mst,system%nspin*info%numo,cg%rgk_ob,cg%rgk_ob,sum_ob0,system%hvol)
    if ( iter==1 ) then
      do iob=1,system%nspin*info%numo
        call calc_allob(iob,iob_allob,itotmst,mst,system%nspin*info%numo)
  !$OMP parallel do private(iz,iy,ix) collapse(2)
        do iz=mg%is(3),mg%ie(3)
        do iy=mg%is(2),mg%ie(2)
        do ix=mg%is(1),mg%ie(1)
          cg%rpk_ob(ix,iy,iz,iob) = -cg%rgk_ob(ix,iy,iz,iob)
        end do
        end do
        end do
      end do
    else
      do iob=1,system%nspin*info%numo
        call calc_allob(iob,iob_allob,itotmst,mst,system%nspin*info%numo)
        uk=sum_ob0(iob_allob)/gkgk_ob(iob_allob)
  !$OMP parallel do private(iz,iy,ix)
        do iz=mg%is(3),mg%ie(3)
        do iy=mg%is(2),mg%ie(2)
        do ix=mg%is(1),mg%ie(1)
          cg%rpk_ob(ix,iy,iz,iob) = -cg%rgk_ob(ix,iy,iz,iob) + uk*cg%rpk_ob(ix,iy,iz,iob)
        end do
        end do
        end do
      end do
    end if 
    gkgk_ob(:)=sum_ob0(:)
    call inner_product7(mg,itotmst,mst,system%nspin*info%numo,cg%rxk_ob,cg%rpk_ob,xkpk_ob,system%hvol)
    call inner_product7(mg,itotmst,mst,system%nspin*info%numo,cg%rpk_ob,cg%rpk_ob,pkpk_ob,system%hvol)
    call inner_product7(mg,itotmst,mst,system%nspin*info%numo,cg%rpk_ob,cg%rhxk_ob,pkhxk_ob,system%hvol)
  
    do iob=1,system%nspin*info%numo
      call calc_allob(iob,iob_allob,itotmst,mst,system%nspin*info%numo)
  !$OMP parallel do private(iz,iy,ix) collapse(2)
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        stpsi%rwf(ix,iy,iz,1,1,1,1) = cg%rpk_ob(ix,iy,iz,iob)
      end do
      end do
      end do
      
      if(iob_allob<=mst(1))then
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

      call hpsi(stpsi,shtpsi,info_ob,mg,v,nspin_1,stencil,srg_ob_1,ppg)
      
  !$OMP parallel do private(iz,iy,ix) collapse(2)
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
         cg%rgk_ob(ix,iy,iz,iob)=shtpsi%rwf(ix,iy,iz,1,1,1,1)
      end do
      end do
      end do
    end do
    call inner_product7(mg,itotmst,mst,system%nspin*info%numo,cg%rpk_ob,cg%rgk_ob,pkHpk_ob,system%hvol)
    do iob=1,system%nspin*info%numo
      call calc_allob(iob,iob_allob,itotmst,mst,system%nspin*info%numo)
      ak=pkHpk_ob(iob_allob)*xkpk_ob(iob_allob)-pkhxk_ob(iob_allob)*pkpk_ob(iob_allob)
      bk=pkHpk_ob(iob_allob)*xkxk_ob(iob_allob)-xkhxk_ob(iob_allob)*pkpk_ob(iob_allob)
      ck=pkhxk_ob(iob_allob)*xkxk_ob(iob_allob)-xkhxk_ob(iob_allob)*xkpk_ob(iob_allob)
      alpha=(-bk+sqrt(bk*bk-4*ak*ck))/(2*ak)
  
  !$OMP parallel do private(iz,iy,ix) collapse(2)
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        cg%rxk_ob(ix,iy,iz,iob) = cg%rxk_ob(ix,iy,iz,iob) + alpha*cg%rpk_ob(ix,iy,iz,iob)
        cg%rhxk_ob(ix,iy,iz,iob) = cg%rhxk_ob(ix,iy,iz,iob) + alpha*cg%rgk_ob(ix,iy,iz,iob)
      end do
      end do
      end do
    end do
    call inner_product7(mg,itotmst,mst,system%nspin*info%numo,cg%rxk_ob,cg%rhxk_ob,xkhxk_ob,system%hvol)
    call inner_product7(mg,itotmst,mst,system%nspin*info%numo,cg%rxk_ob,cg%rxk_ob,xkxk_ob,system%hvol)
    rk_ob(:)=xkhxk_ob(:)/xkxk_ob(:)
  end do Iteration
  call timer_end(LOG_GSCG_ITERATION)


  call timer_begin(LOG_GSCG_DEINIT)
  call inner_product7(mg,itotmst,mst,system%nspin*info%numo,cg%rxk_ob,cg%rxk_ob,sum_ob0,system%hvol)
  do iob=1,system%nspin*info%numo
    call calc_allob(iob,iob_allob,itotmst,mst,system%nspin*info%numo)
    if(ilsda==0.or.ilsda==1.and.iob<=info%numo)then
      is=1
    else
      is=2
    end if
  !$OMP parallel do private(iz,iy,ix) collapse(2)
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      spsi%rwf(ix,iy,iz,is,iob-(is-1)*info%numo,1,1)=cg%rxk_ob(ix,iy,iz,iob)/sqrt(sum_ob0(iob_allob))
    end do
    end do
    end do
  end do
  
  if(iflag.eq.1) then
    iflag=0
  end if
  
  deallocate (gk)

  deallocate(stpsi%rwf,shtpsi%rwf)
  deallocate(v(nspin_1)%f)
  deallocate(v)
  call timer_end(LOG_GSCG_DEINIT)

  call timer_end(LOG_GSCG_TOTAL)

  return
  
end subroutine sgscg

end module gscg_sub
