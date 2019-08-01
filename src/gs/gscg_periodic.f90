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
module gscg_periodic_sub
  implicit none

contains

!=======================================================================
!======================================= Conjugate-Gradient minimization

subroutine gscg_periodic(mg,system,info,stencil,srg_ob_1,spsi,iflag,itotmst,mst,ilsda,nproc_ob,cg,  &
                         info_ob,ppg,vlocal)
  use inputoutput, only: ncg
  use structures, only: s_rgrid,s_dft_system,s_orbital_parallel,s_orbital,s_stencil,s_scalar,s_pp_grid,s_cg
  use salmon_parallel, only: nproc_group_kgrid, nproc_group_korbital, nproc_id_korbital, nproc_group_k
  use salmon_communication, only: comm_bcast, comm_summation
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
  type(s_cg),intent(inout)       :: cg
  type(s_orbital_parallel)       :: info_ob
  type(s_scalar),intent(in) :: vlocal(system%nspin)
  integer,parameter :: nd=4
  integer :: j
  integer :: iter,iob,job
  integer :: ik
  integer :: ix,iy,iz
  integer :: is,iobsta(2),iobend(2)
  integer :: nspin_1
  type(s_orbital)  :: stpsi
  type(s_orbital)  :: shtpsi
  type(s_scalar),allocatable :: v(:)
  complex(8) :: sum0,sum1
  complex(8) :: sum_ob1(itotmst)
  complex(8) :: sum_obmat0(itotmst,itotmst),sum_obmat1(itotmst,itotmst)
  complex(8) :: xkxk_ob(itotmst),xkHxk_ob(itotmst),xkHpk_ob(itotmst),pkHpk_ob(itotmst),gkgk_ob(itotmst)
  complex(8) :: uk
  real(8) :: ev
  complex(8) :: cx,cp
  complex(8) :: zs_ob(itotmst)
  complex(8):: zmatbox_m(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
  integer :: iob_myob,job_myob
  integer :: iob_allob
  integer :: icorr_iob,icorr_job
  integer :: iroot
  integer :: is_sta,is_end
  integer :: iter_bak_ob(itotmst)
  complex(8),parameter :: zi=(0.d0,1.d0)

  call timer_begin(LOG_GSCG_TOTAL)

  call timer_begin(LOG_GSCG_INIT)
  allocate(stpsi%zwf(mg%is_array(1):mg%ie_array(1),  &
                     mg%is_array(2):mg%ie_array(2),  &
                     mg%is_array(3):mg%ie_array(3),1,1,1,1))
  allocate(shtpsi%zwf(mg%is_array(1):mg%ie_array(1),  &
                      mg%is_array(2):mg%ie_array(2),  &
                      mg%is_array(3):mg%ie_array(3),1,1,1,1))

  allocate(stencil%vec_kAc(3,1:1))

  nspin_1=1
  allocate(v(nspin_1))
  allocate(v(nspin_1)%f(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))
  

  call set_isstaend(is_sta,is_end,ilsda)
  
  !$OMP parallel do private(iz,iy,ix) collapse(2)
  do iz=mg%is_array(3),mg%ie_array(3)
  do iy=mg%is_array(2),mg%ie_array(2)
  do ix=mg%is_array(1),mg%ie_array(1)
    stpsi%zwf(ix,iy,iz,1,1,1,1)=0.d0
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

  
  do ik=info%ik_s,info%ik_e

    call timer_begin(LOG_GSCG_INIT_ITERATION)
    do j=1,3
      stencil%vec_kAc(j,1) = system%vec_k(j,ik)
    end do
    call update_kvector_nonlocalpt(ppg,stencil%vec_kAc,1,1)

    iter_bak_ob(:)=0 
  
    do iob_myob=1,system%nspin*info%numo
      call calc_allob(iob_myob,iob_allob,itotmst,mst,system%nspin*info%numo)
      if(ilsda==0.or.ilsda==1.and.iob_myob<=info%numo)then
        is=1
      else
        is=2
      end if
    
    !$omp parallel do private(iz,iy,ix) collapse(2) 
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        cg%zxk_ob(ix,iy,iz,iob_myob)=spsi%zwf(ix,iy,iz,is,iob_myob-(is-1)*info%numo,ik,1)
        stpsi%zwf(ix,iy,iz,1,1,1,1)=cg%zxk_ob(ix,iy,iz,iob_myob)
      end do
      end do
      end do
     
      if(iob_allob<=mst(1))then
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

      call hpsi(stpsi,shtpsi,info_ob,mg,v,nspin_1,stencil,srg_ob_1,ppg)
      
    !$omp parallel do private(iz,iy,ix) collapse(2) 
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        cg%zhxk_ob(ix,iy,iz,iob_myob)=shtpsi%zwf(ix,iy,iz,1,1,1,1)
      end do
      end do
      end do
  
    end do
    call inner_product5(mg,itotmst,mst,system%nspin*info%numo,cg%zxk_ob,cg%zhxk_ob,xkHxk_ob,system%hvol)
    call timer_end(LOG_GSCG_INIT_ITERATION)


    call timer_begin(LOG_GSCG_ITERATION)
    Iteration : do iter=1,Ncg
      do iob_myob=1,system%nspin*info%numo
        call calc_allob(iob_myob,iob_allob,itotmst,mst,system%nspin*info%numo)
    
    !$OMP parallel do private(iz,iy,ix) collapse(2)
        do iz=mg%is(3),mg%ie(3)
        do iy=mg%is(2),mg%ie(2)
        do ix=mg%is(1),mg%ie(1)
          cg%zgk_ob(ix,iy,iz,iob_myob) = cg%zhxk_ob(ix,iy,iz,iob_myob) - xkHxk_ob(iob_allob)*cg%zxk_ob(ix,iy,iz,iob_myob) 
        end do
        end do
        end do
      end do

      if(nproc_ob==1)then
        sum_obmat0(:,:)=0.d0

        do is=is_sta,is_end
        do iob=iobsta(is),iobend(is)
          do job=iobsta(is),iob-1
            sum0=0.d0
    !$omp parallel do private(iz,iy,ix) collapse(2) reduction(+ : sum0)
            do iz=mg%is(3),mg%ie(3)
            do iy=mg%is(2),mg%ie(2)
            do ix=mg%is(1),mg%ie(1)
              sum0=sum0+conjg(spsi%zwf(ix,iy,iz,is,job-(is-1)*info%numo,ik,1))*cg%zgk_ob(ix,iy,iz,iob)
            end do
            end do
            end do
            sum_obmat0(iob,job)=sum0*system%hvol
          end do
        end do
        end do


        call timer_begin(LOG_GSCG_ALLREDUCE)
        call comm_summation(sum_obmat0,sum_obmat1,itotmst*itotmst,nproc_group_k)
        call timer_end(LOG_GSCG_ALLREDUCE)


        do is=is_sta,is_end
        do iob=iobsta(is),iobend(is)
          do job=iobsta(is),iob-1
    !$omp parallel do collapse(2)
            do iz=mg%is(3),mg%ie(3)
            do iy=mg%is(2),mg%ie(2)
            do ix=mg%is(1),mg%ie(1)
              cg%zgk_ob(ix,iy,iz,iob)=cg%zgk_ob(ix,iy,iz,iob)-sum_obmat1(iob,job)*spsi%zwf(ix,iy,iz,is,job-(is-1)*info%numo,ik,1)
            end do
            end do
            end do
          end do
        end do
        end do
      else
        do is=is_sta,is_end
        do iob=iobsta(is),iobend(is)
          call calc_myob(iob,iob_myob,ilsda,nproc_ob,itotmst,mst)
          call check_corrkob(iob,ik,icorr_iob,ilsda,nproc_ob,info%ik_s,info%ik_e,mst)
          do job=iobsta(is),iob-1
            call calc_myob(job,job_myob,ilsda,nproc_ob,itotmst,mst)
            call check_corrkob(job,ik,icorr_job,ilsda,nproc_ob,info%ik_s,info%ik_e,mst)
            if(icorr_job==1)then
    !$omp parallel do private(iz,iy,ix) collapse(2)
              do iz=mg%is(3),mg%ie(3)
              do iy=mg%is(2),mg%ie(2)
              do ix=mg%is(1),mg%ie(1)
                zmatbox_m(ix,iy,iz)=spsi%zwf(ix,iy,iz,is,job_myob-(is-1)*info%numo,ik,1)
              end do
              end do
              end do
            end if
            call calc_iroot(job,iroot,ilsda,nproc_ob,itotmst,mst)
            call comm_bcast(zmatbox_m,nproc_group_kgrid,iroot)
            sum0=0.d0
    !$omp parallel do private(iz,iy,ix) collapse(2) reduction(+ : sum0)
            do iz=mg%is(3),mg%ie(3)
            do iy=mg%is(2),mg%ie(2)
            do ix=mg%is(1),mg%ie(1)
              sum0=sum0+conjg(zmatbox_m(ix,iy,iz))*cg%zgk_ob(ix,iy,iz,iob_myob)
            end do
            end do
            end do
            sum0=sum0*system%hvol
            call comm_summation(sum0,sum1,nproc_group_korbital)
            do iz=mg%is(3),mg%ie(3)
            do iy=mg%is(2),mg%ie(2)
            do ix=mg%is(1),mg%ie(1)
              cg%zgk_ob(ix,iy,iz,iob_myob)=cg%zgk_ob(ix,iy,iz,iob_myob)-sum1*zmatbox_m(ix,iy,iz)
            end do
            end do
            end do
          end do
        end do
        end do
      end if
      call inner_product5(mg,itotmst,mst,system%nspin*info%numo,cg%zgk_ob,cg%zgk_ob,sum_ob1,system%hvol)
        
      do iob_myob=1,system%nspin*info%numo
        call calc_allob(iob_myob,iob_allob,itotmst,mst,system%nspin*info%numo)
    
        if(iter==1)then
    !$OMP parallel do private(iz,iy,ix) collapse(2)
          do iz=mg%is(3),mg%ie(3)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
            cg%zpk_ob(ix,iy,iz,iob_myob)=cg%zgk_ob(ix,iy,iz,iob_myob)
          end do
          end do
          end do
        else
          uk=sum_ob1(iob_allob)/gkgk_ob(iob_allob)
    !$OMP parallel do private(iz,iy,ix) collapse(2)
          do iz=mg%is(3),mg%ie(3)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
            cg%zpk_ob(ix,iy,iz,iob_myob)=cg%zgk_ob(ix,iy,iz,iob_myob)+uk*cg%zpk_ob(ix,iy,iz,iob_myob)
          end do
          end do
          end do
        end if
        gkgk_ob(iob_allob)=sum_ob1(iob_allob)
      end do

      call inner_product5(mg,itotmst,mst,system%nspin*info%numo,cg%zxk_ob,cg%zpk_ob,zs_ob,system%hvol)

      do iob_myob=1,system%nspin*info%numo
        call calc_allob(iob_myob,iob_allob,itotmst,mst,system%nspin*info%numo)
    !$OMP parallel do private(iz,iy,ix) collapse(2)
        do iz=mg%is(3),mg%ie(3)
        do iy=mg%is(2),mg%ie(2)
        do ix=mg%is(1),mg%ie(1)
          cg%zpko_ob(ix,iy,iz,iob_myob)=cg%zpk_ob(ix,iy,iz,iob_myob)-zs_ob(iob_allob)*cg%zxk_ob(ix,iy,iz,iob_myob)
        end do
        end do
        end do
      end do
      call inner_product5(mg,itotmst,mst,system%nspin*info%numo,cg%zpko_ob,cg%zpko_ob,sum_ob1,system%hvol)

      do iob_myob=1,system%nspin*info%numo
        call calc_allob(iob_myob,iob_allob,itotmst,mst,system%nspin*info%numo)
    !$OMP parallel do private(iz,iy,ix) collapse(2)
        do iz=mg%is(3),mg%ie(3)
        do iy=mg%is(2),mg%ie(2)
        do ix=mg%is(1),mg%ie(1)
          cg%zpko_ob(ix,iy,iz,iob_myob)=cg%zpko_ob(ix,iy,iz,iob_myob)/sqrt(sum_ob1(iob_allob))
          stpsi%zwf(ix,iy,iz,1,1,1,1)=cg%zpko_ob(ix,iy,iz,iob_myob)
        end do
        end do
        end do

        if(iob_allob<=mst(1))then
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

        call hpsi(stpsi,shtpsi,info_ob,mg,v,nspin_1,stencil,srg_ob_1,ppg)

    !$OMP parallel do private(iz,iy,ix) collapse(2)
        do iz=mg%is(3),mg%ie(3)
        do iy=mg%is(2),mg%ie(2)
        do ix=mg%is(1),mg%ie(1)
          cg%zhwf_ob(ix,iy,iz,iob_myob)=shtpsi%zwf(ix,iy,iz,1,1,1,1)
        end do
        end do
        end do
      end do
      call inner_product5(mg,itotmst,mst,system%nspin*info%numo,cg%zxk_ob,cg%zhwf_ob,xkHpk_ob,system%hvol)
      call inner_product5(mg,itotmst,mst,system%nspin*info%numo,cg%zpko_ob,cg%zhwf_ob,pkHpk_ob,system%hvol)
        
    
      do iob_myob=1,system%nspin*info%numo
        call calc_allob(iob_myob,iob_allob,itotmst,mst,system%nspin*info%numo)

        ev=0.5d0*((xkHxk_ob(iob_allob)+pkHpk_ob(iob_allob))   &
                 -sqrt((xkHxk_ob(iob_allob)-pkHpk_ob(iob_allob))**2+4.d0*abs(xkHpk_ob(iob_allob))**2))
        cx=xkHpk_ob(iob_allob)/(ev-xkHxk_ob(iob_allob))
        cp=1.d0/sqrt(1.d0+abs(cx)**2)
        cx=cx*cp
        
    !$OMP parallel do private(iz,iy,ix) collapse(2)
        do iz=mg%is(3),mg%ie(3)
        do iy=mg%is(2),mg%ie(2)
        do ix=mg%is(1),mg%ie(1)
          cg%zxk_ob(ix,iy,iz,iob_myob)=cx*cg%zxk_ob(ix,iy,iz,iob_myob)+cp*cg%zpko_ob(ix,iy,iz,iob_myob)
          cg%zhxk_ob(ix,iy,iz,iob_myob)=cx*cg%zhxk_ob(ix,iy,iz,iob_myob)+cp*cg%zhwf_ob(ix,iy,iz,iob_myob)
        end do
        end do
        end do
      end do 
    
      call inner_product5(mg,itotmst,mst,system%nspin*info%numo,cg%zxk_ob,cg%zhxk_ob,xkHxk_ob,system%hvol)
      call inner_product5(mg,itotmst,mst,system%nspin*info%numo,cg%zxk_ob,cg%zxk_ob,xkxk_ob,system%hvol)

      do iob_myob=1,system%nspin*info%numo
        call calc_allob(iob_myob,iob_allob,itotmst,mst,system%nspin*info%numo)
        if(ilsda==0.or.ilsda==1.and.iob_myob<=info%numo)then
          is=1
        else
          is=2
        end if

        if(abs(xkxk_ob(iob_allob))<=1.d30)then
    !$OMP parallel do private(iz,iy,ix) collapse(2)
          do iz=mg%is(3),mg%ie(3)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
            spsi%zwf(ix,iy,iz,is,iob_myob-(is-1)*info%numo,ik,1)=cg%zxk_ob(ix,iy,iz,iob_myob)/sqrt(xkxk_ob(iob_allob))
          end do
          end do
          end do
          iter_bak_ob(iob_allob)=iter
        else
          if(nproc_id_korbital==0)then
            write(*,'("CG termination. orbital ",i6," iter ",i3)') iob_allob,iter_bak_ob(iob_allob)+1
          end if
        end if
    
      end do
    end do Iteration  
    call timer_end(LOG_GSCG_ITERATION)

  end do
  
  call timer_begin(LOG_GSCG_DEINIT)
  if(iflag.eq.1) then
    iflag=0
  end if
  
  deallocate(stpsi%zwf,shtpsi%zwf)
  deallocate(stencil%vec_kAc)
  deallocate(v(nspin_1)%f)
  deallocate(v)
  if(allocated(ppg%zekr_uV)) deallocate(ppg%zekr_uV)
  call timer_end(LOG_GSCG_DEINIT)

  call timer_end(LOG_GSCG_TOTAL)

  return
  
end subroutine gscg_periodic

end module gscg_periodic_sub
