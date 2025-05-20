!
!  Copyright 2019-2020 SALMON developers
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
module Conjugate_Gradient
  implicit none

contains

subroutine gscg_rwf(ncg,mg,system,info,stencil,ppg,vlocal,srg,spsi,cg)
  use structures
  use timer
  use hamiltonian, only: hpsi
  use communication, only: comm_summation
  use salmon_global, only: yn_spinorbit,yn_preconditioning
  use Conjugate_Gradient_so, only: gscg_rwf_so
  !$ use omp_lib
  implicit none
  integer           ,intent(in) :: ncg
  type(s_rgrid)     ,intent(in) :: mg
  type(s_dft_system),intent(in) :: system
  type(s_parallel_info),intent(in) :: info
  type(s_stencil)   ,intent(in) :: stencil
  type(s_pp_grid)   ,intent(in) :: ppg
  type(s_scalar)    ,intent(in) :: vlocal(system%nspin)
  type(s_orbital)               :: spsi
  type(s_sendrecv_grid)         :: srg
  type(s_cg)                    :: cg
  !
  integer,parameter :: nd=4
  integer :: nspin,io,ispin,io_s,io_e,is(3),ie(3),iy,iz
  integer :: ierr,icg
  real(8),parameter :: ep0=0.0d0
  real(8),parameter :: ep1=1.0d-15
  real(8),parameter :: c1  = 2.0d0
  real(8) :: rwork(9),W(2),c
  real(8),allocatable :: rb(:,:),E(:,:),E1(:,:),gkgk(:,:),bk(:,:),res(:,:)
  real(8),allocatable :: utmp3(:,:,:),wtmp2(:,:,:)
  real(8) :: utmp2(2,2),btmp2(2,2)

  if ( yn_spinorbit=='y' ) then
    call gscg_rwf_so(ncg,mg,system,info,stencil,ppg,vlocal,srg,spsi,cg)
    return
  end if

  if(info%im_s/=1 .or. info%im_e/=1) stop "error: im/=1 @ gscg"

  call timer_begin(LOG_GSCG_ISOLATED_CALC)
  nspin = system%nspin
  is = mg%is
  ie = mg%ie
  io_s = info%io_s
  io_e = info%io_e

  if(.not. allocated(cg%hxk%rwf)) then
    call allocate_orbital_real(nspin,mg,info,cg%hxk)
    call allocate_orbital_real(nspin,mg,info,cg%pk)
    call allocate_orbital_real(nspin,mg,info,cg%gk)
    call allocate_orbital_real(nspin,mg,info,cg%pre_gk)
    call allocate_orbital_real(nspin,mg,info,cg%hpk)
    !$acc enter data copyin(cg)
  end if
  
  allocate(rb(system%nspin,system%no))
  allocate(E (system%nspin,system%no))
  allocate(E1(system%nspin,system%no))
  allocate(gkgk(system%nspin,system%no))
  allocate(bk(system%nspin,system%no))
  allocate(wtmp2(6,system%nspin,system%no))
  allocate(utmp3(2,system%nspin,system%no))
  allocate(res(system%nspin,system%no))
  res = 0.0d0
  
  call timer_end(LOG_GSCG_ISOLATED_CALC)

  call timer_begin(LOG_GSCG_ISOLATED_HPSI)
  call hpsi(spsi,cg%hxk,info,mg,vlocal,system,stencil,srg,ppg)
  call timer_end(LOG_GSCG_ISOLATED_HPSI)

  call timer_begin(LOG_GSCG_ISOLATED_CALC)

  E1=1.0d10

  call inner_product(mg,system,info,spsi,cg%hxk,E)

#ifdef USE_OPENACC
!$acc parallel loop private(io,ispin,iz,iy) collapse(4)
#else
!$omp parallel do private(io,ispin,iz,iy) collapse(4)
#endif
  do io=io_s,io_e
  do ispin=1,nspin
  do iz=is(3),ie(3)
  do iy=is(2),ie(2)
    cg%gk%rwf(is(1):ie(1),iy,iz,ispin,io,1,1) = -2.0d0*( cg%hxk%rwf(is(1):ie(1),iy,iz,ispin,io,1,1) &
    & - E(ispin,io)*spsi%rwf(is(1):ie(1),iy,iz,ispin,io,1,1) )
  end do
  end do
  end do
  end do
  call inner_product(mg,system,info,cg%gk,cg%gk,rb)

  do icg=1,Ncg+1

#ifdef USE_OPENACC
!$acc parallel loop private(io,ispin,iz,iy) collapse(4)
#else
!$omp parallel do private(io,ispin,iz,iy) collapse(4)
#endif
    do io=io_s,io_e
    do ispin=1,nspin
    do iz=is(3),ie(3)
    do iy=is(2),ie(2)
      cg%pre_gk%rwf(is(1):ie(1),iy,iz,ispin,io,1,1) = cg%gk%rwf(is(1):ie(1),iy,iz,ispin,io,1,1) ! pre_gk==Pgk
    end do
    end do
    end do
    end do

    res = rb /c1**2

! --- Convergence check ---

    if ( all(rb < ep0) ) exit
    if ( all(abs(E-E1)<ep1) ) exit
    if ( icg==Ncg+1 ) exit

! --- Preconditioning ---

    if(yn_preconditioning=='y')then
      call preconditioning_rgk(mg,system,info,cg%gk,cg%pre_gk)
    end if

! --- orthogonalization
    !call gram_schmidt

! ---

    call inner_product(mg,system,info,cg%pre_gk,cg%gk,rb) ! pre_gk==Pgk

    if ( icg==1 ) then
#ifdef USE_OPENACC
!$acc parallel loop private(io,ispin,iz,iy) collapse(4)
#else
!$omp parallel do private(io,ispin,iz,iy) collapse(4)
#endif
      do io=io_s,io_e
      do ispin=1,nspin
      do iz=is(3),ie(3)
      do iy=is(2),ie(2)
        cg%pk%rwf(is(1):ie(1),iy,iz,ispin,io,1,1) = cg%pre_gk%rwf(is(1):ie(1),iy,iz,ispin,io,1,1) ! pre_gk==Pgk
      end do
      end do
      end do
      end do        
    else
      bk = rb/gkgk
#ifdef USE_OPENACC
!$acc parallel loop private(io,ispin,iz,iy) collapse(4)
#else
!$omp parallel do private(io,ispin,iz,iy) collapse(4)
#endif
      do io=io_s,io_e
      do ispin=1,nspin
      do iz=is(3),ie(3)
      do iy=is(2),ie(2)
        cg%pk%rwf(is(1):ie(1),iy,iz,ispin,io,1,1) = cg%pre_gk%rwf(is(1):ie(1),iy,iz,ispin,io,1,1) &
        & + bk(ispin,io)*cg%pk%rwf(is(1):ie(1),iy,iz,ispin,io,1,1)
      end do
      end do
      end do
      end do
    end if
    gkgk = rb
    call timer_end(LOG_GSCG_ISOLATED_CALC)

    call timer_begin(LOG_GSCG_ISOLATED_HPSI)
    call hpsi(cg%pk,cg%hpk,info,mg,vlocal,system,stencil,srg,ppg)
    call timer_end(LOG_GSCG_ISOLATED_HPSI)

    call timer_begin(LOG_GSCG_ISOLATED_CALC)
    call inner_product(mg,system,info,spsi ,spsi ,wtmp2(1,:,:))
    call inner_product(mg,system,info,cg%pk,spsi ,wtmp2(2,:,:))
    call inner_product(mg,system,info,cg%pk,cg%pk,wtmp2(3,:,:))
    call inner_product(mg,system,info,spsi ,cg%hxk,wtmp2(4,:,:))
    call inner_product(mg,system,info,cg%pk,cg%hxk,wtmp2(5,:,:))
    call inner_product(mg,system,info,cg%pk,cg%hpk,wtmp2(6,:,:))

    do io=io_s,io_e
    do ispin=1,nspin
      btmp2(1,1)=wtmp2(1,ispin,io)
      btmp2(2,1)=wtmp2(2,ispin,io)
      btmp2(1,2)=wtmp2(2,ispin,io)
      btmp2(2,2)=wtmp2(3,ispin,io)
      utmp2(1,1)=wtmp2(4,ispin,io)
      utmp2(2,1)=wtmp2(5,ispin,io)
      utmp2(1,2)=wtmp2(5,ispin,io)
      utmp2(2,2)=wtmp2(6,ispin,io)
      call dsygv(1,'V','U',2,utmp2,2,btmp2,2,W,rwork,9,ierr)
      if ( abs(W(1)-E(ispin,io))>1.d-1 .and. abs(W(2)-E(ispin,io))<=1.d-1 ) then
        utmp2(1,1)=utmp2(1,2)
        utmp2(2,1)=utmp2(2,2)
        W(1)=W(2)
      end if
      
      !- Fix the phase -
      c=utmp2(1,1)
      if( c<0.0d0 ) then
        utmp2(1,1)=-utmp2(1,1)
        utmp2(2,1)=-utmp2(2,1)
      end if
      utmp3(1:2,ispin,io) = utmp2(1:2,1)
      E1(ispin,io)=E(ispin,io)
      E(ispin,io) =W(1)

#ifdef USE_OPENACC
!$acc parallel loop private(iz,iy) collapse(2)
#else
!$omp parallel do private(iz,iy) collapse(2)
#endif
      do iz=is(3),ie(3)
      do iy=is(2),ie(2)
        cg%hxk%rwf(is(1):ie(1),iy,iz,ispin,io,1,1) = &
          & utmp2(1,1)* cg%hxk%rwf(is(1):ie(1),iy,iz,ispin,io,1,1) &
          & + utmp2(2,1)* cg%hpk%rwf(is(1):ie(1),iy,iz,ispin,io,1,1)
        cg%gk%rwf(is(1):ie(1),iy,iz,ispin,io,1,1) = -2.0d0*( &
          & cg%hxk%rwf(is(1):ie(1),iy,iz,ispin,io,1,1) &
          & - W(1)*(utmp2(1,1) * spsi%rwf(is(1):ie(1),iy,iz,ispin,io,1,1) &
          &  + utmp2(2,1) * cg%pk%rwf(is(1):ie(1),iy,iz,ispin,io,1,1) ) )
      end do
      end do

    end do ! ispin
    end do ! io

    call inner_product(mg,system,info,cg%gk,cg%gk,rb)

#ifdef USE_OPENACC
!$acc parallel loop private(io,ispin,iz,iy)
#else
!$omp parallel do private(io,ispin,iz,iy)
#endif
    do io=io_s,io_e
    do ispin=1,nspin
      if ( rb(ispin,io)/res(ispin,io)>1.0d8 ) then
        E(ispin,io)=E1(ispin,io)
        cycle
      end if
      do iz=is(3),ie(3)
      do iy=is(2),ie(2)
        spsi%rwf(is(1):ie(1),iy,iz,ispin,io,1,1) = &
        & utmp3(1,ispin,io) * spsi%rwf(is(1):ie(1),iy,iz,ispin,io,1,1) &
        & + utmp3(2,ispin,io) * cg%pk%rwf(is(1):ie(1),iy,iz,ispin,io,1,1)
      end do
      end do
    end do
    end do

  end do ! icg
  
  deallocate( utmp3,wtmp2 )
  deallocate( bk,gkgk,E1,E,rb,res )

  call timer_end(LOG_GSCG_ISOLATED_CALC)

  return
contains

subroutine inner_product(mg,system,info,psi1,psi2,rbox)
  !$ use omp_lib
  implicit none
  type(s_rgrid),intent(in) :: mg
  type(s_dft_system),intent(in) :: system
  type(s_parallel_info),intent(in) :: info
  type(s_orbital),intent(in) :: psi1,psi2
  real(8),intent(out) :: rbox(system%nspin,system%no)
  !
  integer :: io,ispin,nspin
  integer :: ix,iy,iz
  real(8) :: rbox2(system%nspin,system%no)
  real(8) :: sum0
  nspin = system%nspin

  rbox2 = 0.d0
#ifdef USE_OPENACC
!$acc parallel loop collapse(2) private(io,ispin,sum0,iz,iy,ix)
#else
!$OMP parallel do collapse(2) private(io,ispin,sum0,iz,iy,ix)
#endif
  do io=info%io_s,info%io_e
  do ispin=1,nspin
    sum0 = 0d0
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      sum0 = sum0 + psi1%rwf(ix,iy,iz,ispin,io,1,1) * psi2%rwf(ix,iy,iz,ispin,io,1,1)
    end do
    end do
    end do
    rbox2(ispin,io) = sum0 * system%hvol
  end do
  end do
  call timer_end(LOG_GSCG_ISOLATED_CALC)

  call timer_begin(LOG_GSCG_ISOLATED_COMM_COLL)
  call comm_summation(rbox2,rbox,nspin*system%no,info%icomm_r)
  call timer_end(LOG_GSCG_ISOLATED_COMM_COLL)

  call timer_begin(LOG_GSCG_ISOLATED_CALC)
end subroutine inner_product

subroutine preconditioning_rgk(mg,system,info,gk,pre_gk)
  !$ use omp_lib
  use preconditioning_sub, only: dstencil_preconditioning
  use structures
  use sendrecv_grid, only: update_overlap_real8
  use salmon_global, only: alpha_pre
  implicit none
  type(s_rgrid),intent(in) :: mg
  type(s_dft_system),intent(in) :: system
  type(s_parallel_info),intent(in) :: info
  type(s_orbital),intent(inout) :: gk
  type(s_orbital),intent(inout) :: pre_gk
  !
  integer :: io,ik,ispin,nspin
  integer :: ix,iy,iz
  real(8) :: alpha
  integer :: is(3),ie(3)
  logical :: is_enable_overlapping
  nspin = system%nspin

  alpha = alpha_pre

  if(info%if_divide_rspace) then
    call update_overlap_real8(srg, mg, gk%rwf)
  end if

  do ik=info%ik_s,info%ik_e
  do io=info%io_s,info%io_e
  do ispin=1,nspin
    call dstencil_preconditioning(mg%is_array,mg%ie_array,mg%is,  &
                                  mg%ie,mg%idx,mg%idy,mg%idz,system%hgs, &
                                  gk%rwf(:,:,:,ispin,io,ik,1), &
                                  pre_gk%rwf(:,:,:,ispin,io,ik,1),alpha)
  end do
  end do
  end do

end subroutine preconditioning_rgk

end subroutine gscg_rwf

!===================================================================================================================================

subroutine gscg_zwf(ncg,mg,system,info,stencil,ppg,vlocal,srg,spsi,cg)
  use structures
  use timer
  use hamiltonian, only: hpsi
  use communication, only: comm_summation
  use salmon_global, only: yn_spinorbit,yn_preconditioning
  use Conjugate_Gradient_so, only: gscg_zwf_so
  !$ use omp_lib
  implicit none
  integer           ,intent(in) :: ncg
  type(s_rgrid)     ,intent(in) :: mg
  type(s_dft_system),intent(in) :: system
  type(s_parallel_info),intent(in) :: info
  type(s_stencil)   ,intent(in) :: stencil
  type(s_pp_grid)   ,intent(in) :: ppg
  type(s_scalar)    ,intent(in) :: vlocal(system%nspin)
  type(s_orbital)               :: spsi
  type(s_sendrecv_grid)         :: srg
  type(s_cg)                    :: cg
  !
  integer,parameter :: nd=4
  integer :: nspin,ik,io,ispin,ik_s,ik_e,io_s,io_e,is(3),ie(3),iy,iz
  integer :: ierr,icg
  real(8),parameter :: ep0=0.0d0
  real(8),parameter :: ep1=1.0d-15
  real(8),parameter :: c1  = 2.0d0
  real(8) :: rwork(9),W(2),r,c,d
  real(8),allocatable :: rb(:,:,:),E(:,:,:),E1(:,:,:),gkgk(:,:,:),bk(:,:,:),res(:,:,:)
  complex(8),allocatable :: utmp3(:,:,:,:),wtmp2(:,:,:,:),zb(:,:,:)
  complex(8) :: utmp2(2,2),btmp2(2,2)
  complex(8) :: work(9),zphase,ztmp

  if ( yn_spinorbit=='y' ) then
    call gscg_zwf_so(ncg,mg,system,info,stencil,ppg,vlocal,srg,spsi,cg)
    return
  end if

  if(info%im_s/=1 .or. info%im_e/=1) stop "error: im/=1 @ gscg"

  call timer_begin(LOG_GSCG_PERIODIC_CALC)
  nspin = system%nspin
  is = mg%is
  ie = mg%ie
  ik_s = info%ik_s
  ik_e = info%ik_e
  io_s = info%io_s
  io_e = info%io_e

  if(.not. allocated(cg%hxk%zwf)) then
    call allocate_orbital_complex(nspin,mg,info,cg%hxk)
    call allocate_orbital_complex(nspin,mg,info,cg%pk)
    call allocate_orbital_complex(nspin,mg,info,cg%gk)
    call allocate_orbital_complex(nspin,mg,info,cg%pre_gk)
    call allocate_orbital_complex(nspin,mg,info,cg%hpk)
    !$acc enter data copyin(cg)
  end if
  
  allocate(zb     (system%nspin,system%no,system%nk))
  allocate(rb     (system%nspin,system%no,system%nk))
  allocate(E      (system%nspin,system%no,system%nk))
  allocate(E1     (system%nspin,system%no,system%nk))
  allocate(res    (system%nspin,system%no,system%nk))
  allocate(gkgk   (system%nspin,system%no,system%nk))
  allocate(bk     (system%nspin,system%no,system%nk))
  allocate(wtmp2(6,system%nspin,system%no,system%nk))
  allocate(utmp3(2,system%nspin,system%no,system%nk))
  res = 0.0d0

  E1 = 1.0d10

  call timer_end(LOG_GSCG_PERIODIC_CALC)

  call timer_begin(LOG_GSCG_PERIODIC_HPSI)
  call hpsi(spsi,cg%hxk,info,mg,vlocal,system,stencil,srg,ppg)
  call timer_end(LOG_GSCG_PERIODIC_HPSI)

  call timer_begin(LOG_GSCG_PERIODIC_CALC)
  call inner_product(mg,system,info,spsi,cg%hxk,zb)
  E = dble(zb)
 
#ifdef USE_OPENACC
!$acc parallel loop private(ik,io,ispin,iz,iy) collapse(5)
#else
!$omp parallel do private(ik,io,ispin,iz,iy) collapse(5)
#endif
  do ik=ik_s,ik_e
  do io=io_s,io_e
  do ispin=1,nspin
  do iz=is(3),ie(3)
  do iy=is(2),ie(2)
    cg%gk%zwf(is(1):ie(1),iy,iz,ispin,io,ik,1) = -cg%hxk%zwf(is(1):ie(1),iy,iz,ispin,io,ik,1) &
    & + E(ispin,io,ik) * spsi%zwf(is(1):ie(1),iy,iz,ispin,io,ik,1)
  end do
  end do
  end do
  end do
  end do

  call inner_product(mg,system,info,cg%gk,cg%gk,zb)
  rb = dble(zb)

  do icg=1,Ncg+1

    !Pgk(n1:n2,n)=gk(n1:n2,n)
#ifdef USE_OPENACC
!$acc parallel loop private(ik,io,ispin,iz,iy) collapse(5)
#else
!$omp parallel do private(ik,io,ispin,iz,iy) collapse(5)
#endif
    do ik=ik_s,ik_e
    do io=io_s,io_e
    do ispin=1,nspin
    do iz=is(3),ie(3)
    do iy=is(2),ie(2)
      cg%pre_gk%zwf(is(1):ie(1),iy,iz,ispin,io,ik,1) = cg%gk%zwf(is(1):ie(1),iy,iz,ispin,io,ik,1)
    end do
    end do
    end do
    end do
    end do
    
    res = rb

! --- Convergence check ---

    if ( all(rb < ep0) ) exit
    if ( all(abs(E-E1)<ep1) ) exit
    if ( icg==Ncg+1 ) exit

! --- Preconditioning ---

    if(yn_preconditioning=='y')then
      call preconditioning_zgk(mg,system,info,cg%gk,cg%pre_gk)
    end if

! --- orthogonalization

! ---

    call inner_product(mg,system,info,cg%pre_gk,cg%gk,zb)
    rb = dble(zb)

    if ( icg==1 ) then

#ifdef USE_OPENACC
!$acc parallel loop private(ik,io,ispin,iz,iy) collapse(5)
#else
!$omp parallel do private(ik,io,ispin,iz,iy) collapse(5)
#endif
      do ik=ik_s,ik_e
      do io=io_s,io_e
      do ispin=1,nspin
      do iz=is(3),ie(3)
      do iy=is(2),ie(2)
        cg%pk%zwf(is(1):ie(1),iy,iz,ispin,io,ik,1) = cg%pre_gk%zwf(is(1):ie(1),iy,iz,ispin,io,ik,1)
      end do
      end do
      end do
      end do
      end do

    else

      bk = rb/gkgk

#ifdef USE_OPENACC
!$acc parallel loop private(ik,io,ispin,iz,iy) collapse(5)
#else
!$omp parallel do private(ik,io,ispin,iz,iy) collapse(5)
#endif
      do ik=ik_s,ik_e
      do io=io_s,io_e
      do ispin=1,nspin
      do iz=is(3),ie(3)
      do iy=is(2),ie(2)
        cg%pk%zwf(is(1):ie(1),iy,iz,ispin,io,ik,1) = &
        & cg%pre_gk%zwf(is(1):ie(1),iy,iz,ispin,io,ik,1) + &
        & bk(ispin,io,ik) * cg%pk%zwf(is(1):ie(1),iy,iz,ispin,io,ik,1)
      end do
      end do
      end do
      end do
      end do
      
    end if

    gkgk = rb
    call timer_end(LOG_GSCG_PERIODIC_CALC)

    call timer_begin(LOG_GSCG_PERIODIC_HPSI)
    call hpsi(cg%pk,cg%hpk,info,mg,vlocal,system,stencil,srg,ppg)
    call timer_end(LOG_GSCG_PERIODIC_HPSI)

    call timer_begin(LOG_GSCG_PERIODIC_CALC)

    call inner_product(mg,system,info,spsi,  spsi, wtmp2(1,:,:,:))
    call inner_product(mg,system,info,cg%pk, spsi, wtmp2(2,:,:,:))
    call inner_product(mg,system,info,cg%pk,cg%pk, wtmp2(3,:,:,:))
    call inner_product(mg,system,info,spsi, cg%hxk,wtmp2(4,:,:,:))
    call inner_product(mg,system,info,cg%pk,cg%hxk,wtmp2(5,:,:,:))
    call inner_product(mg,system,info,cg%pk,cg%hpk,wtmp2(6,:,:,:))

    do ik=ik_s,ik_e
    do io=io_s,io_e
    do ispin=1,nspin
      btmp2(1,1)=wtmp2(1,ispin,io,ik)
      btmp2(2,1)=wtmp2(2,ispin,io,ik)
      btmp2(1,2)=wtmp2(2,ispin,io,ik)
      btmp2(2,2)=wtmp2(3,ispin,io,ik)
      utmp2(1,1)=wtmp2(4,ispin,io,ik)
      utmp2(2,1)=wtmp2(5,ispin,io,ik)
      utmp2(1,2)=wtmp2(5,ispin,io,ik)
      utmp2(2,2)=wtmp2(6,ispin,io,ik)
      ztmp=btmp2(1,2)
      ztmp=conjg(ztmp)
      btmp2(1,2)=ztmp
      ztmp=utmp2(1,2)
      ztmp=conjg(ztmp)
      utmp2(1,2)=ztmp
      call zhegv(1,'V','U',2,utmp2,2,btmp2,2,W,work,9,rwork,ierr)
      if ( abs(W(1)-E(ispin,io,ik))>1.d-1 .and. abs(W(2)-E(ispin,io,ik))<=1.d-1 ) then
        utmp2(1,1)=utmp2(1,2)
        utmp2(2,1)=utmp2(2,2)
        W(1)=W(2)
      end if
      !- Fix the phase -
      ztmp=utmp2(1,1)
      r=abs(ztmp)
      c=real(ztmp)/r
      d=aimag(ztmp)/r
      zphase=dcmplx(c,-d)
      utmp2(1,1)=utmp2(1,1)*zphase
      utmp2(2,1)=utmp2(2,1)*zphase

      utmp3(1:2,ispin,io,ik) = utmp2(1:2,1)

      E1(ispin,io,ik) = E(ispin,io,ik)
      E (ispin,io,ik) = W(1)

#ifdef USE_OPENACC
!$acc parallel loop private(iz,iy) collapse(2)
#else
!$omp parallel do private(iz,iy) collapse(2)
#endif
      do iz=is(3),ie(3)
      do iy=is(2),ie(2)
        cg%hxk%zwf(is(1):ie(1),iy,iz,ispin,io,ik,1) = &
        & utmp2(1,1) * cg%hxk%zwf(is(1):ie(1),iy,iz,ispin,io,ik,1) &
        & + utmp2(2,1) * cg%hpk%zwf(is(1):ie(1),iy,iz,ispin,io,ik,1)
      end do
      end do

#ifdef USE_OPENACC
!$acc parallel loop private(iz,iy) collapse(2)
#else
!$omp parallel do private(iz,iy) collapse(2)
#endif
      do iz=is(3),ie(3)
      do iy=is(2),ie(2)
        cg%gk%zwf(is(1):ie(1),iy,iz,ispin,io,ik,1) = - cg%hxk%zwf(is(1):ie(1),iy,iz,ispin,io,ik,1) &
        & + W(1)*( utmp2(1,1) * spsi%zwf(is(1):ie(1),iy,iz,ispin,io,ik,1) &
        & + utmp2(2,1) * cg%pk%zwf(is(1):ie(1),iy,iz,ispin,io,ik,1) )
      end do
      end do

    end do ! ispin
    end do ! io
    end do ! ik
    
    call inner_product(mg,system,info,cg%gk,cg%gk,zb)
    rb = dble(zb)

#ifdef USE_OPENACC
!$acc parallel loop private(ik,io,ispin,iz,iy) collapse(2)
#else
!$omp parallel do private(ik,io,ispin,iz,iy) collapse(2)
#endif
    do ik=ik_s,ik_e
    do io=io_s,io_e
    do ispin=1,nspin
      if ( rb(ispin,io,ik)/res(ispin,io,ik)>1.0d8 ) then
        E(ispin,io,ik) = E1(ispin,io,ik)
        cycle
      end if
      do iz=is(3),ie(3)
      do iy=is(2),ie(2)
        spsi%zwf(is(1):ie(1),iy,iz,ispin,io,ik,1) = &
        & utmp3(1,ispin,io,ik) * spsi%zwf(is(1):ie(1),iy,iz,ispin,io,ik,1) &
        & + utmp3(2,ispin,io,ik) * cg%pk%zwf(is(1):ie(1),iy,iz,ispin,io,ik,1)
      end do
      end do
    end do
    end do
    end do

  end do ! icg

  deallocate( utmp3,wtmp2 )
  deallocate( bk,gkgk,E1,E,rb,res )
  deallocate( zb )

  call timer_end(LOG_GSCG_PERIODIC_CALC)

  return
contains

subroutine inner_product(mg,system,info,psi1,psi2,zbox)
  !$ use omp_lib
  implicit none
  type(s_rgrid),intent(in) :: mg
  type(s_dft_system),intent(in) :: system
  type(s_parallel_info),intent(in) :: info
  type(s_orbital),intent(in) :: psi1,psi2
  complex(8),intent(out) :: zbox(system%nspin,system%no,system%nk)
  !
  integer :: io,ik,ispin,nspin
  integer :: ix,iy,iz
  complex(8) :: zbox2(system%nspin,system%no,system%nk)
  complex(8) :: sum0
  nspin = system%nspin

  zbox2(:,:,:) = 0.d0
#ifdef USE_OPENACC
!$acc parallel loop collapse(2) private(ik,io,ispin,sum0,iz,iy,ix)
#else
!$OMP parallel do collapse(2) private(ik,io,ispin,sum0,iz,iy,ix)
#endif
  do ik=info%ik_s,info%ik_e
  do io=info%io_s,info%io_e
  do ispin=1,nspin
    sum0 = 0d0
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      sum0 = sum0 + conjg(psi1%zwf(ix,iy,iz,ispin,io,ik,1))*psi2%zwf(ix,iy,iz,ispin,io,ik,1)
    end do
    end do
    end do
    zbox2(ispin,io,ik) = sum0 * system%hvol
  end do
  end do
  end do
  call timer_end(LOG_GSCG_PERIODIC_CALC)

  call timer_begin(LOG_GSCG_PERIODIC_COMM_COLL)
  call comm_summation(zbox2,zbox,nspin*system%no*system%nk,info%icomm_r)
  call timer_end(LOG_GSCG_PERIODIC_COMM_COLL)

  call timer_begin(LOG_GSCG_PERIODIC_CALC)
end subroutine inner_product

subroutine preconditioning_zgk(mg,system,info,gk,pre_gk)
  !$ use omp_lib
  use preconditioning_sub, only: zstencil_preconditioning,zstencil_nonorthogonal_preconditioning
  use structures
  use sendrecv_grid, only: update_overlap_complex8
  use salmon_global, only: yn_want_communication_overlapping,alpha_pre
  implicit none
  type(s_rgrid),intent(in) :: mg
  type(s_dft_system),intent(in) :: system
  type(s_parallel_info),intent(in) :: info
  type(s_orbital),intent(inout) :: gk
  type(s_orbital),intent(inout) :: pre_gk
  !
  integer :: io,ik,ispin,nspin
  real(8) :: alpha
  integer :: is(3),ie(3)
  logical :: is_enable_overlapping
  nspin = system%nspin

  alpha = alpha_pre

  is_enable_overlapping = (yn_want_communication_overlapping == 'y') .and. &
                          stencil%if_orthogonal .and. &
                          info%if_divide_rspace

  if(info%if_divide_rspace .and. .not. is_enable_overlapping) then
    call update_overlap_complex8(srg, mg, gk%zwf)
  end if

  if(stencil%if_orthogonal) then
    do ik=info%ik_s,info%ik_e
    do io=info%io_s,info%io_e
    do ispin=1,nspin
      call zstencil_preconditioning(mg%is_array,mg%ie_array,mg%is,  &
                                    mg%ie,mg%idx,mg%idy,mg%idz,system%hgs, &
                                    gk%zwf(:,:,:,ispin,io,ik,1), &
                                    pre_gk%zwf(:,:,:,ispin,io,ik,1),alpha)
    end do
    end do
    end do
  else
    do ik=info%ik_s,info%ik_e
    do io=info%io_s,info%io_e
    do ispin=1,nspin
      call zstencil_nonorthogonal_preconditioning(mg%is_array,mg%ie_array,mg%is,  &
                                    mg%ie,mg%idx,mg%idy,mg%idz, &
                                    gk%zwf(:,:,:,ispin,io,ik,1), &
                                    pre_gk%zwf(:,:,:,ispin,io,ik,1), &
                                    stencil%coef_lap0_nd1, &
                                    stencil%coef_lap_nd1,stencil%coef_nab_nd1, &
                                    stencil%coef_F,alpha)
    end do
    end do
    end do
  end if

end subroutine preconditioning_zgk

end subroutine gscg_zwf

end module Conjugate_Gradient
