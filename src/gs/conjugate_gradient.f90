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
module Conjugate_Gradient
  implicit none

contains

subroutine gscg_isolated(mg,system,info,stencil,ppg,vlocal,srg,spsi,cg)
  use inputoutput, only: ncg
  use structures
  use timer
  use hpsi_sub
  use salmon_communication, only: comm_summation
  !$ use omp_lib
  implicit none
  type(s_rgrid)     ,intent(in) :: mg
  type(s_dft_system),intent(in) :: system
  type(s_orbital_parallel),intent(in) :: info
  type(s_stencil)   ,intent(in) :: stencil
  type(s_pp_grid)   ,intent(in) :: ppg
  type(s_scalar)    ,intent(in) :: vlocal(system%nspin)
  type(s_orbital)               :: spsi
  type(s_sendrecv_grid)         :: srg
  type(s_cg)                    :: cg
  !
  integer,parameter :: nd=4
  integer :: nspin,io,io_all,ispin,io_s,io_e,is(3),ie(3),ix,iy,iz
  integer :: iter
  real(8),dimension(system%nspin,system%no) :: sum,xkxk,xkHxk,xkHpk,pkHpk,gkgk,uk,ev,cx,cp,zs

  call timer_begin(LOG_GSCG_TOTAL)
  call timer_begin(LOG_GSCG_INIT)

  if(info%im_s/=1 .or. info%im_e/=1) stop "error: im/=1 @ gscg"
  if(info%if_divide_orbit) stop "error: nproc_ob/=1 @ gscg"

  nspin = system%nspin
  is = mg%is
  ie = mg%ie
  io_s = info%io_s
  io_e = info%io_e

  if(.not. allocated(cg%xk%rwf)) then
    call allocate_orbital_real(nspin,mg,info,cg%xk)
    call allocate_orbital_real(nspin,mg,info,cg%hxk)
    call allocate_orbital_real(nspin,mg,info,cg%pk)
    call allocate_orbital_real(nspin,mg,info,cg%gk)
    call allocate_orbital_real(nspin,mg,info,cg%pko)
    call allocate_orbital_real(nspin,mg,info,cg%hwf)
  end if

  !$omp parallel do private(io,ispin,iz,iy,ix) collapse(5)
  do io=io_s,io_e
  do ispin=1,nspin
  do iz=is(3),ie(3)
  do iy=is(2),ie(2)
  do ix=is(1),ie(1)
    cg%xk%rwf(ix,iy,iz,ispin,io,1,1) = spsi%rwf(ix,iy,iz,ispin,io,1,1)
  end do
  end do
  end do
  end do
  end do

  call hpsi(cg%xk,cg%hxk,info,mg,vlocal,system,stencil,srg,ppg)
  call inner_product(mg,system,info,cg%xk,cg%hxk,xkHxk)

  call timer_end(LOG_GSCG_INIT)

  call timer_begin(LOG_GSCG_ITERATION)
  Iteration : do iter=1,Ncg

    !$omp parallel do private(io,ispin,iz,iy,ix,io_all) collapse(5)
    do io=io_s,io_e
    do ispin=1,nspin
    do iz=is(3),ie(3)
    do iy=is(2),ie(2)
    do ix=is(1),ie(1)
      io_all = info%io_tbl(io)
      cg%gk%rwf(ix,iy,iz,ispin,io,1,1) = &
      & cg%hxk%rwf(ix,iy,iz,ispin,io,1,1) - xkHxk(ispin,io_all)* cg%xk%rwf(ix,iy,iz,ispin,io,1,1)
    end do
    end do
    end do
    end do
    end do

    call orthogonalization(mg,system,info,spsi,cg%gk)
    call inner_product(mg,system,info,cg%gk,cg%gk,sum)

    if(iter==1)then
      uk = 0d0
    else
      uk=sum/gkgk
    end if
    !$omp parallel do private(io,ispin,iz,iy,ix,io_all) collapse(5)
    do io=io_s,io_e
    do ispin=1,nspin
    do iz=is(3),ie(3)
    do iy=is(2),ie(2)
    do ix=is(1),ie(1)
      io_all = info%io_tbl(io)
      cg%pk%rwf(ix,iy,iz,ispin,io,1,1) = &
      & cg%gk%rwf(ix,iy,iz,ispin,io,1,1) + uk(ispin,io_all) * cg%pk%rwf(ix,iy,iz,ispin,io,1,1)
    end do
    end do
    end do
    end do
    end do

    gkgk = sum
    call inner_product(mg,system,info,cg%xk,cg%pk,zs)

    !$omp parallel do private(io,ispin,iz,iy,ix,io_all) collapse(5)
    do io=io_s,io_e
    do ispin=1,nspin
    do iz=is(3),ie(3)
    do iy=is(2),ie(2)
    do ix=is(1),ie(1)
      io_all = info%io_tbl(io)
      cg%pko%rwf(ix,iy,iz,ispin,io,1,1) = &
      & cg%pk%rwf(ix,iy,iz,ispin,io,1,1) - zs(ispin,io_all) *cg%xk%rwf(ix,iy,iz,ispin,io,1,1)
    end do
    end do
    end do
    end do
    end do

    call inner_product(mg,system,info,cg%pko,cg%pko,sum)

    !$omp parallel do private(io,ispin,iz,iy,ix,io_all) collapse(5)
    do io=io_s,io_e
    do ispin=1,nspin
    do iz=is(3),ie(3)
    do iy=is(2),ie(2)
    do ix=is(1),ie(1)
      io_all = info%io_tbl(io)
      cg%pko%rwf(ix,iy,iz,ispin,io,1,1) = cg%pko%rwf(ix,iy,iz,ispin,io,1,1) /sqrt(sum(ispin,io_all))
    end do
    end do
    end do
    end do
    end do

    call hpsi(cg%pko,cg%hwf,info,mg,vlocal,system,stencil,srg,ppg)

    call inner_product(mg,system,info,cg%xk,cg%hwf,xkHpk)
    call inner_product(mg,system,info,cg%pko,cg%hwf,pkHpk)

    ev=0.5d0*((xkHxk+pkHpk)   &
             -sqrt((xkHxk-pkHpk)**2+4.d0*abs(xkHpk)**2))
    cx=xkHpk/(ev-xkHxk)
    cp=1.d0/sqrt(1.d0+abs(cx)**2)
    cx=cx*cp

    !$omp parallel do private(io,ispin,iz,iy,ix,io_all) collapse(5)
    do io=io_s,io_e
    do ispin=1,nspin
    do iz=is(3),ie(3)
    do iy=is(2),ie(2)
    do ix=is(1),ie(1)
      io_all = info%io_tbl(io)
      cg%xk%rwf(ix,iy,iz,ispin,io,1,1) = &
      & cx(ispin,io_all)* cg%xk%rwf(ix,iy,iz,ispin,io,1,1) + cp(ispin,io_all)* cg%pko%rwf(ix,iy,iz,ispin,io,1,1)
      cg%hxk%rwf(ix,iy,iz,ispin,io,1,1) = &
      & cx(ispin,io_all)* cg%hxk%rwf(ix,iy,iz,ispin,io,1,1) + cp(ispin,io_all)* cg%hwf%rwf(ix,iy,iz,ispin,io,1,1)
    end do
    end do
    end do
    end do
    end do

    call inner_product(mg,system,info,cg%xk,cg%hxk,xkHxk)
    call inner_product(mg,system,info,cg%xk,cg%xk,xkxk)

    !$omp parallel do private(io,ispin,iz,iy,ix,io_all) collapse(5)
    do io=io_s,io_e
    do ispin=1,nspin
    do iz=is(3),ie(3)
    do iy=is(2),ie(2)
    do ix=is(1),ie(1)
      io_all = info%io_tbl(io)
      if(abs(xkxk(ispin,io_all))<=1.d30)then
        spsi%rwf(ix,iy,iz,ispin,io,1,1) = cg%xk%rwf(ix,iy,iz,ispin,io,1,1) /sqrt(xkxk(ispin,io_all))
      end if
    end do
    end do
    end do
    end do
    end do

  end do Iteration
  call timer_end(LOG_GSCG_ITERATION)

  call timer_end(LOG_GSCG_TOTAL)

  return
contains

subroutine orthogonalization(mg,system,info,psi,gk)
  implicit none
  type(s_rgrid),intent(in) :: mg
  type(s_dft_system),intent(in) :: system
  type(s_orbital_parallel),intent(in) :: info
  type(s_orbital),intent(in) :: psi
  type(s_orbital)            :: gk
  !
  integer :: nspin,no,ispin,io_s,io_e,is(3),ie(3),ix,iy,iz,io1,io2,io1_all,io2_all
  complex(8) :: sum0
  complex(8),dimension(system%no,system%no,system%nspin) :: sum_obmat0,sum_obmat1

  nspin = system%nspin
  no = system%no
  is = mg%is
  ie = mg%ie
  io_s = info%io_s
  io_e = info%io_e

  sum_obmat0 = 0.d0
!$omp parallel do private(ispin,io1,io2,sum0,iz,iy,ix) collapse(2)
  do ispin=1,nspin
    do io1=io_s,io_e
      do io2=io_s,io1-1
        sum0 = 0.d0
! speed? !$omp parallel do private(iz,iy,ix) collapse(2) reduction(+ : sum0)
        do iz=is(3),ie(3)
        do iy=is(2),ie(2)
        do ix=is(1),ie(1)
          sum0=sum0+psi%rwf(ix,iy,iz,ispin,io2,1,1)*gk%rwf(ix,iy,iz,ispin,io1,1,1)
        end do
        end do
        end do
        sum_obmat0(io1,io2,ispin) = sum0*system%hvol
      end do
    end do
  end do
  call timer_begin(LOG_GSCG_ALLREDUCE)
  call comm_summation(sum_obmat0,sum_obmat1,no**2*nspin,info%icomm_rko)
  call timer_end(LOG_GSCG_ALLREDUCE)
!$omp parallel do private(ispin,io1,io2,iz,iy,ix) collapse(2)
  do ispin=1,nspin
    do io1=io_s,io_e
      do io2=io_s,io1-1
! speed? !$omp parallel do private(iz,iy,ix) collapse(2)
        do iz=is(3),ie(3)
        do iy=is(2),ie(2)
        do ix=is(1),ie(1)
          gk%rwf(ix,iy,iz,ispin,io1,1,1) = gk%rwf(ix,iy,iz,ispin,io1,1,1) &
          & -sum_obmat1(io1,io2,ispin) * psi%rwf(ix,iy,iz,ispin,io2,1,1)
        end do
        end do
        end do
      end do
    end do
  end do

end subroutine orthogonalization

subroutine inner_product(mg,system,info,psi1,psi2,rbox)
  !$ use omp_lib
  implicit none
  type(s_rgrid),intent(in) :: mg
  type(s_dft_system),intent(in) :: system
  type(s_orbital_parallel),intent(in) :: info
  type(s_orbital),intent(in) :: psi1,psi2
  real(8),intent(out) :: rbox(system%nspin,system%no)
  !
  integer :: io,io_all,ispin,nspin
  integer :: ix,iy,iz
  real(8) :: rbox2(system%nspin,system%no)
  real(8) :: sum0
  nspin = system%nspin

  rbox2 = 0.d0
!$OMP parallel do collapse(2) private(io,ispin,io_all,sum0,iz,iy,ix)
  do io=info%io_s,info%io_e
  do ispin=1,nspin
    io_all = info%io_tbl(io)
    sum0 = 0d0
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      sum0 = sum0 + psi1%rwf(ix,iy,iz,ispin,io,1,1) * psi2%rwf(ix,iy,iz,ispin,io,1,1)
    end do
    end do
    end do
    rbox2(ispin,io_all) = sum0 * system%hvol
  end do
  end do

  call comm_summation(rbox2,rbox,nspin*system%no,info%icomm_r)

end subroutine inner_product

end subroutine gscg_isolated

!===================================================================================================================================

subroutine gscg_periodic(mg,system,info,stencil,ppg,vlocal,srg,spsi,cg)
  use inputoutput, only: ncg
  use structures
  use timer
  use hpsi_sub
  use salmon_communication, only: comm_summation
  !$ use omp_lib
  implicit none
  type(s_rgrid)     ,intent(in) :: mg
  type(s_dft_system),intent(in) :: system
  type(s_orbital_parallel),intent(in) :: info
  type(s_stencil)   ,intent(in) :: stencil
  type(s_pp_grid)   ,intent(in) :: ppg
  type(s_scalar)    ,intent(in) :: vlocal(system%nspin)
  type(s_orbital)               :: spsi
  type(s_sendrecv_grid)         :: srg
  type(s_cg)                    :: cg
  !
  integer,parameter :: nd=4
  integer :: nspin,ik,io,io_all,ispin,ik_s,ik_e,io_s,io_e,is(3),ie(3),ix,iy,iz
  integer :: iter
  complex(8),dimension(system%nspin,system%no,system%nk) :: sum,xkxk,xkHxk,xkHpk,pkHpk,gkgk,uk,ev,cx,cp,zs

  call timer_begin(LOG_GSCG_TOTAL)
  call timer_begin(LOG_GSCG_INIT)

  if(info%im_s/=1 .or. info%im_e/=1) stop "error: im/=1 @ gscg"
  if(info%if_divide_orbit) stop "error: nproc_ob/=1 @ gscg"

  nspin = system%nspin
  is = mg%is
  ie = mg%ie
  ik_s = info%ik_s
  ik_e = info%ik_e
  io_s = info%io_s
  io_e = info%io_e

  if(.not. allocated(cg%xk%zwf)) then
    call allocate_orbital_complex(nspin,mg,info,cg%xk)
    call allocate_orbital_complex(nspin,mg,info,cg%hxk)
    call allocate_orbital_complex(nspin,mg,info,cg%pk)
    call allocate_orbital_complex(nspin,mg,info,cg%gk)
    call allocate_orbital_complex(nspin,mg,info,cg%pko)
    call allocate_orbital_complex(nspin,mg,info,cg%hwf)
  end if

  !$omp parallel do private(ik,io,ispin,iz,iy,ix) collapse(6)
  do ik=ik_s,ik_e
  do io=io_s,io_e
  do ispin=1,nspin
  do iz=is(3),ie(3)
  do iy=is(2),ie(2)
  do ix=is(1),ie(1)
    cg%xk%zwf(ix,iy,iz,ispin,io,ik,1) = spsi%zwf(ix,iy,iz,ispin,io,ik,1)
  end do
  end do
  end do
  end do
  end do
  end do

  call hpsi(cg%xk,cg%hxk,info,mg,vlocal,system,stencil,srg,ppg)
  call inner_product(mg,system,info,cg%xk,cg%hxk,xkHxk)

  call timer_end(LOG_GSCG_INIT)

  call timer_begin(LOG_GSCG_ITERATION)
  Iteration : do iter=1,Ncg

    !$omp parallel do private(ik,io,ispin,iz,iy,ix,io_all) collapse(6)
    do ik=ik_s,ik_e
    do io=io_s,io_e
    do ispin=1,nspin
    do iz=is(3),ie(3)
    do iy=is(2),ie(2)
    do ix=is(1),ie(1)
      io_all = info%io_tbl(io)
      cg%gk%zwf(ix,iy,iz,ispin,io,ik,1) = &
      & cg%hxk%zwf(ix,iy,iz,ispin,io,ik,1) - xkHxk(ispin,io_all,ik)* cg%xk%zwf(ix,iy,iz,ispin,io,ik,1)
    end do
    end do
    end do
    end do
    end do
    end do

    call orthogonalization(mg,system,info,spsi,cg%gk)
    call inner_product(mg,system,info,cg%gk,cg%gk,sum)

    if(iter==1)then
      uk = 0d0
    else
      uk=sum/gkgk
    end if
    !$omp parallel do private(ik,io,ispin,iz,iy,ix,io_all) collapse(6)
    do ik=ik_s,ik_e
    do io=io_s,io_e
    do ispin=1,nspin
    do iz=is(3),ie(3)
    do iy=is(2),ie(2)
    do ix=is(1),ie(1)
      io_all = info%io_tbl(io)
      cg%pk%zwf(ix,iy,iz,ispin,io,ik,1) = &
      & cg%gk%zwf(ix,iy,iz,ispin,io,ik,1) + uk(ispin,io_all,ik) * cg%pk%zwf(ix,iy,iz,ispin,io,ik,1)
    end do
    end do
    end do
    end do
    end do
    end do

    gkgk = sum
    call inner_product(mg,system,info,cg%xk,cg%pk,zs)

    !$omp parallel do private(ik,io,ispin,iz,iy,ix,io_all) collapse(6)
    do ik=ik_s,ik_e
    do io=io_s,io_e
    do ispin=1,nspin
    do iz=is(3),ie(3)
    do iy=is(2),ie(2)
    do ix=is(1),ie(1)
      io_all = info%io_tbl(io)
      cg%pko%zwf(ix,iy,iz,ispin,io,ik,1) = &
      & cg%pk%zwf(ix,iy,iz,ispin,io,ik,1) - zs(ispin,io_all,ik) *cg%xk%zwf(ix,iy,iz,ispin,io,ik,1)
    end do
    end do
    end do
    end do
    end do
    end do

    call inner_product(mg,system,info,cg%pko,cg%pko,sum)

    !$omp parallel do private(ik,io,ispin,iz,iy,ix,io_all) collapse(6)
    do ik=ik_s,ik_e
    do io=io_s,io_e
    do ispin=1,nspin
    do iz=is(3),ie(3)
    do iy=is(2),ie(2)
    do ix=is(1),ie(1)
      io_all = info%io_tbl(io)
      cg%pko%zwf(ix,iy,iz,ispin,io,ik,1) = cg%pko%zwf(ix,iy,iz,ispin,io,ik,1) /sqrt(sum(ispin,io_all,ik))
    end do
    end do
    end do
    end do
    end do
    end do

    call hpsi(cg%pko,cg%hwf,info,mg,vlocal,system,stencil,srg,ppg)

    call inner_product(mg,system,info,cg%xk,cg%hwf,xkHpk)
    call inner_product(mg,system,info,cg%pko,cg%hwf,pkHpk)

    ev=0.5d0*((xkHxk+pkHpk)   &
             -sqrt((xkHxk-pkHpk)**2+4.d0*abs(xkHpk)**2))
    cx=xkHpk/(ev-xkHxk)
    cp=1.d0/sqrt(1.d0+abs(cx)**2)
    cx=cx*cp

    !$omp parallel do private(ik,io,ispin,iz,iy,ix,io_all) collapse(6)
    do ik=ik_s,ik_e
    do io=io_s,io_e
    do ispin=1,nspin
    do iz=is(3),ie(3)
    do iy=is(2),ie(2)
    do ix=is(1),ie(1)
      io_all = info%io_tbl(io)
      cg%xk%zwf(ix,iy,iz,ispin,io,ik,1) = &
      & cx(ispin,io_all,ik)* cg%xk%zwf(ix,iy,iz,ispin,io,ik,1) + cp(ispin,io_all,ik)* cg%pko%zwf(ix,iy,iz,ispin,io,ik,1)
      cg%hxk%zwf(ix,iy,iz,ispin,io,ik,1) = &
      & cx(ispin,io_all,ik)* cg%hxk%zwf(ix,iy,iz,ispin,io,ik,1) + cp(ispin,io_all,ik)* cg%hwf%zwf(ix,iy,iz,ispin,io,ik,1)
    end do
    end do
    end do
    end do
    end do
    end do

    call inner_product(mg,system,info,cg%xk,cg%hxk,xkHxk)
    call inner_product(mg,system,info,cg%xk,cg%xk,xkxk)

    !$omp parallel do private(ik,io,ispin,iz,iy,ix,io_all) collapse(6)
    do ik=ik_s,ik_e
    do io=io_s,io_e
    do ispin=1,nspin
    do iz=is(3),ie(3)
    do iy=is(2),ie(2)
    do ix=is(1),ie(1)
      io_all = info%io_tbl(io)
      if(abs(xkxk(ispin,io_all,ik))<=1.d30)then
        spsi%zwf(ix,iy,iz,ispin,io,ik,1) = cg%xk%zwf(ix,iy,iz,ispin,io,ik,1) /sqrt(xkxk(ispin,io_all,ik))
      end if
    end do
    end do
    end do
    end do
    end do
    end do

  end do Iteration
  call timer_end(LOG_GSCG_ITERATION)

  call timer_end(LOG_GSCG_TOTAL)

  return
contains

subroutine orthogonalization(mg,system,info,psi,gk)
  implicit none
  type(s_rgrid),intent(in) :: mg
  type(s_dft_system),intent(in) :: system
  type(s_orbital_parallel),intent(in) :: info
  type(s_orbital),intent(in) :: psi
  type(s_orbital)            :: gk
  !
  integer :: nspin,no,nk,ik,ispin,ik_s,ik_e,io_s,io_e,is(3),ie(3),ix,iy,iz,io1,io2,io1_all,io2_all
  complex(8) :: sum0
  complex(8),dimension(system%no,system%no,system%nspin,system%nk) :: sum_obmat0,sum_obmat1

  nspin = system%nspin
  no = system%no
  nk = system%nk
  is = mg%is
  ie = mg%ie
  ik_s = info%ik_s
  ik_e = info%ik_e
  io_s = info%io_s
  io_e = info%io_e

  sum_obmat0 = 0.d0
!$omp parallel do private(ik,ispin,io1,io2,sum0,iz,iy,ix) collapse(3)
  do ik=ik_s,ik_e
  do ispin=1,nspin
    do io1=io_s,io_e
      do io2=io_s,io1-1
        sum0 = 0.d0
! speed? !$omp parallel do private(iz,iy,ix) collapse(2) reduction(+ : sum0)
        do iz=is(3),ie(3)
        do iy=is(2),ie(2)
        do ix=is(1),ie(1)
          sum0=sum0+conjg(psi%zwf(ix,iy,iz,ispin,io2,ik,1))*gk%zwf(ix,iy,iz,ispin,io1,ik,1)
        end do
        end do
        end do
        sum_obmat0(io1,io2,ispin,ik) = sum0*system%hvol
      end do
    end do
  end do
  end do
  call timer_begin(LOG_GSCG_ALLREDUCE)
  call comm_summation(sum_obmat0,sum_obmat1,no**2*nspin*nk,info%icomm_rko)
  call timer_end(LOG_GSCG_ALLREDUCE)
!$omp parallel do private(ik,ispin,io1,io2,iz,iy,ix) collapse(3)
  do ik=ik_s,ik_e
  do ispin=1,nspin
    do io1=io_s,io_e
      do io2=io_s,io1-1
! speed? !$omp parallel do private(iz,iy,ix) collapse(2)
        do iz=is(3),ie(3)
        do iy=is(2),ie(2)
        do ix=is(1),ie(1)
          gk%zwf(ix,iy,iz,ispin,io1,ik,1) = gk%zwf(ix,iy,iz,ispin,io1,ik,1) &
          & -sum_obmat1(io1,io2,ispin,ik) * psi%zwf(ix,iy,iz,ispin,io2,ik,1)
        end do
        end do
        end do
      end do
    end do
  end do
  end do

end subroutine orthogonalization

subroutine inner_product(mg,system,info,psi1,psi2,zbox)
  !$ use omp_lib
  implicit none
  type(s_rgrid),intent(in) :: mg
  type(s_dft_system),intent(in) :: system
  type(s_orbital_parallel),intent(in) :: info
  type(s_orbital),intent(in) :: psi1,psi2
  complex(8),intent(out) :: zbox(system%nspin,system%no,system%nk)
  !
  integer :: io,io_all,ik,ispin,nspin
  integer :: ix,iy,iz
  complex(8) :: zbox2(system%nspin,system%no,system%nk)
  complex(8) :: sum0
  nspin = system%nspin

  zbox2(:,:,:) = 0.d0
!$OMP parallel do collapse(2) private(ik,io,ispin,io_all,sum0,iz,iy,ix)
  do ik=info%ik_s,info%ik_e
  do io=info%io_s,info%io_e
  do ispin=1,nspin
    io_all = info%io_tbl(io)
    sum0 = 0d0
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      sum0 = sum0 + conjg(psi1%zwf(ix,iy,iz,ispin,io,ik,1))*psi2%zwf(ix,iy,iz,ispin,io,ik,1)
    end do
    end do
    end do
    zbox2(ispin,io_all,ik) = sum0 * system%hvol
  end do
  end do
  end do

  call comm_summation(zbox2,zbox,nspin*system%no*system%nk,info%icomm_r)

end subroutine inner_product
  
end subroutine gscg_periodic

end module Conjugate_Gradient
