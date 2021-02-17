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

  use Conjugate_Gradient_so, only: gscg_rwf_so, gscg_zwf_so, SPIN_ORBIT_ON

  implicit none
  private
  public :: gscg_rwf
  public :: gscg_zwf

contains

subroutine gscg_rwf(ncg,mg,system,info,stencil,ppg,vlocal,srg,spsi,cg)
  use structures
  use timer
  use hamiltonian, only: hpsi
  use communication, only: comm_summation
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
  integer :: iter
  real(8),dimension(system%nspin,system%no) :: sum,xkxk,xkHxk,xkHpk,pkHpk,gkgk,uk,ev,cx,cp,zs

  if ( SPIN_ORBIT_ON ) then
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

  if(.not. allocated(cg%xk%rwf)) then
    call allocate_orbital_real(nspin,mg,info,cg%xk)
    call allocate_orbital_real(nspin,mg,info,cg%hxk)
    call allocate_orbital_real(nspin,mg,info,cg%pk)
    call allocate_orbital_real(nspin,mg,info,cg%gk)
    call allocate_orbital_real(nspin,mg,info,cg%pko)
    call allocate_orbital_real(nspin,mg,info,cg%hwf)
  end if

!$omp parallel do private(io,ispin,iz,iy) collapse(4)
  do io=io_s,io_e
  do ispin=1,nspin
  do iz=is(3),ie(3)
  do iy=is(2),ie(2)
    cg%xk%rwf(is(1):ie(1),iy,iz,ispin,io,1,1) = spsi%rwf(is(1):ie(1),iy,iz,ispin,io,1,1)
  end do
  end do
  end do
  end do
  call timer_end(LOG_GSCG_ISOLATED_CALC)

  call timer_begin(LOG_GSCG_ISOLATED_HPSI)
  call hpsi(cg%xk,cg%hxk,info,mg,vlocal,system,stencil,srg,ppg)
  call timer_end(LOG_GSCG_ISOLATED_HPSI)

  call timer_begin(LOG_GSCG_ISOLATED_CALC)
  call inner_product(mg,system,info,cg%xk,cg%hxk,xkHxk)

  Iteration : do iter=1,Ncg

!$omp parallel do private(io,ispin,iz,iy) collapse(4)
    do io=io_s,io_e
    do ispin=1,nspin
    do iz=is(3),ie(3)
    do iy=is(2),ie(2)
      cg%gk%rwf(is(1):ie(1),iy,iz,ispin,io,1,1) = &
      & cg%hxk%rwf(is(1):ie(1),iy,iz,ispin,io,1,1) &
      & - xkHxk(ispin,io) * cg%xk%rwf(is(1):ie(1),iy,iz,ispin,io,1,1)
    end do
    end do
    end do
    end do

!    call orthogonalization(mg,system,info,spsi,cg%gk)
    call inner_product(mg,system,info,cg%gk,cg%gk,sum)

    if(iter==1)then
      uk = 0d0
    else
!$omp parallel do private(io,ispin)
      do io=io_s,io_e
      do ispin=1,nspin
        if (abs(gkgk(ispin,io)) > 1d-16) then
          uk(ispin,io) = sum(ispin,io) / gkgk(ispin,io)
        else
          uk(ispin,io) = 0d0
        end if
      end do
      end do
    end if

!$omp parallel do private(io,ispin,iz,iy) collapse(4)
    do io=io_s,io_e
    do ispin=1,nspin
    do iz=is(3),ie(3)
    do iy=is(2),ie(2)
      cg%pk%rwf(is(1):ie(1),iy,iz,ispin,io,1,1) = &
      & cg%gk%rwf(is(1):ie(1),iy,iz,ispin,io,1,1) &
      & + uk(ispin,io) * cg%pk%rwf(is(1):ie(1),iy,iz,ispin,io,1,1)
    end do
    end do
    end do
    end do

    gkgk = sum
    call inner_product(mg,system,info,cg%xk,cg%pk,zs)

!$omp parallel do private(io,ispin,iz,iy) collapse(4)
    do io=io_s,io_e
    do ispin=1,nspin
    do iz=is(3),ie(3)
    do iy=is(2),ie(2)
      cg%pko%rwf(is(1):ie(1),iy,iz,ispin,io,1,1) = &
      & cg%pk%rwf(is(1):ie(1),iy,iz,ispin,io,1,1) &
      & - zs(ispin,io) * cg%xk%rwf(is(1):ie(1),iy,iz,ispin,io,1,1)
    end do
    end do
    end do
    end do

    call inner_product(mg,system,info,cg%pko,cg%pko,sum)

!$omp parallel do private(io,ispin,iz,iy) collapse(4)
    do io=io_s,io_e
    do ispin=1,nspin
    do iz=is(3),ie(3)
    do iy=is(2),ie(2)
      cg%pko%rwf(is(1):ie(1),iy,iz,ispin,io,1,1) = &
      & cg%pko%rwf(is(1):ie(1),iy,iz,ispin,io,1,1) / sqrt(sum(ispin,io))
    end do
    end do
    end do
    end do
    call timer_end(LOG_GSCG_ISOLATED_CALC)

    call timer_begin(LOG_GSCG_ISOLATED_HPSI)
    call hpsi(cg%pko,cg%hwf,info,mg,vlocal,system,stencil,srg,ppg)
    call timer_end(LOG_GSCG_ISOLATED_HPSI)

    call timer_begin(LOG_GSCG_ISOLATED_CALC)
    call inner_product(mg,system,info,cg%xk,cg%hwf,xkHpk)
    call inner_product(mg,system,info,cg%pko,cg%hwf,pkHpk)

!$omp parallel do private(io,ispin)
    do io=io_s,io_e
    do ispin=1,nspin
      ev(ispin,io)=0.5d0*((xkHxk(ispin,io)+pkHpk(ispin,io))   &
                     -sqrt((xkHxk(ispin,io)-pkHpk(ispin,io))**2+4.d0*abs(xkHpk(ispin,io))**2))
      if (abs(ev(ispin,io) - xkHxk(ispin,io)) > 1d-16) then
        cx(ispin,io)=xkHpk(ispin,io)/(ev(ispin,io)-xkHxk(ispin,io))
        cp(ispin,io)=1.d0/sqrt(1.d0+abs(cx(ispin,io))**2)
        cx(ispin,io)=cx(ispin,io)*cp(ispin,io)
      else
        cx(ispin,io) = 1d0
        cp(ispin,io) = 0d0
      end if
    end do
    end do

!$omp parallel do private(io,ispin,iz,iy) collapse(4)
    do io=io_s,io_e
    do ispin=1,nspin
    do iz=is(3),ie(3)
    do iy=is(2),ie(2)
      cg%xk%rwf(is(1):ie(1),iy,iz,ispin,io,1,1) = &
      & cx(ispin,io)* cg%xk%rwf(is(1):ie(1),iy,iz,ispin,io,1,1) &
      & + cp(ispin,io) * cg%pko%rwf(is(1):ie(1),iy,iz,ispin,io,1,1)
      cg%hxk%rwf(is(1):ie(1),iy,iz,ispin,io,1,1) = &
      & cx(ispin,io)* cg%hxk%rwf(is(1):ie(1),iy,iz,ispin,io,1,1) &
      & + cp(ispin,io) * cg%hwf%rwf(is(1):ie(1),iy,iz,ispin,io,1,1)
    end do
    end do
    end do
    end do

    call inner_product(mg,system,info,cg%xk,cg%hxk,xkHxk)
    call inner_product(mg,system,info,cg%xk,cg%xk,xkxk)

!$omp parallel do private(io,ispin,iz,iy) collapse(4)
    do io=io_s,io_e
    do ispin=1,nspin
    do iz=is(3),ie(3)
    do iy=is(2),ie(2)
      if(1d-16 < abs(xkxk(ispin,io)) .and. abs(xkxk(ispin,io)) <= 1d30) then
        spsi%rwf(is(1):ie(1),iy,iz,ispin,io,1,1) = &
        & cg%xk%rwf(is(1):ie(1),iy,iz,ispin,io,1,1) / sqrt(xkxk(ispin,io))
      end if
    end do
    end do
    end do
    end do

  end do Iteration
  call timer_end(LOG_GSCG_ISOLATED_CALC)

  return
contains

subroutine orthogonalization(mg,system,info,psi,gk)
  implicit none
  type(s_rgrid),intent(in) :: mg
  type(s_dft_system),intent(in) :: system
  type(s_parallel_info),intent(in) :: info
  type(s_orbital),intent(in) :: psi
  type(s_orbital)            :: gk
  !
  integer :: nspin,no,ispin,io_s,io_e,is(3),ie(3),ix,iy,iz,io1,io2
  complex(8) :: sum0
  complex(8),dimension(system%no,system%no,system%nspin) :: sum_obmat0,sum_obmat1

  if(info%if_divide_orbit) stop "error: nproc_ob/=1 @ gscg"

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
  call timer_end(LOG_GSCG_ISOLATED_CALC)

  call timer_begin(LOG_GSCG_ISOLATED_COMM_COLL)
  call comm_summation(sum_obmat0,sum_obmat1,no**2*nspin,info%icomm_rko)
  call timer_end(LOG_GSCG_ISOLATED_COMM_COLL)

  call timer_begin(LOG_GSCG_ISOLATED_CALC)
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
!$OMP parallel do collapse(2) private(io,ispin,sum0,iz,iy,ix)
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

end subroutine gscg_rwf

!===================================================================================================================================

subroutine gscg_zwf(ncg,mg,system,info,stencil,ppg,vlocal,srg,spsi,cg)
  use structures
  use timer
  use hamiltonian, only: hpsi
  use communication, only: comm_summation
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
  integer :: iter
  complex(8),dimension(system%nspin,system%no,system%nk) :: sum,xkxk,xkHxk,xkHpk,pkHpk,gkgk,uk,ev,cx,cp,zs

  if ( SPIN_ORBIT_ON ) then
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

  if(.not. allocated(cg%xk%zwf)) then
    call allocate_orbital_complex(nspin,mg,info,cg%xk)
    call allocate_orbital_complex(nspin,mg,info,cg%hxk)
    call allocate_orbital_complex(nspin,mg,info,cg%pk)
    call allocate_orbital_complex(nspin,mg,info,cg%gk)
    call allocate_orbital_complex(nspin,mg,info,cg%pko)
    call allocate_orbital_complex(nspin,mg,info,cg%hwf)
    !$acc enter data copyin(cg)
  end if

!$omp parallel do private(ik,io,ispin,iz,iy) collapse(5)
  do ik=ik_s,ik_e
  do io=io_s,io_e
  do ispin=1,nspin
  do iz=is(3),ie(3)
  do iy=is(2),ie(2)
    cg%xk%zwf(is(1):ie(1),iy,iz,ispin,io,ik,1) = spsi%zwf(is(1):ie(1),iy,iz,ispin,io,ik,1)
  end do
  end do
  end do
  end do
  end do
  call timer_end(LOG_GSCG_PERIODIC_CALC)

  call timer_begin(LOG_GSCG_PERIODIC_HPSI)
  call hpsi(cg%xk,cg%hxk,info,mg,vlocal,system,stencil,srg,ppg)
  call timer_end(LOG_GSCG_PERIODIC_HPSI)

  call timer_begin(LOG_GSCG_PERIODIC_CALC)
  call inner_product(mg,system,info,cg%xk,cg%hxk,xkHxk)

  Iteration : do iter=1,Ncg

!$omp parallel do private(ik,io,ispin,iz,iy) collapse(5)
    do ik=ik_s,ik_e
    do io=io_s,io_e
    do ispin=1,nspin
    do iz=is(3),ie(3)
    do iy=is(2),ie(2)
      cg%gk%zwf(is(1):ie(1),iy,iz,ispin,io,ik,1) = &
      & cg%hxk%zwf(is(1):ie(1),iy,iz,ispin,io,ik,1) &
      & - xkHxk(ispin,io,ik) * cg%xk%zwf(is(1):ie(1),iy,iz,ispin,io,ik,1)
    end do
    end do
    end do
    end do
    end do

!    call orthogonalization(mg,system,info,spsi,cg%gk)
    call inner_product(mg,system,info,cg%gk,cg%gk,sum)

    if(iter==1)then
      uk = 0d0
    else
!$omp parallel do private(ik,io,ispin) collapse(2)
      do ik=ik_s,ik_e
      do io=io_s,io_e
      do ispin=1,nspin
        if (abs(gkgk(ispin,io,ik)) > 1d-16) then
          uk(ispin,io,ik) = sum(ispin,io,ik) / gkgk(ispin,io,ik)
        else
          uk(ispin,io,ik) = 0d0
        end if
      end do
      end do
      end do
    end if

!$omp parallel do private(ik,io,ispin,iz,iy) collapse(5)
    do ik=ik_s,ik_e
    do io=io_s,io_e
    do ispin=1,nspin
    do iz=is(3),ie(3)
    do iy=is(2),ie(2)
      cg%pk%zwf(is(1):ie(1),iy,iz,ispin,io,ik,1) = &
      & cg%gk%zwf(is(1):ie(1),iy,iz,ispin,io,ik,1) &
      & + uk(ispin,io,ik) * cg%pk%zwf(is(1):ie(1),iy,iz,ispin,io,ik,1)
    end do
    end do
    end do
    end do
    end do

    gkgk = sum
    call inner_product(mg,system,info,cg%xk,cg%pk,zs)

!$omp parallel do private(ik,io,ispin,iz,iy) collapse(5)
    do ik=ik_s,ik_e
    do io=io_s,io_e
    do ispin=1,nspin
    do iz=is(3),ie(3)
    do iy=is(2),ie(2)
      cg%pko%zwf(is(1):ie(1),iy,iz,ispin,io,ik,1) = &
      & cg%pk%zwf(is(1):ie(1),iy,iz,ispin,io,ik,1) &
      & - zs(ispin,io,ik) *cg%xk%zwf(is(1):ie(1),iy,iz,ispin,io,ik,1)
    end do
    end do
    end do
    end do
    end do

    call inner_product(mg,system,info,cg%pko,cg%pko,sum)

!$omp parallel do private(ik,io,ispin,iz,iy) collapse(5)
    do ik=ik_s,ik_e
    do io=io_s,io_e
    do ispin=1,nspin
    do iz=is(3),ie(3)
    do iy=is(2),ie(2)
      cg%pko%zwf(is(1):ie(1),iy,iz,ispin,io,ik,1) = &
      & cg%pko%zwf(is(1):ie(1),iy,iz,ispin,io,ik,1) / sqrt(sum(ispin,io,ik))
    end do
    end do
    end do
    end do
    end do
    call timer_end(LOG_GSCG_PERIODIC_CALC)

    call timer_begin(LOG_GSCG_PERIODIC_HPSI)
    call hpsi(cg%pko,cg%hwf,info,mg,vlocal,system,stencil,srg,ppg)
    call timer_end(LOG_GSCG_PERIODIC_HPSI)

    call timer_begin(LOG_GSCG_PERIODIC_CALC)
    call inner_product(mg,system,info,cg%xk,cg%hwf,xkHpk)
    call inner_product(mg,system,info,cg%pko,cg%hwf,pkHpk)

!$omp parallel do private(ik,io,ispin) collapse(2)
    do ik=ik_s,ik_e
    do io=io_s,io_e
    do ispin=1,nspin
      ev(ispin,io,ik)=0.5d0*((xkHxk(ispin,io,ik)+pkHpk(ispin,io,ik))   &
                     -sqrt((xkHxk(ispin,io,ik)-pkHpk(ispin,io,ik))**2+4.d0*abs(xkHpk(ispin,io,ik))**2))
      if (abs(ev(ispin,io,ik) - xkHxk(ispin,io,ik)) > 1d-16) then
        cx(ispin,io,ik)=xkHpk(ispin,io,ik)/(ev(ispin,io,ik)-xkHxk(ispin,io,ik))
        cp(ispin,io,ik)=1.d0/sqrt(1.d0+abs(cx(ispin,io,ik))**2)
        cx(ispin,io,ik)=cx(ispin,io,ik)*cp(ispin,io,ik)
      else
        cx(ispin,io,ik) = 1d0
        cp(ispin,io,ik) = 0d0
      end if
    end do
    end do
    end do

!$omp parallel do private(ik,io,ispin,iz,iy) collapse(5)
    do ik=ik_s,ik_e
    do io=io_s,io_e
    do ispin=1,nspin
    do iz=is(3),ie(3)
    do iy=is(2),ie(2)
      cg%xk%zwf(is(1):ie(1),iy,iz,ispin,io,ik,1) = &
      & cx(ispin,io,ik) * cg%xk%zwf(is(1):ie(1),iy,iz,ispin,io,ik,1) &
      & + cp(ispin,io,ik) * cg%pko%zwf(is(1):ie(1),iy,iz,ispin,io,ik,1)
      cg%hxk%zwf(is(1):ie(1),iy,iz,ispin,io,ik,1) = &
      & cx(ispin,io,ik) * cg%hxk%zwf(is(1):ie(1),iy,iz,ispin,io,ik,1) &
      & + cp(ispin,io,ik) * cg%hwf%zwf(is(1):ie(1),iy,iz,ispin,io,ik,1)
    end do
    end do
    end do
    end do
    end do

    call inner_product(mg,system,info,cg%xk,cg%hxk,xkHxk)
    call inner_product(mg,system,info,cg%xk,cg%xk,xkxk)

!$omp parallel do private(ik,io,ispin,iz,iy) collapse(5)
    do ik=ik_s,ik_e
    do io=io_s,io_e
    do ispin=1,nspin
    do iz=is(3),ie(3)
    do iy=is(2),ie(2)
      if(1d-16 < abs(xkxk(ispin,io,ik)) .and. abs(xkxk(ispin,io,ik)) <= 1d30) then
        spsi%zwf(is(1):ie(1),iy,iz,ispin,io,ik,1) = &
        & cg%xk%zwf(is(1):ie(1),iy,iz,ispin,io,ik,1) / sqrt(xkxk(ispin,io,ik))
      end if
    end do
    end do
    end do
    end do
    end do

  end do Iteration
  call timer_end(LOG_GSCG_PERIODIC_CALC)

  return
contains

subroutine orthogonalization(mg,system,info,psi,gk)
  implicit none
  type(s_rgrid),intent(in) :: mg
  type(s_dft_system),intent(in) :: system
  type(s_parallel_info),intent(in) :: info
  type(s_orbital),intent(in) :: psi
  type(s_orbital)            :: gk
  !
  integer :: nspin,no,nk,ik,ispin,ik_s,ik_e,io_s,io_e,is(3),ie(3),ix,iy,iz,io1,io2,io1_all,io2_all
  complex(8) :: sum0
  complex(8),dimension(system%no,system%no,system%nspin,system%nk) :: sum_obmat0,sum_obmat1

  if(info%if_divide_orbit) stop "error: nproc_ob/=1 @ gscg"

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
  call timer_end(LOG_GSCG_PERIODIC_CALC)

  call timer_begin(LOG_GSCG_PERIODIC_COMM_COLL)
  call comm_summation(sum_obmat0,sum_obmat1,no**2*nspin*nk,info%icomm_rko)
  call timer_end(LOG_GSCG_PERIODIC_COMM_COLL)

  call timer_begin(LOG_GSCG_PERIODIC_CALC)
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
!$OMP parallel do collapse(2) private(ik,io,ispin,sum0,iz,iy,ix)
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

end subroutine gscg_zwf

end module Conjugate_Gradient
