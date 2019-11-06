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
!This file is "psi_rho.f90"
!This file contain one subroutine.
!Subroutine psi_rho
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
#ifdef __INTEL_COMPILER
# define OMP_SIMD simd
#else
# define OMP_SIMD
#endif

subroutine psi_rho_GS
  use global_variables, only: zu_GS,NB
  implicit none
  call psi_rho_impl(zu_GS,NB)
end subroutine

subroutine psi_rho_RT(zu)
  use global_variables, only: NL,NBoccmax,NK_s,NK_e
  implicit none
  complex(8),intent(in) :: zu(NL,NBoccmax,NK_s:NK_e)
  call psi_rho_impl(zu,NBoccmax)
end subroutine

subroutine psi_rho_impl(zutmp,zu_NB)
  use global_variables
  use timer
  use opt_variables
  use salmon_parallel, only: nproc_group_tdks
  use salmon_communication, only: comm_summation
  implicit none
  integer,intent(in)    :: zu_NB
  complex(8),intent(in) :: zutmp(0:NL-1,zu_NB,NK_s:NK_e)

  call timer_begin(LOG_CALC_RHO)
  ! write(*,*) "Sym:", Sym
  ! stop
  select case(Sym)
  case(1)
    call sym1(zutmp,zu_NB,rho_l)
  case(4)
    if(crystal_structure == 'diamond')then
      call sym4(zutmp,zu_NB,rho_l,rho_tmp1)
    else
      call err_finalize('Bad crystal structure')
    end if
  case(8)
    if(crystal_structure == 'diamond')then
       call sym8(zutmp,zu_NB,rho_l,rho_tmp1,rho_tmp2)
    else if(crystal_structure == 'diamond2')then
       call sym8_diamond2(zutmp,zu_NB,rho_l,rho_tmp1,rho_tmp2)
    else if(crystal_structure == 'tetragonal')then
       call sym8_tetragonal(zutmp,zu_NB,rho_l,rho_tmp1,rho_tmp2)
    else
       call err_finalize('Bad crystal structure')
    end if
  case default
    call err_finalize('Bad Symmetry')
  end select
  call timer_end(LOG_CALC_RHO)

  call comm_summation(rho_l,rho,NL,nproc_group_tdks)


contains
  subroutine reduce(tid,zfac,zutmp,zu_NB)
    use global_variables
    use opt_variables, only: zrhotmp
    use misc_routines, only: ceiling_pow2
    use omp_lib, only: omp_get_thread_num
    implicit none
    integer,intent(in)    :: tid
    real(8),intent(in)    :: zfac
    integer,intent(in)    :: zu_NB
    complex(8),intent(in) :: zutmp(0:NL-1,zu_NB,NK_s:NK_e)

    integer :: ib,ik,i,mytid

    mytid = tid

    zrhotmp(:,mytid)=0.d0

!$omp do private(ik,ib,i) collapse(2)
    do ik=NK_s,NK_e
    do ib=1,NBoccmax
    do i=0,NL-1
      zrhotmp(i,mytid)=zrhotmp(i,mytid)+(zfac*occ(ib,ik))*abs(zutmp(i,ib,ik))**2
    end do
    end do
    end do
!$omp end do

    i = ceiling_pow2(NUMBER_THREADS)/2
    do while(i > 0)
      if(mytid < i) then
        zrhotmp(0:NL-1,mytid) = zrhotmp(0:NL-1,mytid) + zrhotmp(0:NL-1,mytid + i)
      end if
      i = i/2
!$omp barrier
    end do
  end subroutine

  subroutine sym1(zutmp,zu_NB,zrho_l)
    use global_variables
    use opt_variables, only: zrhotmp
    use omp_lib, only: omp_get_thread_num
    implicit none
    integer,intent(in)    :: zu_NB
    complex(8),intent(in) :: zutmp(0:NL-1,zu_NB,NK_s:NK_e)
    real(8),intent(out)   :: zrho_l(0:NL-1)

    integer :: tid

!$omp parallel private(tid)
!$  tid=omp_get_thread_num()
    call reduce(tid,1.0d0,zutmp,zu_NB)
!$omp end parallel

    zrho_l(:) = zrhotmp(0:NL-1,0)
  end subroutine

  !====== diamond(4) structure =========================!
  subroutine sym4(zutmp,zu_NB,zrho_l,zrhotmp1)
    use global_variables
    use opt_variables, only: zrhotmp
    use omp_lib, only: omp_get_thread_num
    implicit none
    integer,intent(in)    :: zu_NB
    complex(8),intent(in) :: zutmp(0:NL-1,zu_NB,NK_s:NK_e)
    real(8),intent(out)   :: zrho_l(0:NL-1)

    real(8) :: zrhotmp1(0:NL-1)
    integer :: i,tid
    real(8) :: zfac

    zfac=1.0d0/4d0

!$omp parallel private(tid)
!$  tid=omp_get_thread_num()

    call reduce(tid,zfac,zutmp,zu_NB)

! 1.T_3
!$omp do OMP_SIMD
    do i=0,NL-1
      zrhotmp1(i)=zrhotmp(i,0)+zrhotmp(itable_sym(1,i+1)-1,0)
    end do
!$omp end do OMP_SIMD

! 2.T_3*T_3
!$omp do OMP_SIMD
    do i=0,NL-1
      zrho_l(i)=zrhotmp1(i)+zrhotmp1(itable_sym(2,i+1)-1)
    end do
!$omp end do OMP_SIMD
!$omp end parallel
  end subroutine

  !====== diamond(8) structure =========================!
  subroutine sym8(zutmp,zu_NB,zrho_l,zrhotmp1,zrhotmp2)
    use global_variables
    use opt_variables, only: zrhotmp
    use omp_lib, only: omp_get_thread_num
    implicit none
    integer,intent(in)    :: zu_NB
    complex(8),intent(in) :: zutmp(0:NL-1,zu_NB,NK_s:NK_e)
    real(8),intent(out)   :: zrho_l(0:NL-1)

    real(8) :: zrhotmp1(0:NL-1)
    real(8) :: zrhotmp2(0:NL-1)
    integer :: i,tid
    real(8) :: zfac

    ! wk(ik)=8.0,(ikx==iky >. wk(ik)=4.0)
    zfac=1.0d0/32d0

!$omp parallel private(tid)
!$  tid=omp_get_thread_num()

    call reduce(tid,zfac,zutmp,zu_NB)

! 1.T_4
!$omp do OMP_SIMD
    do i=0,NL-1
      zrhotmp2(i)=zrhotmp(i,0)+zrhotmp(itable_sym(4,i+1)-1,0)
    end do
!$omp end do OMP_SIMD

! 2.T_3*T_3
!$omp do OMP_SIMD
    do i=0,NL-1
      zrhotmp1(i)=zrhotmp2(i)+zrhotmp2(itable_sym(5,i+1)-1)
    end do
!$omp end do OMP_SIMD

! 2.T_3
!$omp do OMP_SIMD
    do i=0,NL-1
      zrhotmp2(i)=zrhotmp1(i)+zrhotmp1(itable_sym(3,i+1)-1)
    end do
!$omp end do OMP_SIMD

! 2.T_1
!$omp do OMP_SIMD
    do i=0,NL-1
      zrhotmp1(i)=zrhotmp2(i)+zrhotmp2(itable_sym(1,i+1)-1)
    end do
!$omp end do OMP_SIMD

! 2.T_2
!$omp do OMP_SIMD
    do i=0,NL-1
      zrho_l(i)=zrhotmp1(i)+zrhotmp1(itable_sym(2,i+1)-1)
    end do
!$omp end do OMP_SIMD
!$omp end parallel
  end subroutine

   !====== diamond2(8) structure =========================!
  subroutine sym8_diamond2(zutmp,zu_NB,zrho_l,zrhotmp1,zrhotmp2)
    use global_variables
    use opt_variables, only: zrhotmp
    use omp_lib, only: omp_get_thread_num
    implicit none
    integer,intent(in)    :: zu_NB
    complex(8),intent(in) :: zutmp(0:NL-1,zu_NB,NK_s:NK_e)
    real(8),intent(out)   :: zrho_l(0:NL-1)

    real(8) :: zrhotmp1(0:NL-1)
    real(8) :: zrhotmp2(0:NL-1)
    integer :: i,tid
    real(8) :: zfac

    ! wk(ik)=8.0,(ikx==iky >. wk(ik)=4.0)
    zfac=1.0d0/16d0

!$omp parallel private(tid)
!$  tid=omp_get_thread_num()

    call reduce(tid,zfac,zutmp,zu_NB)

!$omp do OMP_SIMD
    do i=0,NL-1
      zrhotmp1(i)=zrhotmp(i,0)+zrhotmp(itable_sym(4,i+1)-1,0)
    end do
!$omp end do OMP_SIMD

!$omp do OMP_SIMD
    do i=0,NL-1
      zrhotmp2(i)=zrhotmp1(i)+zrhotmp1(itable_sym(3,i+1)-1)
    end do
!$omp end do OMP_SIMD

!$omp do OMP_SIMD
    do i=0,NL-1
      zrhotmp1(i)=zrhotmp2(i)+zrhotmp2(itable_sym(1,i+1)-1)
    end do
!$omp end do OMP_SIMD

!$omp do OMP_SIMD
    do i=0,NL-1
      zrho_l(i)=zrhotmp1(i)+zrhotmp1(itable_sym(2,i+1)-1)
    end do
!$omp end do OMP_SIMD
!$omp end parallel
  end subroutine

  !====== tetragonal(8) structure =========================!
  subroutine sym8_tetragonal(zutmp,zu_NB,zrho_l,zrhotmp1,zrhotmp2)
    use global_variables
    use opt_variables, only: zrhotmp
    use omp_lib, only: omp_get_thread_num
    implicit none
    integer,intent(in)    :: zu_NB
    complex(8),intent(in) :: zutmp(0:NL-1,zu_NB,NK_s:NK_e)
    real(8),intent(out)   :: zrho_l(0:NL-1)

    real(8) :: zrhotmp1(0:NL-1)
    real(8) :: zrhotmp2(0:NL-1)
    integer :: i,tid
    real(8) :: zfac

    ! wk(ik)=8.0,(ikx==iky >. wk(ik)=4.0)
    zfac=1.0d0/16d0

!$omp parallel private(tid)
!$  tid=omp_get_thread_num()

    call reduce(tid,zfac,zutmp,zu_NB)

! 1.T_3
!$omp do OMP_SIMD
    do i=0,NL-1
      zrhotmp2(i)=zrhotmp(i,0)+zrhotmp(itable_sym(3,i+1)-1,0)
    end do
!$omp end do OMP_SIMD

! 2.T_3*T_3
!$omp do OMP_SIMD
    do i=0,NL-1
      zrhotmp1(i)=zrhotmp2(i)+zrhotmp2(itable_sym(4,i+1)-1)
    end do
!$omp end do OMP_SIMD

! 3.T_1
!$omp do OMP_SIMD
    do i=0,NL-1
      zrhotmp2(i)=zrhotmp1(i)+zrhotmp1(itable_sym(1,i+1)-1)
    end do
!$omp end do OMP_SIMD

! 4.T_2
!$omp do OMP_SIMD
    do i=0,NL-1
      zrho_l(i)=zrhotmp2(i)+zrhotmp2(itable_sym(2,i+1)-1)
    end do
!$omp end do OMP_SIMD
!$omp end parallel
  end subroutine
end subroutine psi_rho_impl

