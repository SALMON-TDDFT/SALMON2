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

#include "config.h"

module nonlocal_potential
  implicit none

! WARNING: We must not call these except for hpsi routine.

contains

subroutine dpseudo(tpsi,htpsi,info,nspin,ppg)
  use structures
  use communication, only: comm_summation
  use timer
  implicit none
  integer,intent(in) :: nspin
  type(s_parallel_info),intent(in) :: info
  type(s_pp_grid),intent(in) :: ppg
  type(s_orbital),intent(in) :: tpsi
  type(s_orbital) :: htpsi
  !
  integer :: ispin,io,ik,im,im_s,im_e,ik_s,ik_e,io_s,io_e,norb
  integer :: ilma,ia,j,ix,iy,iz,Nlma,ilocal
  real(8) :: uVpsi,wrk
  real(8),allocatable :: uVpsibox (:,:,:,:,:)
  real(8),allocatable :: uVpsibox2(:,:,:,:,:)

  call timer_begin(LOG_UHPSI_PSEUDO)

  im_s = info%im_s
  im_e = info%im_e
  ik_s = info%ik_s
  ik_e = info%ik_e
  io_s = info%io_s
  io_e = info%io_e
  norb = Nspin* info%numo * info%numk * info%numm

  Nlma = ppg%Nlma

  if(info%if_divide_rspace) then

    allocate(uVpsibox (Nlma,Nspin,io_s:io_e,ik_s:ik_e,im_s:im_e))
    allocate(uVpsibox2(Nlma,Nspin,io_s:io_e,ik_s:ik_e,im_s:im_e))

!$omp parallel do collapse(4) &
!$omp             private(im,ik,io,ispin,ilma,ia,uVpsi,j,ix,iy,iz)
    do im=im_s,im_e
    do ik=ik_s,ik_e
    do io=io_s,io_e
    do ispin=1,Nspin

      do ilma=1,Nlma
        ia = ppg%ia_tbl(ilma)
        uVpsi = 0.d0
        do j=1,ppg%mps(ia)
          ix = ppg%jxyz(1,j,ia)
          iy = ppg%jxyz(2,j,ia)
          iz = ppg%jxyz(3,j,ia)
          uVpsi = uVpsi + ppg%uV(j,ilma) * tpsi%rwf(ix,iy,iz,ispin,io,ik,im)
        end do
        uVpsi = uVpsi * ppg%rinv_uvu(ilma)
        uVpsibox(ilma,ispin,io,ik,im) = uVpsi
      end do

    end do
    end do
    end do
    end do
!$omp end parallel do

    call timer_end(LOG_UHPSI_PSEUDO)

    call timer_begin(LOG_UHPSI_PSEUDO_COMM)
    call comm_summation(uVpsibox,uVpsibox2,Nlma*Norb,info%icomm_r)
    call timer_end(LOG_UHPSI_PSEUDO_COMM)

    call timer_begin(LOG_UHPSI_PSEUDO)

!$omp parallel do collapse(4) &
!$omp             private(im,ik,io,ispin,ilocal,ilma,ia,uVpsi,j,ix,iy,iz,wrk)
    do im=im_s,im_e
    do ik=ik_s,ik_e
    do io=io_s,io_e
    do ispin=1,Nspin

      do ilocal=1,ppg%ilocal_nlma
        ilma=ppg%ilocal_nlma2ilma(ilocal)
        ia  =ppg%ilocal_nlma2ia  (ilocal)
        uVpsi = uVpsibox2(ilma,ispin,io,ik,im)
!OCL norecurrence
        do j=1,ppg%mps(ia)
          ix = ppg%jxyz(1,j,ia)
          iy = ppg%jxyz(2,j,ia)
          iz = ppg%jxyz(3,j,ia)
          wrk = uVpsi * ppg%uV(j,ilma)
          htpsi%rwf(ix,iy,iz,ispin,io,ik,im) = htpsi%rwf(ix,iy,iz,ispin,io,ik,im) + wrk
        end do
      end do

    end do
    end do
    end do
    end do
!$omp end parallel do

    deallocate(uVpsibox,uVpsibox2)

  else

#ifdef USE_OPENACC
!$acc parallel loop collapse(4) private(im,ik,io,ispin,ilma,ia,uVpsi,j,ix,iy,iz,wrk)
#else
!$omp parallel do collapse(4) &
!$omp             private(im,ik,io,ispin,ilma,ia,uVpsi,j,ix,iy,iz,wrk)
#endif
    do im=im_s,im_e
    do ik=ik_s,ik_e
    do io=io_s,io_e
    do ispin=1,Nspin

      do ilma=1,Nlma
        ia = ppg%ia_tbl(ilma)
        uVpsi = 0.d0
        do j=1,ppg%mps(ia)
          ix = ppg%jxyz(1,j,ia)
          iy = ppg%jxyz(2,j,ia)
          iz = ppg%jxyz(3,j,ia)
          uVpsi = uVpsi + ppg%uV(j,ilma) * tpsi%rwf(ix,iy,iz,ispin,io,ik,im)
        end do
        uVpsi = uVpsi * ppg%rinv_uvu(ilma)
!OCL norecurrence
        do j=1,ppg%mps(ia)
          ix = ppg%jxyz(1,j,ia)
          iy = ppg%jxyz(2,j,ia)
          iz = ppg%jxyz(3,j,ia)
          wrk = uVpsi * ppg%uV(j,ilma)
          htpsi%rwf(ix,iy,iz,ispin,io,ik,im) = htpsi%rwf(ix,iy,iz,ispin,io,ik,im) + wrk
        end do
      end do

    end do
    end do
    end do
    end do
#ifdef USE_OPENACC
!$acc end parallel
#else
!$omp end parallel do
#endif

  end if

  call timer_end(LOG_UHPSI_PSEUDO)

  return
end subroutine dpseudo

!-----------------------------------------------------------------------------------------------------------------------------------

subroutine zpseudo(tpsi,htpsi,info,nspin,ppg)
  use structures
  use timer
  use iso_c_binding
  implicit none
  integer,intent(in) :: nspin
  type(s_parallel_info),intent(in) :: info
  type(s_pp_grid),intent(in) :: ppg
  type(s_orbital),intent(in) :: tpsi
  type(s_orbital) :: htpsi
  !
  integer :: ispin,io,ik,im,im_s,im_e,ik_s,ik_e,io_s,io_e
  integer :: ilma,ia,j,ix,iy,iz,Nlma,ilocal,vi,my_nlma,k
  complex(8) :: uVpsi,wrk
  complex(8),allocatable :: uVpsibox (:,:,:,:,:)
  complex(8),allocatable :: uVpsibox2(:,:,:,:,:)
  complex(8) :: IMAGINARY_UNIT = (0, 1)
#if defined(USE_OPENACC) && defined(USE_CUDA)
  integer :: n,natom
  interface
    subroutine zpseudo_cuda(htpsi_zwf,n,im_s,im_e,ik_s,ik_e,io_s,io_e,Nspin,Nlma,ppg_nps,natom,mg_is_array_1,mg_ie_array_1,mg_is_array_2,mg_ie_array_2,mg_is_array_3,mg_ie_array_3,ppg_ia_tbl,ppg_mps,ppg_jxyz,ppg_zekr_uV,ppg_rinv_uvu,tpsi_zwf) bind(c)
      import
      ! Input integer
      integer(c_int), intent(in), value :: n
      integer(c_int), intent(in), value :: im_s
      integer(c_int), intent(in), value :: im_e
      integer(c_int), intent(in), value :: ik_s
      integer(c_int), intent(in), value :: ik_e
      integer(c_int), intent(in), value :: io_s
      integer(c_int), intent(in), value :: io_e
      integer(c_int), intent(in), value :: Nspin
      integer(c_int), intent(in), value :: Nlma
      integer(c_int), intent(in), value :: ppg_nps
      integer(c_int), intent(in), value :: natom
      integer(c_int), intent(in), value :: mg_is_array_1
      integer(c_int), intent(in), value :: mg_ie_array_1
      integer(c_int), intent(in), value :: mg_is_array_2
      integer(c_int), intent(in), value :: mg_ie_array_2
      integer(c_int), intent(in), value :: mg_is_array_3
      integer(c_int), intent(in), value :: mg_ie_array_3
      ! Output
      complex(c_double_complex), intent(inout) :: htpsi_zwf(mg_is_array_1:mg_ie_array_1, mg_is_array_2:mg_ie_array_2, mg_is_array_2:mg_ie_array_3, Nspin, io_e:io_s, ik_s:ik_e, im_s:im_e)
      ! Input
      integer(c_int)           , intent(in) :: ppg_ia_tbl(n * natom)
      integer(c_int)           , intent(in) :: ppg_mps(natom)
      integer(c_int)           , intent(in) :: ppg_jxyz(3, ppg_nps, natom)
      complex(c_double_complex), intent(in) :: ppg_zekr_uV(ppg_nps, ppg_Nlma, ik_s:ik_e)
      real(c_double)           , intent(in) :: ppg_rinv_uvu(n * natom)
      complex(c_double_complex), intent(in) :: tpsi_zwf(mg_is_array_1:mg_ie_array_1, mg_is_array_2:mg_ie_array_2,mg_is_array_2:mg_ie_array_3, Nspin, io_e:io_s, ik_s:ik_e, im_s:im_e)
    end subroutine zpseudo_cuda
  end interface
#endif

#ifdef FORTRAN_COMPILER_HAS_2MB_ALIGNED_ALLOCATION
!dir$ attributes align : 2097152 :: uVpsibox, uVpsibox2
#endif

  call timer_begin(LOG_UHPSI_PSEUDO)

  im_s = info%im_s
  im_e = info%im_e
  ik_s = info%ik_s
  ik_e = info%ik_e
  io_s = info%io_s
  io_e = info%io_e

  Nlma = ppg%Nlma

  if(info%if_divide_rspace) then

    call timer_end(LOG_UHPSI_PSEUDO)

    call calc_uVpsi_rdivided(nspin,info,ppg,tpsi,uVpsibox,uVpsibox2)

    call timer_begin(LOG_UHPSI_PSEUDO)

!$omp parallel do collapse(4) &
!$omp             private(im,ik,io,ispin,ilocal,ilma,ia,uVpsi,j,ix,iy,iz,wrk)
    do im=im_s,im_e
    do ik=ik_s,ik_e
    do io=io_s,io_e
    do ispin=1,Nspin

      do ilocal=1,ppg%ilocal_nlma
        ilma=ppg%ilocal_nlma2ilma(ilocal)
        ia  =ppg%ilocal_nlma2ia  (ilocal)
        uVpsi = uVpsibox2(ispin,io,ik,im,ilma)
!OCL norecurrence
        do j=1,ppg%mps(ia)
          ix = ppg%jxyz(1,j,ia)
          iy = ppg%jxyz(2,j,ia)
          iz = ppg%jxyz(3,j,ia)
          wrk = uVpsi * ppg%zekr_uV(j,ilma,ik)
          htpsi%zwf(ix,iy,iz,ispin,io,ik,im) = htpsi%zwf(ix,iy,iz,ispin,io,ik,im) + wrk
        end do
      end do

    end do
    end do
    end do
    end do
!$omp end parallel do

    deallocate(uVpsibox,uVpsibox2)

  else
#ifdef USE_OPENACC
#ifdef USE_CUDA
    natom=size(ppg%mps)
    n=size(ppg%rinv_uvu)/natom
    call zpseudo_cuda(htpsi%zwf,&
        n,&
        im_s,im_e,&
        ik_s,ik_e,&
        io_s,io_e,&
        Nspin,&
        Nlma,&
        ppg%nps,&
        natom,&
        lbound(htpsi%zwf,1),ubound(htpsi%zwf,1),&
        lbound(htpsi%zwf,2),ubound(htpsi%zwf,2),&
        lbound(htpsi%zwf,3),ubound(htpsi%zwf,3),&
        ppg%ia_tbl,&
        ppg%mps,&
        ppg%jxyz,&
        ppg%zekr_uV,&
        ppg%rinv_uvu,&
        tpsi%zwf)
#else
!$acc kernels present(ppg,tpsi,htpsi)
!$acc loop private(ilocal,ilma,ia,uvpsi,vi,my_nlma,k,j,ix,iy,iz,wrk) collapse(4) gang
    do im=im_s,im_e
    do ik=ik_s,ik_e
    do io=io_s,io_e
    do ispin=1,Nspin

!$acc loop gang independent
      do ilma=1,Nlma
        ia = ppg%ia_tbl(ilma)
        uVpsi = 0.d0
!$acc loop vector reduction(+:uVpsi)
        do j=1,ppg%mps(ia)
          ix = ppg%jxyz(1,j,ia)
          iy = ppg%jxyz(2,j,ia)
          iz = ppg%jxyz(3,j,ia)
          uVpsi = uVpsi + conjg(ppg%zekr_uV(j,ilma,ik)) * tpsi%zwf(ix,iy,iz,ispin,io,ik,im)
        end do
        ppg%uVpsibox(ilma,ispin,io,ik,im) = uVpsi * ppg%rinv_uvu(ilma)
      end do

!$acc loop gang vector independent
      do vi=0,ppg%max_vi-1
        my_nlma = ppg%v2nlma(vi)
        if (my_nlma < 1) cycle

        wrk = 0d0
!$acc loop seq
        do k=1,my_nlma
          ilma = ppg%k2ilma(vi,k)
          j    = ppg%k2j(vi,k)
          wrk  = wrk + ppg%uVpsibox(ilma,ispin,io,ik,im) * ppg%zekr_uV(j,ilma,ik)
        end do

        ix = ppg%v2j(1,vi)
        iy = ppg%v2j(2,vi)
        iz = ppg%v2j(3,vi)
        htpsi%zwf(ix,iy,iz,ispin,io,ik,im) = htpsi%zwf(ix,iy,iz,ispin,io,ik,im) + wrk
      end do

    end do
    end do
    end do
    end do
!$acc end kernels
#endif
#else
!$omp parallel do collapse(4) &
!$omp             private(im,ik,io,ispin,ilma,ia,uVpsi,j,ix,iy,iz,wrk)
    do im=im_s,im_e
    do ik=ik_s,ik_e
    do io=io_s,io_e
    do ispin=1,Nspin

      do ilma=1,Nlma
        ia = ppg%ia_tbl(ilma)
        uVpsi = 0.d0
        do j=1,ppg%mps(ia)
          ix = ppg%jxyz(1,j,ia)
          iy = ppg%jxyz(2,j,ia)
          iz = ppg%jxyz(3,j,ia)
          uVpsi = uVpsi + conjg(ppg%zekr_uV(j,ilma,ik)) * tpsi%zwf(ix,iy,iz,ispin,io,ik,im)
        end do
        uVpsi = uVpsi * ppg%rinv_uvu(ilma)
!OCL norecurrence
        do j=1,ppg%mps(ia)
          ix = ppg%jxyz(1,j,ia)
          iy = ppg%jxyz(2,j,ia)
          iz = ppg%jxyz(3,j,ia)
          wrk = uVpsi * ppg%zekr_uV(j,ilma,ik)
          htpsi%zwf(ix,iy,iz,ispin,io,ik,im) = htpsi%zwf(ix,iy,iz,ispin,io,ik,im) + wrk
        end do
      end do

    end do
    end do
    end do
    end do
!$omp end parallel do
#endif

  end if

  call timer_end(LOG_UHPSI_PSEUDO)

  return
end subroutine zpseudo

subroutine calc_uVpsi(nspin,info,ppg,tpsi,uVpsibox)
  use structures
  use timer
  implicit none
  integer        ,intent(in) :: nspin
  type(s_parallel_info),intent(in) :: info
  type(s_pp_grid),intent(in) :: ppg
  type(s_orbital),intent(in) :: tpsi
  complex(8)    ,allocatable :: uVpsibox (:,:,:,:,:)
  integer :: ispin,io,ik,im,im_s,im_e,ik_s,ik_e,io_s,io_e,norb
  integer :: ilma,ia,j,ix,iy,iz,Nlma
  complex(8) :: uVpsi
  integer :: ilocal

  call timer_begin(LOG_UHPSI_PSEUDO)

  im_s = info%im_s
  im_e = info%im_e
  ik_s = info%ik_s
  ik_e = info%ik_e
  io_s = info%io_s
  io_e = info%io_e
  norb = Nspin* info%numo * info%numk * info%numm

  Nlma = ppg%Nlma

  allocate(uVpsibox(Nspin,io_s:io_e,ik_s:ik_e,im_s:im_e,Nlma))

#ifdef USE_OPENACC
!$acc parallel loop collapse(4) private(im,ik,io,ispin,ilocal,ilma,ia,uVpsi,j,ix,iy,iz)
#else
!$omp parallel do collapse(4) &
!$omp             private(im,ik,io,ispin,ilocal,ilma,ia,uVpsi,j,ix,iy,iz)
#endif
  do im=im_s,im_e
  do ik=ik_s,ik_e
  do io=io_s,io_e
  do ispin=1,Nspin

    do ilocal=1,ppg%ilocal_nlma
      ilma=ppg%ilocal_nlma2ilma(ilocal)
      ia  =ppg%ilocal_nlma2ia  (ilocal)
      uVpsi = 0.d0
      do j=1,ppg%mps(ia)
        ix = ppg%jxyz(1,j,ia)
        iy = ppg%jxyz(2,j,ia)
        iz = ppg%jxyz(3,j,ia)
        uVpsi = uVpsi + conjg(ppg%zekr_uV(j,ilma,ik)) * tpsi%zwf(ix,iy,iz,ispin,io,ik,im)
      end do
      uVpsi = uVpsi * ppg%rinv_uvu(ilma)
      uVpsibox(ispin,io,ik,im,ilma) = uVpsi
    end do

  end do
  end do
  end do
  end do
#ifdef USE_OPENACC
!$acc end parallel
#else
!$omp end parallel do
#endif

  call timer_end(LOG_UHPSI_PSEUDO)

  return
end subroutine calc_uVpsi

subroutine calc_uVpsi_rdivided(nspin,info,ppg,tpsi,uVpsibox,uVpsibox2)
  use structures
  use timer
#ifdef FORTRAN_COMPILER_HAS_MPI_VERSION3
  use salmon_global, only: natom
  use communication, only: comm_wait_all,comm_show_error
  use mpi, only: MPI_SUM,MPI_DOUBLE_COMPLEX
#else
  use communication, only: comm_summation
#endif
  implicit none
  integer        ,intent(in) :: nspin
  type(s_parallel_info),intent(in) :: info
  type(s_pp_grid),intent(in) :: ppg
  type(s_orbital),intent(in) :: tpsi
  complex(8)    ,allocatable :: uVpsibox (:,:,:,:,:)
  complex(8)    ,allocatable :: uVpsibox2(:,:,:,:,:)
  integer :: ispin,io,ik,im,im_s,im_e,ik_s,ik_e,io_s,io_e,norb
  integer :: ilma,ia,j,ix,iy,iz,Nlma,ilocal
  complex(8) :: uVpsi
#ifdef FORTRAN_COMPILER_HAS_MPI_VERSION3
  integer :: ireqs(natom),nreq,ierr,is,ie,ns
#endif

  call timer_begin(LOG_UHPSI_PSEUDO)

  im_s = info%im_s
  im_e = info%im_e
  ik_s = info%ik_s
  ik_e = info%ik_e
  io_s = info%io_s
  io_e = info%io_e
  norb = Nspin* info%numo * info%numk * info%numm

  Nlma = ppg%Nlma

  allocate(uVpsibox (Nspin,io_s:io_e,ik_s:ik_e,im_s:im_e,Nlma))
  allocate(uVpsibox2(Nspin,io_s:io_e,ik_s:ik_e,im_s:im_e,Nlma))

#ifdef FORTRAN_COMPILER_HAS_MPI_VERSION3
  uVpsibox2 = 0d0
#else
  uVpsibox  = 0d0
#endif

!$omp parallel do collapse(4) &
!$omp             private(im,ik,io,ispin,ilocal,ilma,ia,uVpsi,j,ix,iy,iz)
  do im=im_s,im_e
  do ik=ik_s,ik_e
  do io=io_s,io_e
  do ispin=1,Nspin

    do ilocal=1,ppg%ilocal_nlma
      ilma=ppg%ilocal_nlma2ilma(ilocal)
      ia  =ppg%ilocal_nlma2ia  (ilocal)
      uVpsi = 0.d0
      do j=1,ppg%mps(ia)
        ix = ppg%jxyz(1,j,ia)
        iy = ppg%jxyz(2,j,ia)
        iz = ppg%jxyz(3,j,ia)
        uVpsi = uVpsi + conjg(ppg%zekr_uV(j,ilma,ik)) * tpsi%zwf(ix,iy,iz,ispin,io,ik,im)
      end do
      uVpsi = uVpsi * ppg%rinv_uvu(ilma)
      uVpsibox(ispin,io,ik,im,ilma) = uVpsi
    end do

  end do
  end do
  end do
  end do
!$omp end parallel do

  call timer_end(LOG_UHPSI_PSEUDO)

  call timer_begin(LOG_UHPSI_PSEUDO_COMM)
#ifdef FORTRAN_COMPILER_HAS_MPI_VERSION3
! FIXME: This subroutine uses MPI functions directly...
  nreq = 0
  do ia=1,natom
    is = ppg%irange_atom(1,ia)
    ie = ppg%irange_atom(2,ia)
    ns = ie - is + 1

    if (ppg%ireferred_atom(ia)) then
      nreq = nreq + 1
      call MPI_Iallreduce( uvpsibox (1,io_s,ik_s,im_s,is) &
                         , uvpsibox2(1,io_s,ik_s,im_s,is) &
                         , ns*norb, MPI_DOUBLE_COMPLEX, MPI_SUM, ppg%icomm_atom(ia) &
                         , ireqs(nreq), ierr )
      call comm_show_error(ierr)
    !else
      ! uvpsibox2(:,:,:,:,ppg%irange_ia(1:2,ia)) does not use in this process...
      ! We can skip self copy, but zero clear required
    end if
  end do
  call comm_wait_all(ireqs(1:nreq))
#else
  call comm_summation(uVpsibox,uVpsibox2,Nlma*Norb,info%icomm_r)
#endif
  call timer_end(LOG_UHPSI_PSEUDO_COMM)

  return
end subroutine calc_uVpsi_rdivided

end module nonlocal_potential
