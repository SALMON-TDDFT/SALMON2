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
module pseudo_pt_so_sub

  use spin_orbit_global, only: SPIN_ORBIT_ON

  implicit none

! WARNING: We must not call these except for hpsi routine.

  private
  public :: pseudo_so
  public :: SPIN_ORBIT_ON

contains

!-----------------------------------------------------------------------------------------------------------------------------------

  subroutine pseudo_so(tpsi,htpsi,info,nspin,ppg)
    use structures
    use timer
    implicit none
    integer,intent(in) :: nspin
    type(s_parallel_info),intent(in) :: info
    type(s_pp_grid),intent(in) :: ppg
    type(s_orbital),intent(in) :: tpsi
    type(s_orbital) :: htpsi
    !
    integer :: ispin,io,ik,im,im_s,im_e,ik_s,ik_e,io_s,io_e
    integer :: ilma,ia,j,ix,iy,iz,Nlma
    complex(8) :: uVpsi(2,2),wrk
    complex(8),allocatable :: uVpsibox (:,:,:,:,:)
    complex(8),allocatable :: uVpsibox2(:,:,:,:,:)

    !write(*,*) "------------------ pseudo_so"

    call timer_begin(LOG_UHPSI_PSEUDO)

    im_s = info%im_s
    im_e = info%im_e
    ik_s = info%ik_s
    ik_e = info%ik_e
    io_s = info%io_s
    io_e = info%io_e

    Nlma = size(ppg%ia_tbl_so)

    if ( info%if_divide_rspace ) then

      call timer_end(LOG_UHPSI_PSEUDO)

      call calc_uVpsi_rdivided( nspin, info, ppg, tpsi, uVpsibox, uVpsibox2 )

      call timer_begin(LOG_UHPSI_PSEUDO)

      do im=im_s,im_e
      do ik=ik_s,ik_e
      do io=io_s,io_e
      do ispin=1,Nspin

        do ilma=1,Nlma

          ia = ppg%ia_tbl_so(ilma)

          uVpsi(1,1) = uVpsibox2(ispin,io,ik,im,ilma)

          do j=1,ppg%mps(ia)
            ix = ppg%jxyz(1,j,ia)
            iy = ppg%jxyz(2,j,ia)
            iz = ppg%jxyz(3,j,ia)
            wrk = uVpsi(1,1) * ppg%zekr_uV(j,ilma,ik)
            htpsi%zwf(ix,iy,iz,ispin,io,ik,im) = htpsi%zwf(ix,iy,iz,ispin,io,ik,im) + wrk
          end do

        end do !ilma

      end do
      end do
      end do
      end do

      deallocate( uVpsibox, uVpsibox2 )

    else !if ( .not. info%if_divide_rspace ) then

      do im=im_s,im_e
      do ik=ik_s,ik_e
      do io=io_s,io_e

        do ilma=1,Nlma

          ia = ppg%ia_tbl_so(ilma)

          uVpsi = (0.0d0,0.0d0)

          do ispin=1,2

            do j=1,ppg%mps(ia)
              ix = ppg%jxyz(1,j,ia)
              iy = ppg%jxyz(2,j,ia)
              iz = ppg%jxyz(3,j,ia)
              uVpsi(ispin,1) = uVpsi(ispin,1) &
                   + conjg( ppg%zekr_uV_so(j,ilma,ik,ispin,1) ) &
                   * tpsi%zwf(ix,iy,iz,ispin,io,ik,im)
            end do

          end do !ispin

          uVpsi = uVpsi * ppg%rinv_uvu_so(ilma)

          do j=1,ppg%mps(ia)

            ix = ppg%jxyz(1,j,ia)
            iy = ppg%jxyz(2,j,ia)
            iz = ppg%jxyz(3,j,ia)

            wrk = ppg%zekr_uV_so(j,ilma,ik,1,1)*( uVpsi(1,1) + uVpsi(2,1) )

            htpsi%zwf(ix,iy,iz,1,io,ik,im) = htpsi%zwf(ix,iy,iz,1,io,ik,im) + wrk

            wrk = ppg%zekr_uV_so(j,ilma,ik,2,1)*( uVpsi(1,1) + uVpsi(2,1) )

            htpsi%zwf(ix,iy,iz,2,io,ik,im) = htpsi%zwf(ix,iy,iz,2,io,ik,im) + wrk

          end do !j

        end do !ilma

      end do !io
      end do !ik
      end do !im

    end if

    call timer_end(LOG_UHPSI_PSEUDO)

    return
  end subroutine pseudo_so

  subroutine calc_uVpsi_rdivided(nspin,info,ppg,tpsi,uVpsibox,uVpsibox2)
  use structures
  use timer
  use communication, only: comm_summation
#ifdef SALMON_ENABLE_MPI3
  use communication, only: comm_wait_all
  use mpi, only: MPI_SUM,MPI_DOUBLE_COMPLEX
#endif
  implicit none
  integer        ,intent(in) :: nspin
  type(s_parallel_info),intent(in) :: info
  type(s_pp_grid),intent(in) :: ppg
  type(s_orbital),intent(in) :: tpsi
  complex(8)    ,allocatable :: uVpsibox (:,:,:,:,:)
  complex(8)    ,allocatable :: uVpsibox2(:,:,:,:,:)
  integer :: ispin,io,ik,im,im_s,im_e,ik_s,ik_e,io_s,io_e,norb
  integer :: ilma,ia,j,ix,iy,iz,Nlma
  complex(8) :: uVpsi

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
#ifdef SALMON_ENABLE_MPI3
! FIXME: This subroutine uses MPI functions directly...
  nreq = 0
  do ia=1,natom
    if (ppg%ireferred_atom(ia)) then
      is = ppg%irange_atom(1,ia)
      ie = ppg%irange_atom(2,ia)
      ns = ie - is + 1
      nreq = nreq + 1
      call MPI_Iallreduce( uvpsibox (Nspin,io_s,ik_s,im_s,is) &
                         , uvpsibox2(Nspin,io_s,ik_s,im_s,is) &
                         , ns*norb, MPI_DOUBLE_COMPLEX, MPI_SUM, ppg%icomm_atom(ia) &
                         , ireqs(nreq), ierr )
    !else
      ! uvpsibox2(:,:,:,:,ppg%irange_ia(1:2,ia)) does not use in this process...
      ! We can skip self copy.
    end if
  end do
  call comm_wait_all(ireqs(1:nreq))
#else
  call comm_summation(uVpsibox,uVpsibox2,Nlma*Norb,info%icomm_r)
#endif
  call timer_end(LOG_UHPSI_PSEUDO_COMM)

  return
  end subroutine calc_uVpsi_rdivided

end module pseudo_pt_so_sub
