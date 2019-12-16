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
module pseudo_pt_plusU_sub

  use plusU_global, only: PLUS_U_ON, V_eff

  implicit none

! WARNING: We must not call these except for hpsi routine.

  private
  public :: pseudo_plusU
  public :: PLUS_U_ON

contains

!-----------------------------------------------------------------------------------------------------------------------------------

  subroutine pseudo_plusU(tpsi,htpsi,info,nspin,ppg)
    use structures
    use timer
    implicit none
    integer,intent(in) :: nspin
    type(s_orbital_parallel),intent(in) :: info
    type(s_pp_grid),intent(in) :: ppg
    type(s_orbital),intent(in) :: tpsi
    type(s_orbital) :: htpsi
    !
    integer :: ispin,io,ik,im,im_s,im_e,ik_s,ik_e,io_s,io_e
    integer :: ilma,ia,j,ix,iy,iz,Nlma,jlma,iprj,Nproj_pairs
    complex(8),parameter :: zero=(0.0d0,0.0d0)
    complex(8) :: phipsi,wrk
    complex(8),allocatable :: phipsi_lma(:)
    complex(8),allocatable :: phipsibox(:,:,:,:,:)
    complex(8),allocatable :: phipsibox2(:,:,:,:,:)

    if ( .not.allocated(V_eff) ) return

    call timer_begin(LOG_UHPSI_PSEUDO)

    im_s = info%im_s
    im_e = info%im_e
    ik_s = info%ik_s
    ik_e = info%ik_e
    io_s = info%io_s
    io_e = info%io_e

    Nlma = size(ppg%ia_tbl_ao)
    Nproj_pairs = size(ppg%proj_pairs_ao,2)

    if ( info%if_divide_rspace ) then

      call timer_end(LOG_UHPSI_PSEUDO)

      call calc_phipsi_rdivided( nspin, info, ppg, tpsi, phipsibox, phipsibox2 )

      call timer_begin(LOG_UHPSI_PSEUDO)

      do im=im_s,im_e
      do ik=ik_s,ik_e
      do io=io_s,io_e
      do ispin=1,Nspin

        do ilma=1,Nlma

          ia = ppg%ia_tbl_ao(ilma)

          phipsi = phipsibox2(ispin,io,ik,im,ilma)

          do j=1,ppg%mps_ao(ia)
            ix = ppg%jxyz_ao(1,j,ia)
            iy = ppg%jxyz_ao(2,j,ia)
            iz = ppg%jxyz_ao(3,j,ia)
            wrk = phipsi * ppg%zekr_phi_ao(j,ilma,ik)
            htpsi%zwf(ix,iy,iz,ispin,io,ik,im) = htpsi%zwf(ix,iy,iz,ispin,io,ik,im) + wrk
          end do

        end do !ilma

      end do
      end do
      end do
      end do

      deallocate( phipsibox, phipsibox2 )

    else !if ( .not. info%if_divide_rspace ) then

      allocate( phipsi_lma(Nlma) ); phipsi_lma=zero

      do im=im_s,im_e
      do ik=ik_s,ik_e
      do io=io_s,io_e
      do ispin=1,Nspin

        do ilma=1,Nlma
          ia = ppg%ia_tbl_ao(ilma)
          phipsi = zero
          do j=1,ppg%mps_ao(ia)
            ix = ppg%jxyz_ao(1,j,ia)
            iy = ppg%jxyz_ao(2,j,ia)
            iz = ppg%jxyz_ao(3,j,ia)
            phipsi = phipsi + conjg( ppg%zekr_phi_ao(j,ilma,ik) ) &
                            * tpsi%zwf(ix,iy,iz,ispin,io,ik,im)
          end do
          phipsi_lma(ilma) = phipsi * ppg%Hvol
        end do !ilma

        do iprj=1,Nproj_pairs
          ilma = ppg%proj_pairs_ao(1,iprj)
          jlma = ppg%proj_pairs_ao(2,iprj)
          ia   = ppg%ia_tbl_ao(ilma)
          do j=1,ppg%mps_ao(ia)
            ix = ppg%jxyz_ao(1,j,ia)
            iy = ppg%jxyz_ao(2,j,ia)
            iz = ppg%jxyz_ao(3,j,ia)
            wrk = V_eff(iprj,ispin) * ppg%zekr_phi_ao(j,ilma,ik) * phipsi_lma(jlma)
            htpsi%zwf(ix,iy,iz,ispin,io,ik,im) = htpsi%zwf(ix,iy,iz,ispin,io,ik,im) + wrk
          end do !j
        end do !iprj

      end do !ispin
      end do !io
      end do !ik
      end do !im

      deallocate( phipsi_lma )

    end if

    call timer_end(LOG_UHPSI_PSEUDO)

    return
  end subroutine pseudo_plusU

  subroutine calc_phipsi_rdivided(nspin,info,ppg,tpsi,phipsibox,phipsibox2)
    use structures
    use timer
    use communication, only: comm_summation
#ifdef SALMON_ENABLE_MPI3
    use communication, only: comm_wait_all
    use mpi, only: MPI_SUM,MPI_DOUBLE_COMPLEX
#endif
    implicit none
    integer        ,intent(in) :: nspin
    type(s_orbital_parallel),intent(in) :: info
    type(s_pp_grid),intent(in) :: ppg
    type(s_orbital),intent(in) :: tpsi
    complex(8)    ,allocatable :: phipsibox (:,:,:,:,:)
    complex(8)    ,allocatable :: phipsibox2(:,:,:,:,:)
    complex(8),parameter :: zero=(0.0d0,0.0d0)
    integer :: ispin,io,ik,im,im_s,im_e,ik_s,ik_e,io_s,io_e,norb
    integer :: ilma,ia,j,ix,iy,iz,Nlma,iproj
    complex(8) :: phipsi

    call timer_begin(LOG_UHPSI_PSEUDO)

    im_s = info%im_s
    im_e = info%im_e
    ik_s = info%ik_s
    ik_e = info%ik_e
    io_s = info%io_s
    io_e = info%io_e
    norb = Nspin* info%numo * info%numk * info%numm

    Nlma = size(ppg%ia_tbl_ao)

    allocate( phipsibox(Nspin,io_s:io_e,ik_s:ik_e,im_s:im_e,Nlma) )
    allocate( phipsibox2(Nspin,io_s:io_e,ik_s:ik_e,im_s:im_e,Nlma) )

    do im=im_s,im_e
    do ik=ik_s,ik_e
    do io=io_s,io_e
    do ispin=1,Nspin

      do iproj=1,Nlma
        ia = ppg%ia_tbl_ao(ilma)
        phipsi = zero
        do j=1,ppg%mps_ao(ia)
          ix = ppg%jxyz_ao(1,j,ia)
          iy = ppg%jxyz_ao(2,j,ia)
          iz = ppg%jxyz_ao(3,j,ia)
          phipsi = phipsi + conjg( ppg%zekr_phi_ao(j,ilma,ik) ) &
                          * tpsi%zwf(ix,iy,iz,ispin,io,ik,im)
        end do
        phipsi = phipsi * ppg%rinv_uvu(ilma)
        phipsibox(ispin,io,ik,im,ilma) = phipsi
      end do ! iproj

    end do
    end do
    end do
    end do

    call timer_end(LOG_UHPSI_PSEUDO)

    call timer_begin(LOG_UHPSI_PSEUDO_COMM)
#ifdef SALMON_ENABLE_MPI3
! FIXME: This subroutine uses MPI functions directly...
    nreq = 0
    do ia=1,natom
      if ( ppg%ireferred_atom(ia) ) then
        is = ppg%irange_atom(1,ia)
        ie = ppg%irange_atom(2,ia)
        ns = ie - is + 1
        nreq = nreq + 1
        call MPI_Iallreduce( uvpsibox (Nspin,io_s,ik_s,im_s,is) &
                           , uvpsibox2(Nspin,io_s,ik_s,im_s,is) &
                           , ns*norb, MPI_DOUBLE_COMPLEX, MPI_SUM, ppg%icomm_atom(ia) &
                           , ireqs(nreq), ierr )
      end if
    end do
    call comm_wait_all(ireqs(1:nreq))
#else
    call comm_summation(phipsibox,phipsibox2,Nlma*Norb,info%icomm_r)
#endif
    call timer_end(LOG_UHPSI_PSEUDO_COMM)

    return
  end subroutine calc_phipsi_rdivided

end module pseudo_pt_plusU_sub
