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
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130

module density_matrix_plusU_sub

  use plusU_global, only: dm_plusU

  implicit none

  private
  public :: density_matrix_plusU

contains

  subroutine density_matrix_plusU( psi, ppg, occ )
    use structures, only: s_orbital, s_pp_grid
    implicit none
    type(s_orbital),intent(in) :: psi
    type(s_pp_grid),intent(in) :: ppg
    real(8),intent(in) :: occ(:,:,:,:)
    integer :: Nlma,Nspin
    integer :: im,ik,io,ispin,ilma,jlma,iprj,Nproj_pairs
    integer :: ix,iy,iz,j,ia,io_s,io_e,ik_s,ik_e,im_s,im_e
    complex(8) :: phipsi
    complex(8),allocatable :: phipsi_lma(:)
    complex(8),parameter :: zero=(0.0d0,0.0d0)

    Nspin=size(occ,2)
    Nproj_pairs=size(ppg%proj_pairs_ao,2)
    Nlma=size(ppg%ia_tbl_ao)

    if ( .not.allocated(dm_plusU) ) then
      allocate( dm_plusU(Nproj_pairs) )
      dm_plusU=0.0d0
    end if

    allocate( phipsi_lma(Nlma) ); phipsi_lma=zero

    do im=im_s,im_e
    do ik=ik_s,ik_e
    do io=io_s,io_e
    do ispin=1,Nspin

      if ( abs(occ(io,ik,ispin,im)) < 1.d-10 ) cycle

      do ilma=1,Nlma
        ia = ppg%ia_tbl_ao(ilma)
        phipsi = zero
        do j=1,ppg%mps_ao(ia)
          ix = ppg%jxyz_ao(1,j,ia)
          iy = ppg%jxyz_ao(2,j,ia)
          iz = ppg%jxyz_ao(3,j,ia)
          phipsi = phipsi + conjg( ppg%zekr_phi_ao(j,ilma,ik) ) &
                          * psi%zwf(ix,iy,iz,ispin,io,ik,im)
        end do
        phipsi_lma(ilma) = phipsi
      end do !ilma

      do iprj=1,Nproj_pairs
        ilma = ppg%proj_pairs_ao(1,iprj)
        jlma = ppg%proj_pairs_ao(2,iprj)
        dm_plusU(iprj) = dm_plusU(iprj) + occ(io,ik,ispin,im) &
             *phipsi_lma(ilma)*phipsi_lma(jlma)
      end do

    end do !ispin
    end do !io
    end do !ik
    end do !im

    return
  end subroutine density_matrix_plusU

end module density_matrix_plusU_sub
