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
module eigenexa_module
  implicit none

  public :: init_eigenexa

private

contains
  subroutine init_eigenexa(info,n)
    use structures, only: s_parallel_info
    use communication, only: comm_summation
    use eigen_libs_mod
    implicit none
    type(s_parallel_info),intent(in) :: info
    integer,intent(in) :: n

    integer :: npo, i, j, ip
    integer,allocatable :: icount(:)
    integer :: i_loc,j_loc,prow,pcol

    if (pinfo%flag_eigenexa_init) return

    pinfo%icomm_sl = info%icomm_o

    call eigen_init(pinfo%icomm_sl)
    call eigen_get_procs(pinfo%nprocs, pinfo%nprow, pinfo%npcol)
    call eigen_get_id   (pinfo%iam,    pinfo%myrow, pinfo%mycol)

    call eigen_get_matdims(n, pinfo%nrow_local, pinfo%ncol_local)

    pinfo%flag_eigenexa_init = .true.

    ! --- for reduce memory
    ! MEMO: ScaLAPACK is block-cyclic distribution of matrix
    !       EigenExa  is       cyclic distribution of matrix
    npo = pinfo%nprocs
    allocate( pinfo%ndiv(0:npo-1), icount(0:npo-1) )

    icount = 0
    do i = 1,n !row
    do j = 1,n !col
      prow = eigen_owner_node(i, pinfo%nprow, pinfo%myrow)
      pcol = eigen_owner_node(j, pinfo%npcol, pinfo%mycol)
      if (pinfo%myrow == prow .and. pinfo%mycol == pcol) then
        ip = info%irank_io(j)
        icount( ip ) = icount( ip ) + 1
      end if
    end do
    end do
    pinfo%ndiv = icount

    allocate( pinfo%i_tbl( maxval(pinfo%ndiv), 0:npo-1), &
              pinfo%j_tbl( maxval(pinfo%ndiv), 0:npo-1), &
              pinfo%iloc_tbl( maxval(pinfo%ndiv), 0:npo-1), &
              pinfo%jloc_tbl( maxval(pinfo%ndiv), 0:npo-1) )

    icount(:) = 0
    pinfo%i_tbl = 0
    pinfo%j_tbl = 0
    pinfo%iloc_tbl = 0
    pinfo%jloc_tbl = 0

    do i = 1,n !row
    do j = 1,n !col
      i_loc = eigen_translate_g2l(i, pinfo%nprow, pinfo%myrow)
      j_loc = eigen_translate_g2l(j, pinfo%npcol, pinfo%mycol)
      prow = eigen_owner_node(i, pinfo%nprow, pinfo%myrow)
      pcol = eigen_owner_node(j, pinfo%npcol, pinfo%mycol)
      if (pinfo%myrow == prow .and. pinfo%mycol == pcol) then
        ip = info%irank_io(j)
        icount( ip ) = icount( ip ) + 1
        pinfo%i_tbl( icount(ip), ip ) = i
        pinfo%j_tbl( icount(ip), ip ) = j
        pinfo%iloc_tbl( icount(ip), ip ) = i_loc
        pinfo%jloc_tbl( icount(ip), ip ) = j_loc
      end if
    end do
    end do

    deallocate(icount)

    return

  end subroutine init_eigenexa

end module eigenexa_module
