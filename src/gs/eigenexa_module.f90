!
!  Copyright 2020 SALMON developers
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
    type(s_parallel_info),intent(inout) :: info
    integer,intent(in) :: n

    integer :: npo, i, j, ip
    integer,allocatable :: icount(:)
    integer :: i_loc,j_loc,prow,pcol

    if (info%flag_eigenexa_init) return

    info%icomm_sl = info%icomm_o

    call eigen_init(info%icomm_sl)
    call eigen_get_procs(info%nprocs, info%nprow, info%npcol)
    call eigen_get_id   (info%iam,    info%myrow, info%mycol)

    call eigen_get_matdims(n, info%nrow_local, info%ncol_local)

    info%flag_eigenexa_init = .true.

    ! --- for reduce memory
    ! MEMO: ScaLAPACK is block-cyclic distribution of matrix
    !       EigenExa  is       cyclic distribution of matrix
    npo = info%nprocs
    allocate( info%ndiv(0:npo-1), icount(0:npo-1) )

    icount = 0
    do i = 1,n !row
    do j = 1,n !col
      prow = eigen_owner_node(i, info%nprow, info%myrow)
      pcol = eigen_owner_node(j, info%npcol, info%mycol)
      if (info%myrow == prow .and. info%mycol == pcol) then
        ip = info%irank_io(j)
        icount( ip ) = icount( ip ) + 1
      end if
    end do
    end do
    info%ndiv = icount

    allocate( info%i_tbl( maxval(info%ndiv), 0:npo-1), &
              info%j_tbl( maxval(info%ndiv), 0:npo-1), &
              info%iloc_tbl( maxval(info%ndiv), 0:npo-1), &
              info%jloc_tbl( maxval(info%ndiv), 0:npo-1) )

    icount(:) = 0
    info%i_tbl = 0
    info%j_tbl = 0
    info%iloc_tbl = 0
    info%jloc_tbl = 0

    do i = 1,n !row
    do j = 1,n !col
      i_loc = eigen_translate_g2l(i, info%nprow, info%myrow)
      j_loc = eigen_translate_g2l(j, info%npcol, info%mycol)
      prow = eigen_owner_node(i, info%nprow, info%myrow)
      pcol = eigen_owner_node(j, info%npcol, info%mycol)
      if (info%myrow == prow .and. info%mycol == pcol) then
        ip = info%irank_io(j)
        icount( ip ) = icount( ip ) + 1
        info%i_tbl( icount(ip), ip ) = i
        info%j_tbl( icount(ip), ip ) = j
        info%iloc_tbl( icount(ip), ip ) = i_loc
        info%jloc_tbl( icount(ip), ip ) = j_loc
      end if
    end do
    end do

    deallocate(icount)

    return

  end subroutine init_eigenexa

end module eigenexa_module
