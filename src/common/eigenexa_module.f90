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
  subroutine init_eigenexa(pinfo,info,n)
    use structures, only: s_process_info, s_parallel_info
    use communication, only: comm_summation
    use eigen_libs_mod
    implicit none
    type(s_process_info),intent(inout)  :: pinfo
    type(s_parallel_info),intent(in) :: info
    integer,intent(in) :: n

    integer :: npo, i, j, ip
    integer,allocatable :: icount(:)
    integer :: is(2),ie(2)
    integer :: irow,icol

    if (pinfo%flag_eigenexa_init) return

    call eigen_init(info%icomm_o)
    call eigen_get_procs(pinfo%nprocs, pinfo%nprow, pinfo%npcol)
    call eigen_get_id   (pinfo%iam,    pinfo%myrow, pinfo%mycol)

    call eigen_get_matdims(n, pinfo%nrow_local, pinfo%ncol_local)

    pinfo%icomm_sl = info%icomm_o

    ! --- for reduce memory
    npo = pinfo%nprocs
    allocate( pinfo%ndiv(0:npo-1), icount(0:npo-1) )

    is(2) = eigen_loop_start(1, pinfo%npcol, pinfo%mycol)
    ie(2) = eigen_loop_end  (n, pinfo%npcol, pinfo%mycol)
    is(1) = eigen_loop_start(1, pinfo%nprow, pinfo%myrow)
    ie(1) = eigen_loop_end  (n, pinfo%nprow, pinfo%myrow)

    icount = 0
    do j = is(2),ie(2)
    do i = is(1),ie(1)
      icol = eigen_translate_l2g(j, pinfo%npcol, pinfo%mycol)
      irow = eigen_translate_l2g(i, pinfo%nprow, pinfo%myrow)
      ip = info%irank_io(icol)
      icount( ip ) = icount( ip ) + 1
    end do
    end do
    pinfo%ndiv = icount

    allocate( pinfo%i_tbl( max(1,pinfo%ndiv(info%id_o)), 0:npo-1), &
              pinfo%j_tbl( max(1,pinfo%ndiv(info%id_o)), 0:npo-1), &
              pinfo%iloc_tbl( max(1,pinfo%ndiv(info%id_o)), 0:npo-1), &
              pinfo%jloc_tbl( max(1,pinfo%ndiv(info%id_o)), 0:npo-1) )

    icount(:) = 0
    pinfo%i_tbl = 0
    pinfo%j_tbl = 0
    pinfo%iloc_tbl = 0
    pinfo%jloc_tbl = 0

    do j = is(2),ie(2)
    do i = is(1),ie(1)
      icol = eigen_translate_l2g(j, pinfo%npcol, pinfo%mycol)
      irow = eigen_translate_l2g(i, pinfo%nprow, pinfo%myrow)
      ip = info%irank_io(icol)
      icount( ip ) = icount( ip ) + 1
      pinfo%i_tbl( icount(ip), ip ) = irow
      pinfo%j_tbl( icount(ip), ip ) = icol
      pinfo%iloc_tbl( icount(ip), ip ) = i
      pinfo%jloc_tbl( icount(ip), ip ) = j
    end do
    end do

    deallocate(icount)

    pinfo%flag_eigenexa_init = .true.
    return

  end subroutine init_eigenexa

end module eigenexa_module
