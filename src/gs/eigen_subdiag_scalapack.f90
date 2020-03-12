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
module eigen_scalapack
  implicit none

  public :: eigen_pdsyevd

contains

  subroutine eigen_pdsyevd(pinfo,info,h,e,v)
    use structures, only: s_process_info, s_orbital_parallel
    use communication, only: comm_summation, comm_is_root
    implicit none
    type(s_process_info),intent(in)     :: pinfo
    type(s_orbital_parallel),intent(in) :: info
    real(8), intent(in)  :: h(:,:)
    real(8), intent(out) :: e(:)
    real(8), intent(out) :: v(:,:)
    integer :: n
    integer :: len_iwork
    integer :: i, j, ierr
    integer :: i_loc, j_loc, proc_row, proc_col
    integer, allocatable :: iwork(:)
    real(8), allocatable :: work(:), h_div(:,:), v_div(:,:), v_tmp(:,:)

    if (.not. pinfo%flag_blacs_gridinit) &
      stop 'eigen_subdiag_scalapack: ScaLAPACK not initialized.'

    n  = ubound(h,1)
    len_iwork = 2 + 7*n + 8*pinfo%npcol

    allocate( h_div(pinfo%nrow_local,pinfo%ncol_local), &
              v_div(pinfo%nrow_local,pinfo%ncol_local), &
              v_tmp(n,n), iwork(len_iwork), work(pinfo%len_work) )

!$omp parallel do private(i,j,i_loc,j_loc,proc_row,proc_col) collapse(2)
    do i=1,n
    do j=1,n
      call INFOG2L( i, j, pinfo%desca, pinfo%nprow, pinfo%npcol, pinfo%myrow, pinfo%mycol, i_loc, j_loc, proc_row, proc_col )
      if (pinfo%myrow == proc_row .and. pinfo%mycol == proc_col) then
        h_div(i_loc,j_loc) = h(i,j)
      end if
    end do
    end do

    call PDSYEVD( 'V', 'L', n, h_div, 1, 1, pinfo%desca, e, v_div, 1, 1, pinfo%descz, &
                  work, pinfo%len_work, iwork, len_iwork, ierr )

    v_tmp=0d0
!$omp parallel do private(i,j,i_loc,j_loc,proc_row,proc_col) collapse(2)
    do i=1,n
    do j=1,n
      call INFOG2L( i, j, pinfo%descz, pinfo%nprow, pinfo%npcol, pinfo%myrow, pinfo%mycol, i_loc, j_loc, proc_row, proc_col )
      if (pinfo%myrow == proc_row .and. pinfo%mycol == proc_col) then
        v_tmp(i,j) = v_div(i_loc,j_loc)
      end if
    end do
    end do

    call comm_summation(v_tmp, v, n*n, info%icomm_r)

    deallocate( work, iwork, h_div, v_div, v_tmp )

    return
  end subroutine eigen_pdsyevd

end module eigen_scalapack
