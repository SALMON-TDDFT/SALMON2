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

  public :: eigen_pdsyevd, eigen_pdsyevd_red_mem

contains

  subroutine eigen_pdsyevd(pinfo,h,e,v)
    use structures, only: s_process_info
    use communication, only: comm_summation, comm_is_root
    implicit none
    type(s_process_info),intent(in)  :: pinfo
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

    call comm_summation(v_tmp, v, n*n, pinfo%icomm_sl)

    deallocate( work, iwork, h_div, v_div, v_tmp )

    return
  end subroutine eigen_pdsyevd

  subroutine eigen_pdsyevd_red_mem(system,pinfo,info,h,e,v)
    use structures, only: s_process_info, s_parallel_info, s_dft_system
    use communication, only: comm_summation, comm_is_root, comm_bcast
    implicit none
    type(s_dft_system),intent(in) :: system
    type(s_process_info),intent(in)     :: pinfo
    type(s_parallel_info),intent(in) :: info
    real(8), intent(in)  :: h(system%no, info%io_s:info%io_e)
    real(8), intent(out) :: e(system%no)
    real(8), intent(out) :: v(system%no, info%io_s:info%io_e)
    integer :: n
    integer :: len_iwork
    integer :: i, j, ierr, m, k
    integer, allocatable :: iwork(:)
    real(8), allocatable :: work(:), h_div(:,:), v_div(:,:), tmp_mat(:,:), tmp_mat2(:,:)

    if (.not. pinfo%flag_blacs_gridinit) &
      stop 'eigen_subdiag_scalapack: ScaLAPACK not initialized.'

    n  = ubound(h,1)
    len_iwork = 2 + 7*n + 8*pinfo%npcol

    allocate( h_div(pinfo%nrow_local,pinfo%ncol_local), &
              v_div(pinfo%nrow_local,pinfo%ncol_local), &
              tmp_mat(system%no, info%numo_max), &
              tmp_mat2(system%no, info%numo_max), &
              iwork(len_iwork), work(pinfo%len_work) )

    do m = 0, pinfo%nporbital - 1
      if(m == info%id_o) then
!$omp parallel do private(i,j) collapse(2)
        do i = 1, system%no
        do j = info%io_s, info%io_e
          tmp_mat(i, j-info%io_s+1) = h(i,j)
        end do
        end do
      end if
      call comm_bcast( tmp_mat(:,1:info%numo_all(m)), info%icomm_o, info%irank_io(info%io_s_all(m)) )

!$omp parallel do private(k)
      do k = 1, pinfo%ndiv(m)
        h_div(pinfo%iloc_tbl(k,m), pinfo%jloc_tbl(k,m)) = tmp_mat(pinfo%i_tbl(k,m), pinfo%j_tbl(k,m)-info%io_s_all(m)+1)
      enddo
    end do !m
!write(*,*) info%id_o, sum(h_div(:,:))

    call PDSYEVD( 'V', 'L', n, h_div, 1, 1, pinfo%desca, e, v_div, 1, 1, pinfo%descz, &
                  work, pinfo%len_work, iwork, len_iwork, ierr )

    do m = 0, pinfo%nporbital - 1
      tmp_mat = 0d0
!$omp parallel do private(k)
      do k = 1, pinfo%ndiv(m)
        tmp_mat(pinfo%i_tbl(k,m), pinfo%j_tbl(k,m)-info%io_s_all(m)+1) = v_div(pinfo%iloc_tbl(k,m), pinfo%jloc_tbl(k,m))
      end do
      call comm_summation(tmp_mat, tmp_mat2, system%no*info%numo_all(m), info%icomm_o)
      if(m == info%id_o) then
!$omp parallel do private(i,j) collapse(2)
        do i = 1, system%no
        do j = info%io_s, info%io_e
          v(i,j) = tmp_mat2(i, j-info%io_s+1)
        end do
        end do
      end if
    end do !m

    deallocate( work, iwork, h_div, v_div, tmp_mat, tmp_mat2 )

    return
  end subroutine eigen_pdsyevd_red_mem


  subroutine eigen_pzheevd(pinfo,h,e,v)
    use structures, only: s_process_info
    use communication, only: comm_summation, comm_is_root
    implicit none
    type(s_process_info),intent(in) :: pinfo
    complex(8), intent(in)  :: h(:,:)
    real(8), intent(out)    :: e(:)
    complex(8), intent(out) :: v(:,:)

    complex(8), allocatable :: h_div(:,:), v_div(:,:), v_tmp(:,:)
    complex(8), allocatable :: work(:)
    real(8), allocatable    :: rwork(:)
    integer, allocatable    :: iwork(:)
    integer :: n, len_iwork
    integer :: i, j, ierr
    integer :: i_loc, j_loc, proc_row, proc_col

    if (.not. pinfo%flag_blacs_gridinit) &
      stop 'eigen_subdiag_scalapack: ScaLAPACK not initialized.'

    n  = ubound(h,1)
    len_iwork = 2 + 7*n + 8*pinfo%npcol

    allocate( h_div(pinfo%nrow_local,pinfo%ncol_local), &
              v_div(pinfo%nrow_local,pinfo%ncol_local), &
              v_tmp(n,n), work(pinfo%len_work), rwork(pinfo%len_rwork), &
              iwork(len_iwork) )

!$omp parallel do private(i,j,i_loc,j_loc,proc_row,proc_col) collapse(2)
    do i=1,n
    do j=1,n
      call INFOG2L( i, j, pinfo%desca, pinfo%nprow, pinfo%npcol, pinfo%myrow, pinfo%mycol, i_loc, j_loc, proc_row, proc_col )
      if (pinfo%myrow == proc_row .and. pinfo%mycol == proc_col) then
        h_div(i_loc,j_loc) = h(i,j)
      end if
    end do
    end do

    call PZHEEVD( 'V', 'L', n, h_div, 1, 1, pinfo%desca, e, v_div, 1, 1, pinfo%descz, &
                  work, pinfo%len_work, rwork, pinfo%len_rwork, iwork, len_iwork, ierr )

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

    call comm_summation(v_tmp, v, n*n, pinfo%icomm_sl)

    deallocate( iwork, rwork, work, v_tmp, v_div, h_div )

    return
  end subroutine eigen_pzheevd

end module eigen_scalapack
