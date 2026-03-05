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
module eigen_eigenexa
  implicit none

  public :: eigen_pdsyevd_ex, eigen_pdsyevd_ex_red_mem

contains

  subroutine eigen_pdsyevd_ex(info,h,e,v)
    use structures, only: s_parallel_info
    use communication, only: comm_summation, comm_is_root
    use eigen_libs_mod
    implicit none
    type(s_parallel_info),intent(in) :: info
    real(8), intent(in)  :: h(:,:)
    real(8), intent(out) :: e(:)
    real(8), intent(out) :: v(:,:)
    integer :: n
    integer :: i, j, i_loc, j_loc
    integer :: is(2),ie(2)
    real(8), allocatable :: h_div(:,:), v_div(:,:), v_tmp(:,:)

    if (.not. info%flag_eigenexa_init) &
      stop 'eigen_subdiag_eigenexa: EigenExa not initialized.'

    n  = ubound(h,1)
    is(2) = eigen_loop_start(1, info%npcol, info%mycol)
    ie(2) = eigen_loop_end  (n, info%npcol, info%mycol)
    is(1) = eigen_loop_start(1, info%nprow, info%myrow)
    ie(1) = eigen_loop_end  (n, info%nprow, info%myrow)

    allocate( h_div(info%nrow_local,info%ncol_local), &
              v_div(info%nrow_local,info%ncol_local), &
              v_tmp(n,n) )

!$omp parallel do private(i,j,i_loc,j_loc) collapse(2)
    do j_loc=is(2),ie(2)
    do i_loc=is(1),ie(1)
      j = eigen_translate_l2g(j_loc, info%npcol, info%mycol)
      i = eigen_translate_l2g(i_loc, info%nprow, info%myrow)
      h_div(i_loc,j_loc) = h(i,j)
    end do
    end do

    call eigen_sx(n, n, h_div, info%nrow_local, e, v_div, info%nrow_local)

    v_tmp=0d0
!$omp parallel do private(i,j,i_loc,j_loc) collapse(2)
    do j_loc=is(2),ie(2)
    do i_loc=is(1),ie(1)
      j = eigen_translate_l2g(j_loc, info%npcol, info%mycol)
      i = eigen_translate_l2g(i_loc, info%nprow, info%myrow)
      v_tmp(i,j) = v_div(i_loc,j_loc)
    end do
    end do

    call comm_summation(v_tmp, v, size(v_tmp), info%icomm_sl)

    deallocate( h_div, v_div, v_tmp )

    return
  end subroutine eigen_pdsyevd_ex

  subroutine eigen_pdsyevd_ex_red_mem(system,info,h,e,v)
    use structures, only: s_parallel_info, s_dft_system
    use communication, only: comm_summation, comm_is_root, comm_bcast
    use eigen_libs_mod
    implicit none
    type(s_dft_system)   ,intent(in) :: system
    type(s_parallel_info),intent(in) :: info
    real(8), intent(in)  :: h(system%no, info%io_s:info%io_e)
    real(8), intent(out) :: e(system%no)
    real(8), intent(out) :: v(system%no, info%io_s:info%io_e)
    integer :: n
    integer :: i, j, m, k
    real(8), allocatable :: h_div(:,:), v_div(:,:), tmp_mat(:,:), tmp_mat2(:,:)
    integer :: is(2),ie(2)

    if (.not. info%flag_eigenexa_init) &
      stop 'eigen_subdiag_eigenexa: EigenExa not initialized.'

    n  = ubound(h,1)
    is(2) = eigen_loop_start(1, info%npcol, info%mycol)
    ie(2) = eigen_loop_end  (n, info%npcol, info%mycol)
    is(1) = eigen_loop_start(1, info%nprow, info%myrow)
    ie(1) = eigen_loop_end  (n, info%nprow, info%myrow)

    allocate( h_div(info%nrow_local,info%ncol_local), &
              v_div(info%nrow_local,info%ncol_local), &
              tmp_mat(system%no, info%numo_max), &
              tmp_mat2(system%no, info%numo_max) )

    do m = 0, info%nporbital - 1
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
      do k = 1, info%ndiv(m)
        h_div(info%iloc_tbl(k,m), info%jloc_tbl(k,m)) = tmp_mat(info%i_tbl(k,m), info%j_tbl(k,m)-info%io_s_all(m)+1)
      enddo
    end do !m

    call eigen_sx(n, n, h_div, info%nrow_local, e, v_div, info%nrow_local)

    do m = 0, info%nporbital - 1
      tmp_mat = 0d0
!$omp parallel do private(k)
      do k = 1, info%ndiv(m)
        tmp_mat(info%i_tbl(k,m), info%j_tbl(k,m)-info%io_s_all(m)+1) = v_div(info%iloc_tbl(k,m), info%jloc_tbl(k,m))
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

    deallocate( h_div, v_div, tmp_mat, tmp_mat2 )

    return
  end subroutine eigen_pdsyevd_ex_red_mem

end module eigen_eigenexa
