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
!=======================================================================
module gram_schmidt_orth
  use structures, only: s_rgrid, s_wf_info, s_wavefunction
  implicit none

contains

  subroutine gram_schmidt(sys, rg, wfi, wf)
    use timer
    implicit none
    type(s_system),       intent(in)    :: sys
    type(s_rgrid),        intent(in)    :: rg
    type(s_wf_info),      intent(in)    :: wfi
    type(s_wavefunction), intent(inout) :: wf

    call timer_begin(LOG_GRAM_SCHMIDT)

    !if (if_divide_rspace) then:
    if (allocated(wf%rwf)) then
      call gram_schmidt_col_real8(sys, rg, wfi, wf)
    elseif (allocated(wf%zwf)) then
      call gram_schmidt_col_complex8()
    else
      stop "Wavefunctions are not allocated!"
    end if

    call timer_end(LOG_GRAM_SCHMIDT)

    return
  end subroutine

  subroutine gram_schmidt_col_real8(sys, rg, wfi, rwf)
    ! Only for the colinear L(S)DA:
    use timer
    implicit none
    type(s_system), intent(in) :: sys
    type(s_rgrid),  intent(in) :: rg
    type(s_wf_info),intent(in) :: wfi
    real(8), intent(inout) :: rwf( &
      & rg%is_array(1):rg%ie_array(1), &
      & rg%is_array(2):rg%ie_array(2), &
      & rg%is_array(3):rg%ie_array(3), &
      & 1:sys%nspin, &
      & wfi%io_s:wfi%io_e, &
      & sys%ik_s:sys%ik_e, &
      & sys%im_s:sys%im_e)

    real(8), dimension( &
      & rg%is_array(1):rg%ie_array(1), &
      & rg%is_array(2):rg%ie_array(2), &
      & rg%is_array(3):rg%ie_array(3)) &
      & :: rwf_jo1, s_exc, s_exc_tmp
    
    real(8) :: c_ovlp(1:sys%no), c_ovlp_tmp(1:sys%no)
    
    do ispin = 1, sys%nspin
    do im = sys%im_s, sys%im_e
    do ik = sys%ik_s, sys%ik_e
      ! Loop for each orbit #jo1:
      do jo1 = 1, sys%no

        ! Retrieve orbit #jo1 into `rwf_jo1`:
        if (has_orbit(jo1)) then
          io1 = wfi%jo_tbl(jo1)
          call copy_data( &
            & rwf(:, :, :, ispin, io1, ik, im), &
            & rwf_jo1) 
        end if
        call comm_bcast(rwf_jo1, wfi%icomm_o, wfi%irank_jo(jo1))
        
        ! Calculate overlap coefficients "c(1:jo1-1)": 
        c_ovlp_tmp = 0d0
        do jo2 = 1, jo1 - 1
          if (has_orbit(jo2)) then
            io2 = wfi%jo_tbl(jo2)
            c_ovlp_tmp(jo2) = dot_real8( &
              & tmp1_jo1, &
              & wf%rwf(:, :, :, ispin, io2, ik, im))
          end if
        end do 
        call comm_summation(c_ovlp_tmp, c_ovlp, sys%no, wfi%icomm_ro)

        ! Calculate exclusion term "s_exc":
        s_exc_tmp = 0d0
        do jo2 = 1, jo1 - 1
          if (has_orbit(jo2)) then
            io2 = wfi%jo_tbl(jo2)
            call axpy_real8( &
              & c_ovlp(io2), &
              & wf%rwf(:, :, :, ispin, io2, ik, im), &
              & s_exc_tmp)
          end if
        end do 
        call comm_summation(s_exc_tmp, s_exc, &
          & rg%num(1) * rg%num(2) * rg%num(3), wfi%icomm_o)

        if (has_orbit(jo1)) then
          ! Exclude non-orthonormal component:
          call axpy_real8(-1d0, s_exc, rwf_jo1)
          ! Normalization:
          norm2_tmp = dot_real8(rwf_jo1, rwf_jo1)
          call comm_summation(norm2_tmp, norm2, 1, wfi%icomm_r)
          call scal_real8(1d0 / sqrt(norm2), rwf_jo1)
          ! Write back to "rwf":
          io1 = wfi%jo_tbl(jo1)
          call copy_data( &
            & rwf_jo1,
            & rwf(:, :, :, ispin, io1, ik, im)) 
        end if
      end do !jo1
    
    end do !ik
    end do !im
    end do !ispin
    
    return
  contains

  logical function has_orbit(jo) result(f)
    implicit none
    integer, intent(in) :: jo
    f = (0 <= wfi%jo_tbl(jo))
    return
  end function has_orbit

  ! Dot product of two wavefunctions:
  real(8) function dot_real8(x, y) result(p)
    implicit none
    real(8), intent(in) :: x( &
      & rg%is_array(1):rg%ie_array(1), &
      & rg%is_array(2):rg%ie_array(2), &
      & rg%is_array(3):rg%ie_array(3))
    real(8), intent(in) :: y( &
      & rg%is_array(1):rg%ie_array(1), &
      & rg%is_array(2):rg%ie_array(2), &
      & rg%is_array(3):rg%ie_array(3))
    integer :: i1, i2, i3
    p = 0d0
    !$omp parallel do collapse(2) default(shared) private(i1,i2,i3) reduction(+:p)
    do i3 = rg%is(3), rg%ie(3)
      do i2 = rg%is(2), rg%ie(2)
        do i1 = rg%is(1), rg%ie(1)
          p = p + x(i1, i2, i3) * y(i1, i2, i3)
        end do
      end do
    end do
    !$omp end parallel do
    p = p * sys%hvol
    return
  end function dot_real8

  ! Constant times a wavefunction plus a wavefunction:
  subroutine axpy_real8(a, x, y)
    implicit none
    real(8), intent(in) :: a
    real(8), intent(in) :: x( &
      & rg%is_array(1):rg%ie_array(1), &
      & rg%is_array(2):rg%ie_array(2), &
      & rg%is_array(3):rg%ie_array(3))
    real(8), intent(inout) :: y( &
      & rg%is_array(1):rg%ie_array(1), &
      & rg%is_array(2):rg%ie_array(2), &
      & rg%is_array(3):rg%ie_array(3))
    integer :: i1, i2, i3
    !$omp parallel do collapse(2) default(shared) private(i1,i2,i3)
    do i3 = rg%is_overlap(3), rg%ie_overlap(3)
      do i2 = rg%is_overlap(2), rg%ie_overlap(2)
        do i1 = rg%is_overlap(1), rg%ie_overlap(1)
          y(i1, i2, i3) = a * x(i1, i2, i3) + y(i1, i2, i3)
        end do
      end do
    end do
    !$omp end parallel do
    return
  end subroutine axpy_real8


  ! Scales a wavefunction by a constant 
  subroutine scal_real8(a, x)
    implicit none
    real(8), intent(in) :: a
    real(8), intent(inout) :: x( &
      & rg%is_array(1):rg%ie_array(1), &
      & rg%is_array(2):rg%ie_array(2), &
      & rg%is_array(3):rg%ie_array(3))
    integer :: i1, i2, i3
    !$omp parallel do collapse(2) default(shared) private(i1,i2,i3)
    do i3 = rg%is_overlap(3), rg%ie_overlap(3)
      do i2 = rg%is_overlap(2), rg%ie_overlap(2)
        do i1 = rg%is_overlap(1), rg%ie_overlap(1)
          x(i1, i2, i3) = a * x(i1, i2, i3)
        end do
      end do
    end do
    !$omp end parallel do
    return
  end subroutine scal_real8



  end subroutine gram_schmidt_col_real8





end module gram_schmidt_orth
