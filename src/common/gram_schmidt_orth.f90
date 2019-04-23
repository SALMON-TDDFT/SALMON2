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
  use structures, only: s_system, s_rgrid, s_wf_info, s_wavefunction
  use pack_unpack, only: copy_data
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
      call gram_schmidt_col_complex8(sys, rg, wfi, wf)
    else
      stop "Wavefunctions are not allocated!"
    end if

    call timer_end(LOG_GRAM_SCHMIDT)

    return
  end subroutine

  subroutine gram_schmidt_col_real8(sys, rg, wfi, wf)
    ! Only for the colinear L(S)DA:
    use timer
    use salmon_communication, only: comm_bcast, comm_summation
    implicit none
    type(s_system),       intent(in)    :: sys
    type(s_rgrid),        intent(in)    :: rg
    type(s_wf_info),      intent(in)    :: wfi
    type(s_wavefunction), intent(inout) :: wf

    integer :: nsize_rg
    integer :: ik, im, ispin
    integer :: jo1, jo2, io1, io2
    real(8) :: coeff(1:sys%no), coeff_tmp(1:sys%no)
    real(8) :: norm2, norm2_tmp
    real(8), dimension( &
      & rg%is_array(1):rg%ie_array(1), &
      & rg%is_array(2):rg%ie_array(2), &
      & rg%is_array(3):rg%ie_array(3)) &
      & :: wf_jo1, wf_exc, wf_exc_tmp

    nsize_rg =  (rg%ie_array(1) - rg%is_array(1)) &
      & * (rg%ie_array(2) - rg%is_array(2)) &
      & * (rg%ie_array(3) - rg%is_array(3))
      
    do im = wfi%im_s, wfi%im_e
    do ik = wfi%ik_s, wfi%ik_e
    do ispin = 1, sys%nspin

      ! Loop for all orbit #jo1:
      do jo1 = 1, sys%no

        ! Retrieve orbit #jo1 into "wf_jo1":
        if (has_orbit(jo1)) then
          io1 = wfi%jo_tbl(jo1)
          call copy_data( &
            & wf%rwf(:, :, :, ispin, io1, ik, im), &
            & wf_jo1)
        end if
        call comm_bcast(wf_jo1, wfi%icomm_o, wfi%irank_jo(jo1))
        
        ! Calculate overlap coefficients: 
        coeff_tmp = 0d0
        do jo2 = 1, jo1 - 1
          if (has_orbit(jo2)) then
            io2 = wfi%jo_tbl(jo2)
            coeff_tmp(jo2) = dot_wf( &
              & wf%rwf(:, :, :, ispin, io2, ik, im), &
              & wf_jo1)
          end if
        end do 
        call comm_summation(coeff_tmp, coeff, sys%no, wfi%icomm_ro)

        ! Calculate exclusion term "wf_exc":
        wf_exc_tmp = 0d0
        do jo2 = 1, jo1 - 1
          if (has_orbit(jo2)) then
            io2 = wfi%jo_tbl(jo2)
            call axpy_wf_ovlp( &
              & coeff(jo2), wf%rwf(:, :, :, ispin, io2, ik, im), &
              & wf_exc_tmp)
          end if
        end do
        call comm_summation(wf_exc_tmp, wf_exc, nsize_rg, wfi%icomm_o)

        if (has_orbit(jo1)) then
          ! Exclude non-orthonormal component:
          call axpy_wf_ovlp(-1d0, wf_exc, wf_jo1)
          ! Normalization:
          norm2_tmp = dot_wf(wf_jo1, wf_jo1)
          call comm_summation(norm2_tmp, norm2, 1, wfi%icomm_r)
          call scal_wf_ovlp(1d0 / sqrt(norm2), wf_jo1)
          ! Write back to "rwf":
          io1 = wfi%jo_tbl(jo1)
          call copy_data( &
            & wf_jo1, &
            & wf%rwf(:, :, :, ispin, io1, ik, im)) 
        end if
      end do !jo1
    
    end do !ispin
    end do !ik
    end do !im
    
    return
  contains

  logical function has_orbit(jo) result(f)
    implicit none
    integer, intent(in) :: jo
    f = (1 <= wfi%jo_tbl(jo))
    return
  end function has_orbit

  ! Dot product of two wavefunctions:
  real(8) function dot_wf(x, y) result(p)
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
  end function dot_wf

  ! Constant times a wavefunction plus a wavefunction:
  subroutine axpy_wf_ovlp(a, x, y)
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
  end subroutine axpy_wf_ovlp


  ! Scales a wavefunction by a constant 
  subroutine scal_wf_ovlp(a, x)
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
  end subroutine scal_wf_ovlp

  end subroutine gram_schmidt_col_real8


!=======================================================================



  subroutine gram_schmidt_col_complex8(sys, rg, wfi, wf)
    ! Only for the colinear L(S)DA:
    use timer
    use salmon_communication, only: comm_bcast, comm_summation
    implicit none
    type(s_system),       intent(in)    :: sys
    type(s_rgrid),        intent(in)    :: rg
    type(s_wf_info),      intent(in)    :: wfi
    type(s_wavefunction), intent(inout) :: wf

    integer :: nsize_rg
    integer :: ik, im, ispin
    integer :: jo1, jo2, io1, io2
    complex(8) :: coeff(1:sys%no), coeff_tmp(1:sys%no)
    real(8) :: norm2, norm2_tmp
    complex(8), dimension( &
      & rg%is_array(1):rg%ie_array(1), &
      & rg%is_array(2):rg%ie_array(2), &
      & rg%is_array(3):rg%ie_array(3)) &
      & :: wf_jo1, wf_exc, wf_exc_tmp

    complex(8), parameter :: one = 1d0

    nsize_rg =  (rg%ie_array(1) - rg%is_array(1)) &
      & * (rg%ie_array(2) - rg%is_array(2)) &
      & * (rg%ie_array(3) - rg%is_array(3))
      
    do im = wfi%im_s, wfi%im_e
    do ik = wfi%ik_s, wfi%ik_e
    do ispin = 1, sys%nspin

      ! Loop for all orbit #jo1:
      do jo1 = 1, sys%no

        ! Retrieve orbit #jo1 into "wf_jo1":
        if (has_orbit(jo1)) then
          io1 = wfi%jo_tbl(jo1)
          call copy_data( &
            & wf%zwf(:, :, :, ispin, io1, ik, im), &
            & wf_jo1)
        end if
        call comm_bcast(wf_jo1, wfi%icomm_o, wfi%irank_jo(jo1))
        
        ! Calculate overlap coefficients: 
        coeff_tmp = 0d0
        do jo2 = 1, jo1 - 1
          if (has_orbit(jo2)) then
            io2 = wfi%jo_tbl(jo2)
            coeff_tmp(jo2) = dot_wf( &
              & wf%zwf(:, :, :, ispin, io2, ik, im), &
              & wf_jo1)
          end if
        end do 
        call comm_summation(coeff_tmp, coeff, sys%no, wfi%icomm_ro)

        ! Calculate exclusion term "wf_exc":
        wf_exc_tmp = 0d0
        do jo2 = 1, jo1 - 1
          if (has_orbit(jo2)) then
            io2 = wfi%jo_tbl(jo2)
            call axpy_wf_ovlp( &
              & coeff(jo2), wf%zwf(:, :, :, ispin, io2, ik, im), &
              & wf_exc_tmp)
          end if
        end do
        call comm_summation(wf_exc_tmp, wf_exc, nsize_rg, wfi%icomm_o)

        if (has_orbit(jo1)) then
          ! Exclude non-orthonormal component:
          call axpy_wf_ovlp(-one, wf_exc, wf_jo1)
          ! Normalization:
          norm2_tmp = real(dot_wf(wf_jo1, wf_jo1))
          call comm_summation(norm2_tmp, norm2, 1, wfi%icomm_r)
          call scal_wf_ovlp(one / sqrt(norm2), wf_jo1)
          ! Write back to "zwf":
          io1 = wfi%jo_tbl(jo1)
          call copy_data( &
            & wf_jo1, &
            & wf%zwf(:, :, :, ispin, io1, ik, im)) 
        end if
      end do !jo1
    
    end do !ispin
    end do !ik
    end do !im
    
    return
  contains

  logical function has_orbit(jo) result(f)
    implicit none
    integer, intent(in) :: jo
    f = (1 <= wfi%jo_tbl(jo))
    return
  end function has_orbit

  ! Dot product of two wavefunctions:
  complex(8) function dot_wf(x, y) result(p)
    implicit none
    complex(8), intent(in) :: x( &
      & rg%is_array(1):rg%ie_array(1), &
      & rg%is_array(2):rg%ie_array(2), &
      & rg%is_array(3):rg%ie_array(3))
    complex(8), intent(in) :: y( &
      & rg%is_array(1):rg%ie_array(1), &
      & rg%is_array(2):rg%ie_array(2), &
      & rg%is_array(3):rg%ie_array(3))
    integer :: i1, i2, i3
    p = 0d0
    !$omp parallel do collapse(2) default(shared) private(i1,i2,i3) reduction(+:p)
    do i3 = rg%is(3), rg%ie(3)
      do i2 = rg%is(2), rg%ie(2)
        do i1 = rg%is(1), rg%ie(1)
          p = p + conjg(x(i1, i2, i3)) * y(i1, i2, i3)
        end do
      end do
    end do
    !$omp end parallel do
    p = p * sys%hvol
    return
  end function dot_wf

  ! Constant times a wavefunction plus a wavefunction:
  subroutine axpy_wf_ovlp(a, x, y)
    implicit none
    complex(8), intent(in) :: a
    complex(8), intent(in) :: x( &
      & rg%is_array(1):rg%ie_array(1), &
      & rg%is_array(2):rg%ie_array(2), &
      & rg%is_array(3):rg%ie_array(3))
    complex(8), intent(inout) :: y( &
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
  end subroutine axpy_wf_ovlp


  ! Scales a wavefunction by a constant 
  subroutine scal_wf_ovlp(a, x)
    implicit none
    complex(8), intent(in) :: a
    complex(8), intent(inout) :: x( &
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
  end subroutine scal_wf_ovlp

  end subroutine gram_schmidt_col_complex8







  subroutine debug_var_dump(sys, rg, wfi) 
    use salmon_parallel
    implicit none
    type(s_system),       intent(in)    :: sys
    type(s_rgrid),        intent(in)    :: rg
    type(s_wf_info),      intent(in)    :: wfi

    character(100) :: logfile
    write(logfile,'(i3.3, ".log")') nproc_id_global
    open(unit=777, file=trim(logfile))
    write(777,*) 'system----------'
    write(777,*) 'iperiodic', sys%iperiodic
    write(777,*) 'ngrid', sys%ngrid
    write(777,*) 'nspin', sys%nspin
    write(777,*) 'no', sys%no
    write(777,*) 'nk', sys%nk
    write(777,*) 'nion', sys%nion
    write(777,*) 'Hvol', sys%Hvol
    write(777,*) 'Hgs', sys%Hgs
    write(777,*) 'al', sys%al
    write(777,*) 'det_al', sys%det_al
    write(777,*) 'brl', sys%brl
    write(777,*) 'wtk', sys%wtk
    write(777,*) 'Rion', sys%Rion
    write(777,*) 'rocc', sys%rocc
    write(777,*) 'rgrid----------'
    write(777,*) 'ndir', rg%ndir
    write(777,*) 'Nd', rg%Nd
    write(777,*) 'is', rg%is
    write(777,*) 'ie', rg%ie
    write(777,*) 'is_overlap', rg%is_overlap
    write(777,*) 'ie_overlap', rg%ie_overlap
    write(777,*) 'is_array', rg%is_array
    write(777,*) 'ie_array', rg%ie_array
    write(777,*) 'num', rg%num
    write(777,*) 'wf-info----------'
    write(777,*) 'if_divide_rspace', wfi%if_divide_rspace
    write(777,*) 'irank_r', wfi%irank_r
    write(777,*) 'icomm_r', wfi%icomm_r
    write(777,*) 'icomm_o', wfi%icomm_o
    write(777,*) 'icomm_ro', wfi%icomm_ro
    write(777,*) 'icomm_ko', wfi%icomm_ko
    write(777,*) 'icomm_rko', wfi%icomm_rko
    write(777,*) 'im_s', wfi%im_s
    write(777,*) 'im_e', wfi%im_e
    write(777,*) 'numm', wfi%numm
    write(777,*) 'ik_s', wfi%ik_s
    write(777,*) 'ik_e', wfi%ik_e
    write(777,*) 'numk', wfi%numk
    write(777,*) 'io_s', wfi%io_s
    write(777,*) 'io_e', wfi%io_e
    write(777,*) 'numo', wfi%numo
    write(777,*) 'occ', wfi%occ
    write(777,*) 'io_tbl', wfi%io_tbl
    write(777,*) 'jo_tbl', wfi%jo_tbl
    write(777,*) 'irank_jo', wfi%irank_jo

    close(777)
    return
  end subroutine
  
end module gram_schmidt_orth
