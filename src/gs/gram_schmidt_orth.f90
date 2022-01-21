!
!  Copyright 2019-2020 SALMON developers
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
  implicit none

contains

  subroutine gram_schmidt(sys, rg, wfi, wf)
    use structures, only: s_dft_system, s_rgrid, s_parallel_info, s_orbital
    use gram_schmidt_so_sub, only: gram_schmidt_so
    use salmon_global, only: yn_spinorbit, yn_gramschmidt_blas
    use timer
    implicit none
    type(s_dft_system),   intent(in)    :: sys
    type(s_rgrid),        intent(in)    :: rg
    type(s_parallel_info),intent(in)    :: wfi
    type(s_orbital),      intent(inout) :: wf

    if ( yn_spinorbit=='y' ) then
      call gram_schmidt_so(sys, rg, wfi, wf)
      return
    end if

    call timer_begin(LOG_CALC_GRAM_SCHMIDT)

    if (sys%if_real_orbital) then
      if(yn_gramschmidt_blas =='y')then
        call gram_schmidt_col_rblas(sys, rg, wfi, wf)
      else
        call gram_schmidt_col_real8(sys, rg, wfi, wf)
      end if
    else
      if(yn_gramschmidt_blas =='y')then
        call gram_schmidt_col_cblas(sys, rg, wfi, wf)
      else
        call gram_schmidt_col_complex8(sys, rg, wfi, wf)
      end if
    end if

    call timer_end(LOG_CALC_GRAM_SCHMIDT)

    return
  end subroutine

  subroutine gram_schmidt_col_real8(sys, rg, wfi, wf)
    ! Only for the colinear L(S)DA:
    use structures, only: s_dft_system, s_rgrid, s_parallel_info, s_orbital
    use pack_unpack, only: copy_data
    use timer
    use communication, only: comm_bcast, comm_summation
    implicit none
    type(s_dft_system),   intent(in)    :: sys
    type(s_rgrid),        intent(in)    :: rg
    type(s_parallel_info),intent(in)    :: wfi
    type(s_orbital),      intent(inout) :: wf

    integer :: nsize_rg
    integer :: ik, im, ispin
    integer :: io1, io2
    real(8) :: coeff(1:sys%no), coeff_tmp(1:sys%no)
    real(8) :: norm2, norm2_tmp
    real(8), dimension( &
      & rg%is_array(1):rg%ie_array(1), &
      & rg%is_array(2):rg%ie_array(2), &
      & rg%is_array(3):rg%ie_array(3)) &
      & :: wf_io1, wf_exc, wf_exc_tmp

    nsize_rg =  (rg%ie_array(1) - rg%is_array(1) + 1) &
      & * (rg%ie_array(2) - rg%is_array(2) + 1) &
      & * (rg%ie_array(3) - rg%is_array(3) + 1)

    do im = wfi%im_s, wfi%im_e
    do ik = wfi%ik_s, wfi%ik_e
    do ispin = 1, sys%nspin

      ! Loop for all orbit #io1:
      do io1 = 1, sys%no

        ! Retrieve orbit #io1 into "wf_io1":
        if (wfi%io_s<= io1 .and. io1 <= wfi%io_e) then
          call copy_data( &
            & wf%rwf(:, :, :, ispin, io1, ik, im), &
            & wf_io1)
        end if
        if (wfi%if_divide_orbit) &
          & call comm_bcast(wf_io1, wfi%icomm_o, wfi%irank_io(io1))

        ! Calculate overlap coefficients:
        coeff_tmp = 0d0
        do io2 = 1, io1 - 1
          if (wfi%io_s<= io2 .and. io2 <= wfi%io_e) then
            coeff_tmp(io2) = dot_wf( &
              & wf%rwf(:, :, :, ispin, io2, ik, im), &
              & wf_io1)
          end if
        end do
        if (wfi%if_divide_rspace .or. wfi%if_divide_orbit) then
          call comm_summation(coeff_tmp, coeff, sys%no, wfi%icomm_ro)
        else
          coeff = coeff_tmp
        endif

        ! Calculate exclusion term "wf_exc":
        wf_exc_tmp = 0d0
        do io2 = 1, io1 - 1
          if (wfi%io_s<= io2 .and. io2 <= wfi%io_e) then
            call axpy_wf( &
              & coeff(io2), wf%rwf(:, :, :, ispin, io2, ik, im), &
              & wf_exc_tmp)
          end if
        end do
        if (wfi%if_divide_orbit) then
          call comm_summation(wf_exc_tmp, wf_exc, nsize_rg, wfi%icomm_o)
        else
          call copy_data(wf_exc_tmp, wf_exc)
        end if

        if (wfi%io_s<= io1 .and. io1 <= wfi%io_e) then
          ! Exclude non-orthonormal component:
          call axpy_wf(-1d0, wf_exc, wf_io1)
          ! Normalization:
          norm2_tmp = dot_wf(wf_io1, wf_io1)
          if (wfi%if_divide_rspace) then
            call comm_summation(norm2_tmp, norm2, wfi%icomm_r)
          else
            norm2 = norm2_tmp
          endif
          call scal_wf(1d0 / sqrt(norm2), wf_io1)
          ! Write back to "rwf":
          call copy_data( &
            & wf_io1, &
            & wf%rwf(:, :, :, ispin, io1, ik, im))
        end if
      end do !io1

    end do !ispin
    end do !ik
    end do !im

    return
  contains

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
  subroutine axpy_wf(a, x, y)
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
    do i3 = rg%is(3), rg%ie(3)
      do i2 = rg%is(2), rg%ie(2)
        do i1 = rg%is(1), rg%ie(1)
          y(i1, i2, i3) = a * x(i1, i2, i3) + y(i1, i2, i3)
        end do
      end do
    end do
    !$omp end parallel do
    return
  end subroutine axpy_wf


  ! Scales a wavefunction by a constant
  subroutine scal_wf(a, x)
    implicit none
    real(8), intent(in) :: a
    real(8), intent(inout) :: x( &
      & rg%is_array(1):rg%ie_array(1), &
      & rg%is_array(2):rg%ie_array(2), &
      & rg%is_array(3):rg%ie_array(3))
    integer :: i1, i2, i3
    !$omp parallel do collapse(2) default(shared) private(i1,i2,i3)
    do i3 = rg%is(3), rg%ie(3)
      do i2 = rg%is(2), rg%ie(2)
        do i1 = rg%is(1), rg%ie(1)
          x(i1, i2, i3) = a * x(i1, i2, i3)
        end do
      end do
    end do
    !$omp end parallel do
    return
  end subroutine scal_wf

  end subroutine gram_schmidt_col_real8


!=======================================================================

  subroutine gram_schmidt_col_cblas(sys, rg, wfi, wf)
    ! Only for the colinear L(S)DA:
    use structures, only: s_dft_system, s_rgrid, s_parallel_info, s_orbital
    use pack_unpack, only: copy_data
    use timer
    use communication, only: comm_bcast, comm_summation
    implicit none
    type(s_dft_system),   intent(in)    :: sys
    type(s_rgrid),        intent(in)    :: rg
    type(s_parallel_info),intent(in)    :: wfi
    type(s_orbital),      intent(inout) :: wf

    complex(8), parameter :: zero = 0d0
    complex(8), parameter :: one = 1d0
    integer, parameter :: n_one = 1
    character(1),parameter :: TRANSA='C', TRANSB='N'

    integer :: nsize_rg
    integer :: ik, im, ispin, m
    integer :: io, jo, jo1, jo2
    integer :: idiv1, idiv2, idiv3, io0_s, io0_e, io1_s, io1_e, io2_s, io2_e, io3_s, io3_e
    integer :: numo0, numo1, numo1_0, numo1_1, numo2, numo2_0, numo2_1, numo3, numo3_0, numo3_1
    complex(8) :: coeff(wfi%numo), coeff_tmp(wfi%numo)
    real(8) :: norm2, norm2_tmp
    complex(8) :: wf_block(rg%is(1):rg%ie(1), rg%is(2):rg%ie(2), rg%is(3):rg%ie(3), wfi%numo)
    complex(8) :: wf_block1(rg%is(1):rg%ie(1), rg%is(2):rg%ie(2), rg%is(3):rg%ie(3), wfi%numo)
    complex(8) :: wf_block_send(rg%is(1):rg%ie(1), rg%is(2):rg%ie(2), rg%is(3):rg%ie(3), wfi%numo_max)
    complex(8) :: umat(wfi%numo_max,wfi%numo), umat_tmp(wfi%numo_max,wfi%numo)
    complex(8) :: ZDOTC

    nsize_rg = (rg%ie(1)-rg%is(1)+1)*(rg%ie(2)-rg%is(2)+1)*(rg%ie(3)-rg%is(3)+1)

    do im = wfi%im_s, wfi%im_e
    do ik = wfi%ik_s, wfi%ik_e
    do ispin = 1, sys%nspin

      ! Copy wave function
      do io = wfi%io_s, wfi%io_e
        jo = io - wfi%io_s + 1
        call copy_data( &
          & wf%zwf(rg%is(1):rg%ie(1), rg%is(2):rg%ie(2), rg%is(3):rg%ie(3), ispin, io, ik, im), &
          & wf_block(:, :, :, jo))
      end do

      ! Loop for all:
      do m = 0, wfi%nporbital-1
        if(wfi%id_o==m)then

          if(wfi%numo < 32) then

        ! Loop for orbit #io1 in my rank:
            do jo1 = 1, wfi%numo

              ! normalize orbital jo1
              norm2_tmp = real( ZDOTC( &
               &      nsize_rg, wf_block(:,:,:,jo1), n_one, &
               &                wf_block(:,:,:,jo1), n_one )) &
               &     * sys%hvol
              if (wfi%if_divide_rspace) then
                call comm_summation(norm2_tmp, norm2, wfi%icomm_r)
              else
                norm2 = norm2_tmp
              endif
              call ZDSCAL(nsize_rg, one/sqrt(norm2), wf_block(:,:,:,jo1), n_one)

              ! Calculate overlap coefficients:
              coeff_tmp = 0d0
              do jo2 = jo1+1, wfi%numo
                coeff_tmp(jo2) = ZDOTC( &
                 &      nsize_rg, wf_block(:,:,:,jo1), n_one, &
                 &                wf_block(:,:,:,jo2), n_one ) &
                 &     * sys%hvol
              end do
              if (wfi%if_divide_rspace) then
                call comm_summation(coeff_tmp, coeff, wfi%numo, wfi%icomm_r)
              else
                coeff = coeff_tmp
              endif

              ! Exclude non-orthonormal component:
              do jo2 = jo1+1, wfi%numo
                call ZAXPY( &
                  &  nsize_rg, - coeff(jo2), wf_block(:,:,:,jo1), n_one, &
                  &                          wf_block(:,:,:,jo2), n_one)
              end do

            end do !jo1

          else

            io0_s = 1
            io0_e = wfi%numo
            numo0 = wfi%numo

            numo1_0 = (numo0 + 1)/2
            numo1_1 = numo0 - numo1_0
            do idiv1 = 0, 1
              if( idiv1 == 0 ) then
                io1_s = io0_s
                io1_e = io0_s - 1 + numo1_0
              else
                io1_s = io0_s     + numo1_0
                io1_e = io0_e
              end if
              numo1 = io1_e - io1_s + 1

              numo2_0 = (numo1 + 1)/2
              numo2_1 = numo1 - numo2_0
              do idiv2 = 0, 1
                if( idiv2 == 0 ) then
                  io2_s = io1_s
                  io2_e = io1_s - 1 + numo2_0
                else
                  io2_s = io1_s     + numo2_0
                  io2_e = io1_e
                end if
                numo2 = io2_e - io2_s + 1

                numo3_0 = (numo2 + 1)/2
                numo3_1 = numo2 - numo3_0
                do idiv3 = 0, 1
                  if( idiv3 == 0 ) then
                    io3_s = io2_s
                    io3_e = io2_s - 1 + numo3_0
                  else
                    io3_s = io2_s     + numo3_0
                    io3_e = io2_e
                  end if
                  numo3 = io3_e - io3_s + 1

                  do jo1 = io3_s, io3_e
                    ! normaliza the orbital jo1
                    norm2_tmp = real( ZDOTC( &
                     &      nsize_rg, wf_block(:,:,:,jo1), n_one, &
                     &                wf_block(:,:,:,jo1), n_one )) &
                     &     * sys%hvol
                    if (wfi%if_divide_rspace) then
                      call comm_summation(norm2_tmp, norm2, wfi%icomm_r)
                    else
                      norm2 = norm2_tmp
                    endif
                    call ZDSCAL(nsize_rg, one/sqrt(norm2), wf_block(:,:,:,jo1), n_one)

                ! Calculate overlap coefficients:
                    coeff_tmp = 0d0
                    do jo2 = jo1+1, io3_e
                      coeff_tmp(jo2) = ZDOTC( &
                       &      nsize_rg, wf_block(:,:,:,jo1), n_one, &
                       &                wf_block(:,:,:,jo2), n_one ) &
                       &     * sys%hvol
                    end do
                    if (wfi%if_divide_rspace) then
                      call comm_summation(coeff_tmp(io3_s:io3_e), coeff(io3_s:io3_e), numo3, wfi%icomm_r)
                    else
                      coeff(io3_s:io3_e) = coeff_tmp(io3_s:io3_e)
                    endif

                    ! Exclude non-orthonormal component:
                    do jo2 = jo1+1, io3_e
                      call ZAXPY( &
                        &  nsize_rg, - coeff(jo2), wf_block(:,:,:,jo1), n_one, &
                        &                          wf_block(:,:,:,jo2), n_one)
                    end do

                  end do !jo1

                  if(idiv3 == 0) then

                    umat_tmp = 0.d0
                    wf_block_send(:,:,:,1:numo3_0) = wf_block(:,:,:,io2_s:io2_s-1+numo3_0)
                    wf_block1(:,:,:,1:numo3_1) = wf_block(:,:,:,io2_s+numo3_0:io2_e)
                    call zgemm(TRANSA, TRANSB, numo3_0, numo3_1, nsize_rg,  &
                        &  one*sys%hvol, wf_block_send, nsize_rg,  &
                        &                wf_block1, nsize_rg,  &
                        &          zero, umat_tmp, wfi%numo_max)

                   if(wfi%if_divide_rspace) then
                      call comm_summation(umat_tmp, umat, wfi%numo_max*wfi%numo, wfi%icomm_r)
                    else
                      umat = umat_tmp
                    end if

                    call zgemm(TRANSB, TRANSB, nsize_rg, numo3_1, numo3_0,  &
                      &         - one, wf_block_send, nsize_rg,  &
                      &                umat, wfi%numo_max,  &
                      &           one, wf_block1, nsize_rg)
                    wf_block(:,:,:,io2_s+numo3_0:io2_e) = wf_block1(:,:,:,1:numo3_1)
                  end if ! idiv3 == 0

                end do !idiv3

                if(idiv2 == 0) then

                  umat_tmp = 0.d0
                  wf_block_send(:,:,:,1:numo2_0) = wf_block(:,:,:,io1_s:io1_s-1+numo2_0)
                  wf_block1(:,:,:,1:numo2_1) = wf_block(:,:,:,io1_s+numo2_0:io1_e)
                  call zgemm(TRANSA, TRANSB, numo2_0, numo2_1, nsize_rg,  &
                      &  one*sys%hvol, wf_block_send, nsize_rg,  &
                      &                wf_block1, nsize_rg,  &
                      &          zero, umat_tmp, wfi%numo_max)

                  if(wfi%if_divide_rspace) then
                    call comm_summation(umat_tmp, umat, wfi%numo_max*wfi%numo, wfi%icomm_r)
                  else
                    umat = umat_tmp
                  end if

                  call zgemm(TRANSB, TRANSB, nsize_rg, numo2_1, numo2_0,  &
                    &         - one, wf_block_send, nsize_rg,  &
                    &                umat, wfi%numo_max,  &
                    &           one, wf_block1, nsize_rg)

                  wf_block(:,:,:,io1_s+numo2_0:io1_e) = wf_block1(:,:,:,1:numo2_1)
                end if ! idiv2 == 0

              end do !idiv2

              if(idiv1 == 0) then

                umat_tmp = 0.d0
                wf_block_send(:,:,:,1:numo1_0) = wf_block(:,:,:,io0_s:io0_s-1+numo1_0)
                wf_block1(:,:,:,1:numo1_1) = wf_block(:,:,:,io0_s+numo1_0:io0_e)

                call zgemm(TRANSA, TRANSB, numo1_0, numo1_1, nsize_rg,  &
                    &  one*sys%hvol, wf_block_send, nsize_rg,  &
                    &                wf_block1, nsize_rg,  &
                    &          zero, umat_tmp, wfi%numo_max)

                if(wfi%if_divide_rspace) then
                  call comm_summation(umat_tmp, umat, wfi%numo_max*wfi%numo, wfi%icomm_r)
                else
                  umat = umat_tmp
                end if

                call zgemm(TRANSB, TRANSB, nsize_rg, numo1_1, numo1_0,  &
                  &         - one, wf_block_send, nsize_rg,  &
                  &                umat, wfi%numo_max,  &
                  &           one, wf_block1, nsize_rg)
                wf_block(:,:,:,io0_s+numo1_0:io0_e) = wf_block1(:,:,:,1:numo1_1)
              end if ! idiv1 == 0

            end do !idiv1

          end if ! if(wfi%numo < 32) else

        end if !m==wfi%id_o

        if (wfi%if_divide_orbit) then
          call copy_data( wf_block(:,:,:,:), wf_block_send(:,:,:,1:wfi%numo) )
          wf_block_send(:,:,:,wfi%numo+1:wfi%numo_max) = 0d0
          call comm_bcast( wf_block_send, wfi%icomm_o, wfi%irank_io(wfi%io_s_all(m)))

          umat_tmp=0.d0
          if( wfi%id_o > m )then
            call zgemm(TRANSA, TRANSB, wfi%numo_max, wfi%numo, nsize_rg,  &
              &   one*sys%hvol, wf_block_send(:,:,:,1), nsize_rg,  &
              &                wf_block(:,:,:,1), nsize_rg,  &
              &          zero, umat_tmp, wfi%numo_max )

            if(wfi%if_divide_rspace) then
              call comm_summation(umat_tmp, umat, wfi%numo_max*wfi%numo, wfi%icomm_r)
            else
              umat = umat_tmp
            end if

            call zgemm(TRANSB, TRANSB, nsize_rg, wfi%numo, wfi%numo_max,  &
              &         - one, wf_block_send(:,:,:,1), nsize_rg,  &
              &                umat, wfi%numo_max,  &
              &           one, wf_block(:,:,:,1), nsize_rg)
          end if
        end if ! wfi%if_divide_orbit

      end do !m

      do io = wfi%io_s, wfi%io_e
        jo = io - wfi%io_s + 1
        call copy_data( &
          & wf_block(:, :, :, jo), &
          & wf%zwf(rg%is(1):rg%ie(1), rg%is(2):rg%ie(2), rg%is(3):rg%ie(3), ispin, io, ik, im))
      end do

    end do !ispin
    end do !ik
    end do !im

    return

  end subroutine gram_schmidt_col_cblas

  !=======================================================================

    subroutine gram_schmidt_col_cblas_old(sys, rg, wfi, wf)
      ! Only for the colinear L(S)DA:
      use structures, only: s_dft_system, s_rgrid, s_parallel_info, s_orbital
      use pack_unpack, only: copy_data
      use timer
      use communication, only: comm_bcast, comm_summation
      implicit none
      type(s_dft_system),   intent(in)    :: sys
      type(s_rgrid),        intent(in)    :: rg
      type(s_parallel_info),intent(in)    :: wfi
      type(s_orbital),      intent(inout) :: wf

      complex(8), parameter :: zero = 0d0
      complex(8), parameter :: one = 1d0
      integer, parameter :: n_one = 1
      character(1),parameter :: TRANSA='C', TRANSB='N'

      integer :: nsize_rg
      integer :: ik, im, ispin, m
      integer :: io, jo, jo1, jo2
      complex(8) :: coeff(wfi%numo), coeff_tmp(wfi%numo)
      real(8) :: norm2, norm2_tmp
      complex(8) :: wf_block(rg%is(1):rg%ie(1), rg%is(2):rg%ie(2), rg%is(3):rg%ie(3), wfi%numo)
      complex(8) :: wf_block_send(rg%is(1):rg%ie(1), rg%is(2):rg%ie(2), rg%is(3):rg%ie(3), wfi%numo_max)
      complex(8) :: umat(wfi%numo_max,wfi%numo), umat_tmp(wfi%numo_max,wfi%numo)
      complex(8) :: ZDOTC

      nsize_rg = (rg%ie(1)-rg%is(1)+1)*(rg%ie(2)-rg%is(2)+1)*(rg%ie(3)-rg%is(3)+1)

      do im = wfi%im_s, wfi%im_e
      do ik = wfi%ik_s, wfi%ik_e
      do ispin = 1, sys%nspin

        ! Copy wave function
        do io = wfi%io_s, wfi%io_e
          jo = io - wfi%io_s + 1
          call copy_data( &
            & wf%zwf(rg%is(1):rg%ie(1), rg%is(2):rg%ie(2), rg%is(3):rg%ie(3), ispin, io, ik, im), &
            & wf_block(:, :, :, jo))
        end do

        ! Loop for all:
        do m = 0, wfi%nporbital-1
          if(wfi%id_o==m)then
          ! Loop for orbit #io1 in my rank:
            do jo1 = 1, wfi%numo

              ! normalize orbital jo1
              norm2_tmp = real( ZDOTC( &
               &      nsize_rg, wf_block(:,:,:,jo1), n_one, &
               &                wf_block(:,:,:,jo1), n_one )) &
               &     * sys%hvol
              if (wfi%if_divide_rspace) then
                call comm_summation(norm2_tmp, norm2, wfi%icomm_r)
              else
                norm2 = norm2_tmp
              endif
              call ZDSCAL(nsize_rg, one/sqrt(norm2), wf_block(:,:,:,jo1), n_one)

              ! Calculate overlap coefficients:
              coeff_tmp = 0d0
              do jo2 = jo1+1, wfi%numo
                coeff_tmp(jo2) = ZDOTC( &
                 &      nsize_rg, wf_block(:,:,:,jo1), n_one, &
                 &                wf_block(:,:,:,jo2), n_one ) &
                 &     * sys%hvol
              end do
              if (wfi%if_divide_rspace) then
                call comm_summation(coeff_tmp, coeff, wfi%numo, wfi%icomm_r)
              else
                coeff = coeff_tmp
              endif

              ! Exclude non-orthonormal component:
              do jo2 = jo1+1, wfi%numo
                call ZAXPY( &
                  &  nsize_rg, - coeff(jo2), wf_block(:,:,:,jo1), n_one, &
                  &                          wf_block(:,:,:,jo2), n_one)
              end do

            end do !jo1

          end if !m==wfi%id_o

          if (wfi%if_divide_orbit) then
            call copy_data( wf_block(:,:,:,:), wf_block_send(:,:,:,1:wfi%numo) )
            wf_block_send(:,:,:,wfi%numo+1:wfi%numo_max) = 0d0
            call comm_bcast( wf_block_send, wfi%icomm_o, wfi%irank_io(wfi%io_s_all(m)))

            umat_tmp=0.d0
            if( wfi%id_o > m )then
              call zgemm(TRANSA, TRANSB, wfi%numo_max, wfi%numo, nsize_rg,  &
                &   one*sys%hvol, wf_block_send(:,:,:,1), nsize_rg,  &
                &                wf_block(:,:,:,1), nsize_rg,  &
                &          zero, umat_tmp, wfi%numo_max )

              if(wfi%if_divide_rspace) then
                call comm_summation(umat_tmp, umat, wfi%numo_max*wfi%numo, wfi%icomm_r)
              else
                umat = umat_tmp
              end if

              call zgemm(TRANSB, TRANSB, nsize_rg, wfi%numo, wfi%numo_max,  &
                &         - one, wf_block_send(:,:,:,1), nsize_rg,  &
                &                umat, wfi%numo_max,  &
                &           one, wf_block(:,:,:,1), nsize_rg)
            end if
          end if ! wfi%if_divide_orbit

        end do !m

        do io = wfi%io_s, wfi%io_e
          jo = io - wfi%io_s + 1
          call copy_data( &
            & wf_block(:, :, :, jo), &
            & wf%zwf(rg%is(1):rg%ie(1), rg%is(2):rg%ie(2), rg%is(3):rg%ie(3), ispin, io, ik, im))
        end do

      end do !ispin
      end do !ik
      end do !im

      return

    end subroutine gram_schmidt_col_cblas_old

  !=======================================================================


    subroutine gram_schmidt_col_rblas(sys, rg, wfi, wf)
      ! Only for the colinear L(S)DA:
      use structures, only: s_dft_system, s_rgrid, s_parallel_info, s_orbital
      use pack_unpack, only: copy_data
      use timer
      use communication, only: comm_bcast, comm_summation
      implicit none
      type(s_dft_system),   intent(in)    :: sys
      type(s_rgrid),        intent(in)    :: rg
      type(s_parallel_info),intent(in)    :: wfi
      type(s_orbital),      intent(inout) :: wf

      real(8), parameter :: zero = 0d0
      real(8), parameter :: one = 1d0
      integer, parameter :: n_one = 1
      character(1),parameter :: TRANSA='T', TRANSB='N'

      integer :: nsize_rg
      integer :: ik, im, ispin, m
      integer :: io, jo, jo1, jo2
      integer :: idiv1, idiv2, idiv3, io0_s, io0_e, io1_s, io1_e, io2_s, io2_e, io3_s, io3_e
      integer :: numo0, numo1, numo1_0, numo1_1, numo2, numo2_0, numo2_1, numo3, numo3_0, numo3_1
      real(8) :: coeff(wfi%numo), coeff_tmp(wfi%numo)
      real(8) :: norm2, norm2_tmp
      real(8) :: wf_block(rg%is(1):rg%ie(1), rg%is(2):rg%ie(2), rg%is(3):rg%ie(3), wfi%numo)
      real(8) :: wf_block1(rg%is(1):rg%ie(1), rg%is(2):rg%ie(2), rg%is(3):rg%ie(3), wfi%numo)
      real(8) :: wf_block_send(rg%is(1):rg%ie(1), rg%is(2):rg%ie(2), rg%is(3):rg%ie(3), wfi%numo_max)
      real(8) :: umat(wfi%numo_max,wfi%numo), umat_tmp(wfi%numo_max,wfi%numo)
      real(8) :: DDOT

      nsize_rg = (rg%ie(1)-rg%is(1)+1)*(rg%ie(2)-rg%is(2)+1)*(rg%ie(3)-rg%is(3)+1)

      do im = wfi%im_s, wfi%im_e
      do ik = wfi%ik_s, wfi%ik_e
      do ispin = 1, sys%nspin

        ! Copy wave function
        do io = wfi%io_s, wfi%io_e
          jo = io - wfi%io_s + 1
          call copy_data( &
            &  wf%rwf(rg%is(1):rg%ie(1), rg%is(2):rg%ie(2), rg%is(3):rg%ie(3), ispin, io, ik, im) , &
            & wf_block(:, :, :, jo))
        end do

        ! Loop for all:
        do m = 0, wfi%nporbital-1
          if(wfi%id_o==m)then

          if(wfi%numo < 32) then

          ! Loop for orbit #jo1 in my rank:
            do jo1 = 1, wfi%numo

              ! normaliza the orbital jo1
              norm2_tmp = DDOT( &
               &      nsize_rg, wf_block(:,:,:,jo1), n_one, &
               &                wf_block(:,:,:,jo1), n_one ) &
               &     * sys%hvol
              if (wfi%if_divide_rspace) then
                call comm_summation(norm2_tmp, norm2, wfi%icomm_r)
              else
                norm2 = norm2_tmp
              endif
              call DSCAL(nsize_rg, one/sqrt(norm2), wf_block(:,:,:,jo1), n_one)

              ! Calculate overlap coefficients:
              coeff_tmp = 0d0
              do jo2 = jo1+1, wfi%numo
                coeff_tmp(jo2) = DDOT( &
                 &      nsize_rg, wf_block(:,:,:,jo1), n_one, &
                 &                wf_block(:,:,:,jo2), n_one ) &
                 &     * sys%hvol
              end do
              if (wfi%if_divide_rspace) then
                call comm_summation(coeff_tmp, coeff, wfi%numo, wfi%icomm_r)
              else
                coeff = coeff_tmp
              endif

              ! Exclude non-orthonormal component:
              do jo2 = jo1+1, wfi%numo
                call DAXPY( &
                  &  nsize_rg, - coeff(jo2), wf_block(:,:,:,jo1), n_one, &
                  &                          wf_block(:,:,:,jo2), n_one)
              end do

            end do !jo1

          else

            io0_s = 1
            io0_e = wfi%numo
            numo0 = wfi%numo

            numo1_0 = (numo0 + 1)/2
            numo1_1 = numo0 - numo1_0
            do idiv1 = 0, 1
              if( idiv1 == 0 ) then
                io1_s = io0_s
                io1_e = io0_s - 1 + numo1_0
              else
                io1_s = io0_s     + numo1_0
                io1_e = io0_e
              end if
              numo1 = io1_e - io1_s + 1

              numo2_0 = (numo1 + 1)/2
              numo2_1 = numo1 - numo2_0
              do idiv2 = 0, 1
                if( idiv2 == 0 ) then
                  io2_s = io1_s
                  io2_e = io1_s - 1 + numo2_0
                else
                  io2_s = io1_s     + numo2_0
                  io2_e = io1_e
                end if
                numo2 = io2_e - io2_s + 1

                numo3_0 = (numo2 + 1)/2
                numo3_1 = numo2 - numo3_0
                do idiv3 = 0, 1
                  if( idiv3 == 0 ) then
                    io3_s = io2_s
                    io3_e = io2_s - 1 + numo3_0
                  else
                    io3_s = io2_s     + numo3_0
                    io3_e = io2_e
                  end if
                  numo3 = io3_e - io3_s + 1

                  do jo1 = io3_s, io3_e
                    ! normaliza the orbital jo1
                    norm2_tmp = DDOT( &
                     &      nsize_rg, wf_block(:,:,:,jo1), n_one, &
                     &                wf_block(:,:,:,jo1), n_one ) &
                     &     * sys%hvol
                    if (wfi%if_divide_rspace) then
                      call comm_summation(norm2_tmp, norm2, wfi%icomm_r)
                    else
                      norm2 = norm2_tmp
                    endif
                    call DSCAL(nsize_rg, one/sqrt(norm2), wf_block(:,:,:,jo1), n_one)

                ! Calculate overlap coefficients:
                    coeff_tmp = 0d0
                    do jo2 = jo1+1, io3_e
                      coeff_tmp(jo2) = DDOT( &
                       &      nsize_rg, wf_block(:,:,:,jo1), n_one, &
                       &                wf_block(:,:,:,jo2), n_one ) &
                       &     * sys%hvol
                    end do
                    if (wfi%if_divide_rspace) then
                      call comm_summation(coeff_tmp(io3_s:io3_e), coeff(io3_s:io3_e), numo3, wfi%icomm_r)
                    else
                      coeff(io3_s:io3_e) = coeff_tmp(io3_s:io3_e)
                    endif

                    ! Exclude non-orthonormal component:
                    do jo2 = jo1+1, io3_e
                      call DAXPY( &
                        &  nsize_rg, - coeff(jo2), wf_block(:,:,:,jo1), n_one, &
                        &                          wf_block(:,:,:,jo2), n_one)
                    end do

                  end do !jo1

                  if(idiv3 == 0) then
                    umat_tmp = 0.d0
                    wf_block_send(:,:,:,1:numo3_0) = wf_block(:,:,:,io2_s:io2_s-1+numo3_0)
                    wf_block1(:,:,:,1:numo3_1) = wf_block(:,:,:,io2_s+numo3_0:io2_e)
                    call dgemm(TRANSA, TRANSB, numo3_0, numo3_1, nsize_rg,  &
                        &      sys%hvol, wf_block_send, nsize_rg,  &
                        &                wf_block1, nsize_rg,  &
                        &          zero, umat_tmp, wfi%numo_max)
                    if(wfi%if_divide_rspace) then
                      call comm_summation(umat_tmp, umat, wfi%numo_max*wfi%numo, wfi%icomm_r)
                    else
                      umat = umat_tmp
                    end if

                    call dgemm(TRANSB, TRANSB, nsize_rg, numo3_1, numo3_0,  &
                      &         - one, wf_block_send, nsize_rg,  &
                      &                umat, wfi%numo_max,  &
                      &           one, wf_block1, nsize_rg)
                    wf_block(:,:,:,io2_s+numo3_0:io2_e) = wf_block1(:,:,:,1:numo3_1)
                  end if ! idiv3 == 0

                end do !idiv3

                if(idiv2 == 0) then

                  umat_tmp = 0.d0
                  wf_block_send(:,:,:,1:numo2_0) = wf_block(:,:,:,io1_s:io1_s-1+numo2_0)
                  wf_block1(:,:,:,1:numo2_1) = wf_block(:,:,:,io1_s+numo2_0:io1_e)
                  call dgemm(TRANSA, TRANSB, numo2_0, numo2_1, nsize_rg,  &
                      &      sys%hvol, wf_block_send, nsize_rg,  &
                      &                wf_block1, nsize_rg,  &
                      &          zero, umat_tmp, wfi%numo_max)

                  if(wfi%if_divide_rspace) then
                    call comm_summation(umat_tmp, umat, wfi%numo_max*wfi%numo, wfi%icomm_r)
                  else
                    umat = umat_tmp
                  end if

                  call dgemm(TRANSB, TRANSB, nsize_rg, numo2_1, numo2_0,  &
                    &         - one, wf_block_send, nsize_rg,  &
                    &                umat, wfi%numo_max,  &
                    &           one, wf_block1, nsize_rg)
                  wf_block(:,:,:,io1_s+numo2_0:io1_e) = wf_block1(:,:,:,1:numo2_1)
                end if ! idiv2 == 0

              end do !idiv2

              if(idiv1 == 0) then
                umat_tmp = 0.d0
                wf_block_send(:,:,:,1:numo1_0) = wf_block(:,:,:,io0_s:io0_s-1+numo1_0)
                wf_block1(:,:,:,1:numo1_1) = wf_block(:,:,:,io0_s+numo1_0:io0_e)
                call dgemm(TRANSA, TRANSB, numo1_0, numo1_1, nsize_rg,  &
                    &      sys%hvol, wf_block_send, nsize_rg,  &
                    &                wf_block1, nsize_rg,  &
                    &          zero, umat_tmp, wfi%numo_max)

                if(wfi%if_divide_rspace) then
                  call comm_summation(umat_tmp, umat, wfi%numo_max*wfi%numo, wfi%icomm_r)
                else
                  umat = umat_tmp
                end if

                call dgemm(TRANSB, TRANSB, nsize_rg, numo1_1, numo1_0,  &
                  &         - one, wf_block_send, nsize_rg,  &
                  &                umat, wfi%numo_max,  &
                  &           one, wf_block1, nsize_rg)
                wf_block(:,:,:,io0_s+numo1_0:io0_e) = wf_block1(:,:,:,1:numo1_1)
              end if ! idiv1 == 0

            end do !idiv1

          end if ! if(wfi%numo < 32) else

          end if !m==wfi%id_o

          if (wfi%if_divide_orbit) then
            call copy_data( wf_block(:,:,:,:), wf_block_send(:,:,:,1:wfi%numo) )
            wf_block_send(:,:,:,wfi%numo+1:wfi%numo_max) = 0d0
            call comm_bcast( wf_block_send, wfi%icomm_o, wfi%irank_io(wfi%io_s_all(m)))

            umat_tmp=0.d0
            if( wfi%id_o > m )then
              call dgemm(TRANSA, TRANSB, wfi%numo_max, wfi%numo, nsize_rg,  &
                &      sys%hvol, wf_block_send(:,:,:,1), nsize_rg,  &
                &                wf_block(:,:,:,1), nsize_rg,  &
                &          zero, umat_tmp, wfi%numo_max )

              if(wfi%if_divide_rspace) then
                call comm_summation(umat_tmp, umat, wfi%numo_max*wfi%numo, wfi%icomm_r)
              else
                umat = umat_tmp
              end if

              call dgemm(TRANSB, TRANSB, nsize_rg, wfi%numo, wfi%numo_max,  &
                &         - one, wf_block_send(:,:,:,1), nsize_rg,  &
                &                umat, wfi%numo_max,  &
                &           one, wf_block(:,:,:,1), nsize_rg)
            end if
          end if ! wfi%if_divide_orbit

        end do !m

        do io = wfi%io_s, wfi%io_e
          jo = io - wfi%io_s + 1
          call copy_data( &
            & wf_block(:, :, :, jo), &
            & wf%rwf(rg%is(1):rg%ie(1), rg%is(2):rg%ie(2), rg%is(3):rg%ie(3), ispin, io, ik, im))
        end do

      end do !ispin
      end do !ik
      end do !im

      return

    end subroutine gram_schmidt_col_rblas

    !=======================================================================

      subroutine gram_schmidt_col_rblas_old(sys, rg, wfi, wf)
        ! Only for the colinear L(S)DA:
        use structures, only: s_dft_system, s_rgrid, s_parallel_info, s_orbital
        use pack_unpack, only: copy_data
        use timer
        use communication, only: comm_bcast, comm_summation
        implicit none
        type(s_dft_system),   intent(in)    :: sys
        type(s_rgrid),        intent(in)    :: rg
        type(s_parallel_info),intent(in)    :: wfi
        type(s_orbital),      intent(inout) :: wf

        real(8), parameter :: zero = 0d0
        real(8), parameter :: one = 1d0
        integer, parameter :: n_one = 1
        character(1),parameter :: TRANSA='T', TRANSB='N'

        integer :: nsize_rg
        integer :: ik, im, ispin, m
        integer :: io, jo, jo1, jo2
        real(8) :: coeff(wfi%numo), coeff_tmp(wfi%numo)
        real(8) :: norm2, norm2_tmp
        real(8) :: wf_block(rg%is(1):rg%ie(1), rg%is(2):rg%ie(2), rg%is(3):rg%ie(3), wfi%numo)
        real(8) :: wf_block_send(rg%is(1):rg%ie(1), rg%is(2):rg%ie(2), rg%is(3):rg%ie(3), wfi%numo_max)
        real(8) :: umat(wfi%numo_max,wfi%numo), umat_tmp(wfi%numo_max,wfi%numo)
        real(8) :: DDOT

        nsize_rg = (rg%ie(1)-rg%is(1)+1)*(rg%ie(2)-rg%is(2)+1)*(rg%ie(3)-rg%is(3)+1)

        do im = wfi%im_s, wfi%im_e
        do ik = wfi%ik_s, wfi%ik_e
        do ispin = 1, sys%nspin

          ! Copy wave function
          do io = wfi%io_s, wfi%io_e
            jo = io - wfi%io_s + 1
            call copy_data( &
              &  wf%rwf(rg%is(1):rg%ie(1), rg%is(2):rg%ie(2), rg%is(3):rg%ie(3), ispin, io, ik, im) , &
              & wf_block(:, :, :, jo))
          end do

          ! Loop for all:
          do m = 0, wfi%nporbital-1
            if(wfi%id_o==m)then
            ! Loop for orbit #jo1 in my rank:
              do jo1 = 1, wfi%numo

                ! normaliza the orbital jo1
                norm2_tmp = DDOT( &
                 &      nsize_rg, wf_block(:,:,:,jo1), n_one, &
                 &                wf_block(:,:,:,jo1), n_one ) &
                 &     * sys%hvol
                if (wfi%if_divide_rspace) then
                  call comm_summation(norm2_tmp, norm2, wfi%icomm_r)
                else
                  norm2 = norm2_tmp
                endif
                call DSCAL(nsize_rg, one/sqrt(norm2), wf_block(:,:,:,jo1), n_one)

                ! Calculate overlap coefficients:
                coeff_tmp = 0d0
                do jo2 = jo1+1, wfi%numo
                  coeff_tmp(jo2) = DDOT( &
                   &      nsize_rg, wf_block(:,:,:,jo1), n_one, &
                   &                wf_block(:,:,:,jo2), n_one ) &
                   &     * sys%hvol
                end do
                if (wfi%if_divide_rspace) then
                  call comm_summation(coeff_tmp, coeff, wfi%numo, wfi%icomm_r)
                else
                  coeff = coeff_tmp
                endif

                ! Exclude non-orthonormal component:
                do jo2 = jo1+1, wfi%numo
                  call DAXPY( &
                    &  nsize_rg, - coeff(jo2), wf_block(:,:,:,jo1), n_one, &
                    &                          wf_block(:,:,:,jo2), n_one)
                end do

              end do !jo1

            end if !m==wfi%id_o

            if (wfi%if_divide_orbit) then
              call copy_data( wf_block(:,:,:,:), wf_block_send(:,:,:,1:wfi%numo) )
              wf_block_send(:,:,:,wfi%numo+1:wfi%numo_max) = 0d0
              call comm_bcast( wf_block_send, wfi%icomm_o, wfi%irank_io(wfi%io_s_all(m)))

              umat_tmp=0.d0
              if( wfi%id_o > m )then
                call dgemm(TRANSA, TRANSB, wfi%numo_max, wfi%numo, nsize_rg,  &
                  &      sys%hvol, wf_block_send(:,:,:,1), nsize_rg,  &
                  &                wf_block(:,:,:,1), nsize_rg,  &
                  &          zero, umat_tmp, wfi%numo_max )

                if(wfi%if_divide_rspace) then
                  call comm_summation(umat_tmp, umat, wfi%numo_max*wfi%numo, wfi%icomm_r)
                else
                  umat = umat_tmp
                end if

                call dgemm(TRANSB, TRANSB, nsize_rg, wfi%numo, wfi%numo_max,  &
                  &         - one, wf_block_send(:,:,:,1), nsize_rg,  &
                  &                umat, wfi%numo_max,  &
                  &           one, wf_block(:,:,:,1), nsize_rg)
              end if
            end if ! wfi%if_divide_orbit

          end do !m

          do io = wfi%io_s, wfi%io_e
            jo = io - wfi%io_s + 1
            call copy_data( &
              & wf_block(:, :, :, jo), &
              & wf%rwf(rg%is(1):rg%ie(1), rg%is(2):rg%ie(2), rg%is(3):rg%ie(3), ispin, io, ik, im))
          end do

        end do !ispin
        end do !ik
        end do !im

        return

      end subroutine gram_schmidt_col_rblas_old

    !=======================================================================

    subroutine gram_schmidt_col_complex8(sys, rg, wfi, wf)
      ! Only for the colinear L(S)DA:
      use structures, only: s_dft_system, s_rgrid, s_parallel_info, s_orbital
      use pack_unpack, only: copy_data
      use timer
      use communication, only: comm_bcast, comm_summation
      implicit none
      type(s_dft_system),   intent(in)    :: sys
      type(s_rgrid),        intent(in)    :: rg
      type(s_parallel_info),intent(in)    :: wfi
      type(s_orbital),      intent(inout) :: wf

      integer :: nsize_rg
      integer :: ik, im, ispin
      integer :: io1, io2
      complex(8) :: coeff(1:sys%no), coeff_tmp(1:sys%no)
      real(8) :: norm2, norm2_tmp
      complex(8), dimension( &
        & rg%is_array(1):rg%ie_array(1), &
        & rg%is_array(2):rg%ie_array(2), &
        & rg%is_array(3):rg%ie_array(3)) &
        & :: wf_io1, wf_exc, wf_exc_tmp

      complex(8), parameter :: one = 1d0

      nsize_rg =  (rg%ie_array(1) - rg%is_array(1) + 1) &
        & * (rg%ie_array(2) - rg%is_array(2) + 1) &
        & * (rg%ie_array(3) - rg%is_array(3) + 1)

      do im = wfi%im_s, wfi%im_e
      do ik = wfi%ik_s, wfi%ik_e
      do ispin = 1, sys%nspin

        ! Loop for all orbit #io1:
        do io1 = 1, sys%no

          ! Retrieve orbit #io1 into "wf_io1":
          if (wfi%io_s<= io1 .and. io1 <= wfi%io_e) then
            call copy_data( &
              & wf%zwf(:, :, :, ispin, io1, ik, im), &
              & wf_io1)
          end if
          if (wfi%if_divide_orbit) &
            & call comm_bcast(wf_io1, wfi%icomm_o, wfi%irank_io(io1))

          ! Calculate overlap coefficients:
          coeff_tmp = 0d0
          do io2 = 1, io1 - 1
            if (wfi%io_s<= io2 .and. io2 <= wfi%io_e) then
              coeff_tmp(io2) = dot_wf( &
                & wf%zwf(:, :, :, ispin, io2, ik, im), &
                & wf_io1)
            end if
          end do
          if (wfi%if_divide_rspace .or. wfi%if_divide_orbit) then
            call comm_summation(coeff_tmp, coeff, sys%no, wfi%icomm_ro)
          else
            coeff = coeff_tmp
          endif

          ! Calculate exclusion term "wf_exc":
          wf_exc_tmp = 0d0
          do io2 = 1, io1 - 1
            if (wfi%io_s<= io2 .and. io2 <= wfi%io_e) then
              call axpy_wf( &
                & coeff(io2), wf%zwf(:, :, :, ispin, io2, ik, im), &
                & wf_exc_tmp)
            end if
          end do
          if (wfi%if_divide_orbit) then
            call comm_summation(wf_exc_tmp, wf_exc, nsize_rg, wfi%icomm_o)
          else
            call copy_data(wf_exc_tmp, wf_exc)
          endif

          if (wfi%io_s<= io1 .and. io1 <= wfi%io_e) then
            ! Exclude non-orthonormal component:
            call axpy_wf(-one, wf_exc, wf_io1)
            ! Normalization:
            norm2_tmp = real(dot_wf(wf_io1, wf_io1))
            if (wfi%if_divide_rspace) then
              call comm_summation(norm2_tmp, norm2, wfi%icomm_r)
            else
              norm2 = norm2_tmp
            endif
            call scal_wf(one / sqrt(norm2), wf_io1)
            ! Write back to "zwf":
            call copy_data( &
              & wf_io1, &
              & wf%zwf(:, :, :, ispin, io1, ik, im))
          end if
        end do !io1

      end do !ispin
      end do !ik
      end do !im

      return
    contains

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
    subroutine axpy_wf(a, x, y)
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
      do i3 = rg%is(3), rg%ie(3)
        do i2 = rg%is(2), rg%ie(2)
          do i1 = rg%is(1), rg%ie(1)
            y(i1, i2, i3) = a * x(i1, i2, i3) + y(i1, i2, i3)
          end do
        end do
      end do
      !$omp end parallel do
      return
    end subroutine axpy_wf


    ! Scales a wavefunction by a constant
    subroutine scal_wf(a, x)
      implicit none
      complex(8), intent(in) :: a
      complex(8), intent(inout) :: x( &
        & rg%is_array(1):rg%ie_array(1), &
        & rg%is_array(2):rg%ie_array(2), &
        & rg%is_array(3):rg%ie_array(3))
      integer :: i1, i2, i3
      !$omp parallel do collapse(2) default(shared) private(i1,i2,i3)
      do i3 = rg%is(3), rg%ie(3)
        do i2 = rg%is(2), rg%ie(2)
          do i1 = rg%is(1), rg%ie(1)
            x(i1, i2, i3) = a * x(i1, i2, i3)
          end do
        end do
      end do
      !$omp end parallel do
      return
    end subroutine scal_wf

  end subroutine gram_schmidt_col_complex8

end module gram_schmidt_orth
