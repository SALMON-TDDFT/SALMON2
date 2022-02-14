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

module gram_schmidt_so_sub

  use structures, only: s_dft_system, s_rgrid, s_parallel_info, s_orbital
  use pack_unpack, only: copy_data
  use timer
  use salmon_global, only: yn_gramschmidt_blas
  use communication, only: comm_bcast, comm_summation

  implicit none

  private
  public :: gram_schmidt_so

  complex(8),allocatable :: wf_io1(:,:,:,:)
  complex(8),allocatable :: wf_exc(:,:,:,:)
  complex(8),allocatable :: wf_exc_tmp(:,:,:,:)
  complex(8),allocatable :: coeff_tmp(:), coeff(:)

  complex(8),parameter :: zero=(0.0d0,0.0d0)
  complex(8),parameter :: one =(1.0d0,0.0d0)

  real(8),allocatable :: rwf_io1(:,:,:,:)
  real(8),allocatable :: rwf_exc(:,:,:,:)
  real(8),allocatable :: rwf_exc_tmp(:,:,:,:)
  real(8),allocatable :: rcoeff_tmp(:), rcoeff(:)

contains

  subroutine gram_schmidt_so( sys, rg, wfi, wf )
    implicit none
    type(s_dft_system),       intent(in)    :: sys
    type(s_rgrid),            intent(in)    :: rg
    type(s_parallel_info)   , intent(in)    :: wfi
    type(s_orbital),          intent(inout) :: wf

    call timer_begin(LOG_CALC_GRAM_SCHMIDT)

    if(yn_gramschmidt_blas=='y')then
      call gram_schmidt_so_col_cblas(sys, rg, wfi, wf)
    else
      call gram_schmidt_so_col_complex8(sys, rg, wfi, wf)
    end if

    call timer_end(LOG_CALC_GRAM_SCHMIDT)

    return

  end subroutine gram_schmidt_so

  subroutine gram_schmidt_so_col_cblas(sys, rg, wfi, wf)
    ! Only for the colinear L(S)DA:
    use structures, only: s_dft_system, s_rgrid, s_parallel_info, s_orbital
    use pack_unpack, only: copy_data
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

    integer :: nsize_rg_so
    integer :: ik, im, m, m4, n4
    integer :: io, jo, jo1, jo2
    integer :: idiv1, idiv2, idiv3, io0_s, io0_e, io1_s, io1_e, io2_s, io2_e, io3_s, io3_e
    integer :: numo0, numo1, numo1_0, numo1_1, numo2, numo2_0, numo2_1, numo3, numo3_0, numo3_1
    complex(8) :: coeff(wfi%numo), coeff_tmp(wfi%numo)
    real(8) :: norm2, norm2_tmp
    complex(8) :: wf_block(rg%is(1):rg%ie(1), rg%is(2):rg%ie(2), rg%is(3):rg%ie(3), &
                           lbound(wf%zwf,4):ubound(wf%zwf,4), wfi%numo)
    complex(8) :: wf_block1(rg%is(1):rg%ie(1), rg%is(2):rg%ie(2), rg%is(3):rg%ie(3), &
                            lbound(wf%zwf,4):ubound(wf%zwf,4), wfi%numo)
    complex(8) :: wf_block_send(rg%is(1):rg%ie(1), rg%is(2):rg%ie(2), rg%is(3):rg%ie(3), &
                            lbound(wf%zwf,4):ubound(wf%zwf,4), wfi%numo_max)
    complex(8) :: umat(wfi%numo_max,wfi%numo), umat_tmp(wfi%numo_max,wfi%numo)
    complex(8) :: ZDOTC

    m4 = lbound(wf%zwf,4); n4 = ubound(wf%zwf,4)
    nsize_rg_so = (rg%ie(1)-rg%is(1)+1)*(rg%ie(2)-rg%is(2)+1)*(rg%ie(3)-rg%is(3)+1)*(n4-m4+1)

    do im = wfi%im_s, wfi%im_e
    do ik = wfi%ik_s, wfi%ik_e

      ! Copy wave function
      do io = wfi%io_s, wfi%io_e
        jo = io - wfi%io_s + 1
        call copy_data( &
          & wf%zwf(rg%is(1):rg%ie(1), rg%is(2):rg%ie(2), rg%is(3):rg%ie(3), m4:n4, io, ik, im), &
          & wf_block(:, :, :, :, jo))
      end do

      ! Loop for all:
      do m = 0, wfi%nporbital-1
        if(wfi%id_o==m)then

          if(wfi%numo < 32) then

        ! Loop for orbit #io1 in my rank:
            do jo1 = 1, wfi%numo

              ! normalize orbital jo1
              norm2_tmp = real( ZDOTC( &
               &      nsize_rg_so, wf_block(:,:,:,:,jo1), n_one, &
               &                wf_block(:,:,:,:,jo1), n_one )) &
               &     * sys%hvol
              if (wfi%if_divide_rspace) then
                call comm_summation(norm2_tmp, norm2, wfi%icomm_r)
              else
                norm2 = norm2_tmp
              endif
              call ZDSCAL(nsize_rg_so, one/sqrt(norm2), wf_block(:,:,:,:,jo1), n_one)

              ! Calculate overlap coefficients:
              coeff_tmp = 0d0
              do jo2 = jo1+1, wfi%numo
                coeff_tmp(jo2) = ZDOTC( &
                 &      nsize_rg_so, wf_block(:,:,:,:,jo1), n_one, &
                 &                wf_block(:,:,:,:,jo2), n_one ) &
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
                  &  nsize_rg_so, - coeff(jo2), wf_block(:,:,:,:,jo1), n_one, &
                  &                          wf_block(:,:,:,:,jo2), n_one)
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
                     &      nsize_rg_so, wf_block(:,:,:,:,jo1), n_one, &
                     &                wf_block(:,:,:,:,jo1), n_one )) &
                     &     * sys%hvol
                    if (wfi%if_divide_rspace) then
                      call comm_summation(norm2_tmp, norm2, wfi%icomm_r)
                    else
                      norm2 = norm2_tmp
                    endif
                    call ZDSCAL(nsize_rg_so, one/sqrt(norm2), wf_block(:,:,:,:,jo1), n_one)

                ! Calculate overlap coefficients:
                    coeff_tmp = 0d0
                    do jo2 = jo1+1, io3_e
                      coeff_tmp(jo2) = ZDOTC( &
                       &      nsize_rg_so, wf_block(:,:,:,:,jo1), n_one, &
                       &                wf_block(:,:,:,:,jo2), n_one ) &
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
                        &  nsize_rg_so, - coeff(jo2), wf_block(:,:,:,:,jo1), n_one, &
                        &                          wf_block(:,:,:,:,jo2), n_one)
                    end do

                  end do !jo1

                  if(idiv3 == 0) then

                    umat_tmp = 0.d0
                    wf_block_send(:,:,:,:,1:numo3_0) = wf_block(:,:,:,:,io2_s:io2_s-1+numo3_0)
                    wf_block1(:,:,:,:,1:numo3_1) = wf_block(:,:,:,:,io2_s+numo3_0:io2_e)
                    call zgemm(TRANSA, TRANSB, numo3_0, numo3_1, nsize_rg_so,  &
                        &  one*sys%hvol, wf_block_send, nsize_rg_so,  &
                        &                wf_block1, nsize_rg_so,  &
                        &          zero, umat_tmp, wfi%numo_max)

                   if(wfi%if_divide_rspace) then
                      call comm_summation(umat_tmp, umat, wfi%numo_max*wfi%numo, wfi%icomm_r)
                    else
                      umat = umat_tmp
                    end if

                    call zgemm(TRANSB, TRANSB, nsize_rg_so, numo3_1, numo3_0,  &
                      &         - one, wf_block_send, nsize_rg_so,  &
                      &                umat, wfi%numo_max,  &
                      &           one, wf_block1, nsize_rg_so)
                    wf_block(:,:,:,:,io2_s+numo3_0:io2_e) = wf_block1(:,:,:,:,1:numo3_1)
                  end if ! idiv3 == 0

                end do !idiv3

                if(idiv2 == 0) then

                  umat_tmp = 0.d0
                  wf_block_send(:,:,:,:,1:numo2_0) = wf_block(:,:,:,:,io1_s:io1_s-1+numo2_0)
                  wf_block1(:,:,:,:,1:numo2_1) = wf_block(:,:,:,:,io1_s+numo2_0:io1_e)
                  call zgemm(TRANSA, TRANSB, numo2_0, numo2_1, nsize_rg_so,  &
                      &  one*sys%hvol, wf_block_send, nsize_rg_so,  &
                      &                wf_block1, nsize_rg_so,  &
                      &          zero, umat_tmp, wfi%numo_max)

                  if(wfi%if_divide_rspace) then
                    call comm_summation(umat_tmp, umat, wfi%numo_max*wfi%numo, wfi%icomm_r)
                  else
                    umat = umat_tmp
                  end if

                  call zgemm(TRANSB, TRANSB, nsize_rg_so, numo2_1, numo2_0,  &
                    &         - one, wf_block_send, nsize_rg_so,  &
                    &                umat, wfi%numo_max,  &
                    &           one, wf_block1, nsize_rg_so)

                  wf_block(:,:,:,:,io1_s+numo2_0:io1_e) = wf_block1(:,:,:,:,1:numo2_1)
                end if ! idiv2 == 0

              end do !idiv2

              if(idiv1 == 0) then

                umat_tmp = 0.d0
                wf_block_send(:,:,:,:,1:numo1_0) = wf_block(:,:,:,:,io0_s:io0_s-1+numo1_0)
                wf_block1(:,:,:,:,1:numo1_1) = wf_block(:,:,:,:,io0_s+numo1_0:io0_e)

                call zgemm(TRANSA, TRANSB, numo1_0, numo1_1, nsize_rg_so,  &
                    &  one*sys%hvol, wf_block_send, nsize_rg_so,  &
                    &                wf_block1, nsize_rg_so,  &
                    &          zero, umat_tmp, wfi%numo_max)

                if(wfi%if_divide_rspace) then
                  call comm_summation(umat_tmp, umat, wfi%numo_max*wfi%numo, wfi%icomm_r)
                else
                  umat = umat_tmp
                end if

                call zgemm(TRANSB, TRANSB, nsize_rg_so, numo1_1, numo1_0,  &
                  &         - one, wf_block_send, nsize_rg_so,  &
                  &                umat, wfi%numo_max,  &
                  &           one, wf_block1, nsize_rg_so)
                wf_block(:,:,:,:,io0_s+numo1_0:io0_e) = wf_block1(:,:,:,:,1:numo1_1)
              end if ! idiv1 == 0

            end do !idiv1

          end if ! if(wfi%numo < 32) else

        end if !m==wfi%id_o

        if (wfi%if_divide_orbit) then
          call copy_data( wf_block(:,:,:,:,:), wf_block_send(:,:,:,:,1:wfi%numo) )
          wf_block_send(:,:,:,:,wfi%numo+1:wfi%numo_max) = 0d0
          call comm_bcast( wf_block_send, wfi%icomm_o, wfi%irank_io(wfi%io_s_all(m)))

          umat_tmp=0.d0
          if( wfi%id_o > m )then
            call zgemm(TRANSA, TRANSB, wfi%numo_max, wfi%numo, nsize_rg_so,  &
              &   one*sys%hvol, wf_block_send(:,:,:,:,1), nsize_rg_so,  &
              &                wf_block(:,:,:,:,1), nsize_rg_so,  &
              &          zero, umat_tmp, wfi%numo_max )

            if(wfi%if_divide_rspace) then
              call comm_summation(umat_tmp, umat, wfi%numo_max*wfi%numo, wfi%icomm_r)
            else
              umat = umat_tmp
            end if

            call zgemm(TRANSB, TRANSB, nsize_rg_so, wfi%numo, wfi%numo_max,  &
              &         - one, wf_block_send(:,:,:,:,1), nsize_rg_so,  &
              &                umat, wfi%numo_max,  &
              &           one, wf_block(:,:,:,:,1), nsize_rg_so)
          end if
        end if ! wfi%if_divide_orbit

      end do !m

      do io = wfi%io_s, wfi%io_e
        jo = io - wfi%io_s + 1
        call copy_data( &
          & wf_block(:, :, :, :, jo), &
          & wf%zwf(rg%is(1):rg%ie(1), rg%is(2):rg%ie(2), rg%is(3):rg%ie(3), m4:n4, io, ik, im))
      end do

    end do !ik
    end do !im

    return

  end subroutine gram_schmidt_so_col_cblas

!=======================================================================

  subroutine gram_schmidt_so_col_complex8( sys, rg, wfi, wf )
    use structures, only: s_dft_system, s_rgrid, s_parallel_info, s_orbital
    implicit none
    type(s_dft_system),       intent(in)    :: sys
    type(s_rgrid),            intent(in)    :: rg
    type(s_parallel_info)   , intent(in)    :: wfi
    type(s_orbital),          intent(inout) :: wf

    complex(8),allocatable :: wf_io1(:,:,:,:)
    complex(8),allocatable :: wf_exc(:,:,:,:)
    complex(8),allocatable :: wf_exc_tmp(:,:,:,:)
    complex(8),allocatable :: coeff_tmp(:), coeff(:)
  
    complex(8),parameter :: zero=(0.0d0,0.0d0)
    complex(8),parameter :: one =(1.0d0,0.0d0)
  
    integer :: im, ik, io1, io2, m1,m2,m3,m4,n1,n2,n3,n4
    real(8) :: norm2, norm2_tmp

    m1 = rg%is_array(1); n1 = rg%ie_array(1)
    m2 = rg%is_array(2); n2 = rg%ie_array(2)
    m3 = rg%is_array(3); n3 = rg%ie_array(3)
    m4 = lbound(wf%zwf,4); n4 = ubound(wf%zwf,4)

    if ( .not.allocated(wf_io1) ) then
      allocate( wf_io1(m1:n1,m2:n2,m3:n3,m4:n4) )
      allocate( wf_exc(m1:n1,m2:n2,m3:n3,m4:n4) )
      allocate( wf_exc_tmp(m1:n1,m2:n2,m3:n3,m4:n4) )
      wf_io1=zero
      wf_exc=zero
      wf_exc_tmp=zero
      allocate( coeff_tmp(sys%no-1) ); coeff_tmp=zero
      allocate( coeff(sys%no-1) ); coeff=zero
    end if

    do im = wfi%im_s, wfi%im_e
    do ik = wfi%ik_s, wfi%ik_e

      do io1 = 1, sys%no

        ! Retrieve orbit #io1 into "wf_io1":
        if ( wfi%io_s <= io1 .and. io1 <= wfi%io_e ) then
          call copy_data( wf%zwf(:,:,:,:,io1,ik,im), wf_io1 )
        end if
        if ( wfi%if_divide_orbit ) then
          call comm_bcast( wf_io1, wfi%icomm_o, wfi%irank_io(io1) )
        end if

        ! Calculate overlap coefficients: 
        coeff_tmp=zero

        do io2 = 1, io1 - 1
          if ( wfi%io_s <= io2 .and. io2 <= wfi%io_e ) then
            coeff_tmp(io2) = dot_wf( wf%zwf(:,:,:,:,io2,ik,im), wf_io1 )
          end if
        end do 

        if ( wfi%if_divide_rspace .or. wfi%if_divide_orbit ) then
          call comm_summation( coeff_tmp, coeff, size(coeff), wfi%icomm_ro )
        else
          coeff = coeff_tmp
        end if

        ! Calculate exclusion term "wf_exc":
        wf_exc_tmp = 0.0d0
        do io2 = 1, io1 - 1
          if ( wfi%io_s <= io2 .and. io2 <= wfi%io_e ) then
            call axpy_wf( coeff(io2), wf%zwf(:,:,:,:,io2,ik,im), wf_exc_tmp )
          end if
        end do
        if ( wfi%if_divide_orbit ) then
          call comm_summation( wf_exc_tmp, wf_exc, size(wf_exc), wfi%icomm_o )
        else
          call copy_data( wf_exc_tmp, wf_exc )
        end if

        if ( wfi%io_s <= io1 .and. io1 <= wfi%io_e ) then
          ! Exclude non-orthonormal component:
          call axpy_wf( -one, wf_exc, wf_io1 )
          ! Normalization:
          norm2_tmp = real( dot_wf(wf_io1,wf_io1) )
          if ( wfi%if_divide_rspace ) then
            call comm_summation( norm2_tmp, norm2, wfi%icomm_r )
          else
            norm2 = norm2_tmp
          end if
          call scal_wf( one/sqrt(norm2), wf_io1 )
          ! Write back to "zwf":
          call copy_data( wf_io1, wf%zwf(:,:,:,:,io1,ik,im) ) 
        end if
      end do !io1

    end do !ik
    end do !im
    
    return

  contains

  ! Dot product of two wavefunctions:
  complex(8) function dot_wf(x, y) result(p)
    implicit none
    complex(8), intent(in) :: x(rg%is_array(1):,rg%is_array(2):,rg%is_array(3):,:)
    complex(8), intent(in) :: y(rg%is_array(1):,rg%is_array(2):,rg%is_array(3):,:)
    integer :: i1, i2, i3, i4
    p=zero
    !$omp parallel do collapse(2) default(shared) private(i1,i2,i3) reduction(+:p)
    do i4 = 1, size(x,4)
      do i3 = rg%is(3), rg%ie(3)
        do i2 = rg%is(2), rg%ie(2)
          do i1 = rg%is(1), rg%ie(1)
            p = p + conjg( x(i1,i2,i3,i4) )*y(i1,i2,i3,i4)
          end do
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
    complex(8), intent(in) :: x(rg%is_array(1):,rg%is_array(2):,rg%is_array(3):,:)
    complex(8), intent(inout) :: y(rg%is_array(1):,rg%is_array(2):,rg%is_array(3):,:)
    integer :: i1, i2, i3, i4
    !$omp parallel do collapse(2) default(shared) private(i1,i2,i3)
    do i4 = 1, size(x,4)
      do i3 = rg%is(3), rg%ie(3)
        do i2 = rg%is(2), rg%ie(2)
          do i1 = rg%is(1), rg%ie(1)
            y(i1, i2, i3, i4) = a * x(i1, i2, i3, i4) + y(i1, i2, i3, i4)
          end do
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
    complex(8), intent(inout) :: x(rg%is_array(1):,rg%is_array(2):,rg%is_array(3):,:)
    integer :: i1, i2, i3, i4
    !$omp parallel do collapse(2) default(shared) private(i1,i2,i3)
    do i4 = 1, size(x,4)
      do i3 = rg%is(3), rg%ie(3)
        do i2 = rg%is(2), rg%ie(2)
          do i1 = rg%is(1), rg%ie(1)
            x(i1, i2, i3, i4) = a * x(i1, i2, i3, i4)
          end do
        end do
      end do
    end do
    !$omp end parallel do
    return
  end subroutine scal_wf

  end subroutine gram_schmidt_so_col_complex8

end module gram_schmidt_so_sub
