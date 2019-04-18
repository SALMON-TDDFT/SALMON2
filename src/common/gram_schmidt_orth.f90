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
      call gram_schmidt_col_real8()
    elseif (allocated(wf%zwf)) then
      call gram_schmidt_col_complex8()
    else
      stop "Wavefunctions are not allocated!"
    end if

    call timer_end(LOG_GRAM_SCHMIDT)

    return
    continue

    logical function has_orbit(jo)
      implicit none
      integer, intent(in) :: jo
      has_orbit = (wfi%io_s <= wfi%jo_tbl(jo)) &
          & .and. (wfi%jo_tbl(jo) <= wfi%io_e)
      return
    end function has_orbit

    real(8) function prod_real8(tmp, is, io, ik, im)
      implicit none
      real(8), intent(in) :: tmp(rg%is(1):rg%ie(1),rg%is(2):rg%ie(2),rg%is(3):rg%ie(3))
      integer, intent(in) :: is, io, ik, im
      integer :: i1, i2, i3
      prod_real8 = 0d0
!$omp parallel do collapse(2) default(shared) private(i1, i2, i3) reduction(+: prod_real8)
      do i3 = rg%is(3), rg%ie(3)
        do i2 = rg%is(2), rg%ie(2)
          do i1 = rg%is(1), rg%ie(1)
            prod_real8 = prod_real8 + &
            & tmp(i1, i2, i3) * wf%rwf(i1, i2, i3, is, io, ik, im)
          end do
        end do
      end do
  !$omp end parallel do
      return
    end function prod_real8

    subroutine add_real8(tmp, coeff, is, io, ik, im)
      implicit none
      real(8), intent(in) :: tmp(rg%is(1):rg%ie(1),rg%is(2):rg%ie(2),rg%is(3):rg%ie(3))
      real(8), intent(in) :: coeff
      integer, intent(in) :: io
      integer :: i1, i2, i3
      prod_real8 = 0d0
!$omp parallel do collapse(2) default(shared) private(i1, i2, i3)
      do i3 = rg%is(3), rg%ie(3)
        do i2 = rg%is(2), rg%ie(2)
          do i1 = rg%is(1), rg%ie(1)
            tmp(i1, i2, i3) = tmp(i1, i2, i3) &
            & + coeff * wf%rwf(i1, i2, i3, is, io, ik, im)
          end do
        end do
      end do
  !$omp end parallel do
      return
    end subroutine add_real8



    subroutine gram_schmidt_col_real8
      implicit none
       ! Only for the colinear L(S)DA:
      real(8) :: tmp(rg%is(1):rg%ie(1), rg%is(2):rg%ie(2), rg%is(3):rg%ie(3))
      
      do im = sys%im_s, sys%im_e
        do is = 1, sys%nspin
          do ik = sys%ik_s, sys%ik_e
            do jo1 = 1, sys%no
              ! Retrieve #jo1-th orbit into tmp variable:
              if (has_orbit(jo1))
                call copy_data( &
                & wf%rwf( &
                &   rg%is(1):rg%ie(1), &
                &   rg%is(2):rg%ie(2), &
                &   rg%is(3):rg%ie(3), &
                &   is, io, ik, im), &
                & tmp)
              call comm_bcast(tmp, wfi%irank_jo(jo1))
              
              ! Calculate overlap coefficients: c(1:jo1-1)
              c_ovlp_tmp(:) = 0d0
              do jo2 = 1, jo1 - 1
                if (has_orbit(jo2)) then
                  io2 = wfi%jo_tbl(jo2)
                  c_ovlp_tmp(jo2) = prod_real8(tmp, io2) 
                end if
              end do
              call comm_sum_all(c_ovlp_tmp, c_ovlp, wfi%icomm_ro)

              ! Calculate exclusion term:
              tmp2 = 0d0
              do jo2 = 1, jo1 - 1
                if (has_orbit(jo2)) then
                  io2 = wfi%jo_tbl(jo2)
                  call add_real8(tmp2, c_ovlp(jo2), is, io, ik, im)
                end if
              end do

              call comm_sum(tmp2, tmp3, wfi%icomm_o, wfi%irank_jo(jo1)))

              tmp = tmp - tmp3

              ! Calculate norm term:
              norm_tmp = snorm()
              call comm_sum_all(norm_tmp, norm, comm_r)

              tmp = tmp * (1.0 / sqrt(norm2))
            end do !jo1
          end do !ik
        end do !is
      end do !im
    end subroutine gram_schmidt_real8

    

    

  end subroutine

  ! Gram_Schmitd for orbital(b) omp
  subroutine gram_schmidt_omp_b()
    implicit none
    return
  end subroutine gram_schmidt_omp_b

  ! Gram_Schmitd for k-space omp
  subroutine gram_schmidt_omp_k()
    implicit none
    return
  end subroutine gram_schmidt_omp_k

  ! Gram_Schmitd for realspace(r) omp
  subroutine gram_schmidt_omp_r()
    implicit none
    return
  end subroutine gram_schmidt_omp_r

end module gram_schmidt_orth
