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


    subroutine gram_schmidt_col_real8
      implicit none
       ! Only for the colinear L(S)DA:
      real(8) :: tmp(rg%is(1):rg%ie(1), rg%is(2):rg%ie(2), rg%is(3):rg%ie(3))
      
      do im = sys%im_s, sys%im_e
        do is = 1, sys%nspin
          do ik = sys%ik_s, sys%ik_e
            do io = 1, sys%no
              if (wfi%io_s <= io .and. io <= wfi%io_e) then
                ! orbital #io is stored in the present node:
                call copy_data( &
                  wf%rwf(rg%is(1):rg%ie(1), rg%is(2):rg%ie(2), rg%is(3):rg%ie(3), is, io, ik, im), tmp &
                )
              end if
              ! store root_table
              call comm_bcast(tmp, iroot)

              do jo = 1, io - 1
                if (wfi%io_s <= jo .and. jo <= wfi%io_e) then
                  
                end if
              end do
            end do
          end do
        end do

    

    

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
