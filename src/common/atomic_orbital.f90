!
!  Copyright 2019-2026 SALMON developers
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
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130

#include "config.h"

module atomic_orbital

  implicit none

  private

  ! Step 0: basis table and grid expansion
  ! Step 1: overlap integrals S(R)
  ! Step 2: projection coefficients
  ! Step 3: Löwdin transformation matrix
  ! Step 4: orthogonalized projection coefficients

  public :: verify_atomic_orbital_consistency !!DEBUG

contains

  !> Test routine to verify consistency of the current variable implementation
  !! This subroutine validates the data structures and basic properties
  !! and will be removed before release.
  subroutine verify_atomic_orbital_consistency(system, pp, lg, mg, nproc_id) !!DEBUG
    use structures, only: s_dft_system, s_pp_info, s_rgrid !!DEBUG
    use salmon_global, only: nelem !!DEBUG
    use math_constants, only: pi !!DEBUG
    implicit none !!DEBUG
    type(s_dft_system), intent(inout) :: system !!DEBUG
    type(s_pp_info), intent(inout) :: pp !!DEBUG
    type(s_rgrid), intent(inout) :: lg, mg !!DEBUG
    integer, intent(in) :: nproc_id !!DEBUG

    integer :: ik, il, ir !!DEBUG
    real(8) :: c_l, tmp, r_i0, r_i1, upp_i0, upp_i1, dr, r_i, u_l !!DEBUG
    character(64) :: filename !!DEBUG

    if (nproc_id /= 0) return !!DEBUG

    write(*,*) !!DEBUG
    write(*,'(a)') '===== VERIFICATION: atomic_orbital consistency =====' !!DEBUG

    ! Check system data
    write(*,'(a,i0)') 'system%nion = ', system%nion !!DEBUG
    write(*,'(a,i0)') 'system%no = ', system%no !!DEBUG
    write(*,'(a,i0)') 'system%nk = ', system%nk !!DEBUG
    write(*,'(a,i0)') 'system%nspin = ', system%nspin !!DEBUG

    ! Check pseudopotential info
    write(*,'(a,i0)') 'nelem = ', nelem !!DEBUG
    write(*,'(a,i0)') 'pp%nrmax = ', pp%nrmax !!DEBUG
    write(*,'(a,i0)') 'pp%lmax = ', pp%lmax !!DEBUG

    ! Check grid data
    write(*,'(a,i0,a,i0,a,i0)') 'lg%num(1:3) = ', lg%num(1), ', ', lg%num(2), ', ', lg%num(3) !!DEBUG
    write(*,'(a,i0,a,i0,a,i0)') 'mg%num(1:3) = ', mg%num(1), ', ', mg%num(2), ', ', mg%num(3) !!DEBUG

    ! Check that pp allocations are present
    if (.not. allocated(pp%has_rho_pp)) then !!DEBUG
      write(*,'(a)') 'ERROR: pp%has_rho_pp is not allocated' !!DEBUG
      stop !!DEBUG
    end if !!DEBUG
    if (.not. allocated(pp%has_wf_pp)) then !!DEBUG
      write(*,'(a)') 'ERROR: pp%has_wf_pp is not allocated' !!DEBUG
      stop !!DEBUG
    end if !!DEBUG
    if (.not. allocated(pp%has_proj_pp)) then !!DEBUG
      write(*,'(a)') 'ERROR: pp%has_proj_pp is not allocated' !!DEBUG
      stop !!DEBUG
    end if !!DEBUG

    ! Report flags
    do ik = 1, nelem !!DEBUG
      write(*,'(a,i0,a,l1,a,l1,a,l1)') 'Element ', ik, ': has_rho_pp=', pp%has_rho_pp(ik), &
        ', has_wf_pp=', pp%has_wf_pp(ik), ', has_proj_pp=', pp%has_proj_pp(ik) !!DEBUG
    end do !!DEBUG

    write(*,'(a)') '===== VERIFICATION: pseudowavefunction norm =====' !!DEBUG

    ! Verify norm of u_l(r) restored from upp_f: u_l = upp_f*r**(l+1)/c_l, c_l=sqrt((2l+1)/(4pi)) !!DEBUG
    ! (upp(i) refers to rad(i+1)) !!DEBUG
    do ik = 1, nelem !!DEBUG
      do il = 0, pp%mlps(ik) !!DEBUG
        c_l = sqrt((2d0*il+1d0)/(4d0*pi)) !!DEBUG
        tmp = 0.0d0 !!DEBUG
        do ir = 1, pp%nrmax - 1 !!DEBUG
          r_i0 = pp%rad(ir+0, ik) !!DEBUG
          r_i1 = pp%rad(ir+1, ik) !!DEBUG
          upp_i0 = pp%upptbl_ao_l(ir+0, il, ik) * r_i0**(il+1) / c_l !!DEBUG
          upp_i1 = pp%upptbl_ao_l(ir+1, il, ik) * r_i1**(il+1) / c_l !!DEBUG
          dr = r_i1 - r_i0 !!DEBUG
          tmp = tmp + 0.5d0 * (upp_i0**2 + upp_i1**2) * dr !!DEBUG
        end do !!DEBUG
        write(*,'(a,i0,a,i0,a,e15.8)') 'ik=', ik, ' l=', il, ' norm=', tmp !!DEBUG
      end do !!DEBUG
    end do !!DEBUG

    ! Dump restored u_l(r) = upptbl_ao_l*r**(l+1)/c_l to upp<ik>_<l>.txt !!DEBUG
    do ik = 1, nelem !!DEBUG
      do il = 0, pp%mlps(ik) !!DEBUG
        c_l = sqrt((2d0*il+1d0)/(4d0*pi)) !!DEBUG
        write(filename, '(a,i0,a,i0,a)') 'upp', ik, '_', il, '.txt' !!DEBUG
        open(99, file=trim(filename), action='write') !!DEBUG
        do ir = 1, pp%nrmax - 1 !!DEBUG
          r_i = pp%rad(ir, ik) !!DEBUG
          u_l = pp%upptbl_ao_l(ir, il, ik) * r_i**(il+1) / c_l !!DEBUG
          write(99, *) r_i, u_l !!DEBUG
        end do !!DEBUG
        close(99) !!DEBUG
      end do !!DEBUG
    end do !!DEBUG

    write(*,'(a)') '===== VERIFICATION: completed successfully =====' !!DEBUG
    write(*,*) !!DEBUG

  end subroutine verify_atomic_orbital_consistency !!DEBUG

end module atomic_orbital
