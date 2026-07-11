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
!===================================================================
! Length-gauge degeneracy-mode semantic predicates.
!
! The length-gauge SBE degeneracy handling has FOUR distinct machinery layers,
! and the string value 'sbe_lg_degen' selects a subset of them.  Testing the
! raw string ('gicov') at every call site conflates layers that must diverge
! once a second X-full mode ('gicov_int') exists: the finite-difference and the
! integral forms SHARE the X-full representation (single full-band block,
! reduced 3-column prod_dk, qnm==rho) but DIFFER in propagation (FD covariant
! stencil + RK4 vs. bounded Wilson transport + midpoint exact-exponential) and
! in the T2 term (static frozen-gap gate array vs. per-step moving-gap gate).
! These four pure predicates name the layers so each call site selects exactly
! the machinery it needs:
!
!   uses_prod_dk       reads the k-space overlap file file_sbe_prod_dk
!                      (every covariant mode: gi / gifix / gicov / gicov_int)
!   uses_xfull_links   uses the X-full single-block Wilson transport + the
!                      reduced 3-column prod_dk + the qnm==rho representation
!                      (gicov and gicov_int)
!   uses_fd_gicov      uses the finite-difference covariant-derivative RHS
!                      (gicov_rhs) + self-starting RK4 + the STATIC frozen-gap
!                      T2 gate array (gicov only)
!   uses_integral_gicov  uses the integral (covariant-Houston) transport
!                        propagator + the per-step moving-gap T2 gate
!                        (gicov_int only)
!
! Replacing the scattered trim(sbe_lg_degen)=='gicov' tests by these predicates
! is byte-identical for every pre-existing input (none is 'gicov_int'), so the
! default and existing paths are unchanged (guaranteed structurally here; the
! regression byte test guards the shared edited files).
!===================================================================
module sbe_lg_mode_ssbe
  implicit none
  private
  public :: uses_prod_dk
  public :: uses_xfull_links
  public :: uses_fd_gicov
  public :: uses_integral_gicov

contains

  pure logical function uses_prod_dk(mode)
    implicit none
    character(*), intent(in) :: mode
    uses_prod_dk = (trim(mode) == 'gi') .or. (trim(mode) == 'gifix') &
                 & .or. (trim(mode) == 'gicov') .or. (trim(mode) == 'gicov_int')
  end function uses_prod_dk

  pure logical function uses_xfull_links(mode)
    implicit none
    character(*), intent(in) :: mode
    uses_xfull_links = (trim(mode) == 'gicov') .or. (trim(mode) == 'gicov_int')
  end function uses_xfull_links

  pure logical function uses_fd_gicov(mode)
    implicit none
    character(*), intent(in) :: mode
    uses_fd_gicov = (trim(mode) == 'gicov')
  end function uses_fd_gicov

  pure logical function uses_integral_gicov(mode)
    implicit none
    character(*), intent(in) :: mode
    uses_integral_gicov = (trim(mode) == 'gicov_int')
  end function uses_integral_gicov

end module sbe_lg_mode_ssbe
