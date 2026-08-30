!
!  Copyright 2026 SALMON developers
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
!-----------------------------------------------------------------------------------------
module builtin_pw92
  implicit none
  private

  public :: pw92_correlation

  ! Perdew-Wang 1992 unpolarized correlation, eps_c(r_s) = -2A(1+a1 r_s) ln(1 + 1/(2A Q(r_s)))
  real(8),parameter :: pw_a  = 0.0310907d0
  real(8),parameter :: pw_a1 = 0.21370d0
  real(8),parameter :: pw_b1 = 7.59570d0
  real(8),parameter :: pw_b2 = 3.58760d0
  real(8),parameter :: pw_b3 = 1.63820d0
  real(8),parameter :: pw_b4 = 0.49294d0

contains

  ! eps_c at zeta = 0 with its first two r_s derivatives.
  subroutine pw92_correlation(rs, ec, dec_drs, d2ec_drs2)
    implicit none
    real(8),intent(in)  :: rs
    real(8),intent(out) :: ec, dec_drs, d2ec_drs2
    real(8) :: sqrt_rs, qden, dqden, d2qden, mden, logs, dlogs, d2logs, apoly

    sqrt_rs = sqrt(rs)
    qden    = sqrt_rs * (pw_b1 + pw_b3 * rs) + rs * (pw_b2 + pw_b4 * rs)
    dqden   = 0.5d0 * pw_b1 / sqrt_rs + pw_b2 + 1.5d0 * pw_b3 * sqrt_rs + 2d0 * pw_b4 * rs
    d2qden  = -0.25d0 * pw_b1 / (rs * sqrt_rs) + 0.75d0 * pw_b3 / sqrt_rs + 2d0 * pw_b4

    ! d/drs of ln(1 + 1/(2 A Q)) is -Q'/M with M = Q (2 A Q + 1)
    mden   = qden * (2d0 * pw_a * qden + 1d0)
    logs   = log(1d0 + 1d0 / (2d0 * pw_a * qden))
    dlogs  = -dqden / mden
    d2logs = -d2qden / mden + dqden * dqden * (4d0 * pw_a * qden + 1d0) / (mden * mden)

    apoly     = 1d0 + pw_a1 * rs
    ec        = -2d0 * pw_a * apoly * logs
    dec_drs   = -2d0 * pw_a * (pw_a1 * logs + apoly * dlogs)
    d2ec_drs2 = -2d0 * pw_a * (2d0 * pw_a1 * dlogs + apoly * d2logs)

    return
  end subroutine pw92_correlation

end module builtin_pw92
