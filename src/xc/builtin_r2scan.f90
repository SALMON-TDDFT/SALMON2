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
module builtin_r2scan
  implicit none
  private

  public :: exc_cor_r2scan
  public :: grad_floor

  ! interpolation in alpha, exchange and correlation
  real(8),parameter :: cx_poly(0:7) = (/  1d0             , -0.667d0          , -0.4445555d0      , &
                                         -0.663086601049d0,  1.451297044490d0 , -0.887998041597d0 , &
                                          0.234528941479d0, -0.023185843322d0 /)
  real(8),parameter :: cc_poly(0:7) = (/  1d0             , -0.64d0           , -0.4352d0         , &
                                         -1.535685604549d0,  3.061560252175d0 , -1.915710236206d0 , &
                                          0.516884468372d0, -0.051848879792d0 /)
  real(8),parameter :: c1x = 0.667d0, c2x = 0.8d0, dx_switch = 1.24d0
  real(8),parameter :: c1c = 0.64d0 , c2c = 1.5d0, dc_switch = 0.7d0
  real(8),parameter :: alpha_poly_max = 2.5d0

  ! exchange enhancement
  real(8),parameter :: kappa0 = 0.174d0
  real(8),parameter :: kappa1 = 0.065d0
  real(8),parameter :: mu_gea = 10d0 / 81d0
  real(8),parameter :: a_damp = 4.9479d0
  real(8),parameter :: d_damp = 0.361d0
  real(8),parameter :: d_damp4 = d_damp**4
  real(8),parameter :: eta_reg = 1d-3

  ! correlation
  real(8),parameter :: b1c = 0.0285764d0
  real(8),parameter :: b2c = 0.0889d0
  real(8),parameter :: b3c = 0.125541d0
  real(8),parameter :: chi_inf = 0.12802585262625815d0
  real(8),parameter :: beta0 = 0.066725d0
  real(8),parameter :: beta_num = 0.1d0
  real(8),parameter :: beta_den = 0.1778d0
  real(8),parameter :: gamma_c = 0.031090690869655d0      ! (1 - ln 2)/pi^2
  real(8),parameter :: pw_a  = 0.0310907d0
  real(8),parameter :: pw_a1 = 0.21370d0
  real(8),parameter :: pw_b1 = 7.59570d0
  real(8),parameter :: pw_b2 = 3.58760d0
  real(8),parameter :: pw_b3 = 1.63820d0
  real(8),parameter :: pw_b4 = 0.49294d0

  ! geometry of the reduced variables
  real(8),parameter :: c_rs      = 0.62035049089940001665d0   ! (3/(4 pi))^(1/3)
  real(8),parameter :: c_tau_hom = 2.8712340001881918160d0    ! (3/10)(3 pi^2)^(2/3)
  real(8),parameter :: c_s2      = 0.026121172985233599568d0  ! 1/(4 (3 pi^2)^(2/3))
  real(8),parameter :: c_t2grad  = 0.063468206097703704205d0  ! pi/(16 (3 pi^2)^(1/3))
  real(8),parameter :: c_ex      = -0.73855876638202240587d0  ! -(3/4)(3/pi)^(1/3)

  ! constants fixed by the second-order gradient expansion
  real(8),parameter :: fx_one = cx_poly(1) + 2d0*cx_poly(2) + 3d0*cx_poly(3) + 4d0*cx_poly(4) &
                              + 5d0*cx_poly(5) + 6d0*cx_poly(6) + 7d0*cx_poly(7)
  real(8),parameter :: fc_one = cc_poly(1) + 2d0*cc_poly(2) + 3d0*cc_poly(3) + 4d0*cc_poly(4) &
                              + 5d0*cc_poly(5) + 6d0*cc_poly(6) + 7d0*cc_poly(7)
  real(8),parameter :: c_eta   = 20d0 / 27d0 + 5d0 * eta_reg / 3d0
  real(8),parameter :: c_xgrad = c_eta * kappa0 * fx_one

  ! regularizations.  Each one is applied to the value and to the derivatives of the same
  ! expression, so that a clamped point is evaluated as the constant branch it was clamped
  ! onto and not as a value whose derivative belongs to a different branch.
  real(8),parameter :: rho_floor = 1d-18   ! below this the functional returns zero
  real(8),parameter :: grad_floor = 1d-18  ! |grad n| below this has no defined direction
  real(8),parameter :: s2_floor = 1d-16    ! g_x(s2) has an essential singularity at s2 = 0
  real(8),parameter :: w_floor = 1d-30     ! w_1 -> 0 as eps_c^PW92 -> 0 in the density tail

contains

  ! Grid driver: rho includes any partial core; rho_s, tau_s are the spin-resolved
  ! halves, restored to totals here for this spin-unpolarized functional. Returns
  ! f = n*eps_xc and df/dn, df/dtau, df/d|grad n|.
  subroutine exc_cor_r2scan(nl, rho, rho_s, grho_norm, tau_s, eexc, vxc, vtau, vgrad)
    implicit none
    integer,intent(in)  :: nl
    real(8),intent(in)  :: rho(nl), rho_s(nl), grho_norm(nl), tau_s(nl)
    real(8),intent(out) :: eexc(nl), vxc(nl), vtau(nl), vgrad(nl)
    integer :: i

!$omp parallel do private(i)
    do i = 1, nl
      call r2scan_point(max(rho(i), 2d0*max(rho_s(i), 0d0)), max(grho_norm(i), 0d0), &
                        max(2d0*tau_s(i), 0d0), eexc(i), vxc(i), vtau(i), vgrad(i))
    end do
!$omp end parallel do

    return
  end subroutine exc_cor_r2scan


  ! One grid point of r2SCAN: ingredients, exchange, correlation, then the four
  ! outputs. Exchange and correlation share only the ingredients and are
  ! otherwise independent.
  subroutine r2scan_point(rho, grad_norm, tau, eexc, vxc, vtau, vgrad)
    implicit none
    real(8),intent(in)  :: rho, grad_norm, tau
    real(8),intent(out) :: eexc, vxc, vtau, vgrad

    real(8) :: grad2, rs, drs_drho, tau_hom, tau_w, s2, t2
    real(8) :: ds2_drho, ds2_dg2, dt2_drho, dt2_dg2
    real(8) :: alpha, alpha_den, dalpha_drho, dalpha_dg2, dalpha_dtau
    real(8) :: gx_damp, dgx_ds2, s2_damp, xs, dxs_ds2, hx0, hx1, dhx1_dxs
    real(8) :: fx_alpha, dfx_alpha, fx_grad, fx_enh, dq_ds2, dq_dalpha, dfx_enh_ds2
    real(8) :: ex_dens, dex_dens, dex_drho, dex_dg2, dex_dtau
    real(8) :: ec_lda, dec_lda, d2ec_lda, ec_lda0, dec_lda0, d2ec_lda0
    real(8) :: wc0, dwc0_drs, gc0, dgc0_ds2, sc0, hc0, dhc0_dwc0, dhc0_dgc0
    real(8) :: ec_iso, dec_iso_drho, dec_iso_dg2
    real(8) :: wc1, dwc1_drs, beta, dbeta_drs, dif_lda, ddif_lda, d2dif_lda
    real(8) :: kgrad, dkgrad_drs, s2_gauss, ds2_gauss
    real(8) :: argnum, dargnum_drs, dargnum_drho, dargnum_dg2
    real(8) :: hc1_arg, darg_drho, darg_dg2, arg_quart, gc1, sc1, hc1, dhc1_dwc1, dhc1_darg
    real(8) :: ec_slow, dec_slow_drho, dec_slow_dg2
    real(8) :: fc_alpha, dfc_alpha, ec_diff, ec_r2scan, dec_drho, dec_dg2, dec_dtau

    eexc  = 0d0
    vxc   = 0d0
    vtau  = 0d0
    vgrad = 0d0
    if (rho <= rho_floor) return

    ! ---- ingredients ---------------------------------------------------------------
    grad2    = grad_norm * grad_norm
    rs       = c_rs / rho**(1d0/3d0)
    drs_drho = -rs / (3d0 * rho)
    tau_hom  = c_tau_hom * rho**(5d0/3d0)
    tau_w    = 0.125d0 * grad2 / rho
    s2       = c_s2 * grad2 / rho**(8d0/3d0)
    t2       = c_t2grad * grad2 / rho**(7d0/3d0)

    ds2_drho = -8d0 * s2 / (3d0 * rho)
    ds2_dg2  = c_s2 / rho**(8d0/3d0)
    dt2_drho = -7d0 * t2 / (3d0 * rho)
    dt2_dg2  = c_t2grad / rho**(7d0/3d0)

    alpha_den   = tau_hom + eta_reg * tau_w
    alpha       = (tau - tau_w) / alpha_den
    dalpha_dtau = 1d0 / alpha_den
    dalpha_drho = (tau_w / rho - alpha * (5d0 * tau_hom / (3d0 * rho) - eta_reg * tau_w / rho)) / alpha_den
    dalpha_dg2  = -(1d0 + eta_reg * alpha) * 0.125d0 / (rho * alpha_den)

    ! ---- exchange ------------------------------------------------------------------
    if (s2 <= s2_floor) then
      ! g_x -> 1 with every derivative vanishing faster than any power of s2
      gx_damp = 1d0
      dgx_ds2 = 0d0
    else
      gx_damp = 1d0 - exp(-a_damp / s2**0.25d0)
      dgx_ds2 = -0.25d0 * a_damp * (1d0 - gx_damp) / s2**1.25d0
    end if

    s2_damp  = exp(-s2 * s2 / d_damp4)
    xs       = (c_xgrad * s2_damp + mu_gea) * s2
    dxs_ds2  = c_xgrad * s2_damp * (1d0 - 2d0 * s2 * s2 / d_damp4) + mu_gea
    hx0      = 1d0 + kappa0
    hx1      = 1d0 + kappa1 * xs / (kappa1 + xs)
    dhx1_dxs = (kappa1 / (kappa1 + xs))**2

    call interpolate_in_alpha(alpha, cx_poly, c1x, c2x, dx_switch, fx_alpha, dfx_alpha)

    fx_grad     = hx1 + fx_alpha * (hx0 - hx1)
    fx_enh      = fx_grad * gx_damp
    dq_ds2      = (1d0 - fx_alpha) * dhx1_dxs * dxs_ds2
    dq_dalpha   = dfx_alpha * (hx0 - hx1)
    dfx_enh_ds2 = dq_ds2 * gx_damp + fx_grad * dgx_ds2

    ex_dens  = c_ex * rho**(4d0/3d0)
    dex_dens = (4d0/3d0) * c_ex * rho**(1d0/3d0)
    dex_drho = dex_dens * fx_enh + ex_dens * (dfx_enh_ds2 * ds2_drho + gx_damp * dq_dalpha * dalpha_drho)
    dex_dg2  = ex_dens * (dfx_enh_ds2 * ds2_dg2 + gx_damp * dq_dalpha * dalpha_dg2)
    dex_dtau = ex_dens * gx_damp * dq_dalpha * dalpha_dtau

    ! ---- correlation ---------------------------------------------------------------
    call pw92_correlation(rs, ec_lda, dec_lda, d2ec_lda)
    call iso_orbital_lda(rs, ec_lda0, dec_lda0, d2ec_lda0)

    ! alpha = 0 branch: the single-orbital LDA plus its gradient term
    wc0       = exp(-ec_lda0 / b1c) - 1d0
    dwc0_drs  = -(wc0 + 1d0) * dec_lda0 / b1c
    gc0       = (1d0 + 4d0 * chi_inf * s2)**(-0.25d0)
    dgc0_ds2  = -chi_inf * gc0 / (1d0 + 4d0 * chi_inf * s2)
    sc0       = 1d0 + wc0 * (1d0 - gc0)
    hc0       = b1c * log(sc0)
    dhc0_dwc0 = b1c * (1d0 - gc0) / sc0
    dhc0_dgc0 = -b1c * wc0 / sc0

    ec_iso       = ec_lda0 + hc0
    dec_iso_drho = (dec_lda0 + dhc0_dwc0 * dwc0_drs) * drs_drho + dhc0_dgc0 * dgc0_ds2 * ds2_drho
    dec_iso_dg2  = dhc0_dgc0 * dgc0_ds2 * ds2_dg2

    ! alpha = 1 branch: PW92 + PBE-like term. w_1 is strictly positive in practice
    ! (>=2e-5 at rho_floor); the w_floor clamp below is a NaN guard, applied to
    ! both value and derivative so a clamped point stays self-consistent.
    wc1 = exp(-ec_lda / gamma_c) - 1d0
    if (wc1 <= w_floor) then
      wc1      = w_floor
      dwc1_drs = 0d0
    else
      dwc1_drs = -(wc1 + 1d0) * dec_lda / gamma_c
    end if

    beta       = beta0 * (1d0 + beta_num * rs) / (1d0 + beta_den * rs)
    dbeta_drs  = -(beta_den - beta_num) * beta0 / (1d0 + beta_den * rs)**2
    dif_lda    = ec_lda0 - ec_lda
    ddif_lda   = dec_lda0 - dec_lda
    d2dif_lda  = d2ec_lda0 - d2ec_lda
    kgrad      = (fc_one / 27d0) * (20d0 * rs * ddif_lda - 45d0 * eta_reg * dif_lda)
    dkgrad_drs = (fc_one / 27d0) * (20d0 * ddif_lda + 20d0 * rs * d2dif_lda - 45d0 * eta_reg * ddif_lda)

    s2_gauss  = s2 * s2_damp
    ds2_gauss = s2_damp * (1d0 - 2d0 * s2 * s2 / d_damp4)

    argnum       = beta * t2 - kgrad * s2_gauss
    dargnum_drs  = dbeta_drs * t2 - dkgrad_drs * s2_gauss
    dargnum_drho = dargnum_drs * drs_drho + beta * dt2_drho - kgrad * ds2_gauss * ds2_drho
    dargnum_dg2  = beta * dt2_dg2 - kgrad * ds2_gauss * ds2_dg2

    ! A can be negative in the dilute tail (down to -0.0051); do not clamp it to
    ! zero -- that is genuine functional behaviour, not a singularity to guard.
    ! Bound holds: 1+4A >= 0.98, S >= 0.9999 always.
    hc1_arg   = argnum / (gamma_c * wc1)
    darg_drho = (dargnum_drho - hc1_arg * gamma_c * dwc1_drs * drs_drho) / (gamma_c * wc1)
    darg_dg2  = dargnum_dg2 / (gamma_c * wc1)

    arg_quart = 1d0 + 4d0 * hc1_arg
    gc1       = arg_quart**(-0.25d0)
    sc1       = 1d0 + wc1 * (1d0 - gc1)
    hc1       = gamma_c * log(sc1)
    dhc1_dwc1 = gamma_c * (1d0 - gc1) / sc1
    dhc1_darg = gamma_c * wc1 * gc1 / (sc1 * arg_quart)

    ec_slow       = ec_lda + hc1
    dec_slow_drho = (dec_lda + dhc1_dwc1 * dwc1_drs) * drs_drho + dhc1_darg * darg_drho
    dec_slow_dg2  = dhc1_darg * darg_dg2

    call interpolate_in_alpha(alpha, cc_poly, c1c, c2c, dc_switch, fc_alpha, dfc_alpha)

    ec_diff   = ec_iso - ec_slow
    ec_r2scan = ec_slow + fc_alpha * ec_diff
    dec_drho  = dec_slow_drho + fc_alpha * (dec_iso_drho - dec_slow_drho) + dfc_alpha * ec_diff * dalpha_drho
    dec_dg2   = dec_slow_dg2  + fc_alpha * (dec_iso_dg2  - dec_slow_dg2 ) + dfc_alpha * ec_diff * dalpha_dg2
    dec_dtau  = dfc_alpha * ec_diff * dalpha_dtau

    ! ---- the four outputs ----------------------------------------------------------
    eexc  = ex_dens * fx_enh + rho * ec_r2scan
    vxc   = dex_drho + ec_r2scan + rho * dec_drho
    vtau  = dex_dtau + rho * dec_dtau
    vgrad = 2d0 * grad_norm * (dex_dg2 + rho * dec_dg2)

    return
  end subroutine r2scan_point


  ! Switching function in alpha and its derivative; exchange and correlation
  ! differ only in the polynomial and the three exponential-tail constants, so
  ! one routine serves both (polynomial evaluated by Horner's method).
  subroutine interpolate_in_alpha(alpha, cpoly, c_low, c_high, d_high, f, dfda)
    implicit none
    real(8),intent(in)  :: alpha, cpoly(0:7), c_low, c_high, d_high
    real(8),intent(out) :: f, dfda
    integer :: j

    if (alpha < 0d0) then
      f    = exp(-c_low * alpha / (1d0 - alpha))
      dfda = -c_low * f / (1d0 - alpha)**2
    else if (alpha <= alpha_poly_max) then
      f    = cpoly(7)
      dfda = 7d0 * cpoly(7)
      do j = 6, 1, -1
        dfda = dfda * alpha + dble(j) * cpoly(j)
        f    = f * alpha + cpoly(j)
      end do
      f = f * alpha + cpoly(0)
    else
      f    = -d_high * exp(c_high / (1d0 - alpha))
      dfda = f * c_high / (1d0 - alpha)**2
    end if

    return
  end subroutine interpolate_in_alpha


  ! Perdew-Wang 1992 LDA correlation at zeta=0, with the first two derivatives
  ! w.r.t. r_s -- the second is needed because K carries dD/dr_s and is itself
  ! differentiated.
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


  ! The single-orbital (alpha = 0) correlation LDA and its first two r_s
  ! derivatives, in the same form as pw92_correlation so the two can be subtracted.
  subroutine iso_orbital_lda(rs, ec, dec_drs, d2ec_drs2)
    implicit none
    real(8),intent(in)  :: rs
    real(8),intent(out) :: ec, dec_drs, d2ec_drs2
    real(8) :: sqrt_rs, rpoly, drpoly, d2rpoly

    sqrt_rs = sqrt(rs)
    rpoly   = 1d0 + b2c * sqrt_rs + b3c * rs
    drpoly  = 0.5d0 * b2c / sqrt_rs + b3c
    d2rpoly = -0.25d0 * b2c / (rs * sqrt_rs)

    ec        = -b1c / rpoly
    dec_drs   = b1c * drpoly / rpoly**2
    d2ec_drs2 = b1c * d2rpoly / rpoly**2 - 2d0 * b1c * drpoly**2 / rpoly**3

    return
  end subroutine iso_orbital_lda

end module builtin_r2scan
