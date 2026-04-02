module builtin_r2scan

  implicit none

  real(8), parameter :: pi = 3.1415926535897932384626433832795d0
  real(8), parameter :: one = 1d0
  real(8), parameter :: two = 2d0
  real(8), parameter :: three = 3d0
  real(8), parameter :: four = 4d0
  real(8), parameter :: five = 5d0
  real(8), parameter :: six = 6d0
  real(8), parameter :: seven = 7d0
  real(8), parameter :: third = one / three
  real(8), parameter :: two_third = two / three
  real(8), parameter :: four_third = four / three
  real(8), parameter :: five_third = five / three
  real(8), parameter :: eight_third = 8d0 / 3d0

  real(8), parameter :: cx0 = 1d0
  real(8), parameter :: cx1 = -0.667d0
  real(8), parameter :: cx2 = -0.4445555d0
  real(8), parameter :: cx3 = -0.663086601049d0
  real(8), parameter :: cx4 = 1.451297044490d0
  real(8), parameter :: cx5 = -0.887998041597d0
  real(8), parameter :: cx6 = 0.234528941479d0
  real(8), parameter :: cx7 = -0.023185843322d0

  real(8), parameter :: cc0 = 1d0
  real(8), parameter :: cc1 = -0.64d0
  real(8), parameter :: cc2 = -0.4352d0
  real(8), parameter :: cc3 = -1.535685604549d0
  real(8), parameter :: cc4 = 3.061560252175d0
  real(8), parameter :: cc5 = -1.915710236206d0
  real(8), parameter :: cc6 = 0.516884468372d0
  real(8), parameter :: cc7 = -0.051848879792d0

  real(8), parameter :: scan_a1 = 4.9479d0
  real(8), parameter :: scan_c1x = 0.667d0
  real(8), parameter :: scan_c2x = 0.8d0
  real(8), parameter :: scan_dx = 1.24d0
  real(8), parameter :: scan_c1c = 0.64d0
  real(8), parameter :: scan_c2c = 1.5d0
  real(8), parameter :: scan_dc = 0.7d0

  real(8), parameter :: k0 = 0.174d0
  real(8), parameter :: k1 = 0.065d0
  real(8), parameter :: mu = 10d0 / 81d0
  real(8), parameter :: eta = 1d-3
  real(8), parameter :: dp2 = 0.361d0
  real(8), parameter :: dp24 = dp2**4

  real(8), parameter :: b1c = 0.0285764d0
  real(8), parameter :: b2c = 0.0889d0
  real(8), parameter :: b3c = 0.125541d0
  real(8), parameter :: beta_mb = 0.066725d0
  real(8), parameter :: chi_infinity = 0.12802585262625815d0
  real(8), parameter :: sgam = 0.031090690869655d0

  real(8), parameter :: lda_a = 0.03109070d0
  real(8), parameter :: lda_a1 = 0.21370d0
  real(8), parameter :: lda_b1 = 7.59570d0
  real(8), parameter :: lda_b2 = 3.58760d0
  real(8), parameter :: lda_b3 = 1.63820d0
  real(8), parameter :: lda_b4 = 0.49294d0

  real(8), parameter :: c_p = 0.026121172985233599568d0
  real(8), parameter :: c_t2 = 0.063468206097703704205d0
  real(8), parameter :: c_tau_unif = 2.8712340001881918160d0
  real(8), parameter :: ax_lda = -0.73855876638202240587d0
  real(8), parameter :: vx_lda_pref = -0.98474502184269654116d0

  real(8), parameter :: c_exchange_interp = (20d0 / 27d0 + 5d0 * eta / 3d0) * &
       (cx1 + two * cx2 + three * cx3 + four * cx4 + five * cx5 + six * cx6 + seven * cx7) * k0
  real(8), parameter :: c_correlation_interp = cc1 + two * cc2 + three * cc3 + four * cc4 + &
       five * cc5 + six * cc6 + seven * cc7

  real(8), parameter :: rho_floor = 1d-18
  real(8), parameter :: grad_floor = 1d-18
  real(8), parameter :: tiny_positive = 1d-30

contains

  subroutine exc_cor_r2scan(nl, rho, rho_s, grho_s, tau_s, eexc, vxc, vtau, vgrad)
    implicit none
    integer, intent(in) :: nl
    real(8), intent(in) :: rho(nl)
    real(8), intent(in) :: rho_s(nl)
    real(8), intent(in) :: grho_s(nl,3)
    real(8), intent(in) :: tau_s(nl)
    real(8), intent(out) :: eexc(nl)
    real(8), intent(out) :: vxc(nl)
    real(8), intent(out) :: vtau(nl)
    real(8), intent(out) :: vgrad(nl)
    integer :: i

    do i = 1, nl
      call eval_r2scan_point(rho(i), rho_s(i), grho_s(i,:), tau_s(i), &
           eexc(i), vxc(i), vtau(i), vgrad(i))
    end do
  end subroutine exc_cor_r2scan

  subroutine eval_r2scan_point(rho, rho_s, grho_s, tau_s, eexc, vxc, vtau, vgrad)
    implicit none
    real(8), intent(in) :: rho
    real(8), intent(in) :: rho_s
    real(8), intent(in) :: grho_s(3)
    real(8), intent(in) :: tau_s
    real(8), intent(out) :: eexc
    real(8), intent(out) :: vxc
    real(8), intent(out) :: vtau
    real(8), intent(out) :: vgrad

    real(8) :: rr, grad, sigma, tau, rs, srs, u, du_dn
    real(8) :: p, dp_dn, dp_dsigma_weighted
    real(8) :: t2, dt2_dn, dt2_dsigma_weighted
    real(8) :: w, dw_dn, dw_dsigma_weighted
    real(8) :: alpha, alpha_den_sq, dalpha_dn, dalpha_dsigma_weighted
    real(8) :: gx, gxp, dgx_dn, dgx_dsigma_weighted
    real(8) :: x, xp, h0x, h1x, h1xx, dh1x_dn, dh1x_dsigma_weighted
    real(8) :: fx, dfx_da, ffx
    real(8) :: ex, ex_density, vx_density, exarg, vxtau, vxn, vxs
    real(8) :: elda, dlda_dn, elda0, dlda0_dn, ddlda, beta, dbeta_dn
    real(8) :: fc, ddfc
    real(8) :: w1, dw1_dn, w0, dw0_dn, ginfinity, dgin_dn, dgin_dsigma_weighted
    real(8) :: h0c, dh0c_dn, dh0c_dsigma_weighted, ec0, dec0_dn, dec0_dsigma_weighted
    real(8) :: difddlda, dif2ddlda, y, dy_dn, dy_dsigma_weighted
    real(8) :: ddy, pp, dpp_dn, dpp_dsigma_weighted
    real(8) :: dddy_dn, dddy_dsigma_weighted
    real(8) :: term1, term4, dh1c_wn, dh1c_dn, dh1c_dsigma_weighted
    real(8) :: ec1, dec1_dn, dec1_dsigma_weighted
    real(8) :: corarg, vctau, vcn, vcs
    real(8) :: stf1, stf2, stf3, stf4, stf5, stf6, stf7
    real(8) :: bigdenom, exp_p, h1_arg, h0_arg, h1_den
    real(8) :: rho_total_from_spin

    eexc = 0d0
    vxc = 0d0
    vtau = 0d0
    vgrad = 0d0

    rho_total_from_spin = two * max(rho_s, 0d0)
    rr = max(rho, rho_total_from_spin)
    if (rr <= rho_floor) return

    grad = two * sqrt(max(0d0, dot_product(grho_s, grho_s)))
    sigma = grad * grad
    tau = max(two * tau_s, 0d0)

    rr = max(rr, rho_floor)
    rs = 0.62035049089940001665d0 / rr**third
    srs = sqrt(rs)

    u = c_tau_unif * rr**five_third
    du_dn = five_third * c_tau_unif * rr**two_third

    p = c_p * sigma / rr**eight_third
    dp_dn = -eight_third * p / rr
    dp_dsigma_weighted = c_p * grad / rr**eight_third

    t2 = c_t2 * sigma / rr**(seven / three)
    dt2_dn = -(seven / three) * t2 / rr
    dt2_dsigma_weighted = c_t2 * grad / rr**(seven / three)

    w = 0.125d0 * sigma / rr
    dw_dn = -w / rr
    dw_dsigma_weighted = 0.125d0 * grad / rr

    alpha = (tau - w) / max(u + eta * w, tiny_positive)
    alpha_den_sq = max((u + eta * w)**2, tiny_positive)
    dalpha_dn = -((u + eta * tau) * dw_dn + (tau - w) * du_dn) / alpha_den_sq
    dalpha_dsigma_weighted = -(u + eta * tau) * dw_dsigma_weighted / alpha_den_sq

    if (p <= 1d-16) then
      gx = 1d0
      gxp = 0d0
    else
      gx = 1d0 - exp(-scan_a1 / p**0.25d0)
      gxp = -0.25d0 * scan_a1 * exp(-scan_a1 / p**0.25d0) / p**1.25d0
    end if
    dgx_dn = gxp * dp_dn
    dgx_dsigma_weighted = gxp * dp_dsigma_weighted

    exp_p = exp(-(p * p) / dp24)
    x = (c_exchange_interp * exp_p + mu) * p
    xp = c_exchange_interp * exp_p * (1d0 - two * p * p / dp24) + mu

    h0x = 1d0 + k0
    h1x = 1d0 + k1 - k1 / (1d0 + x / k1)
    h1xx = 1d0 / (1d0 + x / k1)**2
    dh1x_dn = h1xx * xp * dp_dn
    dh1x_dsigma_weighted = h1xx * xp * dp_dsigma_weighted

    call eval_interp_exchange(alpha, fx, dfx_da)
    ffx = (h1x + fx * (h0x - h1x)) * gx

    ex = ax_lda * rr**four_third
    ex_density = ex
    vx_density = vx_lda_pref * rr**third
    exarg = ex_density * ffx
    vxtau = ex_density * dfx_da * (h0x - h1x) * gx / max(u + eta * w, tiny_positive)
    vxn = vx_density * ffx + ex_density * ( &
         dfx_da * dalpha_dn * gx * (h0x - h1x) + &
         dgx_dn * (h1x + fx * (h0x - h1x)) + &
         gx * dh1x_dn * (1d0 - fx))
    vxs = ex_density * ( &
         dfx_da * dalpha_dsigma_weighted * gx * (h0x - h1x) + &
         dgx_dsigma_weighted * (h1x + fx * (h0x - h1x)) + &
         gx * dh1x_dsigma_weighted * (1d0 - fx))

    bigdenom = srs * (lda_b3 * rs + lda_b1) + rs * (lda_b4 * rs + lda_b2)
    elda = -two * lda_a * (1d0 + lda_a1 * rs) * log(1d0 + 0.5d0 / (lda_a * bigdenom))
    dlda_dn = ( &
         -two * lda_a * lda_a1 * log(1d0 + 1d0 / (two * lda_a * bigdenom)) + &
         (lda_a1 * rs + 1d0) * &
         ((lda_b3 * rs + lda_b1) / (two * srs) + srs * lda_b3 + two * lda_b4 * rs + lda_b2) / &
         (bigdenom**2 + bigdenom / (two * lda_a)) ) * (-rs / (three * rr))

    elda0 = -b1c / (1d0 + b2c * srs + b3c * rs)
    dlda0_dn = b1c * (0.5d0 * b2c / srs + b3c) / (1d0 + b2c * srs + b3c * rs)**2 * (-rs / (three * rr))
    ddlda = elda0 - elda

    beta = beta_mb * (1d0 + 0.1d0 * rs) / (1d0 + 0.1778d0 * rs)
    dbeta_dn = -(0.0778d0 * beta_mb / (1d0 + 0.1778d0 * rs)**2) * (-rs / (three * rr))

    call eval_interp_correlation(alpha, fc, ddfc)

    w1 = exp(-elda / sgam) - 1d0
    dw1_dn = -(exp(-elda / sgam) / sgam) * dlda_dn
    w0 = exp(-elda0 / b1c) - 1d0
    dw0_dn = -(exp(-elda0 / b1c) / b1c) * dlda0_dn

    ginfinity = 1d0 / (1d0 + four * chi_infinity * p)**0.25d0
    dgin_dn = -(chi_infinity * ginfinity / (1d0 + four * chi_infinity * p)) * dp_dn
    dgin_dsigma_weighted = -(chi_infinity * ginfinity / (1d0 + four * chi_infinity * p)) * dp_dsigma_weighted

    h0_arg = max(1d0 + w0 * (1d0 - ginfinity), tiny_positive)
    h0c = b1c * log(h0_arg)
    dh0c_dn = (b1c / h0_arg) * ((1d0 - ginfinity) * dw0_dn - w0 * dgin_dn)
    dh0c_dsigma_weighted = -b1c * w0 * dgin_dsigma_weighted / h0_arg
    ec0 = elda0 + h0c
    dec0_dn = dlda0_dn + dh0c_dn
    dec0_dsigma_weighted = dh0c_dsigma_weighted

    difddlda = b1c * (b2c / (two * srs) + b3c) / (1d0 + b2c * srs + b3c * rs)**2 + &
         two * lda_a * lda_a1 * log(1d0 + 1d0 / (two * lda_a * bigdenom)) - &
         (lda_a1 * rs + 1d0) * ((lda_b3 * rs + lda_b1) / (two * srs) + srs * lda_b3 + &
         two * lda_b4 * rs + lda_b2) / (bigdenom**2 * (1d0 + 1d0 / (two * lda_a * bigdenom)))

    y = beta * t2 / max(sgam * w1, tiny_positive)
    dy_dsigma_weighted = beta * dt2_dsigma_weighted / max(sgam * w1, tiny_positive)
    dy_dn = (dbeta_dn * t2 + beta * dt2_dn - beta * t2 * dw1_dn / max(w1, tiny_positive)) / max(sgam * w1, tiny_positive)

    ddy = (c_correlation_interp / max(27d0 * sgam * w1, tiny_positive)) * &
         (20d0 * rs * difddlda - 45d0 * eta * ddlda) * p * exp_p

    pp = p * exp_p
    dpp_dsigma_weighted = exp_p * (dp24 - two * p * p) / dp24 * dp_dsigma_weighted
    dpp_dn = exp_p * (dp24 - two * p * p) / dp24 * dp_dn

    dddy_dsigma_weighted = (c_correlation_interp / max(27d0 * sgam * w1, tiny_positive)) * &
         (20d0 * rs * difddlda - 45d0 * eta * ddlda) * dpp_dsigma_weighted

    stf1 = -b1c * (8d0 * srs * rs * b3c**2 + 9d0 * b2c * b3c * rs + three * srs * b2c**2 + b2c) / &
         (four * rs * srs * (1d0 + b2c * srs + b3c * rs)**3)
    stf2 = bigdenom
    stf3 = 1d0 + 1d0 / (two * lda_a * stf2)
    stf4 = (two * lda_a1 * ((lda_b3 * rs + lda_b1) / (two * srs) + srs * lda_b3 + two * lda_b4 * rs + lda_b2)) / &
         (stf2**2 * stf3)
    stf5 = (two * (lda_a1 * rs + 1d0) * ((lda_b3 * rs + lda_b1) / (two * srs) + srs * lda_b3 + &
         two * lda_b4 * rs + lda_b2)**2) / (stf2**3 * stf3)
    stf6 = (lda_a1 * rs + 1d0) * (-(lda_b3 * rs + lda_b1) / (four * srs * rs) + lda_b3 / srs + two * lda_b4) / &
         (stf2**2 * stf3)
    stf7 = (lda_a1 * rs + 1d0) * ((lda_b3 * rs + lda_b1) / (two * srs) + srs * lda_b3 + &
         two * lda_b4 * rs + lda_b2)**2 / (two * lda_a * stf2**4 * stf3**2)
    dif2ddlda = stf1 - stf4 + stf5 - stf6 - stf7

    dddy_dn = (c_correlation_interp / max(27d0 * sgam * w1, tiny_positive)) * &
         (20d0 * rs * difddlda - 45d0 * eta * ddlda) * dpp_dn - &
         ddy * dw1_dn / max(w1, tiny_positive) + &
         (c_correlation_interp / max(27d0 * sgam * w1, tiny_positive)) * &
         (20d0 * difddlda + (20d0 * rs * dif2ddlda - 45d0 * eta * difddlda)) * pp * (-rs / (three * rr))

    term1 = max(1d0 + four * (y - ddy), tiny_positive)
    term4 = term1**0.25d0
    h1_arg = max(1d0 + w1 * (1d0 - 1d0 / term4), tiny_positive)
    h1_den = max(term4 * (1d0 + w1) - w1, tiny_positive)
    dh1c_dsigma_weighted = (sgam * w1 / (term1 * h1_den)) * (dy_dsigma_weighted - dddy_dsigma_weighted)
    dh1c_wn = sgam * ((term4 - 1d0) / h1_den) * dw1_dn
    dh1c_dn = dh1c_wn + (sgam * w1 / (term1 * h1_den)) * (dy_dn - dddy_dn)

    ec1 = elda + sgam * log(h1_arg)
    dec1_dn = dlda_dn + dh1c_dn
    dec1_dsigma_weighted = dh1c_dsigma_weighted

    corarg = rr * (ec1 + fc * (ec0 - ec1))
    vctau = rr * ddfc * (ec0 - ec1) / max(u + eta * w, tiny_positive)
    vcn = ec1 + fc * (ec0 - ec1) + ddfc * dalpha_dn * rr * (ec0 - ec1) + &
         dec0_dn * rr * fc + dec1_dn * rr * (1d0 - fc)
    vcs = ddfc * dalpha_dsigma_weighted * rr * (ec0 - ec1) + &
         dec0_dsigma_weighted * rr * fc + dec1_dsigma_weighted * rr * (1d0 - fc)

    eexc = exarg + corarg
    vxc = vxn + vcn
    vtau = vxtau + vctau
    vgrad = two * (vxs + vcs)
  end subroutine eval_r2scan_point

  subroutine eval_interp_exchange(alpha, fx, dfx_da)
    implicit none
    real(8), intent(in) :: alpha
    real(8), intent(out) :: fx
    real(8), intent(out) :: dfx_da

    if (alpha < 0d0) then
      fx = exp(-scan_c1x * alpha / (1d0 - alpha))
      dfx_da = -scan_c1x * fx / (1d0 - alpha)**2
    else if (alpha < 2.5d0) then
      fx = cx0 + alpha * (cx1 + alpha * (cx2 + alpha * (cx3 + alpha * (cx4 + alpha * (cx5 + alpha * (cx6 + alpha * cx7))))))
      dfx_da = cx1 + alpha * (two * cx2 + alpha * (three * cx3 + alpha * (four * cx4 + alpha * (five * cx5 + alpha * (six * cx6 + alpha * seven * cx7)))))
    else
      fx = -scan_dx * exp(scan_c2x / (1d0 - alpha))
      dfx_da = -scan_dx * scan_c2x * exp(scan_c2x / (1d0 - alpha)) / (1d0 - alpha)**2
    end if
  end subroutine eval_interp_exchange

  subroutine eval_interp_correlation(alpha, fc, ddfc)
    implicit none
    real(8), intent(in) :: alpha
    real(8), intent(out) :: fc
    real(8), intent(out) :: ddfc

    if (alpha < 0d0) then
      fc = exp(-scan_c1c * alpha / (1d0 - alpha))
      ddfc = -scan_c1c * fc / (1d0 - alpha)**2
    else if (alpha < 2.5d0) then
      fc = cc0 + alpha * (cc1 + alpha * (cc2 + alpha * (cc3 + alpha * (cc4 + alpha * (cc5 + alpha * (cc6 + alpha * cc7))))))
      ddfc = cc1 + alpha * (two * cc2 + alpha * (three * cc3 + alpha * (four * cc4 + alpha * (five * cc5 + alpha * (six * cc6 + alpha * seven * cc7)))))
    else
      fc = -scan_dc * exp(scan_c2c / (1d0 - alpha))
      ddfc = -scan_dc * scan_c2c * exp(scan_c2c / (1d0 - alpha)) / (1d0 - alpha)**2
    end if
  end subroutine eval_interp_correlation
end module builtin_r2scan
