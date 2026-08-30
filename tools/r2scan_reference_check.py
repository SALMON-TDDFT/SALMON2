#!/usr/bin/env python3
#
#  Copyright 2026 SALMON developers
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
#
"""Reference check for the built-in r2SCAN functional (src/xc/builtin_r2scan.f90).

This is a standalone transcription of the point evaluator r2scan_point, written in the
same variable names and the same order as the Fortran, together with the checks that fix
the functional: the uniform-gas limit, the second-order gradient expansion that fixes the
two constants C_x and K, the continuity of the alpha interpolation at its branch points,
and the three analytic derivatives against high-order finite differences of the value.

It exists because the Fortran evaluator is short, dense and entirely made of derivatives:
a sign error in one of them changes the self-consistent field rather than crashing it, and
a plane-wave regression test would not localize it.  Running this script after editing
builtin_r2scan.f90 checks the algebra directly, in a few seconds and without a build.

    python3 tools/r2scan_reference_check.py

Exit status 0 if every check passes.  It has no dependencies beyond the standard library.

What it does and does not establish.  Because it follows the Fortran statement for
statement, it is not an independent second implementation, and a mistake made in both at
the same place would survive it.  What it does test independently of that transcription is
everything the checks compare the transcription against: the closed-form limits, the two
gradient-expansion constraints that fix C_x and K, the branch joins, and the derivatives,
which are compared with finite differences of the value and so do not depend on how the
derivative expressions were written.

References
  [1] J. W. Furness, A. D. Kaplan, J. Ning, J. P. Perdew and J. Sun,
      J. Phys. Chem. Lett. 11, 8208 (2020).
  [2] J. Sun, A. Ruzsinszky and J. P. Perdew, Phys. Rev. Lett. 115, 036402 (2015).
  [3] J. P. Perdew and Y. Wang, Phys. Rev. B 45, 13244 (1992).
"""

import math
import random
import sys

# ---- the parameters of builtin_r2scan.f90, in the same order ---------------------------
cx_poly = (1.0, -0.667, -0.4445555, -0.663086601049,
           1.451297044490, -0.887998041597, 0.234528941479, -0.023185843322)
cc_poly = (1.0, -0.64, -0.4352, -1.535685604549,
           3.061560252175, -1.915710236206, 0.516884468372, -0.051848879792)
c1x, c2x, dx_switch = 0.667, 0.8, 1.24
c1c, c2c, dc_switch = 0.64, 1.5, 0.7
alpha_poly_max = 2.5

kappa0, kappa1 = 0.174, 0.065
mu_gea = 10.0 / 81.0
a_damp, d_damp = 4.9479, 0.361
d_damp4 = d_damp ** 4
eta_reg = 1.0e-3

b1c, b2c, b3c = 0.0285764, 0.0889, 0.125541
chi_inf = 0.12802585262625815
beta0, beta_num, beta_den = 0.066725, 0.1, 0.1778
gamma_c = 0.031090690869655
pw_a, pw_a1 = 0.0310907, 0.21370
pw_b1, pw_b2, pw_b3, pw_b4 = 7.59570, 3.58760, 1.63820, 0.49294

c_rs = 0.62035049089940001665
c_tau_hom = 2.8712340001881918160
c_s2 = 0.026121172985233599568
c_t2grad = 0.063468206097703704205
c_ex = -0.73855876638202240587

fx_one = sum((i + 0.0) * cx_poly[i] for i in range(1, 8))
fc_one = sum((i + 0.0) * cc_poly[i] for i in range(1, 8))
c_eta = 20.0 / 27.0 + 5.0 * eta_reg / 3.0
c_xgrad = c_eta * kappa0 * fx_one

rho_floor = 1.0e-18
s2_floor = 1.0e-16
w_floor = 1.0e-30


def interpolate_in_alpha(alpha, cpoly, c_low, c_high, d_high):
    if alpha < 0.0:
        f = math.exp(-c_low * alpha / (1.0 - alpha))
        return f, -c_low * f / (1.0 - alpha) ** 2
    if alpha <= alpha_poly_max:
        f = cpoly[7]
        dfda = 7.0 * cpoly[7]
        for j in range(6, 0, -1):
            dfda = dfda * alpha + j * cpoly[j]
            f = f * alpha + cpoly[j]
        return f * alpha + cpoly[0], dfda
    f = -d_high * math.exp(c_high / (1.0 - alpha))
    return f, f * c_high / (1.0 - alpha) ** 2


def pw92_correlation(rs):
    sqrt_rs = math.sqrt(rs)
    qden = sqrt_rs * (pw_b1 + pw_b3 * rs) + rs * (pw_b2 + pw_b4 * rs)
    dqden = 0.5 * pw_b1 / sqrt_rs + pw_b2 + 1.5 * pw_b3 * sqrt_rs + 2.0 * pw_b4 * rs
    d2qden = -0.25 * pw_b1 / (rs * sqrt_rs) + 0.75 * pw_b3 / sqrt_rs + 2.0 * pw_b4
    mden = qden * (2.0 * pw_a * qden + 1.0)
    logs = math.log(1.0 + 1.0 / (2.0 * pw_a * qden))
    dlogs = -dqden / mden
    d2logs = -d2qden / mden + dqden * dqden * (4.0 * pw_a * qden + 1.0) / (mden * mden)
    apoly = 1.0 + pw_a1 * rs
    return (-2.0 * pw_a * apoly * logs,
            -2.0 * pw_a * (pw_a1 * logs + apoly * dlogs),
            -2.0 * pw_a * (2.0 * pw_a1 * dlogs + apoly * d2logs))


def iso_orbital_lda(rs):
    sqrt_rs = math.sqrt(rs)
    rpoly = 1.0 + b2c * sqrt_rs + b3c * rs
    drpoly = 0.5 * b2c / sqrt_rs + b3c
    d2rpoly = -0.25 * b2c / (rs * sqrt_rs)
    return (-b1c / rpoly,
            b1c * drpoly / rpoly ** 2,
            b1c * d2rpoly / rpoly ** 2 - 2.0 * b1c * drpoly ** 2 / rpoly ** 3)


def r2scan_point(rho, grad_norm, tau, channel="both"):
    """f = n eps_xc and df/dn, df/dtau, df/d|grad n|; channel selects exchange or
    correlation alone (checks need both separately, the Fortran only needs their sum)."""
    if rho <= rho_floor:
        return 0.0, 0.0, 0.0, 0.0

    grad2 = grad_norm * grad_norm
    rs = c_rs / rho ** (1.0 / 3.0)
    drs_drho = -rs / (3.0 * rho)
    tau_hom = c_tau_hom * rho ** (5.0 / 3.0)
    tau_w = 0.125 * grad2 / rho
    s2 = c_s2 * grad2 / rho ** (8.0 / 3.0)
    t2 = c_t2grad * grad2 / rho ** (7.0 / 3.0)

    ds2_drho = -8.0 * s2 / (3.0 * rho)
    ds2_dg2 = c_s2 / rho ** (8.0 / 3.0)
    dt2_drho = -7.0 * t2 / (3.0 * rho)
    dt2_dg2 = c_t2grad / rho ** (7.0 / 3.0)

    alpha_den = tau_hom + eta_reg * tau_w
    alpha = (tau - tau_w) / alpha_den
    dalpha_dtau = 1.0 / alpha_den
    dalpha_drho = (tau_w / rho
                   - alpha * (5.0 * tau_hom / (3.0 * rho) - eta_reg * tau_w / rho)) / alpha_den
    dalpha_dg2 = -(1.0 + eta_reg * alpha) * 0.125 / (rho * alpha_den)

    # exchange
    if s2 <= s2_floor:
        gx_damp, dgx_ds2 = 1.0, 0.0
    else:
        gx_damp = 1.0 - math.exp(-a_damp / s2 ** 0.25)
        dgx_ds2 = -0.25 * a_damp * (1.0 - gx_damp) / s2 ** 1.25

    s2_damp = math.exp(-s2 * s2 / d_damp4) if s2 * s2 / d_damp4 < 700.0 else 0.0
    xs = (c_xgrad * s2_damp + mu_gea) * s2
    dxs_ds2 = c_xgrad * s2_damp * (1.0 - 2.0 * s2 * s2 / d_damp4) + mu_gea
    hx0 = 1.0 + kappa0
    hx1 = 1.0 + kappa1 * xs / (kappa1 + xs)
    dhx1_dxs = (kappa1 / (kappa1 + xs)) ** 2

    fx_alpha, dfx_alpha = interpolate_in_alpha(alpha, cx_poly, c1x, c2x, dx_switch)

    fx_grad = hx1 + fx_alpha * (hx0 - hx1)
    fx_enh = fx_grad * gx_damp
    dq_ds2 = (1.0 - fx_alpha) * dhx1_dxs * dxs_ds2
    dq_dalpha = dfx_alpha * (hx0 - hx1)
    dfx_enh_ds2 = dq_ds2 * gx_damp + fx_grad * dgx_ds2

    ex_dens = c_ex * rho ** (4.0 / 3.0)
    dex_dens = (4.0 / 3.0) * c_ex * rho ** (1.0 / 3.0)
    dex_drho = dex_dens * fx_enh + ex_dens * (dfx_enh_ds2 * ds2_drho
                                              + gx_damp * dq_dalpha * dalpha_drho)
    dex_dg2 = ex_dens * (dfx_enh_ds2 * ds2_dg2 + gx_damp * dq_dalpha * dalpha_dg2)
    dex_dtau = ex_dens * gx_damp * dq_dalpha * dalpha_dtau

    # correlation
    ec_lda, dec_lda, d2ec_lda = pw92_correlation(rs)
    ec_lda0, dec_lda0, d2ec_lda0 = iso_orbital_lda(rs)

    wc0 = math.exp(-ec_lda0 / b1c) - 1.0
    dwc0_drs = -(wc0 + 1.0) * dec_lda0 / b1c
    gc0 = (1.0 + 4.0 * chi_inf * s2) ** (-0.25)
    dgc0_ds2 = -chi_inf * gc0 / (1.0 + 4.0 * chi_inf * s2)
    sc0 = 1.0 + wc0 * (1.0 - gc0)
    hc0 = b1c * math.log(sc0)
    dhc0_dwc0 = b1c * (1.0 - gc0) / sc0
    dhc0_dgc0 = -b1c * wc0 / sc0

    ec_iso = ec_lda0 + hc0
    dec_iso_drho = (dec_lda0 + dhc0_dwc0 * dwc0_drs) * drs_drho + dhc0_dgc0 * dgc0_ds2 * ds2_drho
    dec_iso_dg2 = dhc0_dgc0 * dgc0_ds2 * ds2_dg2

    wc1 = math.exp(-ec_lda / gamma_c) - 1.0
    if wc1 <= w_floor:
        wc1, dwc1_drs = w_floor, 0.0
    else:
        dwc1_drs = -(wc1 + 1.0) * dec_lda / gamma_c

    beta = beta0 * (1.0 + beta_num * rs) / (1.0 + beta_den * rs)
    dbeta_drs = -(beta_den - beta_num) * beta0 / (1.0 + beta_den * rs) ** 2
    dif_lda = ec_lda0 - ec_lda
    ddif_lda = dec_lda0 - dec_lda
    d2dif_lda = d2ec_lda0 - d2ec_lda
    kgrad = (fc_one / 27.0) * (20.0 * rs * ddif_lda - 45.0 * eta_reg * dif_lda)
    dkgrad_drs = (fc_one / 27.0) * (20.0 * ddif_lda + 20.0 * rs * d2dif_lda
                                    - 45.0 * eta_reg * ddif_lda)

    s2_gauss = s2 * s2_damp
    ds2_gauss = s2_damp * (1.0 - 2.0 * s2 * s2 / d_damp4)

    argnum = beta * t2 - kgrad * s2_gauss
    dargnum_drs = dbeta_drs * t2 - dkgrad_drs * s2_gauss
    dargnum_drho = dargnum_drs * drs_drho + beta * dt2_drho - kgrad * ds2_gauss * ds2_drho
    dargnum_dg2 = beta * dt2_dg2 - kgrad * ds2_gauss * ds2_dg2

    hc1_arg = argnum / (gamma_c * wc1)
    darg_drho = (dargnum_drho - hc1_arg * gamma_c * dwc1_drs * drs_drho) / (gamma_c * wc1)
    darg_dg2 = dargnum_dg2 / (gamma_c * wc1)

    arg_quart = 1.0 + 4.0 * hc1_arg
    gc1 = arg_quart ** (-0.25)
    sc1 = 1.0 + wc1 * (1.0 - gc1)
    hc1 = gamma_c * math.log(sc1)
    dhc1_dwc1 = gamma_c * (1.0 - gc1) / sc1
    dhc1_darg = gamma_c * wc1 * gc1 / (sc1 * arg_quart)

    ec_slow = ec_lda + hc1
    dec_slow_drho = (dec_lda + dhc1_dwc1 * dwc1_drs) * drs_drho + dhc1_darg * darg_drho
    dec_slow_dg2 = dhc1_darg * darg_dg2

    fc_alpha, dfc_alpha = interpolate_in_alpha(alpha, cc_poly, c1c, c2c, dc_switch)

    ec_diff = ec_iso - ec_slow
    ec_r2scan = ec_slow + fc_alpha * ec_diff
    dec_drho = (dec_slow_drho + fc_alpha * (dec_iso_drho - dec_slow_drho)
                + dfc_alpha * ec_diff * dalpha_drho)
    dec_dg2 = (dec_slow_dg2 + fc_alpha * (dec_iso_dg2 - dec_slow_dg2)
               + dfc_alpha * ec_diff * dalpha_dg2)
    dec_dtau = dfc_alpha * ec_diff * dalpha_dtau

    if channel == "x":
        return ex_dens * fx_enh, dex_drho, dex_dtau, 2.0 * grad_norm * dex_dg2
    if channel == "c":
        return rho * ec_r2scan, ec_r2scan + rho * dec_drho, rho * dec_dtau, 2.0 * grad_norm * rho * dec_dg2
    return (ex_dens * fx_enh + rho * ec_r2scan,
            dex_drho + ec_r2scan + rho * dec_drho,
            dex_dtau + rho * dec_dtau,
            2.0 * grad_norm * (dex_dg2 + rho * dec_dg2))


# ---- checks ----------------------------------------------------------------------------

FAIL = []


def report(name, ok, detail):
    print(("  PASS  " if ok else "  FAIL  ") + name + "   " + detail)
    if not ok:
        FAIL.append(name)


def check_uniform_gas():
    """Uniform gas (grad n = 0, tau = tau_hom): F_x must be 1, eps_c must equal PW92,
    to the 1.7e-13 floor set by how well the published coefficients give f_x(1) = 0."""
    print("uniform electron gas")
    for rho in (1.0e-3, 1.0e-2, 0.1, 1.0, 5.0):
        tau = c_tau_hom * rho ** (5.0 / 3.0)
        ex = r2scan_point(rho, 0.0, tau, "x")[0]
        ec = r2scan_point(rho, 0.0, tau, "c")[0]
        fx = ex / (c_ex * rho ** (4.0 / 3.0))
        ec_pw = pw92_correlation(c_rs / rho ** (1.0 / 3.0))[0]
        dfx = abs(fx - 1.0)
        dec = abs(ec / rho - ec_pw) / abs(ec_pw)
        report("F_x = 1, eps_c = PW92 at n = %-6g" % rho, dfx < 1.0e-12 and dec < 1.0e-14,
               "|F_x - 1| = %.2e   |eps_c/eps_c^PW92 - 1| = %.2e" % (dfx, dec))


def check_gradient_expansion():
    """Checks the two constants (C_x, K) fixed by the second-order gradient expansion:
    dF_x/ds2 must equal mu=10/81 on the slowly-varying path, and eps_c's leading gradient
    term must be exactly beta(r_s)*t2 with no linear-in-s2 remainder."""
    print("second-order gradient expansion")
    rho = 1.0
    tau_hom = c_tau_hom * rho ** (5.0 / 3.0)
    grad = 1.0e-3
    s2 = c_s2 * grad * grad / rho ** (8.0 / 3.0)
    tau_w = grad * grad / (8.0 * rho)
    tau = (1.0 - c_eta * s2) * (tau_hom + eta_reg * tau_w) + tau_w
    fx = r2scan_point(rho, grad, tau, "x")[0] / r2scan_point(rho, 0.0, tau_hom, "x")[0]
    slope = (fx - 1.0) / s2
    report("dF_x/ds2 -> mu = 10/81", abs(slope - mu_gea) < 1.0e-7,
           "slope = %.12f   mu = %.12f" % (slope, mu_gea))

    for rho in (1.0e-2, 0.1, 1.0, 5.0):
        tau_hom = c_tau_hom * rho ** (5.0 / 3.0)
        rs = c_rs / rho ** (1.0 / 3.0)
        beta = beta0 * (1.0 + beta_num * rs) / (1.0 + beta_den * rs)
        ec_pw, dec_pw, _ = pw92_correlation(rs)
        ec_0, dec_0, _ = iso_orbital_lda(rs)
        c_corr = -5.0 * eta_reg / 3.0 + (20.0 / 27.0) * rs * (dec_0 - dec_pw) / (ec_0 - ec_pw)
        grad = rho * 1.0e-3
        s2 = c_s2 * grad * grad / rho ** (8.0 / 3.0)
        t2 = c_t2grad * grad * grad / rho ** (7.0 / 3.0)
        tau_w = grad * grad / (8.0 * rho)
        tau = (1.0 + c_corr * s2) * (tau_hom + eta_reg * tau_w) + tau_w
        ec = r2scan_point(rho, grad, tau, "c")[0] / rho
        ratio = (ec - ec_pw) / (beta * t2)
        report("eps_c - eps_c^PW92 -> beta t2 at n = %-6g" % rho, abs(ratio - 1.0) < 1.0e-5,
               "ratio = %.10f" % ratio)


def check_switching_functions():
    """f(1) = 0, and the polynomial must join the exponential tails at alpha = 0 and 2.5."""
    print("alpha interpolation")
    for name, cpoly, cl, ch, dh in (("f_x", cx_poly, c1x, c2x, dx_switch),
                                    ("f_c", cc_poly, c1c, c2c, dc_switch)):
        f1 = interpolate_in_alpha(1.0, cpoly, cl, ch, dh)[0]
        report("%s(1) = 0" % name, abs(f1) < 1.0e-11, "f(1) = %.3e" % f1)
        for a0 in (0.0, alpha_poly_max):
            lo = interpolate_in_alpha(a0 - 1.0e-12, cpoly, cl, ch, dh)
            hi = interpolate_in_alpha(a0 + 1.0e-12, cpoly, cl, ch, dh)
            jump = abs(lo[0] - hi[0])
            djump = abs(lo[1] - hi[1])
            report("%s continuous at alpha = %.1f" % (name, a0), jump < 1.0e-9 and djump < 1.0e-6,
                   "|df| = %.2e   |df'| = %.2e" % (jump, djump))


def check_logarithm_domain():
    """Why the correlation logarithms need no guard: over r_s in [1e-4,1e8] and s2 in
    [1e-8,1e6], A never drops below about -0.005, keeping 1+4A and S within 1e-2 of
    unity. A goes slightly negative in the dilute tail by design (not clamped)."""
    print("domain of the correlation logarithms")
    min_quart, min_s = 1.0e9, 1.0e9
    for i in range(600):
        rs = 10.0 ** (-4.0 + 12.0 * i / 599.0)
        ec_pw, dec_pw, _ = pw92_correlation(rs)
        ec_0, dec_0, _ = iso_orbital_lda(rs)
        wc1 = max(math.exp(-ec_pw / gamma_c) - 1.0, w_floor)
        beta = beta0 * (1.0 + beta_num * rs) / (1.0 + beta_den * rs)
        kgrad = (fc_one / 27.0) * (20.0 * rs * (dec_0 - dec_pw) - 45.0 * eta_reg * (ec_0 - ec_pw))
        for j in range(600):
            s2 = 10.0 ** (-8.0 + 14.0 * j / 599.0)
            t2 = s2 * (c_t2grad / c_s2) * c_rs / rs
            expo = s2 * s2 / d_damp4
            damp = math.exp(-expo) if expo < 700.0 else 0.0
            quart = 1.0 + 4.0 * (beta * t2 - kgrad * s2 * damp) / (gamma_c * wc1)
            min_quart = min(min_quart, quart)
            if quart > 0.0:
                min_s = min(min_s, 1.0 + wc1 * (1.0 - quart ** -0.25))
    report("1 + 4A stays positive", min_quart > 0.5, "min = %.6f" % min_quart)
    report("S_1 stays positive", min_s > 0.5, "min = %.6f" % min_s)


def check_derivatives():
    """Checks the three analytic derivatives against five-point finite differences,
    sweeping the step size (the functional is stiff at large s2) and scaling the
    comparison by max(|analytic|, |f|/x) so cancellation doesn't hide errors."""
    print("analytic derivatives vs five-point finite differences")
    steps = (1.0e-4, 1.0e-5, 1.0e-6)
    random.seed(20260830)
    worst = [0.0, 0.0, 0.0]
    for _ in range(200):
        rho = 10.0 ** random.uniform(-5.0, 0.7)
        grad = rho * 10.0 ** random.uniform(-2.0, 1.0)
        tau_w = grad * grad / (8.0 * rho)
        tau_hom = c_tau_hom * rho ** (5.0 / 3.0)
        alpha = 10.0 ** random.uniform(-3.0, 0.8)
        tau = tau_w + alpha * (tau_hom + eta_reg * tau_w)
        ana = r2scan_point(rho, grad, tau)

        def five(f, x, h):
            return (-f(x + 2 * h) + 8 * f(x + h) - 8 * f(x - h) + f(x - 2 * h)) / (12 * h)

        along = (lambda v: r2scan_point(v, grad, tau)[0],
                 lambda v: r2scan_point(rho, grad, v)[0],
                 lambda v: r2scan_point(rho, v, tau)[0])
        var = (rho, tau, grad)
        for j, k in enumerate((1, 2, 3)):
            scale = max(abs(ana[k]), abs(ana[0]) / var[j])
            err = min(abs(five(along[j], var[j], var[j] * h) - ana[k]) / scale for h in steps)
            worst[j] = max(worst[j], err)
    for j, name in enumerate(("df/dn", "df/dtau", "df/d|grad n|")):
        report("%-14s reproduced by finite differences" % name, worst[j] < 1.0e-8,
               "max scaled error over 200 points = %.2e" % worst[j])


def main():
    print("r2SCAN reference check for src/xc/builtin_r2scan.f90")
    print("  C_x = C_eta kappa0 f_x'(1) = %.15f" % c_xgrad)
    print("  f_x'(1) = %.12f   f_c'(1) = %.12f   C_eta = %.12f" % (fx_one, fc_one, c_eta))
    print()
    check_uniform_gas()
    check_gradient_expansion()
    check_switching_functions()
    check_logarithm_domain()
    check_derivatives()
    print()
    if FAIL:
        print("FAILED: " + ", ".join(FAIL))
        return 1
    print("all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
