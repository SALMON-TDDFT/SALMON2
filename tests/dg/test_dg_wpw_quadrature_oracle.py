#!/usr/bin/env python3
import numpy as np


def volume(u, gu, v, gv, potential, weight):
    s = weight * np.conj(u) * v
    h = weight * (0.5 * np.vdot(gu, gv) + potential * np.conj(u) * v)
    return s, h


def face(u_m, u_p, gu_m, gu_p, v_m, v_p, gv_m, gv_p, normal, h, weight):
    ju = u_m - u_p
    jv = v_m - v_p
    au = 0.5 * np.dot(gu_m + gu_p, normal)
    av = 0.5 * np.dot(gv_m + gv_p, normal)
    return weight * (-0.5 * (np.conj(au) * jv + np.conj(ju) * av)
                     + 10.0 / h * np.conj(ju) * jv)


w = 1.0 + 0.2j
gw = np.array([0.3 + 0.1j, -0.2j, 0.4])
p = 0.7 - 0.3j
gp = np.array([0.1j, 0.2, -0.3j])
s_wp, h_wp = volume(w, gw, p, gp, 1.4, 0.25)
s_pp, h_pp = volume(p, gp, p, gp, 1.4, 0.25)
assert np.isfinite([s_wp.real, h_wp.real, s_pp.real, h_pp.real]).all()

n = np.array([1.0, 0.0, 0.0])
wm, wp = 1.0 + 0.1j, 0.2 - 0.3j
gwm = np.array([0.4, 0.1j, 0.0])
gwp = np.array([-0.2, 0.1j, 0.0])
pm = pp = 0.8 + 0.2j
gpm = gpp = np.array([0.3j, 0.2, 0.0])
f = face(wm, wp, gwm, gwp, pm, pp, gpm, gpp, n, 0.5, 0.4)
fr = face(wp, wm, gwp, gwm, pp, pm, gpp, gpm, -n, 0.5, 0.4)
np.testing.assert_allclose(f, fr, atol=1e-14)
np.testing.assert_allclose(face(wm, wm, gwm, gwm, pm, pp, gpm, gpp, n, 0.5, 0.4), 0.0, atol=0.0)
print("PASS independent one-point WP/PP volume and canonical-face oracle")
