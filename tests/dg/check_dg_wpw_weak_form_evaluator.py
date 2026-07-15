#!/usr/bin/env python3
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
src = ROOT / "src/rt/dg/rt_dg_wpw_weak_form_evaluator.f90"
assert src.exists(), "missing WPW SIPG weak-form evaluator"
text = src.read_text().lower()
for token in ("wpw_volume_weak_pair", "wpw_sipg_face_pair", "jump_u", "jump_v",
              "average_dn_u", "average_dn_v", "sipg_penalty_factor", "penalty_over_h"):
    assert token in text, f"missing SIPG evaluator token: {token}"

def face(u, up, gu, gup, v, vp, gv, gvp, n, h_normal):
    ju, jv = u-up, v-vp
    adu = 0.5*(np.dot(gu,n)+np.dot(gup,n))
    adv = 0.5*(np.dot(gv,n)+np.dot(gvp,n))
    return -0.5*(np.conj(adu)*jv + np.conj(ju)*adv) + (10.0/h_normal)*np.conj(ju)*jv

rng = np.random.default_rng(20260715)
n = np.array([1.,0.,0.])
z = lambda: rng.normal()+1j*rng.normal()
g = lambda: rng.normal(size=3)+1j*rng.normal(size=3)
u, up, v, vp, gu, gup, gv, gvp = z(),z(),z(),z(),g(),g(),g(),g()
auv = face(u,up,gu,gup,v,vp,gv,gvp,n,2.5)
avu = face(v,vp,gv,gvp,u,up,gu,gup,n,2.5)
np.testing.assert_allclose(auv, np.conj(avu), rtol=2e-15, atol=2e-15)
np.testing.assert_allclose(face(u,u,gu,gu,v,v,gv,gv,n,2.5), 0.0, atol=0.0)
print("PASS WPW volume/SIPG face weak-form evaluator contract")
