#!/usr/bin/env python3
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
source = (ROOT / "src/gs/dc/lcfo_flux.f90").read_text()


face_body = source[source.index("subroutine prepare_wpw_canonical_face_trace_provider"):source.index("end subroutine prepare_wpw_canonical_face_trace_provider")]
assert "local_grid<mg%is" in face_body.replace(" ", ""), \
    "the rank owning the sampled core/buffer coordinate must publish the face trace"

for label in ("face_reduce", "face_exchange", "face_bind"):
    assert f"[DG-WPW-LOCAL-FAIL] {label}" in source, \
        f"missing rank-local diagnostic for {label} failure"

for field in ("owned=", "w_fail=", "p_fail=", "mg_is=", "mg_ie="):
    assert field in source[source.index("[DG-WPW-LOCAL-FAIL] face_reduce"):], \
        f"face reduction diagnostic does not report {field}"

print("PASS canonical-face rank ownership and local diagnostics contract")
