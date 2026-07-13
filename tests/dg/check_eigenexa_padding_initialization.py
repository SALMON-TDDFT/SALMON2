#!/usr/bin/env python3
from pathlib import Path

root = Path(__file__).resolve().parents[2]
source = (root / "src/gs/eigen_subdiag_eigenexa.f90").read_text().lower()

for name in ("eigen_pdsyevd_ex", "eigen_pdsyevd_ex_red_mem"):
    body = source.split(f"subroutine {name}", 1)[1].split(f"end subroutine {name}", 1)[0]
    allocate = body.index("allocate( h_div")
    solve = body.index("call eigen_sx")
    between = body[allocate:solve]
    assert "h_div = 0d0" in between or "h_div=0d0" in between, \
        f"{name} leaves EigenExa input padding uninitialized"
    assert "v_div = 0d0" in between or "v_div=0d0" in between, \
        f"{name} leaves EigenExa output padding uninitialized"

print("PASS EigenExa distributed padding is initialized")
