#!/usr/bin/env python3
"""Contract that CheFSI solves and exports the full requested LCFO window."""

from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[2]
caller = (ROOT / "src/gs/dc/lcfo.f90").read_text().lower().replace("&", "")
solver = (ROOT / "src/gs/dc/lcfo_diag_chefsi.f90").read_text().lower().replace("&", "")
runner = (ROOT / "tests/dg/run_lcfo_diag_chefsi_requested_count_mpi.py").read_text().lower()
compact_caller = re.sub(r"\s+", "", caller)
compact_solver = re.sub(r"\s+", "", solver)

assert "lcfo_diag_chefsi_residual_tolerance,coefficient_count,n_basis" in compact_caller, (
    "LCFO does not pass its requested coefficient count to CheFSI"
)
assert re.search(
    r"subroutinediag_chefsi\([^)]*residual_tolerance,requested_state_count,n_basis",
    compact_solver,
), "CheFSI interface does not accept the requested state count"
assert "integer,intent(in)::requested_state_count" in compact_solver
assert "requested_state_count<dc%nstate_tot" in compact_solver
assert "requested_state_count>minval(n_mat)" in compact_solver
assert "ceiling(0.05d0*real(requested_state_count,8))" in compact_solver
assert "minval(n_mat)-requested_state_count" in compact_solver
assert "nsub=requested_state_count+nbuffer" in compact_solver
assert '"chefsisubspace:",nsub,"target:",requested_state_count' in compact_solver
assert "nstate=requested_state_count" in compact_solver
assert "nstate=dc%nstate_tot" not in compact_solver
assert "esp_tot(1:nstate,s)=eigenvalue(1:nstate)" in compact_solver
assert "coef_wf(:,1:nstate,s)=fragment_coef(:,1:nstate)" in compact_solver
assert "pdormtr_work=2_8*int(nsub,8)" in compact_solver
assert "int(dense_block_size-1,8)/2_8" in compact_solver
assert "layout_g%nrow_local+layout_g%ncol_local" in compact_solver
assert "lwork_min=max(pdsyevd_work,pdormtr_work,lwork_query)" in compact_solver
assert "lwork_min>int(huge(lwork),8)/10_8" in compact_solver
assert "stderr=subprocess.stdout" in runner
for diagnostic in ("parameter\\s+number", "illegal\\s+value", "\\bpdormtr\\b", "\\bxerbla\\b"):
    assert diagnostic in runner

# Conventional calls remain exactly the old problem when the explicit request
# equals dc%nstate_tot; expanded retained calls increase only target/buffer size.
def dimensions(nstate_tot, requested, basis_rank):
    assert requested >= nstate_tot
    assert requested <= basis_rank
    buffer = min((requested + 19) // 20, basis_rank - requested)
    return requested, requested + max(0, buffer)


assert dimensions(100, 100, 120) == (100, 105)
assert dimensions(100, 112, 120) == (112, 118)
assert dimensions(100, 120, 120) == (120, 120)

print("PASS CheFSI conventional and expanded retained-state dimensions")
