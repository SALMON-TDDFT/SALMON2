#!/usr/bin/env python3
from pathlib import Path


source = (Path(__file__).resolve().parents[2] / "src/rt/dg/rt_dg_fragment.f90").read_text()
start = source.index("    subroutine build_full_h_seed_global_wannier_position_operator")
end = source.index("    end subroutine build_full_h_seed_global_wannier_position_operator", start)
body = source[start:end]

assert "principal_gram = matmul(conjg(transpose(eig_w" in body
assert "[DG-FULL-H-SEED-XI-PRINCIPAL]" in body
assert "rank_sigma_gt_099" in body
print("PASS Full-H/Wannier principal-angle diagnostics are present")
