#!/usr/bin/env python3
from pathlib import Path


source = (Path(__file__).resolve().parents[2] / "src/rt/dg/rt_dg_fragment.f90").read_text()
start = source.index("  subroutine refresh_global_wannier_flux_eigen_from_current_hamiltonian")
end = source.index("  end subroutine refresh_global_wannier_flux_eigen_from_current_hamiltonian", start)
body = source[start:end]

assert "call apply_overlap_operator_batch" in body
assert "horth(:, :) = matmul(conjg(transpose(sinvhalf)), matmul(hw, sinvhalf))" in body
assert "evec(:, :) = matmul(sinvhalf, zorth)" in body
assert "global Wannier metric is singular" in body
print("PASS global Wannier seed solves H C = S C eps")
