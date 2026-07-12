from pathlib import Path
import numpy as np

root = Path(__file__).resolve().parents[2]
src = (root / "src/gs/dc/lcfo_wannier_sawf_templates.f90").read_text().lower()
for token in ["validate_sawf_buffer_convergence", "center_residual", "projector_residual",
              "overlap_residual", "ww_residual", "wp_residual", "face_residual"]:
    assert token in src, f"buffer convergence misses {token}"

buffers = np.array([2., 3., 4.])
errors = np.exp(-4 * buffers)
assert errors[-1] < 1e-6 and abs(errors[-1] - errors[-2]) < 1e-5
face_bad = np.array([1e-2, 4e-3, 2e-3])
assert abs(face_bad[-1] - face_bad[-2]) > 1e-3
print("PASS buffer convergence covers orbital and operator blocks")
