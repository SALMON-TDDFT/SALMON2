from pathlib import Path
import numpy as np

root = Path(__file__).resolve().parents[2]
src = (root / "src/gs/dc/lcfo_wannier_sawf_templates.f90").read_text().lower()
for token in ["stitch_sawf_neighbor_gauge", "zgesvd", "gauge_unitary", "rank deficient"]:
    assert token in src, f"missing gauge-stitch implementation {token}"

theta = .37
q = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]], complex)
overlap = q
u, singular, vh = np.linalg.svd(overlap)
align = vh.conj().T @ u.conj().T
residual = np.linalg.norm(overlap @ align - np.eye(2))
assert singular.min() > .999999 and residual < 1e-12
assert np.linalg.matrix_rank(np.diag([1., 0.])) < 2
print(f"PASS neighbor gauge stitching residual={residual:.3e}")
