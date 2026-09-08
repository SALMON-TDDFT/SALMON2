from pathlib import Path
import numpy as np

root = Path(__file__).resolve().parents[2]
src = (root / "src/gs/dc/lcfo_wannier_sawf_templates.f90").read_text().lower()
for token in ["validate_sawf_global_local_equivalence", "occupied_projector_residual",
              "overlap_residual", "fixed_hamiltonian_residual", "face_residual"]:
    assert token in src, f"global/local gate misses {token}"

theta = .21
u = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]], complex)
s = np.array([[1., .1], [.1, .9]], complex)
h = np.array([[.2, .03], [.03, .8]], complex)
for a in (s, h):
    local = u.conj().T @ a @ u
    recovered = u @ local @ u.conj().T
    assert np.linalg.norm(recovered - a) < 1e-12
print("PASS small-system global/local SAWF operator equivalence fixture")
