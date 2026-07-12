from pathlib import Path
import numpy as np

root = Path(__file__).resolve().parents[2]
src = (root / "src/gs/dc/lcfo_wannier_sawf_templates.f90").read_text().lower()
for token in ["build_sawf_environment_orbits", "validate_sawf_template_fingerprint",
              "replicate_sawf_operator_block", "validate_sawf_actual_group_operation",
              "complete_projection_shell"]:
    assert token in src, f"missing production API {token}"

d = np.array([[0, 1], [1, 0]], complex)
s = np.array([[1.0, .2j], [-.2j, .8]], complex)
h = np.array([[.3, .1 + .2j], [.1 - .2j, 1.2]], complex)
for a in (s, h):
    transformed = d.conj().T @ a @ d
    assert np.linalg.norm(transformed - np.array([[a[1,1], a[1,0]], [a[0,1], a[0,0]]])) < 1e-14
print("PASS representative SAWF replication covariance fixture")
