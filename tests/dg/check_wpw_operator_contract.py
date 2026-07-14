from pathlib import Path


root = Path(__file__).resolve().parents[2]
contract_path = root / "docs/plans/2026-07-12-wannier-pw-dg-operator-contract.md"

assert contract_path.exists(), f"missing operator contract: {contract_path}"
contract = contract_path.read_text().lower()

required_decisions = {
    "unique fragment-core ownership": "unique grid ownership",
    "fragment-restricted Wannier support": "fragment-restricted wannier",
    "PW support": "plane-wave enrichment support",
    "PW normalization": "plane-wave enrichment normalization",
    "jump convention": "[[u]]",
    "average convention": "{{grad u}}",
    "face orientation": "canonical face normal",
    "penalty scaling": "sigma_f = 10 / h_f",
    "WW block": "h_ww",
    "WP block": "h_wp",
    "PW block": "h_pp",
    "face-neighbor block": "h_face",
    "ownership separation": "physical coupling is independent of mpi ownership",
    "periodic length gauge": "periodic length-gauge",
    "polarization branch": "polarization branch",
    "periodic phase matrix": "u_s = <b|exp(i 2 pi z/l_z)|b>",
    "unitary polar factor": "unitary polar factor",
    "projected position logarithm": "z_tilde = (l_z / 2 pi) (-i log_b q)",
    "position-interface decision": "no independent interface position correction",
    "S-metric velocity relation": "v_z = i s^-1 (h s^-1 z_s - z_s s^-1 h)",
    "DG eigenstate seed": "h_dg c = s c epsilon",
    "namelist policy": "scientific controls belong in the namelist",
    "scope": "gapped, integer-occupation lda",
}

for decision, token in required_decisions.items():
    assert token in contract, f"operator contract is missing {decision}: {token}"

assert "nonzero interface position correction" in contract
assert "rejected" in contract
assert "z_s = x^-h z_tilde x^-1" in contract
assert "q = u_tilde (u_tilde^h u_tilde)^(-1/2)" in contract
assert "z_s = s x z_tilde x^h s" in contract
assert "x is square and full rank" in contract
assert "minimum singular value" in contract and "fatal" in contract
assert "z_s -> z_s + z0 s" in contract
assert "every support-overlap block" in contract
assert "projector-mediated block" in contract
assert "independent of direct support intersection" in contract
assert "not face neighbors" in contract
assert "<b_a|v_nl|b_b>" in contract
assert "finite-difference dpz/dt residual" in contract
assert "face-omitted velocity residual" in contract


import numpy as np


def hermitian_inverse_sqrt(matrix):
    values, vectors = np.linalg.eigh(matrix)
    assert values.min() > 1.0e-12
    return (vectors * values**-0.5) @ vectors.conj().T


# Dense, noncommuting projected periodic-position construction.
s_seed = np.array([[1.7, 0.22 + 0.11j], [0.22 - 0.11j, 1.15]], complex)
x = hermitian_inverse_sqrt(s_seed)
np.testing.assert_allclose(x.conj().T @ s_seed @ x, np.eye(2), atol=2.0e-14)
rotation = np.array([[np.cos(0.37), -np.sin(0.37)], [np.sin(0.37), np.cos(0.37)]], complex)
q_seed = rotation @ np.diag(np.exp(1j * np.array([-0.43, 0.28]))) @ rotation.conj().T
r_positive = np.array([[1.13, 0.17j], [-0.17j, 0.82]], complex)
u_tilde_seed = q_seed @ r_positive
assert np.linalg.norm(u_tilde_seed.conj().T @ u_tilde_seed - u_tilde_seed @ u_tilde_seed.conj().T) > 1.0e-3
x_inv = np.linalg.inv(x)
u_s = x_inv.conj().T @ u_tilde_seed @ x_inv
u_tilde = x.conj().T @ u_s @ x
polar_metric_inv_sqrt = hermitian_inverse_sqrt(u_tilde.conj().T @ u_tilde)
q = u_tilde @ polar_metric_inv_sqrt
np.testing.assert_allclose(q.conj().T @ q, np.eye(2), atol=3.0e-14)

phase_values, phase_vectors = np.linalg.eig(q)
order = np.argsort(np.angle(phase_values))
phase_values = phase_values[order]
phase_vectors = phase_vectors[:, order]
phase_vectors, _ = np.linalg.qr(phase_vectors)
phases = np.angle(phase_values)
l_z = 12.0
z_tilde = phase_vectors @ np.diag(l_z * phases / (2.0 * np.pi)) @ phase_vectors.conj().T
np.testing.assert_allclose(z_tilde, z_tilde.conj().T, atol=2.0e-14)
q_reconstructed = phase_vectors @ np.diag(np.exp(1j * phases)) @ phase_vectors.conj().T
np.testing.assert_allclose(q_reconstructed, q, atol=4.0e-14)
z_s = s_seed @ x @ z_tilde @ x.conj().T @ s_seed
np.testing.assert_allclose(z_s, z_s.conj().T, atol=3.0e-14)
z0 = 0.31
theta0 = 2.0 * np.pi * z0 / l_z
u_s_shift = np.exp(1j * theta0) * u_s
u_tilde_shift = x.conj().T @ u_s_shift @ x
q_shift = u_tilde_shift @ hermitian_inverse_sqrt(u_tilde_shift.conj().T @ u_tilde_shift)
shift_values, shift_vectors = np.linalg.eig(q_shift)
shift_order = np.argsort(np.angle(shift_values))
shift_values = shift_values[shift_order]
shift_vectors = shift_vectors[:, shift_order]
shift_vectors, _ = np.linalg.qr(shift_vectors)
shift_phases = np.angle(shift_values)
for i in range(2):
    while shift_phases[i] - phases[i] - theta0 > np.pi:
        shift_phases[i] -= 2.0 * np.pi
    while shift_phases[i] - phases[i] - theta0 < -np.pi:
        shift_phases[i] += 2.0 * np.pi
z_tilde_shift = shift_vectors @ np.diag(l_z * shift_phases / (2.0 * np.pi)) @ shift_vectors.conj().T
z_s_shift = s_seed @ x @ z_tilde_shift @ x.conj().T @ s_seed
np.testing.assert_allclose(z_s_shift, z_s + z0 * s_seed, atol=8.0e-14)


# Independent Heisenberg check: finite-difference d<Pz>/dt from generalized
# propagation must equal the S-metric commutator expectation.  H_face is an
# explicit DG face block and is required in both propagation and velocity.
h_volume = np.array([[0.82, 0.09 - 0.04j], [0.09 + 0.04j, 1.31]], complex)
h_face = np.array([[0.17, -0.21 + 0.07j], [-0.21 - 0.07j, -0.06]], complex)
h_full = h_volume + h_face
s_inv = np.linalg.inv(s_seed)
v_cov = 1j * (h_full @ s_inv @ z_s - z_s @ s_inv @ h_full)
d0 = np.array([0.77 + 0.13j, -0.31 + 0.52j], complex)
d0 /= np.linalg.norm(d0)
c0 = x @ d0
h_tilde = x.conj().T @ h_full @ x
eval_h, vec_h = np.linalg.eigh(h_tilde)


def propagated_coeff(dt):
    d_t = vec_h @ (np.exp(-1j * eval_h * dt) * (vec_h.conj().T @ d0))
    return x @ d_t


velocity_expectation = np.vdot(c0, v_cov @ c0).real
fd_residuals = []
for dt in (2.0e-2, 1.0e-2, 5.0e-3):
    c_plus, c_minus = propagated_coeff(dt), propagated_coeff(-dt)
    p_plus = np.vdot(c_plus, z_s @ c_plus).real
    p_minus = np.vdot(c_minus, z_s @ c_minus).real
    dpdt_fd = (p_plus - p_minus) / (2.0 * dt)
    fd_residuals.append(abs(dpdt_fd - velocity_expectation))
assert fd_residuals[0] / fd_residuals[1] > 3.8
assert fd_residuals[1] / fd_residuals[2] > 3.8
fd_residual = fd_residuals[-1]
v_face_omitted = 1j * (h_volume @ s_inv @ z_s - z_s @ s_inv @ h_volume)
face_omitted_residual = abs(dpdt_fd - np.vdot(c0, v_face_omitted @ c0).real)
assert fd_residual < 2.0e-6
assert face_omitted_residual > 1.0e-3
print(f"DG velocity residuals: finite-difference={fd_residual:.3e} face-omitted={face_omitted_residual:.3e}")

print("PASS Wannier+PW DG mathematical operator contract")
