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


def matmul(a, b):
    return [
        [sum(a[i][k] * b[k][j] for k in range(len(b))) for j in range(len(b[0]))]
        for i in range(len(a))
    ]


def matsub(a, b):
    return [[a[i][j] - b[i][j] for j in range(len(a[0]))] for i in range(len(a))]


def matadd(a, b):
    return [[a[i][j] + b[i][j] for j in range(len(a[0]))] for i in range(len(a))]


def scale(a, factor):
    return [[factor * value for value in row] for row in a]


# A deterministic non-orthogonal two-state model checks the expanded S-metric
# commutator, its generalized-eigen matrix-element identity, and origin shift.
s = [[2.0, 0.0], [0.0, 3.0]]
s_inv = [[0.5, 0.0], [0.0, 1.0 / 3.0]]
h = [[4.0, 0.0], [0.0, 9.0]]
z_s = [[0.0, 1.25], [1.25, 0.0]]
v_cov = scale(matsub(matmul(matmul(h, s_inv), z_s), matmul(matmul(z_s, s_inv), h)), 1j)
v_action = matmul(s_inv, v_cov)
assert abs(v_action[0][1] + 0.625j) < 1.0e-14
assert abs(v_action[1][0] - (1.25 / 3.0) * 1j) < 1.0e-14

c1 = [1.0 / (2.0**0.5), 0.0]
c2 = [0.0, 1.0 / (3.0**0.5)]
lhs = sum(c1[i].conjugate() * v_cov[i][j] * c2[j] for i in range(2) for j in range(2))
rhs = 1j * (2.0 - 3.0) * sum(
    c1[i].conjugate() * z_s[i][j] * c2[j] for i in range(2) for j in range(2)
)
assert abs(lhs - rhs) < 1.0e-14

z_shift = matadd(z_s, scale(s, 0.37))
v_shift = scale(matsub(matmul(matmul(h, s_inv), z_shift), matmul(matmul(z_shift, s_inv), h)), 1j)
assert max(abs(v_shift[i][j] - v_cov[i][j]) for i in range(2) for j in range(2)) < 1.0e-14

# Exercise the complete periodic projected-position chain on a diagonal but
# non-orthogonal two-state model.  Deliberate radial errors in U_tilde test
# the scalar polar factors, while tracked eigenphases test Log_b and covariance.
import cmath
import math

l_z = 12.0
x = [[1.0 / math.sqrt(2.0), 0.0], [0.0, 1.0 / math.sqrt(3.0)]]
phases = [0.23, -0.41]
radial_errors = [0.91, 1.08]
u_tilde_seed = [radial_errors[i] * cmath.exp(1j * phases[i]) for i in range(2)]
# U_S=X^-H U_tilde X^-1, followed by the normative recovery X^H U_S X.
u_s = [[2.0 * u_tilde_seed[0], 0.0], [0.0, 3.0 * u_tilde_seed[1]]]
u_tilde_matrix = matmul(matmul(x, u_s), x)
u_tilde = [u_tilde_matrix[0][0], u_tilde_matrix[1][1]]
assert abs(u_tilde_matrix[0][1]) < 1.0e-14 and abs(u_tilde_matrix[1][0]) < 1.0e-14
assert min(abs(value) for value in u_tilde) > 0.5
q = [value / abs(value) for value in u_tilde]
assert max(abs(value.conjugate() * value - 1.0) for value in q) < 1.0e-14
tracked_angles = [math.atan2(value.imag, value.real) for value in q]
z_tilde_diag = [l_z * angle / (2.0 * math.pi) for angle in tracked_angles]
assert all(abs(value.imag) < 1.0e-14 for value in [complex(z) for z in z_tilde_diag])
q_reconstructed = [cmath.exp(1j * 2.0 * math.pi * value / l_z) for value in z_tilde_diag]
assert max(abs(q_reconstructed[i] - q[i]) for i in range(2)) < 1.0e-14

# Z_S=S X Z_tilde X^H S; for diagonal S and X this is sqrt(S) Z sqrt(S).
z_from_phase = [[2.0 * z_tilde_diag[0], 0.0], [0.0, 3.0 * z_tilde_diag[1]]]
assert abs(z_from_phase[0][1] - z_from_phase[1][0].conjugate()) < 1.0e-14
z0 = 0.37
theta0 = 2.0 * math.pi * z0 / l_z
q_shift = [value * cmath.exp(1j * theta0) for value in q]
shifted_angles = [math.atan2(value.imag, value.real) for value in q_shift]
z_tilde_shift = [l_z * angle / (2.0 * math.pi) for angle in shifted_angles]
z_from_phase_shift = [[2.0 * z_tilde_shift[0], 0.0], [0.0, 3.0 * z_tilde_shift[1]]]
expected_shift = matadd(z_from_phase, scale(s, z0))
assert max(
    abs(z_from_phase_shift[i][j] - expected_shift[i][j]) for i in range(2) for j in range(2)
) < 1.0e-13

# Expand the declared SIPG face form rather than comparing literal constants.
def sipg_face(v_m, v_p, dnv_m, dnv_p, u_m, u_p, dnu_m, dnu_p, sigma):
    jump_v = v_m.conjugate() - v_p.conjugate()
    jump_u = u_m - u_p
    avg_dv = 0.5 * (dnv_m.conjugate() + dnv_p.conjugate())
    avg_du = 0.5 * (dnu_m + dnu_p)
    return -0.5 * jump_v * avg_du - 0.5 * avg_dv * jump_u + sigma * jump_v * jump_u


sigma = 10.0 / 0.4
assert sipg_face(1, 0, 0, 0, 1, 0, 1, 0, sigma) == -0.25 + sigma
assert sipg_face(1, 0, 0, 0, 0, 1, 0, 1, sigma) == -0.25 - sigma

print("PASS Wannier+PW DG mathematical operator contract")
