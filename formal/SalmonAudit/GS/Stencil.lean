import Mathlib

/-!
# GS: finite-difference kinetic operator (SALMON v2.3.0)

Abstractions: grid functions are replaced by the 1-D Taylor moment conditions of the
stencil weights; the Bloch-transformed kinetic operator is checked as a polynomial identity
in a commuting symbol `d` (= ∂, acting on the periodic part u_k), which is exact because
`k` is constant.
-/

noncomputable section
namespace SalmonAudit.GS

/-- Second-derivative weights `cNmat(·,4)` (8th order), `src/common/initialization.f90:1391-1395`. -/
def c4 : Fin 5 → ℚ := ![-205/72, 8/5, -1/5, 8/315, -1/560]

/-- First-derivative weights `bNmat(·,4)` (8th order), `src/common/initialization.f90:1360-1363`. -/
def b4 : Fin 4 → ℚ := ![4/5, -1/5, 4/105, -1/280]

/-- Even moment `c0·0^m + 2 Σ_{k=1..4} c_k k^m` of the symmetric 2nd-derivative stencil. -/
def c4moment (m : ℕ) : ℚ :=
  (if m = 0 then c4 0 else 0) + 2 * ∑ k : Fin 4, c4 k.succ * ((k : ℚ) + 1) ^ m

/-- Odd moment `2 Σ_k b_k k^m` of the antisymmetric 1st-derivative stencil. -/
def b4moment (m : ℕ) : ℚ := 2 * ∑ k : Fin 4, b4 k * ((k : ℚ) + 1) ^ m

/-- GS-01 (OK): the 8th-order Laplacian weights reproduce f'' exactly for polynomials up to
degree 9: moments 0,2,4,6,8 equal (0, 2, 0, 0, 0), i.e. `Σ c_k f(x+kh) = h² f''(x) + O(h^10)`. -/
theorem c4_moments_ok :
    c4moment 0 = 0 ∧ c4moment 2 = 2 ∧ c4moment 4 = 0 ∧ c4moment 6 = 0 ∧ c4moment 8 = 0 := by
  simp [c4moment, c4, Fin.sum_univ_four]
  norm_num

/-- GS-01 (OK): the 8th-order gradient weights: odd moments 1,3,5,7 equal (1,0,0,0). -/
theorem b4_moments_ok :
    b4moment 1 = 1 ∧ b4moment 3 = 0 ∧ b4moment 5 = 0 ∧ b4moment 7 = 0 := by
  simp [b4moment, b4, Fin.sum_univ_four]
  norm_num

/-- Code (orthogonal, periodic): `htpsi = (V + lap0 + ½|k+A|²) u − ½ Σ lapt·u − i Σ (k+A)_j nab_j u`,
`src/common/hamiltonian.f90:157-161`, `src/common/stencil.f90` (zstencil); here in symbols
(`d` = ∂ acting on u, `k` = k+A, one Cartesian component; d, k commuting complex symbols). -/
def kinCode (k d : ℂ) : ℂ := -(1/2) * d ^ 2 - Complex.I * k * d + (1/2) * k ^ 2

/-- Correct: `e^{-ikx} (−½ ∂²) e^{ikx} = −½ (∂ + i k)²`. -/
def kinCorrect (k d : ℂ) : ℂ := -(1/2) * (d + Complex.I * k) ^ 2

/-- GS-02 (OK): Bloch kinetic operator. -/
theorem kin_ok (k d : ℂ) : kinCode k d = kinCorrect k d := by
  unfold kinCode kinCorrect
  have : (d + Complex.I * k) ^ 2 = d ^ 2 + 2 * Complex.I * k * d + Complex.I ^ 2 * k ^ 2 := by ring
  rw [this, Complex.I_sq]
  ring

/-- Code: nonorthogonal Laplacian `Σ_i F_ii ∂_i² + F(4)∂_2∂_3 + F(5)∂_1∂_3 + F(6)∂_1∂_2`
with `F = B Bᵀ`, `F(4)=F23+F32` etc., `src/common/lattice.f90:58-78`,
`src/common/stencil.f90:171-205`. `x` is the vector of commuting reduced-direction derivatives. -/
def lapNonorthCode (F : Matrix (Fin 3) (Fin 3) ℝ) (x : Fin 3 → ℝ) : ℝ :=
  F 0 0 * x 0 ^ 2 + F 1 1 * x 1 ^ 2 + F 2 2 * x 2 ^ 2
  + (F 1 2 + F 2 1) * (x 1 * x 2) + (F 0 2 + F 2 0) * (x 0 * x 2) + (F 0 1 + F 1 0) * (x 0 * x 1)

/-- Correct: `∇_r² = (Bᵀ∇_u)·(Bᵀ∇_u) = Σ_ij (B Bᵀ)_ij ∂_i ∂_j`. -/
def lapNonorthCorrect (B : Matrix (Fin 3) (Fin 3) ℝ) (x : Fin 3 → ℝ) : ℝ :=
  ∑ a : Fin 3, (∑ i : Fin 3, B i a * x i) ^ 2

/-- GS-03 (OK): nonorthogonal Laplacian coefficients `F = B Bᵀ`, cross terms doubled. -/
theorem lapNonorth_ok (B : Matrix (Fin 3) (Fin 3) ℝ) (x : Fin 3 → ℝ) :
    lapNonorthCode (B * B.transpose) x = lapNonorthCorrect B x := by
  simp [lapNonorthCode, lapNonorthCorrect, Matrix.mul_apply, Fin.sum_univ_three]
  ring

/-- Code: k·∇ term in reduced coords `(B k)·∇_u` (`kAc = matmul(rmatrix_B,kAc)`,
`src/common/hamiltonian.f90:352-356`). Correct: `k·∇_r = k·(Bᵀ ∇_u)`. GS-03 (OK). -/
theorem kdotgrad_nonorth_ok (B : Matrix (Fin 3) (Fin 3) ℝ) (k g : Fin 3 → ℝ) :
    (B.mulVec k) ⬝ᵥ g = k ⬝ᵥ (B.transpose.mulVec g) := by
  simp [Matrix.mulVec, dotProduct, Fin.sum_univ_three]
  ring

end SalmonAudit.GS
end
