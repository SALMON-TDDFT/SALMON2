import Mathlib

/-!
# SBE velocity gauge (SALMON v2.3.0, `src/ssbe/bloch_solver_ssbe.f90`)

Abstractions: bands are `Fin n`, directions `Fin 3`, the k-point sum is a `Finset` sum over an
arbitrary finite type; momentum matrices are arbitrary complex matrices (hermiticity not needed
for the identities proved here).
-/

namespace SalmonAudit.SBE

open Matrix BigOperators Complex

noncomputable section

variable {n : ℕ}

/-- SBE-01 code: `calc_hrho_bloch_k` (`src/ssbe/bloch_solver_ssbe.f90:387-406`):
`hrho_ij = (ε_i-ε_j) ρ_ij + Σ_d A_d ((P_d ρ)_ij - (ρ P_d)_ij)` with `delta_omega(i,j)=ε_i-ε_j`
(`src/ssbe/gs_info_ssbe.f90:289`). -/
def codeHrho (ε : Fin n → ℝ) (A : Fin 3 → ℝ) (P : Fin 3 → Matrix (Fin n) (Fin n) ℂ)
    (ρ : Matrix (Fin n) (Fin n) ℂ) : Matrix (Fin n) (Fin n) ℂ :=
  fun i j => ((ε i - ε j : ℝ) : ℂ) * ρ i j + ∑ d, (A d : ℂ) * ((P d * ρ) i j - (ρ * P d) i j)

/-- Velocity-gauge Bloch Hamiltonian `H = diag ε + A·P + (A²/2)·1`. -/
def hamVG (ε : Fin n → ℝ) (A : Fin 3 → ℝ) (P : Fin 3 → Matrix (Fin n) (Fin n) ℂ) :
    Matrix (Fin n) (Fin n) ℂ :=
  diagonal (fun i => (ε i : ℂ)) + ∑ d, (A d : ℂ) • P d +
    (((∑ d, A d ^ 2) / 2 : ℝ) : ℂ) • (1 : Matrix (Fin n) (Fin n) ℂ)

/-- Correct right-hand side `[H, ρ]` of `i dρ/dt = [H, ρ]`. -/
def correctComm (ε : Fin n → ℝ) (A : Fin 3 → ℝ) (P : Fin 3 → Matrix (Fin n) (Fin n) ℂ)
    (ρ : Matrix (Fin n) (Fin n) ℂ) : Matrix (Fin n) (Fin n) ℂ :=
  hamVG ε A P * ρ - ρ * hamVG ε A P

theorem sbe01_hrho_ok (ε : Fin n → ℝ) (A : Fin 3 → ℝ) (P : Fin 3 → Matrix (Fin n) (Fin n) ℂ)
    (ρ : Matrix (Fin n) (Fin n) ℂ) : codeHrho ε A P ρ = correctComm ε A P ρ := by
  ext i j
  simp only [codeHrho, correctComm, hamVG, add_mul, mul_add, Matrix.sub_apply, Matrix.add_apply,
    Matrix.diagonal_mul, Matrix.mul_diagonal, Finset.sum_mul, Finset.mul_sum, Matrix.smul_mul,
    Matrix.mul_smul, Matrix.one_mul, Matrix.mul_one, Matrix.sum_apply, Matrix.smul_apply,
    smul_eq_mul]
  push_cast
  simp only [mul_sub, Finset.sum_sub_distrib]
  ring

/-- SBE-03 code: `calc_current_bloch` with `norder_correction=0`
(`src/ssbe/bloch_solver_ssbe.f90:91-93,327,330`, `calc_trace` 428-444), one Cartesian direction:
`J = Re(Σ_k w_k Σ_{ib,jb} ρ_{jb,ib} p_{ib,jb}) / W / V + A · (Σ_k w_k Σ_i Re ρ_ii / W) / V`. -/
def codeJvg {K : Type*} [Fintype K] (w : K → ℝ) (V A : ℝ) (p ρ : K → Matrix (Fin n) (Fin n) ℂ) : ℝ :=
  (∑ k, (w k : ℂ) * ∑ ib, ∑ jb, ρ k jb ib * p k ib jb).re / (∑ k, w k) / V +
    A * ((∑ k, ∑ i, (ρ k i i).re * w k) / (∑ k, w k)) / V

/-- Correct (SALMON TDDFT convention, `src/common/density_matrix.f90:422`: matter current of
electrons, no charge sign): `J = (1/(W V)) Σ_k w_k Re Tr[ρ_k (p_k + A·1)]`. -/
def correctJvg {K : Type*} [Fintype K] (w : K → ℝ) (V A : ℝ) (p ρ : K → Matrix (Fin n) (Fin n) ℂ) : ℝ :=
  (∑ k, w k * (Matrix.trace (ρ k * (p k + (A : ℂ) • (1 : Matrix (Fin n) (Fin n) ℂ)))).re) /
    (∑ k, w k) / V

lemma trace_term_re (A : ℝ) (p ρ : Matrix (Fin n) (Fin n) ℂ) :
    (Matrix.trace (ρ * (p + (A : ℂ) • (1 : Matrix (Fin n) (Fin n) ℂ)))).re =
      (∑ ib, ∑ jb, ρ jb ib * p ib jb).re + A * ∑ i, (ρ i i).re := by
  rw [Matrix.mul_add, Matrix.trace_add, Matrix.mul_smul, Matrix.mul_one, Matrix.trace_smul]
  simp only [Matrix.trace, Matrix.diag, Matrix.mul_apply, smul_eq_mul, Complex.add_re,
    Complex.re_sum, Complex.mul_re, Complex.ofReal_re, Complex.ofReal_im, zero_mul, sub_zero,
    Finset.mul_sum]
  rw [Finset.sum_comm]

theorem sbe03_current_vg_ok {K : Type*} [Fintype K] (w : K → ℝ) (V A : ℝ)
    (p ρ : K → Matrix (Fin n) (Fin n) ℂ) : codeJvg w V A p ρ = correctJvg w V A p ρ := by
  unfold codeJvg correctJvg
  simp only [trace_term_re]
  simp only [mul_add, Finset.sum_add_distrib, Complex.re_sum, Complex.re_ofReal_mul]
  have h2 : (∑ k, w k * (A * ∑ i, (ρ k i i).re)) = A * ∑ k, ∑ i, (ρ k i i).re * w k := by
    rw [Finset.mul_sum]
    refine Finset.sum_congr rfl (fun k _ => ?_)
    rw [Finset.mul_sum, Finset.mul_sum, Finset.mul_sum]
    refine Finset.sum_congr rfl (fun i _ => ?_)
    ring
  rw [h2]
  ring

end

end SalmonAudit.SBE
