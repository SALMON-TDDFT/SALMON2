import Mathlib

/-!
# DFT+U (src/plusu) — Dudarev energy, potential, spin factor, total-energy bookkeeping

Abstractions: the occupation matrix of one (atom, l, n, spin) block is an arbitrary
square matrix `n : Matrix (Fin k) (Fin k) 𝕜` (`k = 2l+1`); the projections
`<φ_m|ψ>` are not modelled — only the algebra that SALMON applies to `dm_mms_nla`.
-/

namespace SalmonAudit.MISC.PlusU

open Matrix BigOperators

noncomputable section

variable {k : ℕ}

/-- Code: `src/plusu/density_matrix_and_energy_plusu.f90:136-146`
`ztmp = dm(m2,m2) - Σ_{m1} dm(m2,m1)*dm(m1,m2)`, summed over m2, times `0.5*U_eff`. -/
def codeEU (U : ℝ) (n : Matrix (Fin k) (Fin k) ℝ) : ℝ :=
  0.5 * U * ∑ m2, (n m2 m2 - ∑ m1, n m2 m1 * n m1 m2)

/-- Correct Dudarev energy per spin block: `U/2 · Tr(n - n n)`. -/
def correctEU (U : ℝ) (n : Matrix (Fin k) (Fin k) ℝ) : ℝ :=
  U / 2 * Matrix.trace (n - n * n)

/-- MISC-01 (OK): the loop computes `U/2 Tr(n - n²)`. -/
theorem EU_ok (U : ℝ) (n : Matrix (Fin k) (Fin k) ℝ) : codeEU U n = correctEU U n := by
  unfold codeEU correctEU
  simp only [Matrix.trace, Matrix.diag, Matrix.sub_apply, Matrix.mul_apply]
  ring

/-- Code: `src/plusu/density_matrix_and_energy_plusu.f90:164` diagonal potential element
`V = U_eff*(0.5 - dm(m,m))`, where for `nspin=1` `dm` is built with `rocc=2`
(`src/common/initialization.f90:229`), i.e. `dm = N = n↑ + n↓ = 2 n^σ`. -/
def codeVdiag (U Ntot : ℝ) : ℝ := U * (0.5 - Ntot)

/-- Correct Dudarev potential acting on a spin-σ orbital: `U (1/2 - n^σ)`; for an
unpolarized system `n^σ = N/2`. -/
def correctVdiag (U Ntot : ℝ) : ℝ := U * (1 / 2 - Ntot / 2)

/-- MISC-02 (BUG, nspin=1): the code uses the spin-summed occupation. -/
theorem Vdiag_nspin1_bug (U Ntot : ℝ) :
    codeVdiag U Ntot - correctVdiag U Ntot = - U * Ntot / 2 := by
  unfold codeVdiag correctVdiag; ring

/-- Concrete: U = 1 Ha, half-filled orbital in a closed-shell calc (N = 1, n^σ = 1/2):
correct shift 0, code shift -1/2 Ha. -/
theorem Vdiag_nspin1_counterexample : codeVdiag 1 1 ≠ correctVdiag 1 1 := by
  unfold codeVdiag correctVdiag; norm_num

/-- Code total energy (periodic): `src/common/total_energy.f90:757-768,788-789` puts
`Σ f <ψ|H - T - V_loc|ψ>` into `E_ion_nloc`; `H` contains `pseudo_plusU`
(`src/common/hamiltonian.f90:456-458`), so the +U part entering `E_tot`
(`total_energy.f90:503`) is `Σ_{ab} V_ab n_ba` with `V = U(δ/2 - n)`
(`density_matrix_and_energy_plusu.f90:163-167`); `E_U` itself is never added. -/
def codeEtotU (U : ℝ) (n : Matrix (Fin k) (Fin k) ℝ) : ℝ :=
  ∑ a, ∑ b, (U * ((if a = b then (1 : ℝ) / 2 else 0) - n a b)) * n b a

/-- MISC-03 (BUG): what enters `E_tot` differs from the Dudarev energy by `-U/2 Tr(n²)`. -/
theorem EtotU_bug (U : ℝ) (n : Matrix (Fin k) (Fin k) ℝ) :
    codeEtotU U n - correctEU U n = - U / 2 * Matrix.trace (n * n) := by
  unfold codeEtotU correctEU
  have h1 : ∀ a : Fin k, ∑ b, (U * ((if a = b then (1:ℝ) / 2 else 0) - n a b)) * n b a
      = U / 2 * n a a - U * ∑ b, n a b * n b a := by
    intro a
    simp only [mul_sub, sub_mul, Finset.sum_sub_distrib, Finset.mul_sum]
    congr 1
    rw [Finset.sum_eq_single a]
    · rw [if_pos rfl]; ring
    · intro b _ hb; simp [Ne.symm hb]
    · simp
    · apply Finset.sum_congr rfl; intro b _; ring
  simp only [h1, Finset.sum_sub_distrib, Matrix.trace, Matrix.diag, Matrix.sub_apply,
    Matrix.mul_apply]
  simp only [← Finset.mul_sum]
  ring

/-- Concrete: single orbital (k = 1) fully occupied, U = 1: code 1/2 - 1 = -1/2,
Dudarev energy 0. -/
theorem EtotU_counterexample :
    codeEtotU (k := 1) 1 (Matrix.of fun _ _ => 1) ≠
      correctEU (k := 1) 1 (Matrix.of fun _ _ => 1) := by
  simp [codeEtotU, correctEU, Matrix.trace, Matrix.mul_apply, Matrix.sub_apply]
  norm_num

/-- Code off-diagonal/diagonal potential coefficient multiplying `|φ_a><φ_b|`
(`pseudo_pt_plusu.f90:130` with `V_eff(iprj)`, iprj = (a,b), from
`density_matrix_and_energy_plusu.f90:164-166`): `U(δ_ab/2 - dm_ab)`, where
`dm_ab = Σ f <ψ|φ_a><φ_b|ψ>` (`:121`). -/
def codeVab (U : ℝ) (dm : Matrix (Fin k) (Fin k) ℂ) (a b : Fin k) : ℂ :=
  U * ((if a = b then (1 : ℂ) / 2 else 0) - dm a b)

/-- Correct coefficient `∂E_U/∂(<ψ|φ_a>)` multiplying `|φ_a><φ_b|ψ>`: `U(δ_ab/2 - dm_ba)`. -/
def correctVab (U : ℝ) (dm : Matrix (Fin k) (Fin k) ℂ) (a b : Fin k) : ℂ :=
  U * ((if a = b then (1 : ℂ) / 2 else 0) - dm b a)

/-- MISC-04 (NOTE): equal whenever `dm` is symmetric (real Hermitian, the time-reversal
symmetric GS case) ... -/
theorem Vab_ok_of_symm (U : ℝ) (dm : Matrix (Fin k) (Fin k) ℂ) (h : dm.transpose = dm)
    (a b : Fin k) : codeVab U dm a b = correctVab U dm a b := by
  unfold codeVab correctVab
  have : dm b a = dm a b := by
    have := congrFun (congrFun h a) b
    simpa [Matrix.transpose_apply] using this
  rw [this]

/-- ... but transposed (= complex-conjugated) for a complex Hermitian occupation matrix. -/
theorem Vab_differs_complex :
    codeVab (k := 2) 1 !![0, Complex.I; -Complex.I, 0] 0 1 ≠
    correctVab (k := 2) 1 !![0, Complex.I; -Complex.I, 0] 0 1 := by
  unfold codeVab correctVab
  simp
  intro h
  have := congrArg Complex.im h
  norm_num at this

/-- MISC-05 (OK): `U_eff = U - J` (Dudarev), `plusu_global.f90:95`. -/
def codeUeff (U J : ℝ) : ℝ := U - J
theorem Ueff_ok (U J : ℝ) : codeUeff U J = U - J := rfl

end

end SalmonAudit.MISC.PlusU
