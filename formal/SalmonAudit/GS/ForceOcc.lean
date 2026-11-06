import Mathlib

/-!
# GS: forces (Fourier/nonlocal), occupations, k-weights, mixing, KB normalisation,
# multipole boundary (SALMON v2.3.0)

Abstractions: G-sums are reduced to a single G (the identities are termwise); the derivative
w.r.t. R of a structure-factor phase `e^{-iG·R}` is `-iG e^{-iG·R}` (stated as an algebraic
input `dph = -I*G*ph`).
-/

noncomputable section
namespace SalmonAudit.GS

open Complex

/-- GS-16 (OK). Ewald reciprocal force. Energy term `A |ρ_i|²` with `A = (V/2)·coef·e^{-G²/4a}`,
`∂ρ_i/∂R_a = −iG Z ph/V` (ph = e^{−iGR_a}). `−∂E/∂R_a = −2A·Re(conj ρ_i · ∂ρ_i)`.
Code (`src/common/force.f90:112-117`): `G·Z·coef·e·Im(ρ_i·e^{iGR_a})` with `A = V coef e/2`. -/
theorem ewaldRecForce_ok (V c e G Z : ℝ) (rho : ℂ) (θ : ℝ) (hV : V ≠ 0) :
    -(2 * (V / 2 * c * e) * (star rho * (-I * G * Z * exp (-(θ * I)) / V)).re)
      = G * (Z * c * e * (rho * exp (θ * I)).im) := by
  have hV' : (V : ℂ) ≠ 0 := by exact_mod_cast hV
  have h1 : exp (-(θ * I)) = star (exp (θ * I)) := by
    rw [Complex.star_def, ← Complex.exp_conj]; simp [Complex.conj_ofReal]
  rw [h1]
  have h2 : star rho * (-I * G * Z * star (exp (θ * I)) / V)
      = (-I) * (G * Z / V) * star (rho * exp (θ * I)) := by
    simp only [star_mul']; field_simp
  rw [h2]
  have h3 : (G * Z / V : ℂ) = ((G * Z / V : ℝ) : ℂ) := by push_cast; rfl
  rw [h3]
  simp only [mul_re, mul_im, neg_re, neg_im, I_re, I_im, Complex.star_def, conj_re, conj_im,
    ofReal_re, ofReal_im]
  field_simp
  ring

/-- GS-17 (OK). Local-pseudopotential Fourier force (periodic). `E_loc ∋ Re(conj ρe · VG · ph)`,
`−∂/∂R = −Re(conj ρe · VG · (−iG) ph) `. Code (`src/common/force.f90:116`):
`G·Im(e^{iGR} ρe conj VG)` with VG real. -/
theorem locForce_ok (G VG θ : ℝ) (rhoe : ℂ) :
    -((star rhoe * VG * (-I * G * exp (-(θ * I)))).re)
      = G * (exp (θ * I) * rhoe * star (VG : ℂ)).im := by
  have h1 : exp (-(θ * I)) = star (exp (θ * I)) := by
    rw [Complex.star_def, ← Complex.exp_conj]; simp [Complex.conj_ofReal]
  rw [h1]
  set z := exp (θ * I)
  simp only [mul_re, mul_im, neg_re, neg_im, I_re, I_im, Complex.star_def, conj_re, conj_im,
    ofReal_re, ofReal_im]
  ring

/-- GS-18 (OK). Nonlocal (KB) force. `E = c |P|²` ⇒ `−∂E = −2c Re(conj P · dP)`.
Code (`src/common/force.f90:327-329`): `−2 f w Hvol Re(conj(duVpsi) · uVpsibox2)` with
`uVpsibox2 = P·inorm·Hvol` (`src/common/nonlocal_potential.f90:550`), `duVpsi·Hvol = dP`. -/
theorem nlForce_ok (c : ℝ) (P dP : ℂ) :
    -(2 * c * (star dP * P).re) = -(2 * c * (star P * dP).re) := by
  simp only [Complex.star_def, mul_re, conj_re, conj_im]; ring

/-- Fermi–Dirac occupation (`src/gs/occupation.f90:156`): `wspin/(1+exp((ε−μ)/T))`.
GS-19 (OK): at ε = μ the occupation is wspin/2 and it never exceeds wspin. -/
theorem fermi_ok (w T μ ε : ℝ) (hw : 0 ≤ w) :
    w / (1 + Real.exp ((μ - μ) / T)) = w / 2 ∧ w / (1 + Real.exp ((ε - μ) / T)) ≤ w := by
  constructor
  · simp; norm_num
  · apply div_le_self hw
    linarith [Real.exp_pos ((ε - μ) / T)]

/-- Γ-centred even-N k-weights along one axis (`src/common/lattice.f90:189-201`):
`N+1` points, interior weight `1/N`, both ends `1/(2N)`. -/
def gammaW (N : ℕ) (i : ℕ) : ℝ := if i = 0 ∨ i = N then 1 / (2 * N) else 1 / N

/-- GS-20 (OK). The weights sum to 1 (N ≥ 1). -/
theorem gammaW_sum (N : ℕ) (hN : 1 ≤ N) : ∑ i ∈ Finset.range (N + 1), gammaW N i = 1 := by
  obtain ⟨M, rfl⟩ : ∃ M, N = M + 1 := ⟨N - 1, by omega⟩
  rw [Finset.sum_range_succ, Finset.sum_range_succ']
  have hmid : ∑ i ∈ Finset.range M, gammaW (M + 1) (i + 1) = ∑ i ∈ Finset.range M, (1 / ((M:ℝ) + 1)) := by
    refine Finset.sum_congr rfl ?_
    intro i hi
    have h1 : ¬ (i + 1 = 0 ∨ i + 1 = M + 1) := by simp at hi; omega
    unfold gammaW; rw [ite_eq_right_of_eq_false _ _ (eq_false h1)]
    push_cast; ring_nf
  rw [hmid]
  simp [gammaW]
  field_simp
  ring

/-- Simple mixing (`src/gs/mixing.f90:59-60`, called with c1=1−α, c2=α). GS-21 (OK). -/
theorem simpleMix_ok (α rin rout : ℝ) :
    (1 - α) * rin + α * rout = rin + α * (rout - rin) := by ring

/-- Pulay coefficients normalised by `rc = 1/Σb` (`src/gs/mixing.f90:405-406`). GS-22 (OK):
the normalised coefficients sum to one (constraint of DIIS). -/
theorem pulay_norm_ok {n : ℕ} (b : Fin n → ℝ) (h : ∑ i, b i ≠ 0) :
    ∑ i, (1 / ∑ j, b j) * b i = 1 := by
  rw [← Finset.mul_sum]; field_simp

/-- KB projector normalisation (`src/atom/pp/input_pp.f90:147-155,985-987`,
`src/atom/pp/prep_pp.f90:522`): `uV = β·√|D|`, `V_nl = inorm·|uV⟩⟨uV|`, inorm = sign D.
GS-23 (OK): `inorm·(√|D|)² = D`. -/
theorem kb_ok (D : ℝ) : (SignType.sign D : ℝ) * (Real.sqrt |D|) ^ 2 = D := by
  rw [Real.sq_sqrt (abs_nonneg D)]
  rcases lt_trichotomy D 0 with h | h | h
  · simp [h, abs_of_neg h]
  · simp [h]
  · simp [h, abs_of_pos h]

/-- Multipole boundary (`src/poisson/poisson_isolated.f90:559,962`): with SALMON's
`Ylm(x) = √(4π/(2l+1)) r^l Y⁰_lm` both in `rholm` and in the boundary sum, the product
of the two normalisations is `4π/(2l+1)`, the correct multipole-expansion prefactor
`V = Σ 4π/(2l+1) q_lm Y⁰_lm/r^{l+1}`. GS-24 (OK). -/
theorem multipole_ok (l : ℕ) :
    Real.sqrt (4 * Real.pi / (2 * l + 1)) * Real.sqrt (4 * Real.pi / (2 * l + 1))
      = 4 * Real.pi / (2 * l + 1) := by
  rw [Real.mul_self_sqrt]; positivity

end SalmonAudit.GS
end
