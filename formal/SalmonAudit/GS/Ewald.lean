import Mathlib

/-!
# GS: ion–ion energies and forces (SALMON v2.3.0)

`aEwald` in SALMON is η² (real-space kernel `erfc(sqrt(aEwald*rr))/sqrt(rr)`).
erfc itself is abstracted: we only use `erfc' (x) = -(2/√π) e^{-x²}` as a hypothesis
(`HasDerivAt E (-(2/√π) exp(-x²)) x` for all x).
Structure factors are abstracted as a complex number `S = Σ_a Z_a e^{-iG·R_a}`.
-/

noncomputable section
namespace SalmonAudit.GS

open Real

/-- GS-04 (OK). Self term. Code: `- sqrt(aEwald/Pi)*zps2` (`src/common/total_energy.f90:487`).
Correct (η = √a): `-(η/√π) Σ Z²`. -/
theorem ewald_self_ok (a z2 : ℝ) (ha : 0 ≤ a) :
    -(Real.sqrt (a / π) * z2) = -((Real.sqrt a / Real.sqrt π) * z2) := by
  rw [Real.sqrt_div ha]

/-- GS-05 (OK). Neutralising-background term. Code: `- Pi*zps1**2/(2*aEwald*sysvol)`
(`src/common/total_energy.f90:487`). Correct: `-π Q² / (2 η² V)` with η = √a. -/
theorem ewald_background_ok (a Q V : ℝ) (ha : 0 ≤ a) :
    -(π * Q ^ 2 / (2 * a * V)) = -(π * Q ^ 2 / (2 * (Real.sqrt a) ^ 2 * V)) := by
  rw [Real.sq_sqrt ha]

/-- Code reciprocal Ewald term per G: `sysvol*coef*|rho_i|^2*exp_ewald*0.5` with
`coef = 4π/G²`, `rho_i = S/det_A` (`src/common/total_energy.f90:376`,
`src/atom/pp/prep_pp.f90:877-888`, `src/common/initialization.f90:720-722`). `s2 = |S|²`. -/
def ewaldRecCode (V G2 s2 e : ℝ) : ℝ := V * (4 * π / G2) * (s2 / V ^ 2) * e * 0.5

/-- Correct: `(2π/V) |S(G)|² e^{-G²/4η²} / G²`. -/
def ewaldRecCorrect (V G2 s2 e : ℝ) : ℝ := 2 * π / V * s2 * e / G2

/-- GS-06 (OK). -/
theorem ewaldRec_ok (V G2 s2 e : ℝ) (hV : V ≠ 0) (hG : G2 ≠ 0) :
    ewaldRecCode V G2 s2 e = ewaldRecCorrect V G2 s2 e := by
  unfold ewaldRecCode ewaldRecCorrect
  field_simp
  ring

/-- Code real-space Ewald force magnitude along r̂ (`src/common/force.f90:532-534`):
`F = -ZaZb*(-erfc(√(a rr))/rr - 2*sqrt(a/(rr*Pi))*exp(-a*rr))`, rr = r². -/
def ewaldForceCode (ZZ a r erfcv : ℝ) : ℝ :=
  -ZZ * (-erfcv / r ^ 2 - 2 * Real.sqrt (a / (r ^ 2 * π)) * Real.exp (-a * r ^ 2))

/-- GS-07 (OK). The code force equals `-d/dr [ZZ·erfc(√a r)/r]` (pair energy, the ½ in the
energy is cancelled because each pair appears twice in the bookkeeping list). -/
theorem ewaldForce_ok (ZZ a r : ℝ) (E : ℝ → ℝ) (ha : 0 < a) (hr : 0 < r)
    (hE : ∀ x, HasDerivAt E (-(2 / Real.sqrt π) * Real.exp (-x ^ 2)) x) :
    HasDerivAt (fun x => ZZ * (E (Real.sqrt a * x) / x))
      (-(ewaldForceCode ZZ a r (E (Real.sqrt a * r)))) r := by
  have h1 : HasDerivAt (fun x => E (Real.sqrt a * x))
      ((-(2 / Real.sqrt π) * Real.exp (-(Real.sqrt a * r) ^ 2)) * Real.sqrt a) r := by
    have h0 : HasDerivAt (fun x => Real.sqrt a * x) (Real.sqrt a) r := by
      simpa using (hasDerivAt_id r).const_mul (Real.sqrt a)
    exact (hE (Real.sqrt a * r)).comp r h0
  have h2 := (h1.div (hasDerivAt_id r) hr.ne').const_mul ZZ
  have hfun : (fun x => ZZ * (E (Real.sqrt a * x) / x))
      = (fun y => ZZ * ((fun x => E (Real.sqrt a * x)) y / id y)) := rfl
  rw [hfun]
  convert h2 using 1
  unfold ewaldForceCode
  have hs : (Real.sqrt a * r) ^ 2 = a * r ^ 2 := by
    rw [mul_pow, Real.sq_sqrt ha.le]
  have hsq : Real.sqrt (a / (r ^ 2 * π)) = Real.sqrt a / (r * Real.sqrt π) := by
    rw [Real.sqrt_div ha.le, Real.sqrt_mul (sq_nonneg r), Real.sqrt_sq hr.le]
  rw [hsq, hs]
  have hpi : 0 < Real.sqrt π := Real.sqrt_pos.mpr Real.pi_pos
  simp only [id]
  field_simp
  ring_nf

/-- GS-08 (OK). Isolated ion–ion force. Code (`src/common/force.f90:494-500`):
`F_a += Za Zb (Ra-Rb)/|Ra-Rb|³`; correct: `-∂/∂x (ZZ/x)` at separation x>0 (1-D reduction). -/
theorem ionForce_isolated_ok (ZZ x : ℝ) (hx : 0 < x) :
    HasDerivAt (fun y => ZZ / y) (-(ZZ * x / x ^ 3)) x := by
  have := (hasDerivAt_inv hx.ne').const_mul ZZ
  convert this using 1
  · ext y; simp [div_eq_mul_inv]
  · field_simp

end SalmonAudit.GS
end
