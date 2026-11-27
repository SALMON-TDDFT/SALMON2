import Mathlib

/-!
# SBE velocity-gauge `norder_correction` (SALMON v2.3.0, `src/ssbe/bloch_solver_ssbe.f90:105-323`)

Abstraction: a two-level k-point (band 1 occupied with density `r`, band 2 empty), one Cartesian
direction, real symmetric momentum matrix `p11 = a, p22 = c, p12 = p21 = b`,
`Δ = delta_omega(2,1) = ε₂ - ε₁ ≠ 0`, vector potential `A`.  In the code `nb = 1`,
`ib = jb = kb = 2` is the only term surviving the `nb /= ib` guards.

The perturbative energy of band 1 under `V = A p` (Rayleigh–Schrödinger) is
`E2 = -A²b²/Δ`, `E3 = A³ b²(c-a)/Δ²`, `E4 = A⁴ (b⁴ - b²(c-a)²)/Δ³`
(checked against the series of the exact 2×2 eigenvalue).

Physical requirement (adiabatic theorem + `Σ_k ε_n(k+A)` independent of `A` for filled bands):
the static current of an insulator vanishes order by order, so the correction for the
truncated basis must equal `-r ∂_A E^{(m+1)}_inside` at every order `m`.
Order 1 in the code does this (theorem `sbe04_order1_ok`); orders 2 and 3 have the opposite sign.
-/

namespace SalmonAudit.SBE.Norder

noncomputable section

/-- SBE-04 code, order 1 (`bloch_solver_ssbe.f90:126-128`):
`2·occup(nb)·Re(pni_Ac·pin)/delta_omega(ib,nb)` with `pni_Ac = bA`, `pin = b`. -/
def code1 (f _a b _c Δ A : ℝ) : ℝ := 2 * f * ((b * A) * b) / Δ

/-- `∂_A E2` for `E2 = -A²b²/Δ`. -/
def dE2 (b Δ A : ℝ) : ℝ := -2 * A * b ^ 2 / Δ

theorem sbe04_order1_ok (f a b c Δ A : ℝ) : code1 f a b c Δ A = -f * dE2 b Δ A := by
  unfold code1 dE2; ring

/-- SBE-05 code, order 2 (`bloch_solver_ssbe.f90:165-170` plus `176-180`), two-level transcription:
`r/Δ/Δ·(pjn·pij_Ac·pni_Ac + pjn_Ac·(pij·pni_Ac + pni·pij_Ac))
  - r·(2·pnn_Ac·Re(pni·pin_Ac) + pnn·|pni_Ac|²)/Δ²`. -/
def code2 (r a b c Δ A : ℝ) : ℝ :=
  r / Δ / Δ * (b * (c * A) * (b * A) + (b * A) * (c * (b * A) + b * (c * A)))
    - r * (2 * (a * A) * (b * (b * A)) + a * (b * A) ^ 2) / Δ ^ 2

/-- `∂_A E3`, `E3 = A³ b²(c-a)/Δ²`. -/
def dE3 (a b c Δ A : ℝ) : ℝ := 3 * A ^ 2 * (b ^ 2 * (c - a) / Δ ^ 2)

theorem dE3_hasDerivAt (a b c Δ A : ℝ) :
    HasDerivAt (fun x : ℝ => x ^ 3 * (b ^ 2 * (c - a) / Δ ^ 2)) (dE3 a b c Δ A) A := by
  have h := (hasDerivAt_pow 3 A).mul_const (b ^ 2 * (c - a) / Δ ^ 2)
  unfold dE3
  convert h using 1

/-- What the code computes at order 2 is `+ r ∂_A E3` (sign opposite to order 1). -/
theorem sbe05_code2_eq_plus_dE3 (r a b c Δ A : ℝ) (hΔ : Δ ≠ 0) :
    code2 r a b c Δ A = r * dE3 a b c Δ A := by
  unfold code2 dE3; field_simp; ring

/-- Correct order-2 correction: `-r ∂_A E3` (cancels the in-basis adiabatic current). -/
def correct2 (r a b c Δ A : ℝ) : ℝ := -r * dE3 a b c Δ A

theorem sbe05_order2_bug : code2 2 0 1 1 1 1 ≠ correct2 2 0 1 1 1 1 := by
  unfold code2 correct2 dE3; norm_num

/-- SBE-06 code, order 3 (`bloch_solver_ssbe.f90:252-260`, `267-284`, `311-319`), two-level
transcription with `pji = pij = pkj = pjk = c` (all `p22`), `pnk = pin = pni = pnj = pjn = pkn = b`,
`pnn = a`, `δ = Δ`. -/
def code3 (r a b c Δ A : ℝ) : ℝ :=
  -- triple loop (252-260)
  (- r / Δ / Δ / Δ *
      ((c * A) * ((b * A) * (b * (c * A) + c * (b * A)) + b * (b * A) * (c * A))
        + c * (b * A) * (c * A) * (b * A)))
  -- double loop (267-284)
  + r * (Δ + Δ) / 2 / Δ ^ 2 / Δ ^ 2 *
      (c * (a * A) * (b * A) * (b * A)
        + b * (b * A) * ((b * A) * (b * A) + (a * A) * (c * A))
        + (c * A) * ((a * A) * (b * (b * A) + b * (b * A)) + a * (b * A) * (b * A))
        + (b * A) * (b * ((b * A) * (b * A) + (a * A) * (c * A))
            + (b * A) * (c * (a * A) + b * (b * A) + b * (b * A) + a * (c * A))))
  -- single loop (311-319)
  - r * (a * A) * ((b * (a * A) * (b * A) + (b * A) * (b * (a * A) + 2 * a * (b * A))) / Δ ^ 3)

/-- `∂_A E4`, `E4 = A⁴ (b⁴ - b²(c-a)²)/Δ³`. -/
def dE4 (a b c Δ A : ℝ) : ℝ := 4 * A ^ 3 * ((b ^ 4 - b ^ 2 * (c - a) ^ 2) / Δ ^ 3)

theorem sbe06_code3_eq_plus_dE4 (r a b c Δ A : ℝ) (hΔ : Δ ≠ 0) :
    code3 r a b c Δ A = r * dE4 a b c Δ A := by
  unfold code3 dE4; field_simp; ring

def correct3 (r a b c Δ A : ℝ) : ℝ := -r * dE4 a b c Δ A

theorem sbe06_order3_bug : code3 2 0 1 0 1 1 ≠ correct3 2 0 1 0 1 1 := by
  unfold code3 correct3 dE4; norm_num

end

end SalmonAudit.SBE.Norder
