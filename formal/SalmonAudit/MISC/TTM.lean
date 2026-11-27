import Mathlib

/-!
# Two-temperature model (src/ttm/ttm.f90, coupled in src/maxwell/fdtd_eh.f90)

Unit model: a physical unit is a vector of exponents of (energy, length, time, temperature).
Converting a value from SI to atomic units multiplies by
`E_h^{-e} · a0^{-l} · t_au^{-t} · T_h^{-T}` (SI value of each atomic unit), i.e. the exponent
vector of the conversion factor is the negated unit vector. A code factor
`(1/E_h)·a0^3·T_h` is recorded as exponents `(-1, 3, 0, 1)`.
-/

namespace SalmonAudit.MISC.TTM

noncomputable section

/-- (energy, length, time, temperature) exponents. -/
abbrev Dim := ℤ × ℤ × ℤ × ℤ

/-- Exponents of the correct SI→au factor for a quantity of dimension `d`. -/
def correctFactor (d : Dim) : Dim := (-d.1, -d.2.1, -d.2.2.1, -d.2.2.2)

/-- `Cl` : J m⁻³ K⁻¹. Code `ttm.f90:116-118`: `/E_h · a0³ · T_h`. -/
def dimCl : Dim := (1, -3, 0, -1)
def codeFactorCl : Dim := (-1, 3, 0, 1)
/-- `kappa_e0` : J m⁻¹ s⁻¹ K⁻¹. Code `:119-122`: `/E_h · a0 · T_h · t_au`. -/
def dimKappa : Dim := (1, -1, -1, -1)
def codeFactorKappa : Dim := (-1, 1, 1, 1)
/-- `Ce'` : J m⁻³ K⁻². Code `:123-124`: `/E_h · a0³ · T_h²`. -/
def dimCe : Dim := (1, -3, 0, -2)
def codeFactorCe : Dim := (-1, 3, 0, 2)
/-- `g` (electron–phonon coupling) : W m⁻³ K⁻¹ = J m⁻³ s⁻¹ K⁻¹ — forced by
`Cl dTl/dt = g (Te - Tl)` (`ttm.f90:260,281`). Code `:114-115`: `/E_h · a0³` only. -/
def dimG : Dim := (1, -3, -1, -1)
def codeFactorG : Dim := (-1, 3, 0, 0)

/-- MISC-19a (OK): Cl, kappa_e0, Ce' conversions are dimensionally right. -/
theorem conv_ok :
    codeFactorCl = correctFactor dimCl ∧ codeFactorKappa = correctFactor dimKappa ∧
    codeFactorCe = correctFactor dimCe := by decide

/-- MISC-19 (BUG): the `g` conversion misses `t_au · T_h`. -/
theorem convG_bug : codeFactorG ≠ correctFactor dimG := by decide

/-- Size of the missing factor with the code's own constants (`ttm.f90:52-54`). -/
theorem convG_missing_factor :
    (2.4188843265857e-17 : ℝ) * 3.1577502480407e5 ≠ 1 := by norm_num

/-- Direction and size: code g_au / correct g_au = `1/(t_au T_h)` ≈ 1.31×10¹¹ (too large). -/
theorem convG_too_large :
    1.30e11 < 1 / ((2.4188843265857e-17 : ℝ) * 3.1577502480407e5) ∧
    1 / ((2.4188843265857e-17 : ℝ) * 3.1577502480407e5) < 1.32e11 := by
  constructor
  · rw [lt_div_iff₀ (by norm_num)]; norm_num
  · rw [div_lt_iff₀ (by norm_num)]; norm_num


/-- Code source term: `fdtd_eh.f90:3318-3323` accumulates
`u_energy ← u_energy - divS·dt` (never reset; `allocate_poynting` zeroes once,
`fdtd_eh.f90:6451`) and passes it as `source` to `ttm_main`, where it is added to the
power balance `rhs_e` (`ttm.f90:279`). With `p n = -divS` at step n: -/
def codeSource (p : ℕ → ℝ) (dt : ℝ) (n : ℕ) : ℝ := ∑ i ∈ Finset.range (n + 1), p i * dt

/-- Correct `S(r,t)` in `Ce dTe/dt = ∇κ∇Te - g(Te-Tl) + S`: the instantaneous absorbed
power density. -/
def correctSource (p : ℕ → ℝ) (n : ℕ) : ℝ := p n

/-- MISC-18 (BUG): a pulse absorbed only at step 0 keeps heating electrons forever.
(Modelled only the accumulation; the missing `c/4π` in the Poynting vector is EM-18, and a
complete fix must also remove stored field / reversible material energy.) -/
theorem source_bug :
    codeSource (fun i => if i = 0 then 1 else 0) 1 5 ≠
      correctSource (fun i => if i = 0 then 1 else 0) 5 := by
  simp [codeSource, correctSource, Finset.sum_range_succ]

/-- Code electron update `ttm.f90:280,288`: `Te += rhs·dt/(Ce'·Te)`. -/
def codeDTe (rhs dt Cep Te : ℝ) : ℝ := rhs * (dt / Cep) / Te
/-- Correct explicit Euler with `Ce(Te) = Ce'·Te`: `ΔTe = dt·rhs/Ce(Te)`. -/
def correctDTe (rhs dt Cep Te : ℝ) : ℝ := dt * rhs / (Cep * Te)

/-- MISC-20 (OK). -/
theorem dTe_ok (rhs dt Cep Te : ℝ) : codeDTe rhs dt Cep Te = correctDTe rhs dt Cep Te := by
  unfold codeDTe correctDTe; simp only [div_eq_mul_inv, mul_inv]; ring

end

end SalmonAudit.MISC.TTM
