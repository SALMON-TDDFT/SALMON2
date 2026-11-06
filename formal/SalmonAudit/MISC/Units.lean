import Mathlib

/-!
# Output unit conversions for `unit_system='A_eV_fs'` (src/io/inputoutput.f90:1936-2103)

`L = ulength_from_au = 0.52917721067` (Å per bohr), `t = utime_from_au = 0.02418884326505`
(fs per a.u.), `Ef = 51.42206707` (V/Å per a.u. field). All quantities are converted as
value_out = value_au · conv.
-/

namespace SalmonAudit.MISC.Units

noncomputable section

def L : ℝ := 0.52917721067
def t : ℝ := 0.02418884326505
def Ef : ℝ := 51.42206707

/-- Code conductivity factor `inputoutput.f90:2099`: `1/t/Ef/L`. -/
def codeSigma : ℝ := 1 / t / Ef / L
/-- Correct: `σ = j/E`; `j` converts with `1/(t L²)` (`t_unit_current`, `:1997`),
`E` with `Ef` ⇒ `1/(t L² Ef)` [e/(fs·V·Å)]. -/
def correctSigma : ℝ := (1 / (t * L ^ 2)) / Ef

/-- MISC-29 (BUG): printed σ is too small by the factor `L ≈ 0.529`. -/
theorem sigma_bug : codeSigma = L * correctSigma ∧ codeSigma ≠ correctSigma := by
  have hL : L ≠ 0 := by unfold L; norm_num
  have ht : t ≠ 0 := by unfold t; norm_num
  have hE : Ef ≠ 0 := by unfold Ef; norm_num
  constructor
  · unfold codeSigma correctSigma; field_simp
  · unfold codeSigma correctSigma L t Ef; norm_num

/-- Code polarizability factor `inputoutput.f90:2086`: `L²/Ef`. -/
def codeAlpha : ℝ := L ^ 2 / Ef
/-- Correct: `α = p/E`, `p` in e·Å converts with `L`, `E` with `Ef` ⇒ `L/Ef` [Å²/V]. -/
def correctAlpha : ℝ := L / Ef

/-- MISC-30 (BUG): printed α is too small by the factor `L`. -/
theorem alpha_bug : codeAlpha = L * correctAlpha ∧ codeAlpha ≠ correctAlpha := by
  constructor
  · unfold codeAlpha correctAlpha; ring
  · unfold codeAlpha correctAlpha L Ef; norm_num

/-- Current density `:1997`: `(L/t)/L³`; correct `e/(t_au a0²)` ⇒ `1/(t L²)`. MISC-31 (OK). -/
theorem current_ok (hL : L ≠ 0) : (L / t) / L ^ 3 = 1 / (t * L ^ 2) := by
  field_simp

/-- Vector potential `:2028`: `t·E/(L·q)` with `E` in eV: `A = E·t` in V·fs/Å. MISC-31b (OK):
`t·(eV)/L = t·Ef` because `Ef = (eV per Ha)/L`. Checked to 1e-7 relative. -/
theorem ac_consistent :
    |t * 27.21138505 / L - t * Ef| / (t * Ef) < 1e-7 := by
  unfold t L Ef; rw [div_lt_iff₀ (by norm_num), abs_lt]; constructor <;> norm_num

/-- MISC-32 (NOTE): `Ef = 51.42206707` vs `E_h/(e a0) = 27.21138505/0.52917721067`
from the same file's constants: relative mismatch < 1e-7. -/
theorem efield_note : |27.21138505 / L - Ef| / Ef < 1e-7 := by
  unfold L Ef; rw [div_lt_iff₀ (by norm_num), abs_lt]; constructor <;> norm_num

/-- Force output `write.f90:1547`: `F·(eV/Ha)/(Å/bohr)` = eV/Å. MISC-33 (OK). -/
theorem force_ok (F : ℝ) : F * 27.21138505 / L = F * (27.21138505 / L) := by ring

end

end SalmonAudit.MISC.Units
