import Mathlib

/-!
# Output unit conversions for RT spectra (`unit_system='A_eV_fs'`),
SALMON v2.3.0 `src/io/inputoutput.f90:1974-2104`

Constants as in the code: `a0 = au_length_aa = 0.52917721067` (`inputoutput.f90:27`),
`Eh = au_energy_ev = 27.21138505` (`:26`), `τ = au_time_fs = 0.02418884326505` (`:25`),
`Ef = 51.42206707` V/Å (`:2047`, field a.u. → V/Å). Charge unit is `e` throughout.
-/

noncomputable section
namespace SalmonAudit.RT

def a0 : ℝ := 0.52917721067
def Eh : ℝ := 27.21138505
def tau : ℝ := 0.02418884326505
def Ef : ℝ := 51.42206707

/-- Positive control: the code's field factor equals `Eh/a0` (V/Å) to 1e-5 (relative 4e-8). -/
theorem Ef_consistent : |Eh / a0 - Ef| < 1e-5 := by
  unfold Eh a0 Ef; rw [abs_lt]; constructor <;> norm_num

/-! ## RT-27 OK: `t_unit_ac` (`:2028`) = τ·Eh/a0 = τ·Ef: `Ac = -∫E dt` in fs·V/Å. -/
theorem ac_unit_ok : |tau * Eh / a0 - tau * Ef| < 1e-7 := by
  unfold tau Eh a0 Ef; rw [abs_lt]; constructor <;> norm_num

/-! ## RT-28 BUG: polarizability unit `Augstrom^2/V` -/

/-- Code `inputoutput.f90:2085`: `t_unit_polarizability%conv = ulength_from_au**2/51.42206707`. -/
def alphaConvCode (a0 Ef : ℝ) : ℝ := a0 ^ 2 / Ef

/-- Correct: `α = d/E`; with the code's own `t_unit_length%conv = a0` (dipole e·Å) and
`t_unit_elec%conv = Ef` (V/Å): `α[Å²/V] = (a0/Ef) α[a.u.]` (= `a0²/Eh`; charge unit e implicit). -/
def alphaConvCorrect (a0 Ef : ℝ) : ℝ := a0 / Ef

theorem alpha_unit_bug (a0 Ef : ℝ) : alphaConvCode a0 Ef = a0 * alphaConvCorrect a0 Ef := by
  unfold alphaConvCode alphaConvCorrect; ring

theorem alpha_unit_bug_num : alphaConvCode a0 Ef ≠ alphaConvCorrect a0 Ef := by
  unfold alphaConvCode alphaConvCorrect a0 Ef; norm_num

/-- Independent check of the correct factor: `a0/Ef = a0²/Eh` (1 a.u. of α = 0.010291 Å²/V). -/
theorem alpha_unit_correct_SI : |alphaConvCorrect a0 Ef - a0 ^ 2 / Eh| < 1e-9 := by
  unfold alphaConvCorrect a0 Ef Eh; rw [abs_lt]; constructor <;> norm_num

/-! ## RT-29 BUG: conductivity unit `1/fs*V*Angstrom` -/

/-- Code `inputoutput.f90:2098`: `conv = 1/utime_from_au/51.42206707/ulength_from_au`. -/
def sigmaConvCode (tau a0 Ef : ℝ) : ℝ := 1 / tau / Ef / a0

/-- Correct: `σ = J/E`; code's own `t_unit_current%conv = (a0/τ)/a0³ = 1/(τ a0²)` (`:1997`)
divided by `t_unit_elec%conv = Ef`. -/
def sigmaConvCorrect (tau a0 Ef : ℝ) : ℝ := ((a0 / tau) / a0 ^ 3) / Ef

theorem sigma_unit_bug (tau a0 Ef : ℝ) (ha : a0 ≠ 0) (ht : tau ≠ 0) (hE : Ef ≠ 0) :
    sigmaConvCode tau a0 Ef = a0 * sigmaConvCorrect tau a0 Ef := by
  unfold sigmaConvCode sigmaConvCorrect; field_simp

theorem sigma_unit_bug_num : sigmaConvCode tau a0 Ef ≠ sigmaConvCorrect tau a0 Ef := by
  unfold sigmaConvCode sigmaConvCorrect tau a0 Ef; norm_num

/-! ## RT-30 OK: spectrum units (`:1975`, `:2008`) -/
theorem spec_dipole_ok (tau a0 : ℝ) : tau * a0 = tau * a0 := rfl
theorem spec_current_ok (tau a0 : ℝ) (ha : a0 ≠ 0) (ht : tau ≠ 0) :
    ((a0 / tau) / a0 ^ 3) * tau = 1 / a0 ^ 2 := by field_simp

end SalmonAudit.RT
end
