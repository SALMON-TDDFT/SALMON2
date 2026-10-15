import Mathlib

/-!
# Linear-response post-processing (SALMON v2.3.0, `src/io/write.f90:1001-1130`, projection
`src/io/write.f90:2122-2387`)

Abstractions: the windowed discrete Fourier sums are treated as given complex numbers
(`Dω = dt Σ e^{iωt} w(t) ΔD(t)` etc.); the theorems check how the code turns them into
α, σ, ε. Kick convention established in `Fields.lean` (RT-12): `E(t) = -κ δ(t)`, hence the
Fourier transform of the applied field is `E(ω) = -κ`.
-/

noncomputable section
namespace SalmonAudit.RT

open Complex

/-! ## RT-20 isolated polarizability -/

/-- Code `write.f90:1047`: `zalpha = Σ ... / (-e_impulse) * dt`. -/
def alphaCode (Dω : ℂ) (κ : ℝ) : ℂ := Dω / (-(κ : ℂ))

/-- Correct: `α(ω) = d(ω)/E(ω)` with `E(ω) = -κ` (electron charge −1: `Δp = -∫E dt = κ`). -/
def alphaCorrect (Dω : ℂ) (κ : ℝ) : ℂ :=
  let Eω : ℂ := -(κ : ℂ)
  Dω / Eω

theorem alpha_ok (Dω : ℂ) (κ : ℝ) : alphaCode Dω κ = alphaCorrect Dω κ := rfl

/-- Sign sanity: free electrons, kick `+κ` moves electrons by `+κ t`, electronic dipole
`d(t) = -Nκt`, `d(ω) = Nκ/ω²` (damped FT of `t`: `-1/ω²`); code gives `α = -N/ω²` (Drude). -/
theorem alpha_free_electron (N κ ω : ℝ) (hκ : κ ≠ 0) (hω : ω ≠ 0) :
    alphaCode ((N * κ / ω ^ 2 : ℝ) : ℂ) κ = ((-N / ω ^ 2 : ℝ) : ℂ) := by
  unfold alphaCode
  have : (κ : ℂ) ≠ 0 := by exact_mod_cast hκ
  have : (ω : ℂ) ≠ 0 := by exact_mod_cast hω
  push_cast
  field_simp

/-! ## RT-21 strength function -/

/-- Code `write.f90:1048`: `sf = 2ω/π Im α`. Correct: for `Im α = (π f/(2ω0)) δ(ω-ω0)` (one
oscillator of strength `f`) the weight of `sf` at `ω0` is `f` ⇒ `∫ sf dω = Σ f = N` (TRK). -/
theorem sf_weight_ok (f ω0 : ℝ) (hω : ω0 ≠ 0) :
    (2 * ω0 / Real.pi) * (Real.pi * f / (2 * ω0)) = f := by
  have := Real.pi_ne_zero
  field_simp

/-! ## RT-22 window function -/

/-- Code `write.f90:1044,1119,1131`: `w(x) = 1 - 3x² + 2x³`, `x = t/T`. -/
def win (x : ℝ) : ℝ := 1 - 3 * x ^ 2 + 2 * x ^ 3
def winDeriv (x : ℝ) : ℝ := -6 * x + 6 * x ^ 2

theorem win_ok : win 0 = 1 ∧ win 1 = 0 ∧ winDeriv 0 = 0 ∧ winDeriv 1 = 0 ∧
    ∀ x, HasDerivAt win (winDeriv x) x := by
  refine ⟨by norm_num [win], by norm_num [win], by norm_num [winDeriv], by norm_num [winDeriv], ?_⟩
  intro x
  have := (((hasDerivAt_pow 2 x).const_mul (-3)).add ((hasDerivAt_pow 3 x).const_mul 2)).const_add 1
  unfold win winDeriv
  convert this using 1
  · funext y; simp only [Pi.add_apply]; ring
  · push_cast; ring

/-! ## RT-23 periodic, transverse (`trans_longi='tr'`): σ and ε -/

/-- Code `write.f90:1136`: `zsigma = Σ e^{iωt} w curr dt / e_impulse`, `curr` = electron matter current. -/
def sigmaTrCode (Cω : ℂ) (κ : ℝ) : ℂ := Cω / (κ : ℂ)

/-- Correct: `σ = J_el(ω)/E(ω)`, `J_el = -curr`, `E(ω) = -κ`. -/
def sigmaTrCorrect (Cω : ℂ) (κ : ℝ) : ℂ := (-Cω) / (-(κ : ℂ))

theorem sigma_tr_ok (Cω : ℂ) (κ : ℝ) : sigmaTrCode Cω κ = sigmaTrCorrect Cω κ := by
  unfold sigmaTrCode sigmaTrCorrect; rw [neg_div_neg_eq]

/-- Drude check: after the step `Ac = κ θ(t)`, `curr = nκ`; `∫_0^∞ e^{iωt} dt = i/ω` ⇒
`σ = i n/ω` and `ε = 1 + 4πiσ/ω = 1 - 4πn/ω²` (`write.f90:1138`). -/
theorem drude_ok (n κ ω : ℝ) (hκ : κ ≠ 0) (hω : ω ≠ 0) :
    let σ := sigmaTrCode ((n * κ : ℝ) * (I / (ω : ℂ))) κ
    1 + 4 * (Real.pi : ℂ) * I * σ / (ω : ℂ) = ((1 - 4 * Real.pi * n / ω ^ 2 : ℝ) : ℂ) := by
  simp only [sigmaTrCode]
  have : (κ : ℂ) ≠ 0 := by exact_mod_cast hκ
  have : (ω : ℂ) ≠ 0 := by exact_mod_cast hω
  push_cast
  field_simp
  ring_nf
  rw [I_sq]
  ring

/-! ## RT-24 periodic, longitudinal (`trans_longi='lo'`): ε from E_tot -/

/-- Code `write.f90:1109,1136,1140`: `rtdata = E_tot(1:nt)` (the δ-kick at `n=0` is not in
the sum since `E_ext(n≥1)=0`), `zsigma = Σ E_tot ... / e_impulse`, `ε = 1/(1 - zsigma)`.
Here `Eind` = FT of the induced field. -/
def epsLoCode (Eind : ℂ) (κ : ℝ) : ℂ := 1 / (1 - Eind / (κ : ℂ))

/-- Correct: `1/ε = E_tot(ω)/E_ext(ω)` with `E_ext(ω) = -κ`, `E_tot = E_ext + E_ind`. -/
def epsLoCorrect (Eind : ℂ) (κ : ℝ) : ℂ := (-(κ : ℂ)) / (-(κ : ℂ) + Eind)

theorem eps_lo_ok (Eind : ℂ) (κ : ℝ) (hκ : κ ≠ 0) (hd : -(κ : ℂ) + Eind ≠ 0) :
    epsLoCode Eind κ = epsLoCorrect Eind κ := by
  unfold epsLoCode epsLoCorrect
  have hk : (κ : ℂ) ≠ 0 := by exact_mod_cast hκ
  have h2 : (1 : ℂ) - Eind / (κ : ℂ) ≠ 0 := by
    intro h
    apply hd
    field_simp at h
    linear_combination -h
  have h3 : (κ : ℂ) - Eind ≠ 0 := by
    intro h; apply hd; linear_combination -h
  rw [div_eq_div_iff h2 hd]
  field_simp
  ring

/-- `σ = (ε-1)/(4πi/ω)` (`write.f90:1141`) is the inverse of `ε = 1 + 4πiσ/ω`. -/
theorem sigma_from_eps_ok (σ : ℂ) (ω : ℝ) (hω : ω ≠ 0) :
    ((1 + 4 * (Real.pi : ℂ) * I * σ / (ω : ℂ)) - 1) / (4 * (Real.pi : ℂ) * I / (ω : ℂ)) = σ := by
  have : (ω : ℂ) ≠ 0 := by exact_mod_cast hω
  have hpi : (Real.pi : ℂ) ≠ 0 := by exact_mod_cast Real.pi_ne_zero
  field_simp
  ring

/-! ## RT-25 isolated dipole -/

/-- Code `dip.f90:65`: `Dp_e = -Σ (ix - s) h ρ Hvol` with `s = 0.5` for even `N`
(`dip.f90:186`) — identical to `lg%coordinate(ix) = (ix-0.5) h` (`initialization.f90:1170`);
correct electronic dipole `-∫ r ρ`. -/
theorem dip_ok (ix h ρ Hvol : ℝ) :
    -((ix - 1 / 2) * ρ) * h * Hvol = -(((ix - 1 / 2) * h) * ρ * Hvol) := by ring

/-! ## RT-26 excitation numbers from projection -/

/-- Code `write.f90:2259-2260`: `nee += (w - r0)/w · coef`, `neh -= r0/w · coef` with
`r0 ∈ {0, w}` the GS occupation, `coef = Σ_m rocc_m wk |⟨w_n|u_m⟩|²`. For an occupied GS orbital
the `nee` weight is 0 and `neh` weight 1; for an empty orbital `nee` weight is 1. -/
theorem proj_weights_ok (w : ℝ) (hw : w ≠ 0) :
    ((w - w) / w = 0 ∧ w / w = 1) ∧ ((w - 0) / w = 1 ∧ 0 / w = 0) := by
  refine ⟨⟨by simp, div_self hw⟩, ⟨by simp [div_self hw], by simp⟩⟩

end SalmonAudit.RT
end
