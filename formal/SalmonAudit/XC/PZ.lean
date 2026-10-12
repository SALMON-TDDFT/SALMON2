import Mathlib

/-!
# XC audit: Perdew–Zunger LDA (builtin `pz`, `pzm`, spin-polarized `pz`)

Source: SALMON v2.3.0 `src/xc/builtin_pz.f90`, `builtin_pzm.f90`, `builtin_pz_sp.f90`.

Abstractions: the regularisation `trho+1d-10` / `+1d-20` is dropped (domain `ρ > 0`);
the PZ parameters (γ, β₁, β₂, A, B, C, D) are arbitrary reals (only positivity of the
denominator is assumed), so the same theorems cover `pz` (C=0.002, D=-0.0116) and `pzm`.
`rs ρ = (3/(4πρ))^(1/3)` is the Wigner–Seitz radius.
-/

namespace SalmonAudit.XC

open Real

/-- Wigner–Seitz radius as in `builtin_pz.f90:70`. -/
noncomputable def rs (ρ : ℝ) : ℝ := (3 / (4 * π * ρ)) ^ ((1 : ℝ) / 3)

lemma rs_eq {ρ : ℝ} (hρ : 0 < ρ) :
    rs ρ = (3 / (4 * π)) ^ ((1 : ℝ) / 3) * ρ ^ (-((1 : ℝ) / 3)) := by
  unfold rs
  have h1 : (0 : ℝ) ≤ 3 / (4 * π) := by positivity
  rw [show 3 / (4 * π * ρ) = 3 / (4 * π) * ρ⁻¹ by field_simp]
  rw [Real.mul_rpow h1 (inv_nonneg.2 hρ.le), Real.inv_rpow hρ.le, Real.rpow_neg hρ.le]

/-- `d rs / d ρ = - rs / (3 ρ)`. -/
lemma rs_hasDerivAt {ρ : ℝ} (hρ : 0 < ρ) :
    HasDerivAt rs (-(rs ρ) / (3 * ρ)) ρ := by
  set C : ℝ := (3 / (4 * π)) ^ ((1 : ℝ) / 3)
  have h := ((Real.hasDerivAt_rpow_const (p := -((1 : ℝ) / 3)) (Or.inl hρ.ne')).const_mul C)
  have hev : (fun x => C * x ^ (-((1 : ℝ) / 3))) =ᶠ[nhds ρ] rs := by
    filter_upwards [Ioi_mem_nhds hρ] with x hx
    exact (rs_eq hx).symm
  have h2 := h.congr_of_eventuallyEq hev.symm
  refine h2.congr_deriv ?_
  rw [rs_eq hρ, Real.rpow_sub hρ, Real.rpow_one]
  field_simp
  ring

/-- chain rule through `rs`. -/
lemma comp_rs {f : ℝ → ℝ} {f' ρ : ℝ} (hρ : 0 < ρ) (hf : HasDerivAt f f' (rs ρ)) :
    HasDerivAt (fun x => f (rs x)) (f' * (-(rs ρ) / (3 * ρ))) ρ :=
  hf.comp ρ (rs_hasDerivAt hρ)

lemma rs_pos {ρ : ℝ} (hρ : 0 < ρ) : 0 < rs ρ := by
  unfold rs; positivity

/-! ## XC-01  Slater exchange constant -/

/-- `const = 0.75*(3/(2π))^(2/3)` (`builtin_pz.f90:62`). -/
noncomputable def pzConstCode : ℝ := 0.75 * (3 / (2 * π)) ^ ((2 : ℝ) / 3)

/-- textbook: `ε_x = -(3/(4π)) (9π/4)^{1/3} / r_s`. -/
noncomputable def pzConstCorrect : ℝ := 3 / (4 * π) * (9 * π / 4) ^ ((1 : ℝ) / 3)

theorem pz_exchange_constant_ok : pzConstCode = pzConstCorrect := by
  have hπ := Real.pi_pos
  have ha : (0 : ℝ) ≤ pzConstCode := by unfold pzConstCode; positivity
  have hb : (0 : ℝ) ≤ pzConstCorrect := by unfold pzConstCorrect; positivity
  refine (pow_left_inj₀ ha hb (n := 3) (by norm_num)).1 ?_
  unfold pzConstCode pzConstCorrect
  have e1 : ((3 / (2 * π)) ^ ((2 : ℝ) / 3)) ^ (3 : ℕ) = (3 / (2 * π)) ^ (2 : ℕ) := by
    rw [← Real.rpow_natCast, ← Real.rpow_mul (by positivity)]
    norm_num
  have e2 : ((9 * π / 4) ^ ((1 : ℝ) / 3)) ^ (3 : ℕ) = 9 * π / 4 := by
    rw [← Real.rpow_natCast, ← Real.rpow_mul (by positivity)]
    norm_num
  rw [mul_pow, mul_pow, e1, e2]
  field_simp
  ring

/-! ## XC-02  unpolarized PZ: dε/dρ in the three pieces (`builtin_pz.f90:72,76,80`) -/

/-- exchange part `-K/rs`. -/
noncomputable def epsX (K r : ℝ) : ℝ := -K / r

/-- code `dexc_drho = exc/(3ρ)` (`builtin_pz.f90:72`). -/
noncomputable def dEpsXCode (K ρ : ℝ) : ℝ := epsX K (rs ρ) / (3 * ρ)

theorem pz_exchange_deriv_ok (K : ℝ) {ρ : ℝ} (hρ : 0 < ρ) :
    HasDerivAt (fun x => epsX K (rs x)) (dEpsXCode K ρ) ρ := by
  have hr := rs_pos hρ
  have hf : HasDerivAt (fun r => -K * r⁻¹) (-K * (-(rs ρ ^ 2)⁻¹)) (rs ρ) :=
    (hasDerivAt_inv hr.ne').const_mul (-K)
  have h := comp_rs hρ hf
  have hfun : (fun x => -K * (rs x)⁻¹) = fun x => epsX K (rs x) := by
    funext x; unfold epsX; ring
  rw [hfun] at h
  refine h.congr_deriv ?_
  unfold dEpsXCode epsX
  field_simp

/-- PZ correlation, `rs > 1` branch: `γ/(1+β₁√rs+β₂ rs)`. -/
noncomputable def ecHigh (g b1 b2 r : ℝ) : ℝ := g / (1 + b1 * √r + b2 * r)

/-- code derivative `γ(½β₁√rs+β₂rs)/(3ρ)/(1+β₁√rs+β₂rs)^2` (`builtin_pz.f90:76`). -/
noncomputable def dEcHighCode (g b1 b2 ρ : ℝ) : ℝ :=
  g * (0.5 * b1 * √(rs ρ) + b2 * rs ρ) / (3 * ρ) / (1 + b1 * √(rs ρ) + b2 * rs ρ) ^ 2

lemma ecHigh_hasDerivAt (g b1 b2 : ℝ) {r : ℝ} (hr : 0 < r) (hD : 1 + b1 * √r + b2 * r ≠ 0) :
    HasDerivAt (ecHigh g b1 b2)
      (g * (-(b1 * (1 / (2 * √r)) + b2) / (1 + b1 * √r + b2 * r) ^ 2)) r := by
  have hD' : HasDerivAt (fun x => 1 + b1 * √x + b2 * x) (b1 * (1 / (2 * √r)) + b2) r := by
    have := ((hasDerivAt_const r (1:ℝ)).add ((Real.hasDerivAt_sqrt hr.ne').const_mul b1)).add
      ((hasDerivAt_id r).const_mul b2)
    exact this.congr_deriv (by simp)
  have := (hD'.inv hD).const_mul g
  refine (this.congr_deriv (by ring)).congr_of_eventuallyEq ?_
  filter_upwards with x
  unfold ecHigh; simp [div_eq_mul_inv]

theorem pz_corr_high_deriv_ok (g b1 b2 : ℝ) (hb1 : 0 ≤ b1) (hb2 : 0 ≤ b2) {ρ : ℝ} (hρ : 0 < ρ) :
    HasDerivAt (fun x => ecHigh g b1 b2 (rs x)) (dEcHighCode g b1 b2 ρ) ρ := by
  have hr := rs_pos hρ
  have hD : 0 < 1 + b1 * √(rs ρ) + b2 * rs ρ := by positivity
  have h := comp_rs hρ (ecHigh_hasDerivAt g b1 b2 hr hD.ne')
  refine h.congr_deriv ?_
  unfold dEcHighCode
  obtain ⟨s, hs, hrs⟩ : ∃ s, 0 < s ∧ rs ρ = s ^ 2 :=
    ⟨√(rs ρ), Real.sqrt_pos.2 hr, (Real.sq_sqrt hr.le).symm⟩
  rw [hrs, Real.sqrt_sq hs.le] at *
  field_simp
  ring

/-- PZ correlation, `rs ≤ 1` branch: `A ln rs + B + C rs ln rs + D rs`. -/
noncomputable def ecLow (A B C D r : ℝ) : ℝ := A * Real.log r + B + C * r * Real.log r + D * r

/-- code derivative `-rs/(3ρ)*(A/rs + C(ln rs+1) + D)` (`builtin_pz.f90:80`). -/
noncomputable def dEcLowCode (A C D ρ : ℝ) : ℝ :=
  -(rs ρ) / (3 * ρ) * (A / rs ρ + C * (Real.log (rs ρ) + 1) + D)

lemma ecLow_hasDerivAt (A B C D : ℝ) {r : ℝ} (hr : 0 < r) :
    HasDerivAt (ecLow A B C D) (A / r + C * (Real.log r + 1) + D) r := by
  have hl := Real.hasDerivAt_log hr.ne'
  have := ((((hl.const_mul A).add (hasDerivAt_const r B)).add
    (((hasDerivAt_id r).const_mul C).mul hl)).add ((hasDerivAt_id r).const_mul D))
  refine (this.congr_deriv ?_).congr_of_eventuallyEq ?_
  · simp; field_simp
  · filter_upwards with x; unfold ecLow; simp

theorem pz_corr_low_deriv_ok (A B C D : ℝ) {ρ : ℝ} (hρ : 0 < ρ) :
    HasDerivAt (fun x => ecLow A B C D (rs x)) (dEcLowCode A C D ρ) ρ := by
  have h := comp_rs hρ (ecLow_hasDerivAt A B C D (rs_pos hρ))
  refine h.congr_deriv ?_
  unfold dEcLowCode; ring

/-- code `Vexc = e_xc + trho*de_xc_drho` (`builtin_pz.f90:45`) is `d(ρ ε)/dρ`. -/
theorem pz_potential_assembly_ok {ε : ℝ → ℝ} {ε' ρ : ℝ} (h : HasDerivAt ε ε' ρ) :
    HasDerivAt (fun x => x * ε x) (ε ρ + ρ * ε') ρ := by
  exact ((hasDerivAt_id ρ).mul h).congr_deriv (by simp)

/-! ## XC-03  spin-polarized PZ (`builtin_pz_sp.f90`) -/

/-- code closed form `vcu` for `rs>1` (`builtin_pz_sp.f90:84`). -/
noncomputable def vcHighCode (g b1 b2 r : ℝ) : ℝ :=
  g * (1 + 7 / 6 * b1 * √r + 4 / 3 * b2 * r) / (1 + b1 * √r + b2 * r) ^ 2

/-- correct: `ε_c - (rs/3) dε_c/drs`. -/
theorem pz_sp_vc_high_ok (g b1 b2 : ℝ) {r : ℝ} (hr : 0 < r) (hD : 0 < 1 + b1 * √r + b2 * r) :
    vcHighCode g b1 b2 r =
      ecHigh g b1 b2 r - r / 3 * (g * (-(b1 * (1 / (2 * √r)) + b2) / (1 + b1 * √r + b2 * r) ^ 2)) := by
  unfold vcHighCode ecHigh
  obtain ⟨s, hs, rfl⟩ : ∃ s, 0 < s ∧ r = s ^ 2 :=
    ⟨√r, Real.sqrt_pos.2 hr, (Real.sq_sqrt hr.le).symm⟩
  rw [Real.sqrt_sq hs.le] at *
  field_simp
  ring

/-- code closed form `Vcu` for `rs≤1` (`builtin_pz_sp.f90:89`). -/
noncomputable def vcLowCode (A B C D r : ℝ) : ℝ :=
  A * Real.log r + (B - A / 3) + 2 / 3 * C * r * Real.log r + (2 * D - C) / 3 * r

theorem pz_sp_vc_low_ok (A B C D : ℝ) {r : ℝ} (hr : 0 < r) :
    vcLowCode A B C D r = ecLow A B C D r - r / 3 * (A / r + C * (Real.log r + 1) + D) := by
  unfold vcLowCode ecLow
  field_simp
  ring

/-- ρ ∂ζ/∂ρ_↑ = 1 - ζ and ρ ∂ζ/∂ρ_↓ = -1 - ζ, i.e. `(sgnsigma - zeta)` (`builtin_pz_sp.f90:78,95`). -/
theorem zeta_deriv_up_ok {a b : ℝ} (hab : 0 < a + b) :
    HasDerivAt (fun x => (x - b) / (x + b)) ((1 - (a - b) / (a + b)) / (a + b)) a := by
  have := ((hasDerivAt_id a).sub (hasDerivAt_const a b)).div
    ((hasDerivAt_id a).add (hasDerivAt_const a b)) hab.ne'
  refine this.congr_deriv ?_
  simp; field_simp; ring

theorem zeta_deriv_dn_ok {a b : ℝ} (hab : 0 < a + b) :
    HasDerivAt (fun x => (a - x) / (a + x)) ((-1 - (a - b) / (a + b)) / (a + b)) b := by
  have := ((hasDerivAt_const b a).sub (hasDerivAt_id b)).div
    ((hasDerivAt_const b a).add (hasDerivAt_id b)) hab.ne'
  refine this.congr_deriv ?_
  simp; field_simp; ring

/-- `f(ζ)` spin interpolation and its derivative `dsf` (`builtin_pz_sp.f90:69,71`), |ζ|<1. -/
noncomputable def fz (z : ℝ) : ℝ :=
  ((1 + z) ^ ((4 : ℝ) / 3) + (1 - z) ^ ((4 : ℝ) / 3) - 2) / (2 ^ ((4 : ℝ) / 3) - 2)

noncomputable def dfzCode (z : ℝ) : ℝ :=
  4 / 3 * ((1 + z) ^ ((1 : ℝ) / 3) - (1 - z) ^ ((1 : ℝ) / 3)) / (2 ^ ((4 : ℝ) / 3) - 2)

theorem fz_deriv_ok {z : ℝ} (h1 : -1 < z) (h2 : z < 1) : HasDerivAt fz (dfzCode z) z := by
  have ha : HasDerivAt (fun x => 1 + x) 1 z := by simpa using (hasDerivAt_id z).const_add 1
  have hb : HasDerivAt (fun x => 1 - x) (-1) z := by simpa using (hasDerivAt_id z).const_sub 1
  have hA := ha.rpow_const (p := (4 : ℝ) / 3) (Or.inl (by linarith))
  have hB := hb.rpow_const (p := (4 : ℝ) / 3) (Or.inl (by linarith))
  have := ((hA.add hB).sub (hasDerivAt_const z (2:ℝ))).div_const (2 ^ ((4 : ℝ) / 3) - 2)
  refine this.congr_deriv ?_
  unfold dfzCode
  norm_num
  ring

/-- ρ-scaling: with `ε(ρ,ζ) = E(rs(ρ),ζ)`, `v_σ = ε - (rs/3) ∂E/∂rs + (sgn_σ - ζ) ∂E/∂ζ`.
Code exchange spin potential (`builtin_pz_sp.f90:75-78`) with `k = 2^{1/3}` arbitrary. -/
noncomputable def vxSpinCode (cx k r sf dsf zeta sgn : ℝ) : ℝ :=
  (-cx / r + (-k + 1) * cx / r * sf) - 1 / 3 * cx / r + 1 / 3 * (-k + 1) * cx / r * sf
    + (-k + 1) * cx / r * (sgn - zeta) * dsf

/-- correct: `E_x = -(cx/rs)(1+(k-1)f(ζ))`, `∂E/∂rs = -E/rs`, `∂E/∂ζ = -(cx/rs)(k-1)f'`. -/
noncomputable def vxSpinCorrect (cx k r sf dsf zeta sgn : ℝ) : ℝ :=
  let E := -(cx / r) * (1 + (k - 1) * sf)
  E - r / 3 * (-E / r) + (sgn - zeta) * (-(cx / r) * (k - 1) * dsf)

theorem pz_sp_exchange_potential_ok (cx k r sf dsf zeta sgn : ℝ) (hr : r ≠ 0) :
    vxSpinCode cx k r sf dsf zeta sgn = vxSpinCorrect cx k r sf dsf zeta sgn := by
  unfold vxSpinCode vxSpinCorrect
  field_simp
  ring

/-- correlation spin assembly (`builtin_pz_sp.f90:93-96`):
`v_σ = vcu + f(vcp-vcu) + (ecp-ecu)(sgn-ζ)f'` equals `ε - rs/3 ∂ε/∂rs + (sgn-ζ)∂ε/∂ζ`
when `vcu = ecu - rs/3 ecu'`, `vcp = ecp - rs/3 ecp'` (proved above). -/
theorem pz_sp_corr_assembly_ok (ecu ecp decu decp r f df zeta sgn : ℝ) :
    (ecu - r / 3 * decu) + f * ((ecp - r / 3 * decp) - (ecu - r / 3 * decu))
      + (ecp - ecu) * (sgn - zeta) * df
    = (ecu + f * (ecp - ecu)) - r / 3 * (decu + f * (decp - decu))
      + (sgn - zeta) * (df * (ecp - ecu)) := by ring

end SalmonAudit.XC
