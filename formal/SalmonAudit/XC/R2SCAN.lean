import Mathlib

/-!
# XC audit: builtin r2SCAN (`src/xc/builtin_r2scan.f90`, driver `salmon_xc.f90:1453-1552`)

Reference: Furness, Kaplan, Ning, Perdew, Sun, JPCL 11, 8208 (2020) (+ erratum).
Abstractions: constants are arbitrary reals unless the statement is about the constant;
chain-rule assembly lines are checked as derivative statements of the individual pieces.
-/

set_option linter.unusedSimpArgs false

namespace SalmonAudit.XC

open Real

/-- unfold pointwise `Pi` operations produced by `HasDerivAt.add/mul/...`. -/
macro "pis" : tactic => `(tactic| try simp only [Pi.add_apply, Pi.sub_apply, Pi.mul_apply,
  Pi.div_apply, Pi.neg_apply, Pi.inv_apply, id, id_eq] at *)

/-! ## XC-22  ᾱ and its partial derivatives (`builtin_r2scan.f90:133-147`) -/

/-- `ᾱ(ρ) = (τ - g/(8ρ)) / (c ρ^{5/3} + η g/(8ρ))` at fixed `g = |∇n|²`, `τ`. -/
noncomputable def alphaBar (c η g τ ρ : ℝ) : ℝ :=
  (τ - 0.125 * g / ρ) / (c * ρ ^ ((5 : ℝ) / 3) + η * (0.125 * g / ρ))

/-- code `dalpha_drho` (`builtin_r2scan.f90:146`). -/
noncomputable def dAlphaDrhoCode (c η g τ ρ : ℝ) : ℝ :=
  let th := c * ρ ^ ((5 : ℝ) / 3)
  let tw := 0.125 * g / ρ
  let al := alphaBar c η g τ ρ
  (tw / ρ - al * (5 * th / (3 * ρ) - η * tw / ρ)) / (th + η * tw)

theorem alpha_drho_ok (c η g τ : ℝ) {ρ : ℝ} (hρ : 0 < ρ)
    (hden : c * ρ ^ ((5 : ℝ) / 3) + η * (0.125 * g / ρ) ≠ 0) :
    HasDerivAt (alphaBar c η g τ) (dAlphaDrhoCode c η g τ ρ) ρ := by
  have hinv := hasDerivAt_inv hρ.ne'
  have htw : HasDerivAt (fun x => 0.125 * g / x) (0.125 * g * (-(ρ ^ 2)⁻¹)) ρ := by
    have := hinv.const_mul (0.125 * g)
    exact this.congr_of_eventuallyEq (by filter_upwards with x; simp [div_eq_mul_inv])
  have hp := (Real.hasDerivAt_rpow_const (p := (5 : ℝ) / 3) (Or.inl hρ.ne')).const_mul c
  have hnum := (hasDerivAt_const ρ τ).sub htw
  have hd := hp.add (htw.const_mul η)
  have := hnum.div hd hden
  refine (this.congr_deriv ?_).congr_of_eventuallyEq ?_
  · unfold dAlphaDrhoCode alphaBar
    simp only
    pis
    rw [Real.rpow_sub hρ, Real.rpow_one]
    generalize ρ ^ ((5 : ℝ) / 3) = P at hden ⊢
    set D := c * P + η * (0.125 * g / ρ) with hD
    clear_value D
    field_simp
    ring
  · filter_upwards with x; rfl

/-- code `dalpha_dg2 = -(1+ηᾱ)·0.125/(ρ·den)` (`:147`). -/
theorem alpha_dg2_ok (c η τ ρ : ℝ) {g : ℝ} (hρ : 0 < ρ)
    (hden : c * ρ ^ ((5 : ℝ) / 3) + η * (0.125 * g / ρ) ≠ 0) :
    HasDerivAt (fun x => alphaBar c η x τ ρ)
      (-(1 + η * alphaBar c η g τ ρ) * 0.125 / (ρ * (c * ρ ^ ((5 : ℝ) / 3) + η * (0.125 * g / ρ)))) g := by
  have htw : HasDerivAt (fun x => 0.125 * x / ρ) (0.125 / ρ) g := by
    simpa [mul_div_assoc] using ((hasDerivAt_id g).const_mul 0.125).div_const ρ
  have := ((hasDerivAt_const g τ).sub htw).div
    ((hasDerivAt_const g (c * ρ ^ ((5 : ℝ) / 3))).add (htw.const_mul η)) hden
  refine (this.congr_deriv ?_).congr_of_eventuallyEq ?_
  · unfold alphaBar
    pis
    generalize ρ ^ ((5 : ℝ) / 3) = P at hden ⊢
    set D := c * P + η * (0.125 * g / ρ) with hD
    clear_value D
    field_simp
    ring
  · filter_upwards with x; rfl

/-! ## XC-23  switching function `f(ᾱ)` (`interpolate_in_alpha`, `:260-283`) -/

/-- Horner value and derivative as coded (`:270-276`), unrolled. -/
def hornerF (c : Fin 8 → ℝ) (a : ℝ) : ℝ :=
  ((((((c 7 * a + c 6) * a + c 5) * a + c 4) * a + c 3) * a + c 2) * a + c 1) * a + c 0
def hornerD (c : Fin 8 → ℝ) (a : ℝ) : ℝ :=
  (((((7 * c 7 * a + 6 * c 6) * a + 5 * c 5) * a + 4 * c 4) * a + 3 * c 3) * a + 2 * c 2) * a + 1 * c 1

theorem horner_ok (c : Fin 8 → ℝ) (a : ℝ) : HasDerivAt (hornerF c) (hornerD c a) a := by
  have P : ∀ n : ℕ, ∀ k : ℝ, HasDerivAt (fun x : ℝ => k * x ^ n) (k * (n * a ^ (n - 1))) a :=
    fun n k => (hasDerivAt_pow n a).const_mul k
  have := ((((((((P 7 (c 7)).add (P 6 (c 6))).add (P 5 (c 5))).add (P 4 (c 4))).add
    (P 3 (c 3))).add (P 2 (c 2))).add (P 1 (c 1))).add (hasDerivAt_const a (c 0)))
  refine (this.congr_deriv ?_).congr_of_eventuallyEq ?_
  · unfold hornerD; push_cast; ring
  · filter_upwards with x; unfold hornerF; simp only [Pi.add_apply]; ring

/-- `ᾱ < 0` branch: `f = exp(-c₁ᾱ/(1-ᾱ))`, code `dfda = -c₁ f/(1-ᾱ)²`. -/
theorem switch_low_ok (c1 : ℝ) {a : ℝ} (ha : a < 1) :
    HasDerivAt (fun x => Real.exp (-c1 * x / (1 - x)))
      (-c1 * Real.exp (-c1 * a / (1 - a)) / (1 - a) ^ 2) a := by
  have h1 : (1 - a) ≠ 0 := by linarith
  have hq := ((hasDerivAt_id a).const_mul (-c1)).div ((hasDerivAt_const a (1:ℝ)).sub (hasDerivAt_id a)) h1
  have := hq.exp
  refine (this.congr_deriv ?_).congr_of_eventuallyEq ?_
  · pis
    generalize Real.exp (-c1 * a / (1 - a)) = E
    field_simp; ring
  · filter_upwards with x; simp

/-- `ᾱ > 2.5` branch: `f = -d exp(c₂/(1-ᾱ))`, code `dfda = f c₂/(1-ᾱ)²`. -/
theorem switch_high_ok (c2 d : ℝ) {a : ℝ} (ha : 1 < a) :
    HasDerivAt (fun x => -d * Real.exp (c2 / (1 - x)))
      ((-d * Real.exp (c2 / (1 - a))) * c2 / (1 - a) ^ 2) a := by
  have h1 : (1 - a) ≠ 0 := by linarith
  have hq := (hasDerivAt_const a c2).div ((hasDerivAt_const a (1:ℝ)).sub (hasDerivAt_id a)) h1
  have := hq.exp.const_mul (-d)
  refine this.congr_deriv ?_
  pis
  generalize Real.exp (c2 / (1 - a)) = E
  field_simp; ring

/-- XC-24: `C₂ₓ = κ₀ Σ i cₓᵢ ≈ -0.162742` and `Δf_c2 = Σ i c_cᵢ ≈ -0.711402` (paper values). -/
noncomputable def fxOne : ℝ := (-0.667 : ℝ) + 2 * (-0.4445555) + 3 * (-0.663086601049)
  + 4 * 1.451297044490 + 5 * (-0.887998041597) + 6 * 0.234528941479 + 7 * (-0.023185843322)

theorem r2scan_C2x_ok : -0.1627425 < 0.174 * fxOne ∧ 0.174 * fxOne < -0.1627415 := by
  unfold fxOne; norm_num

noncomputable def fcOne : ℝ := (-0.64 : ℝ) + 2 * (-0.4352) + 3 * (-1.535685604549)
  + 4 * 3.061560252175 + 5 * (-1.915710236206) + 6 * 0.516884468372 + 7 * (-0.051848879792)

theorem r2scan_dfc2_ok : -0.7114025 < fcOne ∧ fcOne < -0.7114015 := by
  unfold fcOne; norm_num

/-! ## XC-25  exchange pieces -/

/-- `g_x(p) = 1 - exp(-a / p^{1/4})`, code `dgx = -0.25 a (1-g)/p^{1.25}` (`:155-156`). -/
theorem gx_deriv_ok (a : ℝ) {p : ℝ} (hp : 0 < p) :
    HasDerivAt (fun x => 1 - Real.exp (-a / x ^ (0.25 : ℝ)))
      (-0.25 * a * (1 - (1 - Real.exp (-a / p ^ (0.25 : ℝ)))) / p ^ (1.25 : ℝ)) p := by
  have hq := Real.hasDerivAt_rpow_const (p := (0.25 : ℝ)) (Or.inl hp.ne')
  have hP : 0 < p ^ (0.25 : ℝ) := Real.rpow_pos_of_pos hp _
  have hinv := (hq.inv hP.ne').const_mul (-a)
  have := ((hinv.exp).const_sub 1)
  refine (this.congr_deriv ?_).congr_of_eventuallyEq ?_
  · have e1 : p ^ ((0.25 : ℝ) - 1) = p ^ (0.25 : ℝ) / p := Real.rpow_sub_one hp.ne' _
    have e2 : p ^ (1.25 : ℝ) = p ^ (0.25 : ℝ) * p := by
      rw [show (1.25 : ℝ) = 0.25 + 1 by norm_num, Real.rpow_add_one hp.ne']
    rw [e1, e2]
    pis
    rw [div_eq_mul_inv (-a)]
    generalize p ^ (0.25 : ℝ) = P at hP ⊢
    generalize Real.exp (-a * P⁻¹) = E
    field_simp
    ring
  · filter_upwards with x; simp [div_eq_mul_inv]

/-- `x(p) = (C e^{-p²/d⁴} + μ) p`, code `dxs_ds2 = C e^{..}(1 - 2p²/d⁴) + μ` (`:160-161`). -/
theorem xs_deriv_ok (C μ d4 : ℝ) (hd : d4 ≠ 0) (p : ℝ) :
    HasDerivAt (fun x => (C * Real.exp (-x * x / d4) + μ) * x)
      (C * Real.exp (-p * p / d4) * (1 - 2 * p * p / d4) + μ) p := by
  have he := (((hasDerivAt_id p).neg.mul (hasDerivAt_id p)).div_const d4).exp
  have := ((he.const_mul C).add_const μ).mul (hasDerivAt_id p)
  refine (this.congr_deriv ?_).congr_of_eventuallyEq ?_
  · pis
    generalize Real.exp (-p * p / d4) = E
    field_simp; ring
  · filter_upwards with x; simp

/-- `h_x¹ = 1 + κ₁ x/(κ₁+x)`, code `dhx1_dxs = (κ₁/(κ₁+x))²` (`:163-164`). -/
theorem hx1_deriv_ok (k1 : ℝ) {x : ℝ} (h : k1 + x ≠ 0) :
    HasDerivAt (fun y => 1 + k1 * y / (k1 + y)) ((k1 / (k1 + x)) ^ 2) x := by
  have := (((hasDerivAt_id x).const_mul k1).div ((hasDerivAt_const x k1).add (hasDerivAt_id x)) h).const_add 1
  refine (this.congr_deriv ?_).congr_of_eventuallyEq ?_
  · pis; field_simp; ring
  · filter_upwards with y; simp

/-! ## XC-26  correlation pieces -/

/-- `g∞(p) = (1+4χp)^{-1/4}`, code `dgc0 = -χ g/(1+4χp)` (`:187-188`). -/
theorem gc0_deriv_ok (χ : ℝ) {p : ℝ} (hz : 0 < 1 + 4 * χ * p) :
    HasDerivAt (fun x => (1 + 4 * χ * x) ^ (-(0.25 : ℝ)))
      (-χ * (1 + 4 * χ * p) ^ (-(0.25 : ℝ)) / (1 + 4 * χ * p)) p := by
  have hl : HasDerivAt (fun x => 1 + 4 * χ * x) (4 * χ) p := by
    simpa using ((hasDerivAt_id p).const_mul (4 * χ)).const_add 1
  have := hl.rpow_const (p := -(0.25 : ℝ)) (Or.inl hz.ne')
  refine this.congr_deriv ?_
  rw [Real.rpow_sub_one hz.ne']
  field_simp
  ring

/-- `H = γ ln(1 + w(1 - (1+4y)^{-1/4}))`, code `dH/dy = γ w g/(S(1+4y))` (`:229-233`, `:190-192`
is the same shape with `b1c`). -/
theorem hc_darg_ok (γ w : ℝ) {y : ℝ} (hz : 0 < 1 + 4 * y)
    (hS : 0 < 1 + w * (1 - (1 + 4 * y) ^ (-(0.25 : ℝ)))) :
    HasDerivAt (fun x => γ * Real.log (1 + w * (1 - (1 + 4 * x) ^ (-(0.25 : ℝ)))))
      (γ * w * (1 + 4 * y) ^ (-(0.25 : ℝ)) / ((1 + w * (1 - (1 + 4 * y) ^ (-(0.25 : ℝ)))) * (1 + 4 * y))) y := by
  have hl : HasDerivAt (fun x => 1 + 4 * x) 4 y := by
    simpa using ((hasDerivAt_id y).const_mul 4).const_add 1
  have hg := hl.rpow_const (p := -(0.25 : ℝ)) (Or.inl hz.ne')
  have hS' := (((hasDerivAt_const y (1:ℝ)).sub hg).const_mul w).const_add 1
  have := (hS'.log hS.ne').const_mul γ
  refine (this.congr_deriv ?_).congr_of_eventuallyEq ?_
  · pis
    rw [Real.rpow_sub_one hz.ne']
    generalize (1 + 4 * y) ^ (-(0.25 : ℝ)) = G at hS ⊢
    field_simp
    ring
  · filter_upwards with x; simp

/-- `∂H/∂w = γ(1-g)/S` (`:232`, `:191`). -/
theorem hc_dw_ok (γ g : ℝ) {w : ℝ} (hS : 0 < 1 + w * (1 - g)) :
    HasDerivAt (fun x => γ * Real.log (1 + x * (1 - g))) (γ * (1 - g) / (1 + w * (1 - g))) w := by
  have h : HasDerivAt (fun x => 1 + x * (1 - g)) (1 - g) w := by
    simpa using ((hasDerivAt_id w).mul_const (1 - g)).const_add 1
  exact ((h.log hS.ne').const_mul γ).congr_deriv (by ring)

/-- `w = exp(-ε/γ) - 1`, code `dw/drs = -(w+1) ε'/γ` (`:186`, `:204`). -/
theorem w_deriv_ok {ε : ℝ → ℝ} {de r γ : ℝ} (h : HasDerivAt ε de r) :
    HasDerivAt (fun x => Real.exp (-ε x / γ) - 1) (-((Real.exp (-ε r / γ) - 1) + 1) * de / γ) r := by
  have := ((h.neg.div_const γ).exp).sub_const 1
  refine this.congr_deriv ?_
  pis
  ring

/-- `β(rs) = β₀(1+0.1rs)/(1+0.1778rs)`, code `dβ = -(0.1778-0.1)β₀/(1+0.1778rs)²` (`:207-208`). -/
theorem beta_deriv_ok (b0 bn bd : ℝ) {r : ℝ} (h : 1 + bd * r ≠ 0) :
    HasDerivAt (fun x => b0 * (1 + bn * x) / (1 + bd * x)) (-(bd - bn) * b0 / (1 + bd * r) ^ 2) r := by
  have hn : HasDerivAt (fun x => b0 * (1 + bn * x)) (b0 * bn) r := by
    simpa using (((hasDerivAt_id r).const_mul bn).const_add 1).const_mul b0
  have hd : HasDerivAt (fun x => 1 + bd * x) bd r := by
    simpa using ((hasDerivAt_id r).const_mul bd).const_add 1
  exact (hn.div hd h).congr_deriv (by field_simp; ring)

/-- `K(rs) = (fc/27)(20 rs D' - 45η D)` with `D = ε⁰-ε¹`; code `dkgrad_drs` (`:212-213`). -/
theorem kgrad_deriv_ok {D D' : ℝ → ℝ} {r d1 d2 fc η : ℝ} (hD : HasDerivAt D d1 r)
    (hD' : HasDerivAt D' d2 r) (h1 : D' r = d1) :
    HasDerivAt (fun x => fc / 27 * (20 * x * D' x - 45 * η * D x))
      (fc / 27 * (20 * D' r + 20 * r * d2 - 45 * η * d1)) r := by
  have := ((((hasDerivAt_id r).const_mul 20).mul hD').sub (hD.const_mul (45 * η))).const_mul (fc / 27)
  refine (this.congr_deriv ?_).congr_of_eventuallyEq ?_
  · pis; rw [h1]; ring
  · filter_upwards with x; simp [mul_assoc]

/-! ## XC-27  outputs (`:248-251`, `salmon_xc.f90:1515-1531`) -/

/-- `vgrad = 2|∇n| ∂f/∂g2` and `rdedd_i = -vgrad ∂_i n/|∇n| = -2 (∂f/∂g2) ∂_i n = -∂f/∂(∂_i n)`. -/
theorem r2scan_rdedd_ok (fg2 gi gnorm : ℝ) (h : gnorm ≠ 0) :
    -(2 * gnorm * fg2) * gi / gnorm = -(2 * fg2 * gi) := by
  field_simp

/-- `e_tau` stored is `Σ vtau·(2 τ_s) = Σ vtau τ_GI` (`salmon_xc.f90:1544`, `τ_s = τ_GI/2`). -/
theorem r2scan_etau_ok (vtau tauGI : ℝ) : vtau * (2 * (tauGI * 0.5)) = vtau * tauGI := by ring

end SalmonAudit.XC
