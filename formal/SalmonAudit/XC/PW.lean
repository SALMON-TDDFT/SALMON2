import Mathlib

/-!
# XC audit: Perdew–Wang 92 correlation

Sources: `src/xc/builtin_pw.f90` (`PWxc`, used by `xc='pw'`), the identical copy `PWc` in
`src/xc/builtin_tbmbj.f90:158-207`, and `src/xc/builtin_pw92.f90` (used by r2SCAN).

Abstractions: PW parameters are arbitrary positive reals; the Fortran real power
`rs**(p+1)` with `p=1` is modelled as `r^2`; `|ζ|<1` (the `eps` clamp is dropped).
-/

set_option linter.unusedSimpArgs false

namespace SalmonAudit.XC

open Real

/-- `Q1 = 2A(β₁√r+β₂r+β₃r√r+β₄r²)` (`builtin_pw.f90:88`). -/
noncomputable def pwQ1 (A b1 b2 b3 b4 r : ℝ) : ℝ := 2 * A * (b1 * √r + b2 * r + b3 * r * √r + b4 * r ^ 2)

/-- `G = -2A(1+α₁r) ln(1+1/Q1)` (`builtin_pw.f90:85,94`). -/
noncomputable def pwG (A a1 b1 b2 b3 b4 r : ℝ) : ℝ :=
  -2 * A * (1 + a1 * r) * Real.log (1 + 1 / pwQ1 A b1 b2 b3 b4 r)

/-- code `dG = -2Aα₁ ln(1+1/Q1) - Q0·Qd/(Q1²+Q1)` (`builtin_pw.f90:91,97`). -/
noncomputable def pwdGCode (A a1 b1 b2 b3 b4 r : ℝ) : ℝ :=
  -2 * A * a1 * Real.log (1 + 1 / pwQ1 A b1 b2 b3 b4 r)
    - (-2 * A * (1 + a1 * r)) * (A * (b1 / √r + 2 * b2 + 3 * b3 * √r + 2 * (1 + 1) * b4 * r))
      / (pwQ1 A b1 b2 b3 b4 r ^ 2 + pwQ1 A b1 b2 b3 b4 r)

theorem pw_G_deriv_ok (A a1 b1 b2 b3 b4 : ℝ) (hA : 0 < A) (h1 : 0 < b1) (h2 : 0 ≤ b2)
    (h3 : 0 ≤ b3) (h4 : 0 ≤ b4) {r : ℝ} (hr : 0 < r) :
    HasDerivAt (pwG A a1 b1 b2 b3 b4) (pwdGCode A a1 b1 b2 b3 b4 r) r := by
  obtain ⟨s, hs, rfl⟩ : ∃ s, 0 < s ∧ r = s ^ 2 :=
    ⟨√r, Real.sqrt_pos.2 hr, (Real.sq_sqrt hr.le).symm⟩
  have hsq : √(s ^ 2) = s := Real.sqrt_sq hs.le
  have hsqrt := Real.hasDerivAt_sqrt (x := s ^ 2) (by positivity)
  rw [hsq] at hsqrt
  have hid := hasDerivAt_id (s ^ 2)
  have hQ : HasDerivAt (pwQ1 A b1 b2 b3 b4)
      (2 * A * (b1 * (1 / (2 * s)) + b2 * 1 + b3 * (1 * s + s ^ 2 * (1 / (2 * s)))
        + b4 * ((2:ℕ) * (s ^ 2) ^ (2 - 1)))) (s ^ 2) := by
    have := ((((hsqrt.const_mul b1).add (hid.const_mul b2)).add
      ((hid.mul hsqrt).const_mul b3)).add ((hasDerivAt_pow 2 (s ^ 2)).const_mul b4)).const_mul (2 * A)
    refine (this.congr_deriv (by simp [hsq])).congr_of_eventuallyEq ?_
    filter_upwards with x
    unfold pwQ1; simp [mul_assoc]
  have hQpos : 0 < pwQ1 A b1 b2 b3 b4 (s ^ 2) := by unfold pwQ1; rw [hsq]; positivity
  have hin : HasDerivAt (fun x => 1 + 1 / pwQ1 A b1 b2 b3 b4 x)
      (-(2 * A * (b1 * (1 / (2 * s)) + b2 * 1 + b3 * (1 * s + s ^ 2 * (1 / (2 * s)))
        + b4 * ((2:ℕ) * (s ^ 2) ^ (2 - 1)))) / (pwQ1 A b1 b2 b3 b4 (s ^ 2)) ^ 2) (s ^ 2) := by
    have := (hQ.inv hQpos.ne').const_add 1
    refine this.congr_of_eventuallyEq ?_
    filter_upwards with x; simp
  have hpos : 0 < 1 + 1 / pwQ1 A b1 b2 b3 b4 (s ^ 2) := by positivity
  have hlog := hin.log hpos.ne'
  have hpre : HasDerivAt (fun x => -2 * A * (1 + a1 * x)) (-2 * A * (a1 * 1)) (s ^ 2) := by
    have := ((hid.const_mul a1).const_add 1).const_mul (-2 * A)
    exact this
  have := hpre.mul hlog
  refine (this.congr_deriv ?_).congr_of_eventuallyEq ?_
  · unfold pwdGCode
    simp only [hsq]
    generalize pwQ1 A b1 b2 b3 b4 (s ^ 2) = Q at hQpos ⊢
    generalize Real.log (1 + 1 / Q) = Lg
    field_simp
    ring
  · filter_upwards with x; unfold pwG; simp only [Pi.mul_apply, Pi.add_apply, Pi.sub_apply, Pi.neg_apply, Pi.inv_apply, Pi.div_apply, id]; try ring

/-- PW ζ-interpolation `ε = G0 - Gα f/f''(0) (1-ζ⁴) + (G1-G0) f ζ⁴` (`builtin_pw.f90:104`). -/
noncomputable def pwEz (G0 G1 Ga f2 : ℝ) (f : ℝ → ℝ) (z : ℝ) : ℝ :=
  G0 - Ga * f z / f2 * (1 - z ^ 4) + (G1 - G0) * f z * z ^ 4

/-- code `decdzeta` (`builtin_pw.f90:101-102`). -/
noncomputable def pwdEdzCode (G0 G1 Ga f2 fz dfz z : ℝ) : ℝ :=
  dfz * (-Ga / f2 * (1 - z ^ 4) + (G1 - G0) * z ^ 4) + 4 * z ^ 3 * (fz * (G1 - G0) + Ga * fz / f2)

theorem pw_dEdzeta_ok (G0 G1 Ga f2 : ℝ) (f : ℝ → ℝ) {z dfz : ℝ} (hf : HasDerivAt f dfz z) :
    HasDerivAt (pwEz G0 G1 Ga f2 f) (pwdEdzCode G0 G1 Ga f2 (f z) dfz z) z := by
  have h4 := hasDerivAt_pow 4 z
  have := ((hasDerivAt_const z G0).sub (((hf.const_mul Ga).div_const f2).mul
    ((hasDerivAt_const z (1:ℝ)).sub h4))).add ((hf.const_mul (G1 - G0)).mul h4)
  refine (this.congr_deriv ?_).congr_of_eventuallyEq ?_
  · unfold pwdEdzCode; simp only [Pi.mul_apply, Pi.add_apply, Pi.sub_apply, Pi.neg_apply, Pi.inv_apply, Pi.div_apply, id]; push_cast; ring
  · filter_upwards with x; unfold pwEz; simp only [Pi.mul_apply, Pi.add_apply, Pi.sub_apply, Pi.neg_apply, Pi.inv_apply, Pi.div_apply, id]; try ring

/-- code `decdrs` (`builtin_pw.f90:100`) is `∂ε/∂rs` given `G0',G1',Gα'`. -/
theorem pw_dEdrs_ok {G0 G1 Ga : ℝ → ℝ} {d0 d1 da r : ℝ} (f2 fz z : ℝ)
    (h0 : HasDerivAt G0 d0 r) (h1 : HasDerivAt G1 d1 r) (ha : HasDerivAt Ga da r) :
    HasDerivAt (fun x => G0 x - Ga x * fz / f2 * (1 - z ^ 4) + (G1 x - G0 x) * fz * z ^ 4)
      (d0 * (1 - fz * z ^ 4) + d1 * fz * z ^ 4 - da * fz / f2 * (1 - z ^ 4)) r := by
  have := (h0.sub (((ha.mul_const fz).div_const f2).mul_const (1 - z ^ 4))).add
    (((h1.sub h0).mul_const fz).mul_const (z ^ 4))
  exact this.congr_deriv (by ring)

/-- chain rule to `ρ_a` (`builtin_pw.f90:105`): with `∂rs/∂ρa = -rs/(3ρ)` and
`∂ζ/∂ρa = (1-ζ)/ρ` the code `(-rs/3 decdrs - (ζ-1) decdzeta)/ρ` is the total derivative. -/
theorem pw_drhoa_ok (r z ρ Er Ez : ℝ) (hρ : ρ ≠ 0) :
    (-r / 3 * Er - (z - 1) * Ez) / ρ = Er * (-r / (3 * ρ)) + Ez * ((1 - z) / ρ) := by
  field_simp; ring

theorem pw_drhob_ok (r z ρ Er Ez : ℝ) (hρ : ρ ≠ 0) :
    (-r / 3 * Er - (z + 1) * Ez) / ρ = Er * (-r / (3 * ρ)) + Ez * ((-1 - z) / ρ) := by
  field_simp; ring

/-- XC-07: unpolarized `exc_cor_pw` (`builtin_pw.f90:41-45`) calls `PWxc(ρ_s,ρ_s)`, so ζ=0,
`decdzeta = f'(0)(...) + 0 = 0`, and `Vexc = ε + 2ρ_s·∂ε/∂ρa` equals `d(ρε)/dρ = ε - (rs/3)∂ε/∂rs`. -/
theorem pw_unpolarized_potential_ok (eps r rhos Er A B : ℝ) (hρ : rhos ≠ 0) :
    let z : ℝ := 0
    let dfz0 : ℝ := 4 / 3 * ((1 + z) ^ ((1:ℝ)/3) - (1 - z) ^ ((1:ℝ)/3)) / A
    let Ez := dfz0 * B + 4 * z ^ 3 * B
    eps + (2 * rhos) * ((-r / 3 * Er - (z - 1) * Ez) / (2 * rhos))
      = eps - r / 3 * Er := by
  simp only
  norm_num
  field_simp
  ring

/-! ## builtin_pw92.f90 (r2SCAN's PW92) : first and second derivatives -/

noncomputable def p92Q (b1 b2 b3 b4 r : ℝ) : ℝ := √r * (b1 + b3 * r) + r * (b2 + b4 * r)
noncomputable def p92dQ (b1 b2 b3 b4 r : ℝ) : ℝ := 0.5 * b1 / √r + b2 + 1.5 * b3 * √r + 2 * b4 * r
noncomputable def p92d2Q (b1 b3 b4 r : ℝ) : ℝ := -0.25 * b1 / (r * √r) + 0.75 * b3 / √r + 2 * b4

lemma sqrt_deriv_s {s : ℝ} (hs : 0 < s) : HasDerivAt Real.sqrt (1 / (2 * s)) (s ^ 2) := by
  have h := Real.hasDerivAt_sqrt (x := s ^ 2) (by positivity)
  rwa [Real.sqrt_sq hs.le] at h

theorem p92_dQ_ok (b1 b2 b3 b4 : ℝ) {r : ℝ} (hr : 0 < r) :
    HasDerivAt (p92Q b1 b2 b3 b4) (p92dQ b1 b2 b3 b4 r) r := by
  obtain ⟨s, hs, rfl⟩ : ∃ s, 0 < s ∧ r = s ^ 2 :=
    ⟨√r, Real.sqrt_pos.2 hr, (Real.sq_sqrt hr.le).symm⟩
  have hq := sqrt_deriv_s hs
  have hid := hasDerivAt_id (s ^ 2)
  have := (hq.mul ((hid.const_mul b3).const_add b1)).add (hid.mul ((hid.const_mul b4).const_add b2))
  refine (this.congr_deriv ?_).congr_of_eventuallyEq ?_
  · unfold p92dQ; simp only [Pi.mul_apply, Pi.add_apply, Pi.sub_apply, Pi.neg_apply, Pi.inv_apply, Pi.div_apply, id]; rw [Real.sqrt_sq hs.le]; field_simp; ring
  · filter_upwards with x; unfold p92Q; simp

theorem p92_d2Q_ok (b1 b2 b3 b4 : ℝ) {r : ℝ} (hr : 0 < r) :
    HasDerivAt (p92dQ b1 b2 b3 b4) (p92d2Q b1 b3 b4 r) r := by
  obtain ⟨s, hs, rfl⟩ : ∃ s, 0 < s ∧ r = s ^ 2 :=
    ⟨√r, Real.sqrt_pos.2 hr, (Real.sq_sqrt hr.le).symm⟩
  have hq := sqrt_deriv_s hs
  have hqs : √(s ^ 2) ≠ 0 := by rw [Real.sqrt_sq hs.le]; exact hs.ne'
  have hid := hasDerivAt_id (s ^ 2)
  have hinv := (hq.inv hqs)
  have := (((hinv.const_mul (0.5 * b1)).add_const b2).add (hq.const_mul (1.5 * b3))).add
    (hid.const_mul (2 * b4))
  refine (this.congr_deriv ?_).congr_of_eventuallyEq ?_
  · unfold p92d2Q; rw [Real.sqrt_sq hs.le]; field_simp; ring
  · filter_upwards with x; unfold p92dQ; simp [div_eq_mul_inv]

/-- `L = ln(1+1/(2AQ))`; code `dlogs = -Q'/M`, `M = Q(2AQ+1)` (`builtin_pw92.f90:46-48`). -/
theorem p92_dlog_ok {Q : ℝ → ℝ} {Q' r A : ℝ} (hA : 0 < A) (hQ : HasDerivAt Q Q' r) (hQp : 0 < Q r) :
    HasDerivAt (fun x => Real.log (1 + 1 / (2 * A * Q x))) (-Q' / (Q r * (2 * A * Q r + 1))) r := by
  have h2 : 0 < 2 * A * Q r := by positivity
  have hin := ((hQ.const_mul (2 * A)).inv h2.ne').const_add 1
  have hpos : 0 < 1 + (2 * A * Q r)⁻¹ := by positivity
  have := hin.log hpos.ne'
  refine (this.congr_deriv ?_).congr_of_eventuallyEq ?_
  · simp only [Pi.mul_apply, Pi.add_apply, Pi.sub_apply, Pi.neg_apply, Pi.inv_apply, Pi.div_apply, id]; field_simp
  · filter_upwards with x; simp

/-- code `d2logs = -Q''/M + Q'^2 (4AQ+1)/M^2` (`builtin_pw92.f90:49`). -/
theorem p92_d2log_ok {Q dQ : ℝ → ℝ} {r A Q' Q'' : ℝ} (hQ : HasDerivAt Q Q' r)
    (hdQ : HasDerivAt dQ Q'' r) (hQ'eq : dQ r = Q') (hM : Q r * (2 * A * Q r + 1) ≠ 0) :
    HasDerivAt (fun x => -dQ x / (Q x * (2 * A * Q x + 1)))
      (-Q'' / (Q r * (2 * A * Q r + 1))
        + Q' * Q' * (4 * A * Q r + 1) / ((Q r * (2 * A * Q r + 1)) * (Q r * (2 * A * Q r + 1)))) r := by
  have hM' := hQ.mul (((hQ.const_mul (2 * A)).add_const 1))
  have := hdQ.neg.div hM' hM
  refine this.congr_deriv ?_
  simp only [Pi.mul_apply, Pi.add_apply, Pi.sub_apply, Pi.neg_apply, Pi.inv_apply, Pi.div_apply, id]
  rw [hQ'eq]
  field_simp
  ring

/-- `ec = -2A(1+a1 r)L`, code `dec_drs = -2A(a1 L + (1+a1 r) L')` and
`d2ec = -2A(2 a1 L' + (1+a1 r) L'')` (`builtin_pw92.f90:52-54`). -/
theorem p92_ec_derivs_ok {L dL : ℝ → ℝ} {r A a1 L' L'' : ℝ} (hL : HasDerivAt L L' r)
    (hdL : HasDerivAt dL L'' r) (hdLeq : dL r = L') :
    HasDerivAt (fun x => -2 * A * (1 + a1 * x) * L x) (-2 * A * (a1 * L r + (1 + a1 * r) * L')) r ∧
    HasDerivAt (fun x => -2 * A * (a1 * L x + (1 + a1 * x) * dL x))
      (-2 * A * (2 * a1 * L' + (1 + a1 * r) * L'')) r := by
  have hp : HasDerivAt (fun x => 1 + a1 * x) a1 r := by
    simpa using ((hasDerivAt_id r).const_mul a1).const_add 1
  constructor
  · have := (hp.const_mul (-2 * A)).mul hL
    refine (this.congr_deriv (by ring)).congr_of_eventuallyEq ?_
    filter_upwards with x; simp only [Pi.mul_apply, Pi.add_apply, Pi.sub_apply, Pi.neg_apply, Pi.inv_apply, Pi.div_apply, id]; try ring
  · have := ((hL.const_mul a1).add (hp.mul hdL)).const_mul (-2 * A)
    refine this.congr_deriv ?_
    rw [hdLeq]; ring

end SalmonAudit.XC
