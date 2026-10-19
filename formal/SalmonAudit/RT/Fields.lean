import Mathlib

/-!
# RT external fields (SALMON v2.3.0, `src/rt/em_field.f90`, `src/rt/initialization_rt.f90`)

Conventions (SALMON, atomic units, electron charge `-1`): `Ac = A/c`, `E = -dAc/dt`,
velocity-gauge kinetic momentum `p + Ac`. Polarisation is taken linear and real
(`epdir_im = 0`), one Cartesian component; envelope/time origin shifted so that `tt` is the
code variable `tt = t - tw/2 - t1_start`.
-/

noncomputable section
namespace SalmonAudit.RT

open Real

/-! ## RT-08 intensity → field amplitude -/

/-- Code `initialization_rt.f90:146`: `E0 = sqrt(I)*1e2*27.4492/5.14223e11`, and
`em_field.f90:110`: `f0 = 5.338e-9*sqrt(I)`. The two prefactors agree. -/
theorem I_to_E_two_code_constants :
    |(1e2 * 27.4492 / 5.14223e11 : ℝ) - 5.338e-9| < 1e-12 := by
  rw [abs_lt]; constructor <;> norm_num

/-- Correct: peak intensity `I = c ε0 E0²/2` (SI) ⇒ `E0[V/m] = sqrt(2/(c ε0)) sqrt(I[W/m²])`;
`27.4492 = sqrt(2/(c ε0))` to 4e-6 relative (`c = 299792458`, `ε0 = 8.8541878128e-12`). -/
theorem I_to_E_physical :
    |(27.4492 : ℝ) ^ 2 * (299792458 * 8.8541878128e-12) / 2 - 1| < 1e-5 := by
  rw [abs_lt]; constructor <;> norm_num

/-! ## RT-09 Acos^n envelope: `E = -dA/dt` and CEP -/

/-- `Im(e * exp(iθ)) = e sin θ` used in `em_field.f90:147-150` (`aimag((epdir_re+i epdir_im) exp(iθ))`). -/
theorem aimag_lin (e θ : ℝ) :
    ((e : ℂ) * Complex.exp (Complex.I * (θ : ℂ))).im = e * Real.sin θ := by
  rw [mul_comm Complex.I, Complex.im_ofReal_mul, Complex.exp_ofReal_mul_I_im]

/-- Code `em_field.f90:147`: `A(tt) = -f0/ω cos(π tt/tw)^n sin(ω tt + 2π φ)`. -/
def AcosA (f0 ω tw φ : ℝ) (n : ℕ) (tt : ℝ) : ℝ :=
  -f0 / ω * Real.cos (π / tw * tt) ^ n * Real.sin (ω * tt + 2 * π * φ)

/-- Derivative of the code's `A`. -/
def AcosAdot (f0 ω tw φ : ℝ) (n : ℕ) (tt : ℝ) : ℝ :=
  -f0 / ω * ((n * Real.cos (π / tw * tt) ^ (n - 1) * (-Real.sin (π / tw * tt) * (π / tw))) *
      Real.sin (ω * tt + 2 * π * φ) +
    Real.cos (π / tw * tt) ^ n * (Real.cos (ω * tt + 2 * π * φ) * ω))

theorem AcosA_hasDeriv (f0 ω tw φ : ℝ) (n : ℕ) (tt : ℝ) :
    HasDerivAt (AcosA f0 ω tw φ n) (AcosAdot f0 ω tw φ n tt) tt := by
  have h1 : HasDerivAt (fun t => Real.cos (π / tw * t)) (-Real.sin (π / tw * tt) * (π / tw)) tt := by
    have := ((hasDerivAt_id tt).const_mul (π / tw)).cos
    simpa using this
  have h2 := h1.pow n
  have h3 : HasDerivAt (fun t => Real.sin (ω * t + 2 * π * φ)) (Real.cos (ω * tt + 2 * π * φ) * ω) tt := by
    have := (((hasDerivAt_id tt).const_mul ω).add_const (2 * π * φ)).sin
    simpa using this
  have h4 := (h2.mul h3).const_mul (-f0 / ω)
  unfold AcosA AcosAdot
  convert h4 using 1
  funext t; simp only [Pi.mul_apply, Pi.pow_apply]; ring

/-- RT-09 OK: at the envelope peak `tt = 0`, `E = -dA/dt = f0 cos(2π φ)`: amplitude `f0`, and
`phi_cep` is the carrier-envelope phase in units of 2π of the cosine-like electric field. -/
theorem Acos_peak_field (f0 ω tw φ : ℝ) (n : ℕ) (hω : ω ≠ 0) :
    -AcosAdot f0 ω tw φ n 0 = f0 * Real.cos (2 * π * φ) := by
  unfold AcosAdot
  simp
  field_simp

/-! ## RT-10 Ecos2: the closed-form `A(t)` really integrates `E = f0 cos²(π t/tw) cos(ωt + 2π·0.75)` -/

/-- Code `em_field.f90:167-171` (one Cartesian component, `epdir_re = 1`). -/
def Ecos2A (f0 ω tw : ℝ) (tt : ℝ) : ℝ :=
  -f0 / (8 * π ^ 2 * ω - 2 * tw ^ 2 * ω ^ 3) *
    ((-4 * π ^ 2 + tw ^ 2 * ω ^ 2 + tw ^ 2 * ω ^ 2 * Real.cos (2 * (π / tw * tt))) * Real.cos (ω * tt) +
      2 * π * (2 * π * Real.cos (tw * ω / 2) +
        tw * ω * Real.sin (2 * (π / tw * tt)) * Real.sin (ω * tt)))

/-- Correct field for `ae_shape='Ecos2'` (documented `E = f0 cos²(π t/T) cos(ωt + 2π φ)`, `φ = 0.75`). -/
def Ecos2E (f0 ω tw : ℝ) (tt : ℝ) : ℝ :=
  f0 * Real.cos (π / tw * tt) ^ 2 * Real.cos (ω * tt + 2 * π * (3 / 4))

theorem Ecos2_hasDeriv (f0 ω tw tt : ℝ) (htw : tw ≠ 0)
    (hD : 8 * π ^ 2 * ω - 2 * tw ^ 2 * ω ^ 3 ≠ 0) :
    HasDerivAt (Ecos2A f0 ω tw) (-Ecos2E f0 ω tw tt) tt := by
  have hu : HasDerivAt (fun t => 2 * (π / tw * t)) (2 * (π / tw)) tt := by
    simpa using ((hasDerivAt_id tt).const_mul (π / tw)).const_mul 2
  have hcu := hu.cos
  have hsu := hu.sin
  have hw : HasDerivAt (fun t => ω * t) ω tt := by simpa using (hasDerivAt_id tt).const_mul ω
  have hcw := hw.cos
  have hsw := hw.sin
  have hF := (((hcu.const_mul (tw ^ 2 * ω ^ 2)).const_add (-4 * π ^ 2 + tw ^ 2 * ω ^ 2)).mul hcw).add
    ((((hsu.const_mul (tw * ω)).mul hsw).const_add (2 * π * Real.cos (tw * ω / 2))).const_mul (2 * π))
  have hA := hF.const_mul (-f0 / (8 * π ^ 2 * ω - 2 * tw ^ 2 * ω ^ 3))
  unfold Ecos2A
  convert hA using 1
  · unfold Ecos2E
    have hc : Real.cos (ω * tt + 2 * π * (3 / 4)) = Real.sin (ω * tt) := by
      have : ω * tt + 2 * π * (3 / 4) = ω * tt + π / 2 + π := by ring
      rw [this, Real.cos_add_pi, Real.cos_add_pi_div_two]; ring
    rw [hc, Real.cos_two_mul, Real.sin_two_mul]
    set D := 8 * π ^ 2 * ω - 2 * tw ^ 2 * ω ^ 3 with hDdef
    field_simp
    rw [hDdef]
    ring

/-- `A` vanishes at the pulse start `tt = -tw/2` (continuity with `A = 0` outside). -/
theorem Ecos2A_start (f0 ω tw : ℝ) (htw : tw ≠ 0) : Ecos2A f0 ω tw (-tw / 2) = 0 := by
  unfold Ecos2A
  have : 2 * (π / tw * (-tw / 2)) = -π := by field_simp
  rw [this, Real.cos_neg, Real.sin_neg, Real.cos_pi, Real.sin_pi]
  have h2 : ω * (-tw / 2) = -(tw * ω / 2) := by ring
  rw [h2, Real.cos_neg]
  ring

/-! ## RT-11 length gauge potential sign (`em_field.f90:486-489`) -/

/-- Code: `Vbox = +E·r` added to the electron potential. Correct: electron (charge `q=-1`) in the
scalar potential `φ = -E·r` has potential energy `qφ`. -/
theorem vbox_sign_ok (E x : ℝ) : E * x = (-1) * (-(E * x)) := by ring

/-! ## RT-12 impulse: isolated kick vs periodic step in `Ac` give the same field -/

/-- Isolated (`initialization_rt.f90:433-443`): `ψ → exp(i κ e·r) ψ` shifts the kinetic momentum by `+κ`.
Periodic (`em_field.f90:122-128`): `Ac(t) = κ θ(t)`, kinetic momentum `k + Ac` jumps by `+κ`.
Both correspond to `∫E dt = -κ` for charge `-1` (`Δp = -∫E dt`). -/
theorem impulse_consistent (k κ : ℝ) :
    ((k + κ) - k = κ) ∧ (let dA := κ - 0; -(dA) = -κ) := by
  constructor <;> simp

/-! ## RT-13 NOTE: dipole-source `calc_E_ext('Acos2')` has a quarter-cycle CEP shift -/

/-- Code `em_field.f90:350-351` (used only for `num_dipole_source ≥ 1`):
`E = -(-α sin 2θ1 cos θ2 - β cos²θ1 sin θ2)/β`. At the envelope peak (`θ1 = 0`) this is `sin θ2`. -/
def dipsrcE (α β θ1 θ2 : ℝ) : ℝ :=
  -(-α * Real.sin (2 * θ1) * Real.cos θ2 - β * Real.cos θ1 ^ 2 * Real.sin θ2) / β

theorem dipsrc_cep_shift :
    dipsrcE 1 1 0 0 ≠ -AcosAdot 1 1 1 0 2 0 := by
  rw [Acos_peak_field 1 1 1 0 2 one_ne_zero]
  unfold dipsrcE
  simp

/-! ## RT-14 two-pulse delay (`em_field.f90:145,153,249`) -/

/-- Pulse-1 centre `t1_start + tw1/2`; pulse-2 centre `tw1/2 + (T1_T2 + t1_start)`;
difference = `T1_T2` (centre-to-centre delay, as documented). -/
theorem delay_ok (tw1 t1s T12 : ℝ) : (tw1 / 2 + (T12 + t1s)) - (t1s + tw1 / 2) = T12 := by ring

end SalmonAudit.RT
end
