import Mathlib

/-!
# SALMON v2.3.0 classical FDTD (src/maxwell/fdtd_eh.f90) — update coefficients

Gaussian atomic units: `ε ∂E/∂t = c ∇×H − 4πσE − 4πJ`, `μ ∂H/∂t = −c ∇×E`.
Grids are abstracted: the discrete curl component is represented by a single real `D`
(the finite difference of the neighbouring H components divided by `h`).
-/

namespace SalmonAudit.EM

noncomputable section

/-- EM-01 code: `c1_e` (src/maxwell/fdtd_eh.f90:2569). -/
def c1e (σ ε dt : ℝ) : ℝ := (1 - 2 * Real.pi * σ / ε * dt) / (1 + 2 * Real.pi * σ / ε * dt)
/-- EM-01 code: `c2_e_*` times `h` (src/maxwell/fdtd_eh.f90:2571-2576); the `/hgs` is folded into `D`. -/
def c2e (c σ ε dt : ℝ) : ℝ := (c / ε * dt) / (1 + 2 * Real.pi * σ / ε * dt)
/-- EM-01 code: `c2_j` (src/maxwell/fdtd_eh.f90:2577); enters as `-c2_j*J` (line 2649 etc.). -/
def c2j (σ ε dt : ℝ) : ℝ := (4 * Real.pi / ε * dt) / (1 + 2 * Real.pi * σ / ε * dt)

/-- Code E update `E1 = c1 E0 + c2 D − c2j J` (eh_fd + eh_add_curr). -/
def eUpdateCode (c σ ε dt E0 D J : ℝ) : ℝ := c1e σ ε dt * E0 + c2e c σ ε dt * D - c2j σ ε dt * J

/-- Correct semi-implicit (Crank–Nicolson in the loss term) discretisation of
`ε ∂E/∂t = c D − 4πσE − 4πJ` : residual must vanish. -/
def eResidual (c σ ε dt E0 E1 D J : ℝ) : ℝ :=
  ε * (E1 - E0) / dt - (c * D - 4 * Real.pi * σ * (E1 + E0) / 2 - 4 * Real.pi * J)

theorem em01_yee_E_ok (c σ ε dt E0 D J : ℝ) (hε : 0 < ε) (hdt : 0 < dt) (hσ : 0 ≤ σ) :
    eResidual c σ ε dt E0 (eUpdateCode c σ ε dt E0 D J) D J = 0 := by
  have hpos : 0 < 1 + 2 * Real.pi * σ / ε * dt := by
    have : 0 ≤ 2 * Real.pi * σ / ε * dt := by positivity
    linarith
  unfold eResidual eUpdateCode c1e c2e c2j
  field_simp
  ring

/-- EM-02 code H coefficient `c2_h = c/μ dt` (×1/h folded) (src/maxwell/fdtd_eh.f90:2580). -/
def hUpdateCode (c μ dt H0 D : ℝ) : ℝ := H0 + (c / μ * dt) * D
/-- Correct `μ (H1−H0)/dt = c D'` where `D' = −(curl E)`; sign is carried by the ± assignment
(c2_hx_y = −c2_h_y, c2_hx_z = +c2_h_z, lines 2610-2611) and verified separately below. -/
theorem em02_yee_H_ok (c μ dt H0 D : ℝ) (hμ : 0 < μ) (hdt : 0 < dt) :
    μ * (hUpdateCode c μ dt H0 D - H0) / dt = c * D := by
  unfold hUpdateCode; field_simp; ring

/-- EM-02 sign table: `∂Hx/∂t ∝ −(∂y Ez − ∂z Ey)` ⇒ coefficient of ∂y Ez is −, of ∂z Ey is +;
`∂Ex/∂t ∝ (∂y Hz − ∂z Hy)` ⇒ + and −. Code signs (lines 2593-2618) as a 4-tuple. -/
def signsCode : ℤ × ℤ × ℤ × ℤ := (1, -1, -1, 1)   -- (ex_y, ex_z, hx_y, hx_z)
def signsCorrect : ℤ × ℤ × ℤ × ℤ := (1, -1, -1, 1)
theorem em02_signs_ok : signsCode = signsCorrect := rfl

/-! ## EM-03 Lorentz–Drude ADE (src/maxwell/fdtd_eh.f90:726-740, 3801-3859)
Model per pole: `dJ/dt + γJ + ω0² P = (f ωp²/4π) E`, `dP/dt = J`, trapezoidal in time. -/

def ldC0 (γ ω0 dt : ℝ) : ℝ := 1 + γ * dt / 2 + (ω0 * dt) ^ 2 / 4
/-- code tmp_c1_ld (line 729-732) -/
def ldC1 (γ ω0 dt : ℝ) : ℝ := (1 - γ * dt / 2 - (ω0 * dt) ^ 2 / 4) / ldC0 γ ω0 dt
/-- code tmp_c2_ld (line 733-734) -/
def ldC2 (γ ω0 f ωp dt : ℝ) : ℝ := f * ωp ^ 2 * dt / (8 * Real.pi * ldC0 γ ω0 dt)
/-- code tmp_c3_ld (line 735) -/
def ldC3 (γ ω0 dt : ℝ) : ℝ := ω0 ^ 2 * dt / ldC0 γ ω0 dt

/-- code J update (line 3821-3824) -/
def ldJ1 (γ ω0 f ωp dt J0 P0 E0 E1 : ℝ) : ℝ :=
  ldC1 γ ω0 dt * J0 + ldC2 γ ω0 f ωp dt * (E1 + E0) - ldC3 γ ω0 dt * P0

/-- correct trapezoidal residual of the ADE with `P1 = P0 + dt/2 (J1+J0)` -/
def ldResidual (γ ω0 f ωp dt J0 J1 P0 E0 E1 : ℝ) : ℝ :=
  let P1 := P0 + dt / 2 * (J1 + J0)
  (J1 - J0) / dt + γ * (J1 + J0) / 2 + ω0 ^ 2 * (P1 + P0) / 2 - f * ωp ^ 2 / (4 * Real.pi) * (E1 + E0) / 2

theorem em03_lorentz_drude_ok (γ ω0 f ωp dt J0 P0 E0 E1 : ℝ) (hdt : 0 < dt) (hγ : 0 ≤ γ) :
    ldResidual γ ω0 f ωp dt J0 (ldJ1 γ ω0 f ωp dt J0 P0 E0 E1) P0 E0 E1 = 0 := by
  have h0 : 0 < ldC0 γ ω0 dt := by unfold ldC0; positivity
  have hpi : 0 < Real.pi := Real.pi_pos
  unfold ldResidual ldJ1 ldC1 ldC2 ldC3
  simp only
  field_simp
  unfold ldC0
  ring

/-- EM-04 code E-update in an LD medium (one pole): `c1 = c5/c4`, `c2 = c dt/c4`, `c2j = 4π dt/c4`,
current fed = `½((1+c1)J0 − c3 P0)` (lines 739-740, 2562-2566, 3842-3850). -/
def ldEUpdateCode (γ ω0 f ωp c dt E0 D J0 P0 : ℝ) : ℝ :=
  let c4 := 1 + 2 * Real.pi * dt * ldC2 γ ω0 f ωp dt
  let c5 := 1 - 2 * Real.pi * dt * ldC2 γ ω0 f ωp dt
  c5 / c4 * E0 + c * dt / c4 * D
    - 4 * Real.pi * dt / c4 * (((1 + ldC1 γ ω0 dt) * J0 - ldC3 γ ω0 dt * P0) / 2)

/-- correct: `(E1−E0)/dt = c D − 4π (J1+J0)/2` with J1 from the ADE (implicit in E1) -/
theorem em04_ld_E_ok (γ ω0 f ωp c dt E0 D J0 P0 : ℝ) (hdt : 0 < dt) (hγ : 0 ≤ γ) (hf : 0 ≤ f) :
    let E1 := ldEUpdateCode γ ω0 f ωp c dt E0 D J0 P0
    (E1 - E0) / dt = c * D - 4 * Real.pi * ((ldJ1 γ ω0 f ωp dt J0 P0 E0 E1) + J0) / 2 := by
  intro E1
  have h0 : 0 < ldC0 γ ω0 dt := by unfold ldC0; positivity
  have hc2 : 0 ≤ ldC2 γ ω0 f ωp dt := by unfold ldC2; positivity
  have hc4 : 0 < 1 + 2 * Real.pi * dt * ldC2 γ ω0 f ωp dt := by
    have : 0 ≤ 2 * Real.pi * dt * ldC2 γ ω0 f ωp dt := by positivity
    linarith
  simp only [E1, ldEUpdateCode, ldJ1]
  field_simp
  ring

/-! ## EM-05 PML (src/maxwell/fdtd_eh.f90:2873-2899)
Graded σ(x)=σmax (x/d)^m; field reflection `R = exp(−2 ∫ α dx)`, α = (4π/c)√(μ/ε) σ. -/

/-- code s_max (line 2873-2874) -/
def pmlSmaxCode (m R d c ε μ : ℝ) : ℝ :=
  -(m + 1) * Real.log R / (2 * d) * c / (4 * Real.pi) * Real.sqrt (ε / μ)

/-- correct: `2 (4π/c) √(μ/ε) σmax d/(m+1) = −ln R` -/
theorem em05_pml_smax_ok (m R d c ε μ : ℝ) (hm : 0 < m + 1) (hd : 0 < d) (hc : 0 < c)
    (hε : 0 < ε) (hμ : 0 < μ) :
    2 * (4 * Real.pi / c) * Real.sqrt (μ / ε) * pmlSmaxCode m R d c ε μ * d / (m + 1)
      = -Real.log R := by
  have hs : Real.sqrt (μ / ε) * Real.sqrt (ε / μ) = 1 := by
    rw [← Real.sqrt_mul (by positivity)]
    rw [show μ / ε * (ε / μ) = 1 by field_simp]
    simp
  have hpi : 0 < Real.pi := Real.pi_pos
  unfold pmlSmaxCode
  have : 2 * (4 * Real.pi / c) * Real.sqrt (μ / ε) *
      (-(m + 1) * Real.log R / (2 * d) * c / (4 * Real.pi) * Real.sqrt (ε / μ)) * d / (m + 1)
      = -Real.log R * (Real.sqrt (μ / ε) * Real.sqrt (ε / μ)) := by
    field_simp
  rw [this, hs, mul_one]

/-- EM-05 matched magnetic conductivity `σ* = (μ/ε) σ` (line 2884) ⇒ `σ*/μ = σ/ε`. -/
theorem em05_pml_matched_ok (σ ε μ : ℝ) (hε : 0 < ε) (hμ : 0 < μ) :
    (μ / ε * σ) / μ = σ / ε := by field_simp

/-! ## EM-06 incident current sheet (src/maxwell/fdtd_eh.f90:1337-1342, 3942-4118)
A sheet current K in a medium (ε, μ, η=√(μ/ε)) radiates `E = −(2π/c) η K` to both sides.
To launch amplitude E0, `K = −c E0/(2π η)`, and the E-update adds `−(4π dt/ε)(K/h)`. -/

def c2incCode (c ε μ dt h : ℝ) : ℝ := (c / ε * dt) / 1 * 2 / (h * Real.sqrt (μ / ε))
def c2incCorrect (c ε μ dt h : ℝ) : ℝ := -(4 * Real.pi * dt / ε) * ((-(c / (2 * Real.pi * Real.sqrt (μ / ε)))) / h)

theorem em06_incident_sheet_ok (c ε μ dt h : ℝ) (hε : 0 < ε) (hμ : 0 < μ) (hh : 0 < h) :
    c2incCode c ε μ dt h = c2incCorrect c ε μ dt h := by
  have hs : 0 < Real.sqrt (μ / ε) := Real.sqrt_pos.mpr (by positivity)
  have hpi : 0 < Real.pi := Real.pi_pos
  unfold c2incCode c2incCorrect
  field_simp
  ring

/-! ## EM-07 plane-wave H for the incident field: `H = k̂ × E` (Gaussian, vacuum)
(src/maxwell/fdtd_eh.f90:5931-5978). Code table: for p_axis and E-component, H-component sign. -/

def cross (a b : Fin 3 → ℤ) : Fin 3 → ℤ :=
  ![a 1 * b 2 - a 2 * b 1, a 2 * b 0 - a 0 * b 2, a 0 * b 1 - a 1 * b 0]

def unit (i : Fin 3) : Fin 3 → ℤ := fun j => if i = j then 1 else 0

/-- code: E along x, propagation z ⇒ hy = +ex ; propagation y ⇒ hz = −ex ;
E along y: prop z ⇒ hx = −ey, prop x ⇒ hz = +ey ; E along z: prop x ⇒ hy = −ez, prop y ⇒ hx = +ez -/
def hCode : List ((Fin 3) × (Fin 3) × (Fin 3 → ℤ)) :=
  [(0, 2, ![0, 1, 0]), (0, 1, ![0, 0, -1]), (1, 2, ![-1, 0, 0]),
   (1, 0, ![0, 0, 1]), (2, 0, ![0, -1, 0]), (2, 1, ![1, 0, 0])]

theorem em07_polarization_ok :
    hCode.all (fun t => decide (cross (unit t.2.1) (unit t.1) = t.2.2)) = true := by
  decide

/-! ## EM-08 LR dielectric function from conductivity (src/maxwell/fdtd_eh.f90:4859-4872)
`ε(ω) = 1 + 4πiσ(ω)/ω`; code prints Re ε = 1 − 4π Im σ/ω, Im ε = 4π Re σ/ω. -/
theorem em08_eps_from_sigma_ok (σ : ℂ) (ω : ℝ) (hω : ω ≠ 0) :
    (1 + 4 * Real.pi * Complex.I * σ / ω) =
      ⟨1 - 4 * Real.pi * σ.im / ω, 4 * Real.pi * σ.re / ω⟩ := by
  apply Complex.ext <;> simp [Complex.div_re, Complex.div_im, Complex.normSq] <;>
    field_simp <;> try ring

/-- EM-09 strength function `df/dE = (2ω/π) Im α(ω)` (line 4750): TRK normalisation check for a
single Lorentz oscillator is not formalised; we record the code expression equals the textbook one. -/
def strengthCode (ω imα : ℝ) : ℝ := 2 * ω / Real.pi * imα
def strengthCorrect (ω imα : ℝ) : ℝ := (2 * ω / Real.pi) * imα
theorem em09_strength_ok (ω imα : ℝ) : strengthCode ω imα = strengthCorrect ω imα := rfl

end

end SalmonAudit.EM
