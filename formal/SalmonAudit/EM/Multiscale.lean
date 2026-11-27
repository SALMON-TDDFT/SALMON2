import Mathlib

/-!
# SALMON v2.3.0 — Weyl-gauge wave equation (src/maxwell/fdtd_weyl_gauge.f90), multiscale
coupling (src/ms/main_ms.f90), single-cell Maxwell coupling (src/rt/em_field.f90), and output
units / Poynting flux in src/maxwell/fdtd_eh.f90.

Convention (Gaussian a.u.): `Ac = A/c`, `E = −∂t Ac`, `H = c ∇×Ac`, physical current
`J = −curr` (SALMON `curr` is the electron-number current; cf. main_ms.f90:469, em_field.f90:54).
Wave equation: `ε ∂²Ac/∂t² = 4πJ − c² ∇×∇×Ac`.
Grids abstracted: `R` stands for the discrete `∇×∇×Ac` component at the cell.
-/

namespace SalmonAudit.EM

noncomputable section

/-- correct leap-frog residual of `ε ∂²Ac/∂t² = 4πJ − c² R` -/
def waveResidual (ε c dt A0 A A1 J R : ℝ) : ℝ :=
  ε * (A1 - 2 * A + A0) / dt ^ 2 - (4 * Real.pi * J - c ^ 2 * R)

/-- EM-10 code 3D update (src/maxwell/fdtd_weyl_gauge.f90:237-240) -/
def weyl3dCode (ε c dt A0 A J R : ℝ) : ℝ :=
  2 * A - A0 + 4 * Real.pi * dt ^ 2 * J / ε - R * (c * dt) ^ 2 / ε

theorem em10_weyl3d_ok (ε c dt A0 A J R : ℝ) (hε : 0 < ε) (hdt : 0 < dt) :
    waveResidual ε c dt A0 A (weyl3dCode ε c dt A0 A J R) J R = 0 := by
  unfold waveResidual weyl3dCode; field_simp; ring

/-- EM-11 code 1D update (src/maxwell/fdtd_weyl_gauge.f90:141-143): the whole bracket,
including `2A − A_old`, is divided by ε a second time. -/
def weyl1dCode (ε c dt A0 A J R : ℝ) : ℝ :=
  (2 * A - A0 + J * 4 * Real.pi * dt ^ 2 / ε - R * (c * dt) ^ 2) / ε

/-- EM-11 BUG: for a uniform static field (J = R = 0, A = A0 = 1) in a dielectric ε = 2
the code does not keep Ac constant (gives 1/2), while the correct update keeps 1. -/
theorem em11_weyl1d_bug :
    weyl1dCode 2 1 1 1 1 0 0 ≠ weyl3dCode 2 1 1 1 1 0 0 := by
  unfold weyl1dCode weyl3dCode; norm_num

/-- EM-11: the 1D code is correct only for ε = 1 (general statement). -/
theorem em11_weyl1d_eps1 (c dt A0 A J R : ℝ) :
    weyl1dCode 1 c dt A0 A J R = weyl3dCode 1 c dt A0 A J R := by
  unfold weyl1dCode weyl3dCode; ring

/-- EM-12 1st-order Mur ABC with incident field, left end
(src/maxwell/fdtd_weyl_gauge.f90:155-159). Code value of the ghost cell `A0^{n+1}`. -/
def murCode (c dt h A0n A1n A1np dInc : ℝ) : ℝ :=
  A1n + (c * dt - h) / (c * dt + h) * A1np - (c * dt - h) / (c * dt + h) * A0n
    + 4 * h / (c * dt + h) * dInc

/-- correct: box-centred discretisation of `(∂t − c∂x)A = 2 ∂t A_inc` at x_{1/2}, t_{n+1/2}. -/
def murResidual (c dt h A0n A0np A1n A1np dInc : ℝ) : ℝ :=
  ((A0np - A0n) + (A1np - A1n)) / (2 * dt) - c * ((A1np + A1n) - (A0np + A0n)) / (2 * h)
    - 2 * dInc / dt

theorem em12_mur_incident_ok (c dt h A0n A1n A1np dInc : ℝ) (hc : 0 < c) (hdt : 0 < dt)
    (hh : 0 < h) :
    murResidual c dt h A0n (murCode c dt h A0n A1n A1np dInc) A1n A1np dInc = 0 := by
  have : 0 < c * dt + h := by positivity
  unfold murResidual murCode; field_simp; ring

/-- EM-13 incident/reflected split (src/ms/main_ms.f90:744-756).
With `A = f(t−x/c) + g(t+x/c)`: `∂tA = f' + g'`, `∂xA = (−f' + g')/c`; `E_inc = −f'`, `E_ref = −g'`. -/
theorem em13_wave_split_ok (c fp gp : ℝ) (hc : c ≠ 0) :
    let dtA := fp + gp
    let dxA := (-fp + gp) / c
    (-(1/2) * (dtA - c * dxA) = -fp) ∧ (-(1/2) * (dtA + c * dxA) = -gp) := by
  simp only
  constructor <;> field_simp <;> ring

/-- EM-14 macro→micro frame (src/ms/main_ms.f90:435-438): code `Ac_ext = R (R v)`;
correct is `R v` (= Ac_tot, crystal frame) or `Rᵀ(R v) = v` (lab frame). For a 90° rotation
about z and v = x̂ the code gives −x̂, neither x̂ nor ŷ. -/
def rot90 : Matrix (Fin 3) (Fin 3) ℝ := !![0, -1, 0; 1, 0, 0; 0, 0, 1]
def vx : Fin 3 → ℝ := ![1, 0, 0]

theorem em14_rmat_squared_bug :
    rot90.mulVec (rot90.mulVec vx) ≠ rot90.mulVec vx ∧
    rot90.mulVec (rot90.mulVec vx) ≠ vx := by
  constructor
  · intro h
    have := congrFun h 0
    simp [rot90, vx, Matrix.mulVec, dotProduct, Fin.sum_univ_succ] at this
  · intro h
    have := congrFun h 0
    simp [rot90, vx, Matrix.mulVec, dotProduct, Fin.sum_univ_succ] at this
    norm_num at this

theorem em14_rmat_orthogonal_check : rot90.transpose * rot90 = 1 := by
  ext i j; fin_cases i <;> fin_cases j <;> simp [rot90, Matrix.mul_apply, Fin.sum_univ_succ]

/-- EM-15 single-cell 'lo' (src/rt/em_field.f90:54): `Ac'' = −4π curr`, with `J = −curr`
this is `Ac'' = 4πJ` (= longitudinal: `E_ind = −4πP`). -/
def loCode (A Aold curr dt : ℝ) : ℝ := 2 * A - Aold - 4 * Real.pi * curr * dt ^ 2
theorem em15_lo_ok (A Aold curr dt : ℝ) (hdt : 0 < dt) :
    (loCode A Aold curr dt - 2 * A + Aold) / dt ^ 2 = 4 * Real.pi * (-curr) := by
  unfold loCode; field_simp; ring

/-- EM-16 thin film '2d' (src/rt/em_field.f90:57-66): code rate gives
`E_r = r E_i − 4π d J/(c(n1+n2))` with `r = (n1−n2)/(n1+n2)`, `J = −curr`.
Correct: `E_t = E_i + E_r` and `n2 E_t − n1 (E_i − E_r) = −(4π/c) K`, `K = d J`. -/
def filmErCode (n1 n2 c d Ei J : ℝ) : ℝ :=
  (n1 - n2) / (n1 + n2) * Ei - 4 * Real.pi / (c * (n1 + n2)) * d * J

theorem em16_thin_film_ok (n1 n2 c d Ei J : ℝ) (hn : 0 < n1 + n2) (hc : 0 < c) :
    let Er := filmErCode n1 n2 c d Ei J
    let Et := Ei + Er
    n2 * Et - n1 * (Ei - Er) = -(4 * Real.pi / c) * (d * J) := by
  simp only [filmErCode]; field_simp; ring

/-- EM-17 field energy density (src/maxwell/fdtd_weyl_gauge.f90:379; fdtd_coulomb_gauge.f90:218)
`u = (E²+H²)/(8π)` and Poynting `S = (c/4π) E×H` in the single-scale path (line 201). -/
def uCode (E2 H2 : ℝ) : ℝ := (1 / (8 * Real.pi)) * (E2 + H2)
theorem em17_energy_density_ok (E2 H2 : ℝ) : uCode E2 H2 = (E2 + H2) / (8 * Real.pi) := by
  unfold uCode; ring

/-- EM-18 Poynting for TTM (src/maxwell/fdtd_eh.f90:6476-6478): code `S = E×H` (no prefactor);
correct Gaussian `S = (c/4π) E×H`. Counterexample: one component with E·H = 1, c = 137.036. -/
def poyntingCode (EH : ℝ) : ℝ := EH
def poyntingCorrect (c EH : ℝ) : ℝ := c / (4 * Real.pi) * EH

theorem em18_poynting_ttm_bug : poyntingCode 1 ≠ poyntingCorrect 137.036 1 := by
  unfold poyntingCode poyntingCorrect
  have h1 : Real.pi < 3.15 := Real.pi_lt_d2
  have h0 : 0 < Real.pi := Real.pi_pos
  intro h
  rw [mul_one, eq_div_iff (by positivity)] at h
  linarith

/-- general: code = (4π/c) × correct -/
theorem em18_poynting_ratio (c EH : ℝ) (hc : c ≠ 0) :
    poyntingCode EH = (4 * Real.pi / c) * poyntingCorrect c EH := by
  unfold poyntingCode poyntingCorrect
  have : Real.pi ≠ 0 := Real.pi_ne_zero
  field_simp

/-- EM-19 H output in A_eV_fs units (src/maxwell/fdtd_eh.f90:379-380, 1236-1238, 3393-3401):
code `H_out = H_au · u_E` (u_E = au of E in V/Å, label "A/Å"); correct SI `H = H_au · u_E / Z0`,
Z0 = 376.730313668 Ω (Gaussian H and E share units; SI H = B/μ0 = E/Z0 for a plane wave). -/
def Z0 : ℝ := 376.730313668
def hOutCode (H uE : ℝ) : ℝ := H * uE
def hOutCorrect (H uE : ℝ) : ℝ := H * uE / Z0

theorem em19_h_unit_bug : hOutCode 1 1 ≠ hOutCorrect 1 1 := by
  unfold hOutCode hOutCorrect Z0; norm_num

/-- EM-20 intensity → peak field (src/rt/em_field.f90:113): `E0 = 5.338e-9 √I`; correct
`E0 = √(I / I_au)`, `I_au = ½ c ε0 E_au² = 3.50945e16 W/cm²`. Agreement to 1e-4 relative. -/
theorem em20_intensity_ok : |(5.338e-9 : ℝ) ^ 2 * 3.50945e16 - 1| < 1e-4 := by
  norm_num [abs_lt]

/-- EM-21 NOTE: Acos2 carrier in classical FDTD (src/maxwell/fdtd_eh.f90:5917-5920) has leading
term `E0 cos²(·) sin θ2`, whereas TDDFT Acos2 (src/rt/em_field.f90:147-149) gives
`E0 cos²(·) cos θ2` for the same phi_cep: a quarter-cycle CEP offset. -/
theorem em21_acos2_cep_note :
    Real.sin (0 : ℝ) ≠ Real.cos 0 ∧ ∀ θ : ℝ, Real.sin θ = Real.cos (θ - Real.pi / 2) := by
  refine ⟨by simp, fun θ => ?_⟩
  rw [Real.cos_sub_pi_div_two]

/-- EM-21 fix: with `φ_FDTD = φ_TDDFT + 1/4` the FDTD carrier `cos(x + 2πφ_FDTD)` equals the
TDDFT carrier `−sin(x + 2πφ_TDDFT)` identically in x, so A (hence E = −∂tA including the
envelope-derivative term) coincide. -/
theorem em21_cep_shift_ok (x φ : ℝ) :
    Real.cos (x + 2 * Real.pi * (φ + 1 / 4)) = -Real.sin (x + 2 * Real.pi * φ) := by
  have : x + 2 * Real.pi * (φ + 1 / 4) = (x + 2 * Real.pi * φ) + Real.pi / 2 := by ring
  rw [this, Real.cos_add_pi_div_two]

/-- EM-23 polarizability unit conversion (src/maxwell/fdtd_eh.f90:4744-4749):
code `α·a0²/uE`, correct `α·a0/uE` (α = p/E, p in e·Å, E in V/Å). -/
def alphaOutCode (α a0 uE : ℝ) : ℝ := α * (a0 ^ 2 / uE)
def alphaOutCorrect (α a0 uE : ℝ) : ℝ := α * (a0 / uE)

theorem em23_alpha_unit_bug : alphaOutCode 1 0.529177 1 ≠ alphaOutCorrect 1 0.529177 1 := by
  unfold alphaOutCode alphaOutCorrect; norm_num

theorem em23_alpha_ratio (α a0 uE : ℝ) : alphaOutCode α a0 uE = a0 * alphaOutCorrect α a0 uE := by
  unfold alphaOutCode alphaOutCorrect; ring

end

end SalmonAudit.EM
