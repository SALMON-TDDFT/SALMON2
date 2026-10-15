import Mathlib

/-!
# SBE length gauge (SALMON v2.3.0, `src/ssbe/bloch_solver_ssbe.f90:484-767`,
`src/ssbe/gs_info_ssbe.f90:279-306`, `src/ssbe/common_ssbe.f90`)

Conventions (a.u., electron charge −1): `H = H0 + E·r`, `r = i∇_k + ξ`,
`ρ̇ = -i[H0,ρ] + E·∇_k ρ - i E·[ξ,ρ]`;  `ξ_ij = -i v_ij/(ε_i-ε_j)` (from `[H0,r] = -i v`).
The code's intraband drift `+E·∇q` (line 655/673) agrees with this.
-/

namespace SalmonAudit.SBE.LG

open Complex

noncomputable section

/-- SBE-10 code `d_matrix` (`gs_info_ssbe.f90:297-299`): `d_ij = i·p_ij/(ε_i-ε_j)`
(`p` = `p_mod_matrix` = velocity matrix element incl. nonlocal commutator). -/
def codeD (p : ℂ) (ω : ℝ) : ℂ := I * p / ω

/-- Interband Berry connection / position matrix element `ξ_ij = -i v_ij/(ε_i-ε_j)`. -/
def xi (p : ℂ) (ω : ℝ) : ℂ := -I * p / ω

theorem sbe10_d_eq_neg_xi (p : ℂ) (ω : ℝ) : codeD p ω = - xi p ω := by
  unfold codeD xi; ring

theorem sbe10_d_ne_xi : codeD 1 1 ≠ xi 1 1 := by
  unfold codeD xi
  intro h
  have := congrArg Complex.im h
  simp at this
  norm_num at this

/-! Two-band equations of motion (gradient and shift-vector terms omitted; both commute with the
constant sign matrix `S = diag(1,-1)` used below).  Interband part of the RHS with a coupling
matrix element `x` (`x = ξ` correct, `x = d` code; `bloch_solver_ssbe.f90:655-687`): -/

/-- off-diagonal: `ρ̇₁₂ = -iω ρ₁₂ + i E x₁₂ (ρ₁₁ - ρ₂₂)` -/
def rhsOff (x12 : ℂ) (ω E : ℝ) (ρ11 ρ22 ρ12 : ℂ) : ℂ :=
  -I * ω * ρ12 + I * E * x12 * (ρ11 - ρ22)

/-- diagonal: `ρ̇₁₁ = -i E (x₁₂ ρ₂₁ - ρ₁₂ x₂₁)` -/
def rhsDiag (x12 x21 : ℂ) (E : ℝ) (ρ12 ρ21 : ℂ) : ℂ :=
  -I * E * (x12 * ρ21 - ρ12 * x21)

/-- With `d = -ξ`, the code's dynamics is the correct dynamics conjugated by `S = diag(1,-1)`:
off-diagonal elements flip sign, populations unchanged (exact for two bands). -/
theorem sbe10_two_band_S_map (ξ12 ξ21 ρ11 ρ22 ρ12 ρ21 : ℂ) (ω E : ℝ) :
    rhsOff (-ξ12) ω E ρ11 ρ22 (-ρ12) = - rhsOff ξ12 ω E ρ11 ρ22 ρ12 ∧
    rhsDiag (-ξ12) (-ξ21) E (-ρ12) (-ρ21) = rhsDiag ξ12 ξ21 E ρ12 ρ21 := by
  unfold rhsOff rhsDiag; constructor <;> ring

/-- SBE-10 code current (`bloch_solver_ssbe.f90:752-765`): `J = -(I_intra + X_inter[ρ_code])/V`
with `X_inter[ρ] = Σ_{i≠j} ω_ij|d_ij| Im q_ji = Tr(ρ p)_inter`.  By the S-map,
`X_inter[ρ_code] = -X_inter[ρ_true]`, so `V·J_code = -I + X` with `I, X` the true intraband and
interband `Tr(ρ v)` parts.  Correct (any global sign convention `s`): `s·(I + X)`. -/
def codeJ (Iintra X : ℝ) : ℝ := -Iintra + X

theorem sbe10_lg_current_relative_sign_bug :
    ¬ ∃ s : ℝ, ∀ Iintra X : ℝ, codeJ Iintra X = s * (Iintra + X) := by
  rintro ⟨s, hs⟩
  have h1 := hs 1 0
  have h2 := hs 0 1
  unfold codeJ at h1 h2
  linarith

/-- SBE-11 code, off-diagonal RHS (`bloch_solver_ssbe.f90:669-687`, three-band term omitted):
the free precession `-iω q` is inside `if (abs_E >= 1e-13)`. -/
def codeOffRHS (E ω T2 xabs : ℝ) (qii qjj q : ℂ) : ℂ :=
  if |E| ≥ 1e-13 then -q / T2 + (-I * ω * q) + I * E * xabs * (qii - qjj) else -q / T2

/-- Correct: `-q/T2 - iωq + iE|x|(q_ii - q_jj)` for every `E`, including `E = 0`. -/
def correctOffRHS (E ω T2 xabs : ℝ) (qii qjj q : ℂ) : ℂ :=
  -q / T2 + (-I * ω * q) + I * E * xabs * (qii - qjj)

theorem sbe11_precession_dropped_at_E0 :
    codeOffRHS 0 1 1 0 0 0 1 ≠ correctOffRHS 0 1 1 0 0 0 1 := by
  unfold codeOffRHS correctOffRHS
  intro h
  have := congrArg Complex.im h
  norm_num at this

/-- SBE-12 code interband current, component `jj` (`bloch_solver_ssbe.f90:545,755-757`):
amplitude `|epdir_jj · d^jj_ij|`; correct amplitude of `d^jj` in the frame of the projected phase
is `|d^jj|` (for phase-aligned components). -/
def codeAmp (e d : ℝ) : ℝ := |e * d|
def correctAmp (d : ℝ) : ℝ := |d|

theorem sbe12_ok_for_axis (d : ℝ) : codeAmp 1 d = correctAmp d := by
  unfold codeAmp correctAmp; simp

theorem sbe12_oblique_bug : codeAmp (1/2) 1 ≠ correctAmp 1 := by
  unfold codeAmp correctAmp; norm_num

/-- SBE-13 Adams–Bashforth coefficients (`bloch_solver_ssbe.f90:782-815`) are the standard
AB4/AB8 ones (reference AB8 numerators over 120960) and are consistent (sum = 1). -/
theorem sbe13_ab4_ok : (55/24 : ℚ) + (-59/24) + 37/24 + (-3/8) = 1 := by norm_num

theorem sbe13_ab4_order2 : (55/24 : ℚ) * 0 + (-59/24) * 1 + 37/24 * 2 + (-3/8) * 3 = -1/2 := by
  norm_num

theorem sbe13_ab8_ok :
    (16083/4480 : ℚ) = 434241/120960 ∧ (-1152169/120960 : ℚ) = -1152169/120960 ∧
    (242653/13440 : ℚ) = 2183877/120960 ∧ (-296053/13440 : ℚ) = -2664477/120960 ∧
    (2102243/120960 : ℚ) = 2102243/120960 ∧ (-115747/13440 : ℚ) = -1041723/120960 ∧
    (32863/13440 : ℚ) = 295767/120960 ∧ (-5257/17280 : ℚ) = -36799/120960 ∧
    (16083/4480 : ℚ) + (-1152169/120960) + 242653/13440 + (-296053/13440) + 2102243/120960
      + (-115747/13440) + 32863/13440 + (-5257/17280) = 1 := by
  norm_num

/-- SBE-14 k-gradient stencil (`common_ssbe.f90:47-85`, `set_bN` 4th column):
`f'(k) ≈ Σ_j c_j (f(k+jh) - f(k-jh))/h` with `c = (4/5,-1/5,4/105,-1/280)`:
exact for k, k³, k⁵, k⁷ (8th order). -/
theorem sbe14_fd_ok :
    2 * ((4/5 : ℚ) * 1 + (-1/5) * 2 + 4/105 * 3 + (-1/280) * 4) = 1 ∧
    (4/5 : ℚ) * 1 + (-1/5) * 2^3 + 4/105 * 3^3 + (-1/280) * 4^3 = 0 ∧
    (4/5 : ℚ) * 1 + (-1/5) * 2^5 + 4/105 * 3^5 + (-1/280) * 4^5 = 0 ∧
    (4/5 : ℚ) * 1 + (-1/5) * 2^7 + 4/105 * 3^7 + (-1/280) * 4^7 = 0 := by
  norm_num

/-- SBE-15 complex k-gradient (`common_ssbe.f90:125,159-175`, `grad_k_array_nb2d_dcomplex`):
the same stencil applied to a complex array `f`, at spacing `h`. -/
def fdStencilC (f : ℤ → ℂ) (h : ℝ) : ℂ :=
  ((4/5 : ℂ) * (f 1 - f (-1)) + (-1/5 : ℂ) * (f 2 - f (-2)) + (4/105 : ℂ) * (f 3 - f (-3))
    + (-1/280 : ℂ) * (f 4 - f (-4))) / h

/-- Code: the complex sum is stored in `real(8) :: w` (Fortran keeps the real part) before it is
written to the complex output. -/
def codeGradC (f : ℤ → ℂ) (h : ℝ) : ℂ := ((fdStencilC f h).re : ℂ)

/-- For `f(k) = i k` the correct gradient is `i`, the code returns `0`. -/
theorem sbe15_correct_val : fdStencilC (fun j => I * (j : ℂ)) 1 = I := by
  unfold fdStencilC; push_cast; ring

theorem sbe15_code_val : codeGradC (fun j => I * (j : ℂ)) 1 = 0 := by
  unfold codeGradC; rw [sbe15_correct_val]; simp

theorem sbe15_real_temp_drops_imag :
    codeGradC (fun j => I * (j : ℂ)) 1 ≠ fdStencilC (fun j => I * (j : ℂ)) 1 := by
  rw [sbe15_code_val, sbe15_correct_val]; exact I_ne_zero.symm

/-- The code agrees with the complex stencil exactly when the stencil value is real. -/
theorem sbe15_ok_iff_real (f : ℤ → ℂ) (h : ℝ) :
    codeGradC f h = fdStencilC f h ↔ (fdStencilC f h).im = 0 := by
  unfold codeGradC
  constructor
  · intro hc; have := congrArg Complex.im hc; simpa using this.symm
  · intro hi; apply Complex.ext <;> simp [hi]

end

end SalmonAudit.SBE.LG
