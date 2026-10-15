import Mathlib

/-!
# Spin–orbit / noncollinear spin (src/so)

Abstractions: a two-component spinor at one grid point `ψ = (ψ1, ψ2) : ℂ × ℂ`;
the spin density matrix is SALMON's `den_mat(is,js) = Σ occ ψ_is conj(ψ_js)`
(`src/so/noncollinear_module.f90:56-57`); the occupation-weighted sum is irrelevant
for the algebra checked here (all formulas are linear in `den_mat`).
-/

namespace SalmonAudit.MISC.SpinOrbit

open Complex ComplexConjugate

noncomputable section

/-- SALMON's density matrix element `den_mat(i,j) = ψ_i * conj ψ_j`. -/
def rho (ψ1 ψ2 : ℂ) : Fin 2 → Fin 2 → ℂ :=
  fun i j => (if i = 0 then ψ1 else ψ2) * conj (if j = 0 then ψ1 else ψ2)

/-- Correct spin expectation `(1/2) ψ† σ_y ψ` with `σ_y = [[0,-i],[i,0]]`. -/
def correctSy (ψ1 ψ2 : ℂ) : ℝ :=
  (1 / 2) * (conj ψ1 * (-I * ψ2) + conj ψ2 * (I * ψ1)).re

/-- Code: `src/so/noncollinear_module.f90:332` (also `:343`, `:384`)
`m(2) = 0.5*dble( -zi*zmat(1,2) + zi*zmat(2,1) )`. -/
def codeSy (ψ1 ψ2 : ℂ) : ℝ :=
  (1 / 2) * (-I * rho ψ1 ψ2 0 1 + I * rho ψ1 ψ2 1 0).re

/-- MISC-09 (BUG): the reported y-magnetization has the opposite sign. -/
theorem Sy_bug (ψ1 ψ2 : ℂ) : codeSy ψ1 ψ2 = - correctSy ψ1 ψ2 := by
  unfold codeSy correctSy rho
  simp [Complex.mul_re, Complex.mul_im]
  ring

/-- Concrete: spinor pointing along +y, `ψ = (1, i)`: correct `S_y = +1`, code `-1`. -/
theorem Sy_counterexample : codeSy 1 I ≠ correctSy 1 I := by
  rw [Sy_bug]
  unfold correctSy
  simp [Complex.mul_re]
  norm_num

/-- Correct `(1/2) ψ† σ_x ψ` and `(1/2) ψ† σ_z ψ`. -/
def correctSx (ψ1 ψ2 : ℂ) : ℝ := (1 / 2) * (conj ψ1 * ψ2 + conj ψ2 * ψ1).re
def correctSz (ψ1 ψ2 : ℂ) : ℝ := (1 / 2) * (conj ψ1 * ψ1 - conj ψ2 * ψ2).re

/-- Code `noncollinear_module.f90:331,333`. -/
def codeSx (ψ1 ψ2 : ℂ) : ℝ := (1 / 2) * (rho ψ1 ψ2 0 1 + rho ψ1 ψ2 1 0).re
def codeSz (ψ1 ψ2 : ℂ) : ℝ := (1 / 2) * (rho ψ1 ψ2 0 0 - rho ψ1 ψ2 1 1).re

/-- MISC-10 (OK): x and z components are right. -/
theorem Sx_ok (ψ1 ψ2 : ℂ) : codeSx ψ1 ψ2 = correctSx ψ1 ψ2 := by
  unfold codeSx correctSx rho; simp [Complex.mul_re]; try ring
theorem Sz_ok (ψ1 ψ2 : ℂ) : codeSz ψ1 ψ2 = correctSz ψ1 ψ2 := by
  unfold codeSz correctSz rho; simp [Complex.mul_re]; try ring

/-- Code spin-current σ_y channel, `noncollinear_module.f90:450`:
`sig(:,2) = -zi*( conjg(p1)*g2 - conjg(p2)*g1 )`. -/
def codeJy (p1 p2 g1 g2 : ℂ) : ℂ := -I * (conj p1 * g2 - conj p2 * g1)
/-- Correct `ψ† σ_y (g)`. -/
def correctJy (p1 p2 g1 g2 : ℂ) : ℂ := conj p1 * (-I * g2) + conj p2 * (I * g1)

/-- MISC-11 (OK): the spin current uses the correct σ_y (so it is inconsistent with MISC-09). -/
theorem Jy_ok (p1 p2 g1 g2 : ℂ) : codeJy p1 p2 g1 g2 = correctJy p1 p2 g1 g2 := by
  unfold codeJy correctJy; ring

/-- Code XC matrix `noncollinear_module.f90:207-210`:
`V11 = v0 + v1 cosθ`, `V12 = v1 e^{-iφ} sinθ`, `V21 = v1 e^{iφ} sinθ`, `V22 = v0 - v1 cosθ`. -/
def codeVxc (v0 v1 θ φ : ℝ) : Fin 2 → Fin 2 → ℂ :=
  ![![v0 + v1 * Real.cos θ, v1 * (Real.cos φ - Real.sin φ * I) * Real.sin θ],
    ![v1 * (Real.cos φ + Real.sin φ * I) * Real.sin θ, v0 - v1 * Real.cos θ]]

/-- Correct `v0 I + v1 (n̂·σ)`, `n̂ = (sinθ cosφ, sinθ sinφ, cosθ)`. -/
def correctVxc (v0 v1 θ φ : ℝ) : Fin 2 → Fin 2 → ℂ :=
  let nx := Real.sin θ * Real.cos φ
  let ny := Real.sin θ * Real.sin φ
  let nz := Real.cos θ
  ![![v0 + v1 * nz, v1 * (nx - ny * I)],
    ![v1 * (nx + ny * I), v0 - v1 * nz]]

/-- MISC-12 (OK): the rotated XC potential is `v0 + v1 n̂·σ`. -/
theorem Vxc_ok (v0 v1 θ φ : ℝ) : codeVxc v0 v1 θ φ = correctVxc v0 v1 θ φ := by
  funext i j
  fin_cases i <;> fin_cases j <;> simp [codeVxc, correctVxc] <;> ring

/-- Code local-frame densities `noncollinear_module.f90:136-143` with `(θ,φ)` given:
`rho1 = n/2 + (ρ11-ρ22)/2 cosθ + Re(ρ12 e^{iφ}) sinθ`. Here `ρ12 = (mx - i my)/2`,
`ρ11-ρ22 = mz` (correct convention `m = Tr σ ρ`). -/
def codeRho1 (n mx my mz θ φ : ℝ) : ℝ :=
  n / 2 + mz / 2 * Real.cos θ + (mx / 2 * Real.cos φ - (-my / 2) * Real.sin φ) * Real.sin θ

/-- Correct: density projected on `+n̂`: `(n + n̂·m)/2`. -/
def correctRho1 (n mx my mz θ φ : ℝ) : ℝ :=
  (n + (Real.sin θ * Real.cos φ * mx + Real.sin θ * Real.sin φ * my + Real.cos θ * mz)) / 2

/-- MISC-12b (OK): `rot_dm_noncollinear` builds `(n ± n̂·m)/2` consistently with
`rot_vxc_noncollinear` (any branch of `atan` gives the same consistent pair). -/
theorem rho1_ok (n mx my mz θ φ : ℝ) :
    codeRho1 n mx my mz θ φ = correctRho1 n mx my mz θ φ := by
  unfold codeRho1 correctRho1; ring

/-- MISC-13 (BUG, `method_mixing='simple_dm'`): `simple_mixing_so`
(`noncollinear_module.f90:292-299`) overwrites the rotated densities by the diagonal
`ρ11, ρ22`, so XC sees the splitting `ρ11 - ρ22 = m_z`, while `rot_vxc` rotates with
the local axis whose splitting should be `|m|`. -/
def codeSplit (mz : ℝ) : ℝ := mz
/-- Correct: XC input = eigenvalues `n/2 ± |m|/2` of the 2×2 density matrix, so the
polarization `ρ₊ - ρ₋` is `|m|` (no claim about how the nonlinear XC then responds). -/
def correctSplit (mx my mz : ℝ) : ℝ := Real.sqrt (mx ^ 2 + my ^ 2 + mz ^ 2)

/-- The eigenvalues of `ρ = (n + m·σ)/2` are `n/2 ± |m|/2`: check `ρ₊+ρ₋ = Tr ρ` and
`ρ₊ρ₋ = det ρ` with `ρ11 = (n+mz)/2`, `ρ22 = (n-mz)/2`, `|ρ12|² = (mx²+my²)/4`. -/
theorem eig_ok (n mx my mz : ℝ) :
    let r := Real.sqrt (mx ^ 2 + my ^ 2 + mz ^ 2)
    (n / 2 + r / 2) + (n / 2 - r / 2) = (n + mz) / 2 + (n - mz) / 2 ∧
    (n / 2 + r / 2) * (n / 2 - r / 2) = (n + mz) / 2 * ((n - mz) / 2) - (mx ^ 2 + my ^ 2) / 4 := by
  intro r
  have hr : r ^ 2 = mx ^ 2 + my ^ 2 + mz ^ 2 := Real.sq_sqrt (by positivity)
  constructor
  · ring
  · nlinarith [hr]

/-- Concrete: `m = (0.8, 0, 0.6)`, `|m| = 1`; code uses 0.6. -/
theorem split_counterexample : codeSplit (3 / 5) ≠ correctSplit (4 / 5) 0 (3 / 5) := by
  unfold codeSplit correctSplit
  have : (4 / 5 : ℝ) ^ 2 + 0 ^ 2 + (3 / 5) ^ 2 = 1 := by norm_num
  rw [this, Real.sqrt_one]; norm_num

/-- Code Clebsch–Gordan coefficients, `src/so/prep_pp_so.f90:279-311`
(`l` = `ll`, `mj` half-integer): `j=l+1/2`: (↑: √((l+mj+1/2)/(2l+1)), ↓: √((l-mj+1/2)/(2l+1)));
`j=l-1/2`: (↑: √((l-mj+1/2)/(2l+1)), ↓: -√((l+mj+1/2)/(2l+1))). -/
def cgPlus (l mj : ℝ) : ℝ × ℝ :=
  (Real.sqrt ((l + mj + 1 / 2) / (2 * l + 1)), Real.sqrt ((l - mj + 1 / 2) / (2 * l + 1)))
def cgMinus (l mj : ℝ) : ℝ × ℝ :=
  (Real.sqrt ((l - mj + 1 / 2) / (2 * l + 1)), - Real.sqrt ((l + mj + 1 / 2) / (2 * l + 1)))

/-- MISC-07 (OK): the coefficients are normalized and the two j-channels orthogonal
(standard Condon–Shortley CG values for coupling `l ⊗ 1/2`). -/
theorem cg_ok (l mj : ℝ) (hl : 0 ≤ l) (h1 : 0 ≤ l + mj + 1 / 2) (h2 : 0 ≤ l - mj + 1 / 2) :
    (cgPlus l mj).1 ^ 2 + (cgPlus l mj).2 ^ 2 = 1 ∧
    (cgMinus l mj).1 ^ 2 + (cgMinus l mj).2 ^ 2 = 1 ∧
    (cgPlus l mj).1 * (cgMinus l mj).1 + (cgPlus l mj).2 * (cgMinus l mj).2 = 0 := by
  have hd : 0 < 2 * l + 1 := by linarith
  have e1 : Real.sqrt ((l + mj + 1 / 2) / (2 * l + 1)) ^ 2 = (l + mj + 1 / 2) / (2 * l + 1) :=
    Real.sq_sqrt (div_nonneg h1 hd.le)
  have e2 : Real.sqrt ((l - mj + 1 / 2) / (2 * l + 1)) ^ 2 = (l - mj + 1 / 2) / (2 * l + 1) :=
    Real.sq_sqrt (div_nonneg h2 hd.le)
  refine ⟨?_, ?_, ?_⟩
  · simp only [cgPlus]; rw [e1, e2]; field_simp; ring
  · simp only [cgMinus]; rw [neg_sq, e1, e2]; field_simp; ring
  · simp only [cgPlus, cgMinus]; ring

/-- Code `zylm` from the real harmonics `Ylm(m)`, `Ylm(-m)` (`prep_pp_so.f90:452-458`),
for `m>0`: `zylm(m) = (Y_m + i Y_{-m})/√2`, `zylm(-m) = (-1)^m (Y_m - i Y_{-m})/√2`. -/
def zylmPos (Ym Ymm : ℝ) : ℂ := (Real.sqrt (1 / 2) : ℂ) * (Ym + Ymm * I)
def zylmNeg (m : ℕ) (Ym Ymm : ℝ) : ℂ := (-1) ^ m * (Real.sqrt (1 / 2) : ℂ) * (Ym - Ymm * I)

/-- MISC-08 (OK): `zylm(-m) = (-1)^m conj(zylm(m))` — the Condon–Shortley relation,
given SALMON's real harmonics (`salmon_math.f90` `Ylm`: e.g. `Y_{1,1} = -x`, `Y_{1,-1} = -y`,
i.e. `√2 Re/Im` of the CS `Y_{l|m|}`), so `zylm` is the CS complex harmonic. -/
theorem zylm_ok (m : ℕ) (Ym Ymm : ℝ) :
    zylmNeg m Ym Ymm = (-1) ^ m * conj (zylmPos Ym Ymm) := by
  unfold zylmNeg zylmPos
  simp [map_mul, Complex.conj_ofReal]
  ring

end

end SalmonAudit.MISC.SpinOrbit
