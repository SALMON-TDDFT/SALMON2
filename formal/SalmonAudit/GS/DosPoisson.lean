import Mathlib

/-!
# GS: DOS / PDOS output and the isolated 'ft' Poisson solver (SALMON v2.3.0)
-/

noncomputable section
namespace SalmonAudit.GS

open Real MeasureTheory

/-- Code Gaussian DOS kernel (`src/io/write.f90:1653-1657`):
`fk·exp(−w²/(2σ²))`, `fk = 2/(√(2π) σ)`, summed with weight wtk for every (io,ik,is). -/
def dosGaussCode (σ w : ℝ) : ℝ := 2 / (Real.sqrt (2 * π) * σ) * Real.exp (-(0.5 / σ ^ 2) * w ^ 2)

/-- The code kernel integrates to 2 for every orbital and every spin channel. -/
theorem dosGaussCode_integral (σ : ℝ) (hσ : 0 < σ) : ∫ w, dosGaussCode σ w = 2 := by
  unfold dosGaussCode
  rw [integral_const_mul, integral_gaussian]
  have hpos : 0 < Real.sqrt (2 * π) * σ := by positivity
  rw [show π / (0.5 / σ ^ 2) = (Real.sqrt (2 * π) * σ) ^ 2 by
        rw [mul_pow, Real.sq_sqrt (by positivity)]; field_simp; ring]
  rw [Real.sqrt_sq hpos.le]
  field_simp

/-- Correct spin weight per orbital: 2 for `nspin=1`, 1 per channel for `nspin=2`
(cf. `wspin` in `src/gs/occupation.f90:38-42`). -/
def wspin (nspin : ℕ) : ℝ := if nspin = 1 then 2 else 1

/-- GS-25 (BUG, nspin=2): each spin-resolved DOS column integrates to 2 per orbital
instead of 1. (OK for nspin=1.) Same prefactor in the Lorentzian branch and in PDOS. -/
theorem dos_spin_factor (σ : ℝ) (hσ : 0 < σ) :
    (∫ w, dosGaussCode σ w) = wspin 1 ∧ (∫ w, dosGaussCode σ w) = 2 * wspin 2 := by
  rw [dosGaussCode_integral σ hσ]; simp [wspin]

/-- GS-25 (BUG, SOC): with `yn_spinorbit='y'` both spin indices carry the same spinor eigenvalue
(`src/common/total_energy.f90:693-695`); DOS sums `is=1,2` each with integral 2, so each spinor state
(which holds 1 electron) contributes 2+2 = 4. -/
theorem dos_soc_factor (σ : ℝ) (hσ : 0 < σ) :
    (∫ w, dosGaussCode σ w) + (∫ w, dosGaussCode σ w) = 4 * 1 := by
  rw [dosGaussCode_integral σ hσ]; norm_num

/-- PDOS angular-radial function, code (`src/io/write.f90:1776-1777`), default `yn_psmask='n'` path.
The stored table is already `upp_f(r,l) = u_l(r)·√((2l+1)/4π)/r^{l+1}` (`src/atom/pp/input_pp.f90:1143-1148`,
copied at `:234`). The code then uses channel `lref` for every L, multiplies **again** by
`√((2 lref+1)/4π)/r^{lref+1}`, and by `Ylm(r̂,L,m) = √(4π/(2L+1)) Y⁰_Lm(r̂)` (unit vector,
`src/math/salmon_math.f90:386-388`). Returned: coefficient of `Y⁰_Lm(r̂)`; `u` = raw radial u_l(r). -/
def pdosRadialCode (u : ℕ → ℝ → ℝ) (lref L : ℕ) (r : ℝ) : ℝ :=
  (u lref r * Real.sqrt ((2 * lref + 1) / (4 * π)) / r ^ (lref + 1))
    / r ^ (lref + 1) * Real.sqrt ((2 * lref + 1) / (4 * π)) * Real.sqrt (4 * π / (2 * L + 1))

/-- Correct: `φ_Lm = (u_L/r) Y⁰_Lm = upp_f(r,L) · Ylm(r⃗,L,m)` (solid harmonic, non-unit vector). -/
def pdosRadialCorrect (u : ℕ → ℝ → ℝ) (L : ℕ) (r : ℝ) : ℝ := u L r / r

/-- GS-26 (BUG): counterexample — lref = 0, L = 0, u ≡ 1, r = 1: code gives `1/(4π)·√(4π)`·…
i.e. `√(1/4π)·√(1/4π)·√(4π) = 1/√(4π) ≠ 1` (double normalisation even when L = lref). -/
theorem pdosRadial_ne :
    pdosRadialCode (fun _ _ => (1:ℝ)) 0 0 1 ≠ pdosRadialCorrect (fun _ _ => (1:ℝ)) 0 1 := by
  unfold pdosRadialCode pdosRadialCorrect
  norm_num
  intro h
  have hpi : 0 < π := Real.pi_pos
  have hq : 0 < Real.sqrt π := Real.sqrt_pos.mpr hpi
  have hsq : Real.sqrt π ^ 2 = π := Real.sq_sqrt hpi.le
  field_simp at h
  nlinarith [Real.pi_gt_three]

/-- GS-27 (BUG): PDOS overlap stored in a `real(8)` array (`rbox_pdos`, `src/io/write.f90:1732,1782`)
so for complex orbitals only `(Re⟨φ|ψ⟩)²` is accumulated instead of `|⟨φ|ψ⟩|²`. Witness ⟨φ|ψ⟩ = i. -/
theorem pdos_realpart_ne : (Complex.I.re) ^ 2 ≠ Complex.normSq Complex.I := by
  simp

/-- GS-28 (BUG): PDOS sums over k without `wtk` (`src/io/write.f90:1795-1806`), DOS uses `wtk`
(`:1651`). Witness: two k-points, weights ½, unit overlaps: code 2, correct 1. -/
theorem pdos_noWtk_ne :
    (∑ _k : Fin 2, (1:ℝ)) ≠ ∑ _k : Fin 2, (1/2 : ℝ) * 1 := by
  simp

/-- GS-29 (SUSPECT). Isolated `method_poisson='ft'` solves the *periodic* problem on the doubled box
with `coef = 4π/G²` and `coef(G=0)=0` (`src/common/initialization.f90:1086-1093`,
`src/poisson/poisson_isolated.f90:1161`): the resulting potential has zero cell average.
The free-space Hartree potential of a non-negative, non-zero density is strictly positive
everywhere, so the two cannot coincide (abstracted to a finite grid; the zero mean is over the doubled cell,
so the deviation is not necessarily a uniform shift inside the physical box). -/
theorem isoFT_zeroMean_ne {ι : Type*} [Fintype ι] [Nonempty ι] (Vft Vtrue : ι → ℝ)
    (hmean : ∑ i, Vft i = 0) (hpos : ∀ i, 0 < Vtrue i) : Vft ≠ Vtrue := by
  intro h
  subst h
  have : 0 < ∑ i, Vft i := Finset.sum_pos (fun i _ => hpos i) Finset.univ_nonempty
  linarith

/-- GS-29 evidence: SALMON's own test references (testsuites/101_C2H2_gs vs 181_C2H2_gs_ft,
identical inputs except `method_poisson`) differ by a uniform shift of the eigenvalues:
io=1: −31.1347 vs −18.4139 eV, io=5: −20.0607 vs −7.3403 eV. -/
theorem isoFT_testsuite_shift :
    (-31.1347 : ℝ) - (-18.4139) = -12.7208 ∧ (-20.0607 : ℝ) - (-7.3403) = -12.7204 := by
  norm_num

end SalmonAudit.GS
end
