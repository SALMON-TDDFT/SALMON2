import Mathlib

/-!
# Macroscopic current and induced field (SALMON v2.3.0, `src/common/density_matrix.f90`,
`src/rt/em_field.f90`)

Abstractions: one Cartesian direction of a periodic grid `ZMod N`; stencil = finite set of
offsets `j` with real coefficients `c j` (the code's `nabt`, already divided by `h`);
grid sums without `Hvol` exactly as the code accumulates them.
-/

noncomputable section
namespace SalmonAudit.RT

open Complex

/-! ## RT-15 local (stencil) part of ⟨p⟩ -/

/-- `X_j = Σ_x conj ψ(x) ψ(x+j)` = the code's `xtmp` contribution of offset `j`
(`density_matrix.f90:474-478`). -/
def Xj {N : ℕ} [NeZero N] (ψ : ZMod N → ℂ) (j : ZMod N) : ℂ :=
  ∑ x, (starRingEnd ℂ) (ψ x) * ψ (x + j)

/-- Reindexing on the periodic grid: `Σ_x conj ψ(x) ψ(x-j) = conj X_j`. -/
lemma sum_back {N : ℕ} [NeZero N] (ψ : ZMod N → ℂ) (j : ZMod N) :
    ∑ x, (starRingEnd ℂ) (ψ x) * ψ (x - j) = (starRingEnd ℂ) (Xj ψ j) := by
  unfold Xj
  rw [map_sum]
  simp only [map_mul, Complex.conj_conj]
  exact Fintype.sum_equiv (Equiv.subRight j) _ _ (fun x => by simp [mul_comm])

/-- Correct: `⟨ψ|p|ψ⟩ = Σ_x conj ψ(x) (-i) ∇ψ(x)`, `∇ψ(x) = Σ_j c_j (ψ(x+j) - ψ(x-j))`
(antisymmetric central-difference stencil). -/
def pCorrect {N : ℕ} [NeZero N] (S : Finset (ZMod N)) (c : ZMod N → ℝ) (ψ : ZMod N → ℂ) : ℂ :=
  ∑ x, (starRingEnd ℂ) (ψ x) * (-I * ∑ j ∈ S, (c j : ℂ) * (ψ (x + j) - ψ (x - j)))

/-- Code `density_matrix.f90:512`: `j2 = aimag(2*xtmp)` with `xtmp = Σ_j c_j X_j`. -/
def pCode {N : ℕ} [NeZero N] (S : Finset (ZMod N)) (c : ZMod N → ℝ) (ψ : ZMod N → ℂ) : ℝ :=
  (2 * ∑ j ∈ S, (c j : ℂ) * Xj ψ j).im

theorem p_stencil_ok {N : ℕ} [NeZero N] (S : Finset (ZMod N)) (c : ZMod N → ℝ) (ψ : ZMod N → ℂ) :
    pCorrect S c ψ = (pCode S c ψ : ℂ) := by
  have h1 : pCorrect S c ψ = -I * ∑ j ∈ S, (c j : ℂ) * (Xj ψ j - (starRingEnd ℂ) (Xj ψ j)) := by
    unfold pCorrect
    simp_rw [Finset.mul_sum, mul_sub]
    rw [Finset.sum_comm]
    refine Finset.sum_congr rfl (fun j _ => ?_)
    rw [← sum_back ψ j]
    simp only [Xj, Finset.mul_sum, ← Finset.sum_sub_distrib]
    refine Finset.sum_congr rfl (fun x _ => ?_)
    ring
  rw [h1, pCode]
  apply Complex.ext
  all_goals simp [Finset.mul_sum, Complex.mul_re, Complex.mul_im, Complex.sub_re, Complex.sub_im]
  all_goals (refine Finset.sum_congr rfl (fun x _ => ?_); ring)

/-! ## RT-16 nonlocal commutator part `⟨ψ| i[V_NL, r] |ψ⟩` -/

/-- Code `density_matrix.f90:518-575`: `jw = 2 Σ aimag(conj(uVpsi_r) * uVpsi)` with
`uVpsi = D⟨u|ψ⟩` (`rinv_uvu = inorm*Hvol`), `uVpsi_r = ⟨u|rψ⟩`. One projector: `a = ⟨u|ψ⟩`, `b = ⟨u|rψ⟩`. -/
def nlCode (D : ℝ) (a b : ℂ) : ℝ := 2 * ((starRingEnd ℂ) b * ((D : ℂ) * a)).im

/-- Correct: velocity `v = i[H,r]` ⊃ `i[V_NL, r]`, `V_NL = |u⟩D⟨u|`:
`⟨ψ|i(V r - r V)|ψ⟩ = i D (conj a · b - conj b · a)`. -/
def nlCorrect (D : ℝ) (a b : ℂ) : ℂ :=
  I * (D : ℂ) * ((starRingEnd ℂ) a * b - (starRingEnd ℂ) b * a)

theorem nl_current_ok (D : ℝ) (a b : ℂ) : nlCorrect D a b = (nlCode D a b : ℂ) := by
  apply Complex.ext <;> simp [nlCorrect, nlCode, Complex.mul_re, Complex.mul_im] <;> ring_nf <;> simp

/-! ## RT-17 normalisation of the macroscopic current -/

/-- Code `density_matrix.f90:426`: `curr = Σ / ngrid` with grid sums lacking `Hvol`.
Correct: `(1/Ω) Hvol Σ` with `Ω = ngrid·Hvol`. -/
theorem curr_norm_ok (S Hvol : ℝ) (ngrid : ℕ) (hH : Hvol ≠ 0) (hn : (ngrid : ℝ) ≠ 0) :
    S / ngrid = (1 / ((ngrid : ℝ) * Hvol)) * (Hvol * S) := by
  field_simp

/-! ## RT-18 longitudinal geometry: `Ac_ind` update -/

/-- Code `em_field.f90:54`: `Ac_ind(t+dt) = 2Ac_ind(t) - Ac_ind(t-dt) - 4π curr dt²`, `curr` =
electron *matter* current (`calc_current` returns `Σ occ ⟨p+k+A⟩/Ω`, header "Matter current"). -/
def loCode (A0 A1 curr dt : ℝ) : ℝ := 2 * A0 - A1 - 4 * Real.pi * curr * dt ^ 2

/-- Correct: depolarisation `dE_ind/dt = -4π J_el`, `J_el = -curr`, `E_ind = -dAc_ind/dt`
⇒ `d²Ac_ind/dt² = -4π curr`; Störmer–Verlet `A(t+h) = 2A(t) - A(t-h) + h² A''`. -/
def loCorrect (A0 A1 curr dt : ℝ) : ℝ :=
  let Jel := -curr
  let Add := 4 * Real.pi * Jel
  2 * A0 - A1 + dt ^ 2 * Add

theorem lo_ok (A0 A1 curr dt : ℝ) : loCode A0 A1 curr dt = loCorrect A0 A1 curr dt := by
  unfold loCode loCorrect; ring

/-! ## RT-19 spin: total current `curr = Σ_σ curr_σ`, occupation already contains degeneracy -/

/-- Code `em_field.f90:46-50` + `rocc` = 2 for `nspin=1`: same physical current either way. -/
theorem spin_sum_ok (j : ℝ) : (2 * j : ℝ) = j + j := by ring

/-! ## RT-32 BUG: stale nonlocal projector phase in `calc_current` -/

/-- Nonlocal current term of one projector on a two-point model (`x = 0, 1`, `u = ψ = 1`, `D = 1`)
with the cached phase factor `z = e^{iqx}` at `x = 1` (`zekr_uV = conj(ekr)·uV`,
`hamiltonian.f90:1273`; used in `density_matrix.f90:554-557`): `a = 1 + z`, `b = z` (`x·ψ` kills x = 0). -/
def nlTwoPoint (z : ℂ) : ℝ := nlCode 1 (1 + z) z

/-- The nonlocal term depends on the phase `q = k + A` (here `q = 0` vs `q = π/2`, `z = 1` vs `z = i`),
so evaluating it with the cached `A(itt-1/2)` (middlepoint, `time_evolution_step.f90:124,130`) or
`A(itt+1)` (aetrs, `:508-510`) instead of `A(itt)` (set at `:285`, no `update_kvector_nonlocalpt`
before `calc_current` at `:286`) changes the result. -/
theorem nl_phase_matters : nlTwoPoint 1 ≠ nlTwoPoint Complex.I := by
  unfold nlTwoPoint nlCode
  norm_num [Complex.mul_im, Complex.mul_re]

/-- Code: phase uses `A_mid = (A_n + A_{n-1})/2` while the local term uses `A_n`; they differ
whenever `A_n ≠ A_{n-1}` (i.e. `E ≠ 0`): `A_n - A_mid = (A_n - A_{n-1})/2 = -E Δt/2`. -/
theorem stale_phase_offset (An An1 : ℝ) : An - (An + An1) / 2 = (An - An1) / 2 := by ring

end SalmonAudit.RT
end
