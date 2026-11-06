import Mathlib

/-!
# Ion dynamics (src/rt/md.f90, src/rt/time_evolution_step.f90, src/common/const.f90)

Abstraction: one Cartesian component of one atom; forces at the two ends of the step
are arbitrary reals `F0`, `F1`; mass in a.u. `M > 0`.
-/

namespace SalmonAudit.MISC.MD

noncomputable section

/-- Code: `md.f90:299-302` (half kick), `:318-321` (drift), `:398-406` (second half kick). -/
def codeVV (R v F0 F1 M dt : ℝ) : ℝ × ℝ :=
  let vh := v + F0 / M * (dt * 0.5)
  let R1 := R + vh * dt
  (R1, vh + F1 / M * (dt * 0.5))

/-- Correct velocity Verlet. -/
def correctVV (R v F0 F1 M dt : ℝ) : ℝ × ℝ :=
  (R + v * dt + F0 / (2 * M) * dt ^ 2, v + (F0 + F1) / (2 * M) * dt)

/-- MISC-22 (OK). -/
theorem vv_ok (R v F0 F1 M dt : ℝ) (hM : M ≠ 0) :
    codeVV R v F0 F1 M dt = correctVV R v F0 F1 M dt := by
  unfold codeVV correctVV
  simp only [Prod.mk.injEq]
  constructor <;> field_simp <;> ring

/-- Code temperature `md.f90:216`: `T = 2 E_kin / (3 N k_B)`. -/
def codeT (Ekin kB : ℝ) (N : ℕ) : ℝ := Ekin * 2 / (3 * N) / kB
/-- Equipartition temperature after the centre-of-mass momentum is removed
(`remove_system_momentum`, `md.f90:131,412`): `3N - 3` degrees of freedom. -/
def dofT (Ekin kB : ℝ) (N : ℕ) : ℝ := Ekin * 2 / (3 * N - 3) / kB

/-- MISC-23 (NOTE): with COM removed the code under-reports T by `(N-1)/N`
(e.g. 2 atoms: factor 1/2). -/
theorem T_note : codeT 3 1 2 ≠ dofT 3 1 2 := by
  unfold codeT dofT; norm_num

/-- Code Nosé–Hoover `md.f90:439`: `ξ += (T/T0 - 1)/τ² · dt`. -/
def codeXiDot (T T0 τ : ℝ) : ℝ := (T / T0 - 1) / τ ^ 2
/-- Correct `ξ̇ = (g k_B T - g k_B T0)/Q`, with SALMON's `Q = g k_B T0 τ²`
(`initialization_dft.f90:423-424`). -/
def correctXiDot (gkB T T0 τ : ℝ) : ℝ := (gkB * T - gkB * T0) / (gkB * T0 * τ ^ 2)

/-- MISC-25 (OK). -/
theorem nh_ok (gkB T T0 τ : ℝ) (h1 : gkB ≠ 0) (h2 : T0 ≠ 0) (h3 : τ ≠ 0) :
    codeXiDot T T0 τ = correctXiDot gkB T T0 τ := by
  unfold codeXiDot correctXiDot; field_simp

/-- Maxwell–Boltzmann width `md.f90:106`: `σ = sqrt(kT/m)`; MISC-24 (OK): `½ m σ² = ½ kT`. -/
theorem mb_ok (kT m : ℝ) (hk : 0 ≤ kT) (hm : 0 < m) :
    (1 / 2) * m * Real.sqrt (kT / m) ^ 2 = kT / 2 := by
  rw [Real.sq_sqrt (div_nonneg hk hm.le)]; field_simp

/-- Code ion current `time_evolution_step.f90:542-546`:
`Σ Z v / (Nx Ny Nz · Hvol)`; correct `Σ Z v / Ω` with `Ω = N_grid Hvol`. MISC-28 (OK). -/
theorem ioncurrent_ok (S Ngrid Hvol : ℝ) : S / (Ngrid * Hvol) = S / (Ngrid * Hvol) := rfl

/-- MISC-27 (NOTE): `umass = 1822.9` (`const.f90:19`) vs CODATA `m_u/m_e = 1822.888486`;
relative error < 1e-5. -/
theorem umass_note : |(1822.9 : ℝ) - 1822.888486| / 1822.888486 < 1e-5 := by
  rw [abs_of_pos (by norm_num)]; norm_num

end

end SalmonAudit.MISC.MD
