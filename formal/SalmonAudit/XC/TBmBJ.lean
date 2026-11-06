import Mathlib

/-!
# XC audit: builtin TB-mBJ (`src/xc/builtin_tbmbj.f90`, driver `salmon_xc.f90:1403-1447`)

Conventions in the code: `tau` from `calc_tau` is `½ Σ f |∇ψ|²` (`salmon_xc.f90:633`);
the driver halves ρ, ∇ρ, ∇²ρ, τ, j to spin densities (`salmon_xc.f90:1417-1427`).
-/

namespace SalmonAudit.XC

/-! ## XC-09  Becke–Roussel curvature `D_σ` (`builtin_tbmbj.f90:80`) -/

/-- code: `D = 2 t - ¼|∇ρ|²/ρ`, with `t = ½Σ|∇ψ|²` (half-convention τ). -/
noncomputable def brDCode (t g2 ρ : ℝ) : ℝ := 2 * t - 0.25 * g2 / ρ

/-- BR89: `D = τ_BR - ¼|∇ρ|²/ρ` with `τ_BR = Σ|∇ψ|² = 2t`; it must vanish for one orbital,
where `t = |∇ρ|²/(8ρ)` (ψ = √ρ). -/
theorem br_D_single_orbital_ok (g2 ρ : ℝ) (hρ : ρ ≠ 0) : brDCode (g2 / (8 * ρ)) g2 ρ = 0 := by
  unfold brDCode; field_simp; ring

/-- BR equation solved by `BR_Newton`: `x e^{-2x/3}/rhs - (x-2) = 0` ⇔ `x e^{-2x/3}/(x-2) = rhs`. -/
theorem br_equation_ok (x rhs : ℝ) (hr : rhs ≠ 0) (hx : x - 2 ≠ 0) (hx0 : x ≠ 0) :
    x * Real.exp (-2 / 3 * x) / rhs - (x - 2) = 0 ↔ x * Real.exp (-2 / 3 * x) / (x - 2) = rhs := by
  constructor
  · intro h; field_simp at h ⊢; linarith
  · intro h; rw [← h]; field_simp; ring

/-- gauge invariance of `t - j²/(2ρ)` (`builtin_tbmbj.f90:79`) under `k → k + A` (one Cartesian
component; `τ' = τ + A j + ½ρA²`, `j' = j + ρ A`).  The kinetic-energy combination TB-mBJ uses
is therefore unaffected by `calc_tau` using `vec_k` without `A` in RT. -/
theorem tbmbj_tau_gauge_invariant_ok (t j ρ A : ℝ) (hρ : ρ ≠ 0) :
    (t + A * j + ρ * A ^ 2 / 2) - (j + ρ * A) ^ 2 / ρ / 2 = t - j ^ 2 / ρ / 2 := by
  field_simp; ring

/-- With NLCC the code divides by `n = ρ_v + ρ_c` while `τ, j` are valence (`salmon_xc.f90:1420`),
so the subtraction is no longer gauge invariant: counterexample `t=0, j=0, ρ_v=1, ρ_c=1, A=1`. -/
theorem tbmbj_tau_gauge_nlcc_counterexample :
    let t : ℝ := 0; let j : ℝ := 0; let ρv : ℝ := 1; let n : ℝ := 2; let A : ℝ := 1
    (t + A * j + ρv * A ^ 2 / 2) - (j + ρv * A) ^ 2 / n / 2 ≠ t - j ^ 2 / n / 2 := by
  norm_num

/-! ## XC-10  BUG: correlation potential with NLCC (`builtin_tbmbj.f90:88-95`)

`rho` passed to `exc_cor_tbmbj` is the valence density (`salmon_xc.f90:1416`), while
`rho_s = ρ/2 + ρ_c/2` includes the partial core.  PW correlation is evaluated at
`n = ρ + ρ_c`, but the potential uses `ε_c(n) + ρ ∂ε_c/∂ρ_a` instead of `ε_c(n) + n ∂ε_c/∂n`. -/

/-- code: `ec + rho*dec_drhoa` with `rho` = valence. (At ζ=0, `dec_drhoa = dε_c/dn`.) -/
def tbVcCode (ec dec ρv _ρc : ℝ) : ℝ := ec + ρv * dec

/-- correct: `d(n ε_c(n))/dn = ε_c + n ε_c'`, `n = ρ_v + ρ_c`. -/
def tbVcCorrect (ec dec ρv ρc : ℝ) : ℝ := ec + (ρv + ρc) * dec

theorem tbmbj_vc_correct_is_derivative {ε : ℝ → ℝ} {d n : ℝ} (h : HasDerivAt ε d n) :
    HasDerivAt (fun x => x * ε x) (ε n + n * d) n := by
  exact ((hasDerivAt_id n).mul h).congr_deriv (by simp)

theorem tbmbj_vc_nlcc_bug (ec dec ρv ρc : ℝ) :
    tbVcCode ec dec ρv ρc = tbVcCorrect ec dec ρv ρc - ρc * dec := by
  unfold tbVcCode tbVcCorrect; ring

theorem tbmbj_vc_nlcc_counterexample : tbVcCode (-1) (1/10) 1 1 ≠ tbVcCorrect (-1) (1/10) 1 1 := by
  unfold tbVcCode tbVcCorrect; norm_num

/-! ## XC-11  BUG: the TB c-parameter is averaged over the MPI-local subdomain only

`builtin_tbmbj.f90:57-64`: `c = α + β sqrt( (1/nl) Σ_{i∈local} |∇ρ|/ρ )`, `nl` = local grid
points (`salmon_xc.f90:1099`, `mg%num`), with no reduction over `info%icomm_r`.
TB09 defines `g = (1/V_cell) ∫_cell |∇ρ|/ρ`.  Model: two equal-size subdomains. -/

/-- code value on rank r: average over its own points. -/
noncomputable def tbCCode (α β : ℝ) (local_vals : List ℝ) : ℝ :=
  α + β * Real.sqrt (local_vals.sum / local_vals.length)

/-- correct: average over all points of the cell. -/
noncomputable def tbCCorrect (α β : ℝ) (all_vals : List ℝ) : ℝ :=
  α + β * Real.sqrt (all_vals.sum / all_vals.length)

theorem tbmbj_c_rank_local_counterexample :
    tbCCode (-0.012) 1.023 [8, 8] ≠ tbCCorrect (-0.012) 1.023 [8, 8, 0, 0] := by
  unfold tbCCode tbCCorrect
  norm_num
  intro h
  have h8 : Real.sqrt 8 = 2 := by linarith
  have := Real.sq_sqrt (show (0:ℝ) ≤ 8 by norm_num)
  rw [h8] at this
  norm_num at this

end SalmonAudit.XC
