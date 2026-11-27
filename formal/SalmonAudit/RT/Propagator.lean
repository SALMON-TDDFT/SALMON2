import Mathlib

/-!
# RT propagators (SALMON v2.3.0, `src/rt/taylor.f90`, `src/rt/initialization_rt.f90`,
`src/rt/time_evolution_step.f90`)

Abstractions: wavefunctions are vectors of an arbitrary complex module `M`; the Hamiltonian
is an arbitrary linear map `H : M →ₗ[ℂ] M`; local potentials are real numbers (one grid point,
the formulas are pointwise).
-/

noncomputable section
namespace SalmonAudit.RT

open Complex

/-! ## RT-01 Taylor coefficients -/

/-- The inner loop `do m=1,n; zc(n)=zc(n)/m; end do` starting from `z`. -/
def divLoop (z : ℂ) : ℕ → ℂ
  | 0 => z
  | m + 1 => divLoop z m / ((m + 1 : ℕ) : ℂ)

/-- Code: `src/rt/initialization_rt.f90:523-527` (middlepoint) and `:516-520` (aetrs, with `dt/2`):
`zc(n) = (-zi*dt)**n` then divided by `m` for `m = 1..n`. -/
def zcCode (dt : ℝ) (n : ℕ) : ℂ := divLoop ((-I * (dt : ℂ)) ^ n) n

/-- Correct: the n-th Taylor coefficient of `exp(-i dt H)` is `(-i dt)^n / n!`. -/
def zcCorrect (dt : ℝ) (n : ℕ) : ℂ := (-I * (dt : ℂ)) ^ n / (n.factorial : ℂ)

theorem divLoop_eq (z : ℂ) (n : ℕ) : divLoop z n = z / (n.factorial : ℂ) := by
  induction n with
  | zero => simp [divLoop]
  | succ m ih =>
    simp only [divLoop, ih, Nat.factorial_succ, Nat.cast_mul]
    have h1 : ((m.factorial : ℕ) : ℂ) ≠ 0 := by exact_mod_cast Nat.factorial_ne_zero m
    have h2 : ((m + 1 : ℕ) : ℂ) ≠ 0 := by exact_mod_cast Nat.succ_ne_zero m
    field_simp

theorem zc_ok (dt : ℝ) (n : ℕ) : zcCode dt n = zcCorrect dt n := by
  simp [zcCode, zcCorrect, divLoop_eq]

/-- ETRS half step uses `dt/2` (`initialization_rt.f90:517`): the same theorem with `dt/2`. -/
theorem zc_etrs_ok (dt : ℝ) (n : ℕ) : zcCode (dt / 2) n = zcCorrect (dt / 2) n := zc_ok _ _

/-! ## RT-02 Taylor accumulation with in-place buffer swapping (`taylor.f90:45-120`) -/

/-- The three buffers `tspsi_in`, `tspsi_out`, `sshtpsi`. -/
structure St (M : Type*) where
  tin : M
  tout : M
  ssh : M

variable {M : Type*} [AddCommGroup M] [Module ℂ M]

/-- One pass `nn` of the loop in `src/rt/taylor.f90:45-120`.
odd `nn`: `sshtpsi = H tspsi_in`; `nn==1`: `out = in + zc*ssh`, else `out += zc*ssh`.
even `nn`: `tspsi_in = H sshtpsi`; `out += zc*tspsi_in`. -/
def step (H : M →ₗ[ℂ] M) (c : ℕ → ℂ) (nn : ℕ) (s : St M) : St M :=
  if nn % 2 = 1 then
    if nn = 1 then ⟨s.tin, s.tin + c nn • H s.tin, H s.tin⟩
    else ⟨s.tin, s.tout + c nn • H s.tin, H s.tin⟩
  else ⟨H s.ssh, s.tout + c nn • H s.ssh, s.ssh⟩

/-- State after `n` passes; `junk` = arbitrary initial content of `tspsi_out`, `sshtpsi`. -/
def run (H : M →ₗ[ℂ] M) (c : ℕ → ℂ) (ψ junk : M) : ℕ → St M
  | 0 => ⟨ψ, junk, junk⟩
  | n + 1 => step H c (n + 1) (run H c ψ junk n)

/-- Correct: truncated Taylor series `ψ + Σ_{k=1}^{n} c_k H^k ψ`. -/
def taylorSum (H : M →ₗ[ℂ] M) (c : ℕ → ℂ) (ψ : M) (n : ℕ) : M :=
  ψ + ∑ k ∈ Finset.range n, c (k + 1) • (H ^ (k + 1)) ψ

lemma pow_succ_apply' (H : M →ₗ[ℂ] M) (n : ℕ) (x : M) : (H ^ (n + 1)) x = H ((H ^ n) x) := by
  rw [pow_succ']; rfl

theorem taylor_invariant (H : M →ₗ[ℂ] M) (c : ℕ → ℂ) (ψ junk : M) (n : ℕ) (hn : 1 ≤ n) :
    (run H c ψ junk n).tout = taylorSum H c ψ n ∧
    (n % 2 = 1 → (run H c ψ junk n).ssh = (H ^ n) ψ ∧ (run H c ψ junk n).tin = (H ^ (n - 1)) ψ) ∧
    (n % 2 = 0 → (run H c ψ junk n).tin = (H ^ n) ψ) := by
  induction n, hn using Nat.le_induction with
  | base =>
    simp [run, step, taylorSum]
  | succ n hn ih =>
    obtain ⟨h1, h2, h3⟩ := ih
    rcases Nat.mod_two_eq_zero_or_one n with hpar | hpar
    · -- n even, n+1 odd
      have hodd : (n + 1) % 2 = 1 := by omega
      have hne : n + 1 ≠ 1 := by omega
      have htin := h3 hpar
      refine ⟨?_, ?_, ?_⟩
      · simp only [run, step, hodd, hne, ite_true, ite_false]
        rw [h1, htin, ← pow_succ_apply', taylorSum, taylorSum, Finset.sum_range_succ, add_assoc]
      · intro _
        simp only [run, step, hodd, hne, ite_true, ite_false]
        rw [htin, ← pow_succ_apply']
        exact ⟨rfl, by simp⟩
      · intro h; omega
    · -- n odd, n+1 even
      have heven : (n + 1) % 2 = 0 := by omega
      have hne : ¬ ((n + 1) % 2 = 1) := by omega
      obtain ⟨hssh, _⟩ := h2 hpar
      refine ⟨?_, ?_, ?_⟩
      · simp only [run, step, hne, ite_false]
        rw [h1, hssh, ← pow_succ_apply', taylorSum, taylorSum, Finset.sum_range_succ, add_assoc]
      · intro h; omega
      · intro _
        simp only [run, step, hne, ite_false]
        rw [hssh, ← pow_succ_apply']

/-- RT-02 OK: the code's buffer juggling yields `Σ_{k=0}^{N} (-i dt)^k/k! H^k ψ` for `N = n_hamil ≥ 1`. -/
theorem taylor_ok (H : M →ₗ[ℂ] M) (dt : ℝ) (ψ junk : M) (N : ℕ) (hN : 1 ≤ N) :
    (run H (zcCode dt) ψ junk N).tout = taylorSum H (zcCorrect dt) ψ N := by
  rw [(taylor_invariant H (zcCode dt) ψ junk N hN).1]
  simp only [taylorSum, zc_ok]

/-! ## RT-03 / RT-04 extrapolations used by the ETRS propagator -/

/-- Code `time_evolution_step.f90:481`: `vloc_new = 3V(t) - 3V(t-dt) + V(t-2dt)`;
exact for any quadratic in time, i.e. it is the quadratic extrapolation to `t+dt`. -/
theorem quad_extrap_ok (a b c t h : ℝ) :
    let f := fun s : ℝ => a + b * s + c * s ^ 2
    3 * f t - 3 * f (t - h) + f (t - 2 * h) = f (t + h) := by
  simp only; ring

/-- Code `time_evolution_step.f90:508`: `Ac_ind(itt+1) ≈ 2 Ac_ind(itt) - Ac_ind(itt-1)`:
exact linear extrapolation. -/
theorem lin_extrap_ok (a b t h : ℝ) :
    let f := fun s : ℝ => a + b * s
    2 * f t - f (t - h) = f (t + h) := by
  simp only; ring

/-! ## RT-05 BUG: predictor-corrector halves the isolated-system laser potential -/

/-- Code (`time_evolution_step.f90:105-118` adds `Vbox` to `V_local`; then
`predictor_corrector`, `:394-465`): `V_tmp = V_KS(t) + Vbox`; after the predictor
`update_vlocal` rebuilds `V_local = V_KS(t+dt)` WITHOUT `Vbox` (`hamiltonian.f90:1213`);
corrector potential `0.5*(V_tmp + V_local)`. -/
def pcCode (V0 V1 Vb : ℝ) : ℝ := (1 / 2) * ((V0 + Vb) + V1)

/-- Correct: time-averaged KS potential plus the full external potential `E·r`. -/
def pcCorrect (V0 V1 Vb : ℝ) : ℝ := (1 / 2) * (V0 + V1) + Vb

theorem pc_bug (V0 V1 Vb : ℝ) : pcCode V0 V1 Vb = pcCorrect V0 V1 Vb - Vb / 2 := by
  unfold pcCode pcCorrect; ring

theorem pc_bug_ex : pcCode 0 0 1 ≠ pcCorrect 0 0 1 := by
  unfold pcCode pcCorrect; norm_num

/-! ## RT-06 BUG: aetrs doubles the isolated-system laser potential in the 2nd half step -/

/-- Code (`time_evolution_step.f90:467-519`): `V_local` already contains `Vbox(itt)` when it is
stored in `vloc_old`/extrapolated, then `Vbox(itt+1)` is added again.
`V0,V1,V2` = KS potential at steps 0,-1,-2; `B0,B1,B2` = `Vbox` added at those steps;
`Bn` = `Vbox(itt+1)`. -/
def etrsCode (V0 V1 V2 B0 B1 B2 Bn : ℝ) : ℝ :=
  3 * (V0 + B0) - 3 * (V1 + B1) + (V2 + B2) + Bn

/-- Correct: extrapolated KS potential plus the external potential once. -/
def etrsCorrect (V0 V1 V2 Bn : ℝ) : ℝ := 3 * V0 - 3 * V1 + V2 + Bn

/-- Exact excess of the code over the correct potential: `3B⁰ − 3B⁻¹ + B⁻²` (no smoothness assumed).
At startup `vloc_old` holds no `Vbox` (`initialization_rt.f90:305-308`), i.e. `B⁻¹ = B⁻² = 0` there. -/
theorem etrs_excess_exact (V0 V1 V2 B0 B1 B2 Bn : ℝ) :
    etrsCode V0 V1 V2 B0 B1 B2 Bn = etrsCorrect V0 V1 V2 Bn + (3 * B0 - 3 * B1 + B2) := by
  unfold etrsCode etrsCorrect; ring

/-- For an external potential linear in time (`B_k = β + γ k`, slowly varying field),
the extrapolated part reproduces `Bn` once more: code = correct + Bn. -/
theorem etrs_bug (V0 V1 V2 β γ : ℝ) :
    etrsCode V0 V1 V2 β (β - γ) (β - 2 * γ) (β + γ) =
      etrsCorrect V0 V1 V2 (β + γ) + (β + γ) := by
  unfold etrsCode etrsCorrect; ring

theorem etrs_bug_ex : etrsCode 0 0 0 1 1 1 1 ≠ etrsCorrect 0 0 0 1 := by
  unfold etrsCode etrsCorrect; norm_num

/-! ## RT-07 NOTE: field time index (isolated: forward difference; aetrs: shifted by +dt) -/

/-- Isolated: `rt%E_ext(itt) = -(Ac(itt+1)-Ac(itt))/dt` (`initialization_rt.f90:184-188`) is the field
at `(itt+1/2)dt`; it is used (`calcVbox(itt)`) for the step `(itt-1)dt → itt dt` whose midpoint is
`(itt-1/2)dt`. The offset is exactly one time step. -/
theorem field_time_offset (itt dt : ℝ) : (itt + 1 / 2) * dt - (itt - 1 / 2) * dt = dt := by ring

/-- RT-31 NOTE: `calc_eigen_energy`/`projection` at step `itt` (`time_evolution_step.f90:247-261`)
run while `system%vec_Ac = (Ac(itt)+Ac(itt-1))/2` (middlepoint, `:124`); `Ac_tot(itt)` is set only
at `:285` before `calc_current`. With `propagator='middlepoint'` the energy/projection labelled `t = itt·dt` uses `A` at `(itt-1/2)dt`
(aetrs differs: `vec_Ac` is then the extrapolated `A(itt+1)`, `:508`). -/
theorem energy_A_offset (itt dt : ℝ) : itt * dt - (itt - 1 / 2) * dt = dt / 2 := by ring

end SalmonAudit.RT
end
