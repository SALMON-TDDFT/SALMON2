import Mathlib

/-!
# SBE observables / output (SALMON v2.3.0, `src/ssbe/realtime_ssbe.f90`,
`src/ssbe/datafile_ssbe.f90`, `src/ssbe/gs_info_ssbe.f90:116-118`)
-/

namespace SalmonAudit.SBE.Out

open BigOperators

noncomputable section

/-- SBE-02 initial condition (`gs_info_ssbe.f90:117-118`, `bloch_solver_ssbe.f90:56-60`):
`ρ_ii = 2` for `i < ne/2`, else 0 (spin-degenerate insulator). Trace = `2·(ne/2)` electrons
per k-point. -/
def occ (nv : ℕ) (i : ℕ) : ℝ := if i < nv then 2 else 0

theorem sbe02_trace_ok (nv nb : ℕ) (h : nv ≤ nb) :
    (∑ i ∈ Finset.range nb, occ nv i) = 2 * nv := by
  have : (∑ i ∈ Finset.range nb, occ nv i) = ∑ i ∈ Finset.range nv, (2 : ℝ) := by
    rw [← Finset.sum_range_add_sum_Ico _ h]
    have h2 : ∑ i ∈ Finset.Ico nv nb, occ nv i = 0 :=
      Finset.sum_eq_zero (fun i hi => by
        simp [occ, not_lt.mpr (Finset.mem_Ico.mp hi).1])
    rw [h2, add_zero]
    exact Finset.sum_congr rfl (fun i hi => by simp [occ, Finset.mem_range.mp hi])
  rw [this]; simp; ring

/-- SBE-07 excited electrons (`realtime_ssbe.f90:106-109`): `nelec = tr(nb) - tr(nv)` where
`tr(m) = Σ_k w_k Σ_{i<m} Re ρ_ii / W`; correct: population of bands `nv ≤ i < nb`. -/
theorem sbe07_nex_ok (pop : ℕ → ℝ) (nv nb : ℕ) (h : nv ≤ nb) :
    (∑ i ∈ Finset.range nb, pop i) - (∑ i ∈ Finset.range nv, pop i)
      = ∑ i ∈ Finset.Ico nv nb, pop i := by
  rw [← Finset.sum_range_add_sum_Ico _ h]; ring

/-- SBE-08 code output of nex (`datafile_ssbe.f90:92,106`): value `N · L⁻³` labelled
`[length^-3]`, where `N` = excited electrons per unit cell and `L = t_unit_length%conv`
(Å per bohr = 0.529… for `unit_system='A_eV_fs'`). -/
def codeNexOut (N L : ℝ) : ℝ := N * L ^ (-3 : ℤ)

/-- Correct options: TDDFT convention (`src/io/write.f90:2051-2053`, unit "none") = `N`; or a
true density `N/(V L³)` (V = cell volume in bohr³). -/
def correctCount (N : ℝ) : ℝ := N
def correctDensity (N V L : ℝ) : ℝ := N / (V * L ^ 3)

theorem sbe08_nex_units_bug :
    codeNexOut 1 (1/2) ≠ correctCount 1 ∧ codeNexOut 1 (1/2) ≠ correctDensity 1 2 (1/2) := by
  unfold codeNexOut correctCount correctDensity; norm_num

/-- SBE-09 absorbed energy (`realtime_ssbe.f90:90`): `E_abs += E·(-J)·V·dt`; with `J` the
electron number current (TDDFT convention) the charge current is `-J`, so power per cell is
`(-J)·E·V`. -/
def codeDE (E J V dt : ℝ) : ℝ := E * (-J) * V * dt
def correctDE (E J V dt : ℝ) : ℝ := ((-1) * J) * E * V * dt

theorem sbe09_energy_ok (E J V dt : ℝ) : codeDE E J V dt = correctDE E J V dt := by
  unfold codeDE correctDE; ring

end

end SalmonAudit.SBE.Out
