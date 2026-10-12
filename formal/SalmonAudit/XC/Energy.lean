import Mathlib

/-!
# XC audit: how E_xc and V_xc enter the total energy; τ / j construction

Sources: `src/xc/salmon_xc.f90:444-464` (E_xc), `576-686` (calc_tau),
`src/common/total_energy.f90:80-268` (isolated), `:503`, `:712-795` (periodic).
Grids abstracted to a `Finset`; `Hvol` an arbitrary real.
-/

namespace SalmonAudit.XC

open Finset

/-! ## XC-18  `E_xc = Hvol Σ_r e_xc(r)` with `e_xc = n ε_xc` (per volume) -/

theorem exc_integral_ok {ι : Type*} (s : Finset ι) (Hvol : ℝ) (n ε : ι → ℝ) :
    (∑ i ∈ s, ε i * n i) * Hvol = ∑ i ∈ s, Hvol * (n i * ε i) := by
  rw [Finset.sum_mul]; refine Finset.sum_congr rfl ?_; intro i _; ring

/-! ## XC-19  isolated total energy double counting (`total_energy.f90:88,105-106,262-263`)

Band energy `Σ f ε = T + E_nl + ∫(V_ion + V_H + V_xc) ρ + e_tau` (KS equation; `e_tau = ∫ vtau τ_GI`
is the meta-GGA τ-operator expectation).  Code: `Etot = Σfε - ½∫V_Hρ - ∫V_xcρ + E_xc + E_ii - e_tau`. -/

theorem isolated_total_energy_ok (T Enl Eloc VHrho Vxcrho etau Exc Eii : ℝ) :
    let band := T + Enl + Eloc + VHrho + Vxcrho + etau
    band - etau + (-(1/2) * VHrho - Vxcrho) + Exc + Eii
      = T + Enl + Eloc + (1/2) * VHrho + Exc + Eii := by
  simp only; ring

/-! ## XC-20  periodic: `E_ion_nloc = Σ f⟨ψ|H - T - V_local|ψ⟩ - e_tau` (`total_energy.f90:766-795`)
where `V_local = V_ion,loc + V_H + V_xc`; then `E_tot` is the direct sum (`:503`). -/

theorem periodic_nloc_ok (T Enl Vlocrho etau : ℝ) :
    ((T + Enl + Vlocrho + etau) - (T + Vlocrho)) - etau = Enl := by ring

/-! ## XC-21  `calc_tau` (`salmon_xc.f90:631-634`): with Bloch `ψ = e^{ik·r} u`, the code's
`zs = ∇u + i k u` equals `e^{-ik·r}∇ψ`; hence `½|zs|²` and `Im(u* zs)` are the true
`½|∇ψ|²` and `Im(ψ*∇ψ)` (one Cartesian component; `e` is the phase, `|e| = 1`). -/

theorem calc_tau_phase_ok (e u du : ℂ) (k : ℝ) (he : Complex.normSq e = 1) :
    Complex.normSq (e * (du + Complex.I * k * u)) = Complex.normSq (du + Complex.I * k * u) := by
  rw [Complex.normSq_mul, he, one_mul]

theorem calc_j_phase_ok (e u du : ℂ) (k : ℝ) (he : Complex.normSq e = 1) :
    (starRingEnd ℂ (e * u) * (e * (du + Complex.I * k * u))).im
      = (starRingEnd ℂ u * (du + Complex.I * k * u)).im := by
  have : starRingEnd ℂ (e * u) * (e * (du + Complex.I * k * u))
      = (starRingEnd ℂ e * e) * (starRingEnd ℂ u * (du + Complex.I * k * u)) := by
    rw [map_mul]; ring
  rw [this, mul_comm (starRingEnd ℂ e) e, Complex.mul_conj, he]
  simp

end SalmonAudit.XC
