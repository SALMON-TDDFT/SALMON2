import Mathlib

/-!
# GS: total-energy assembly, Hartree / local-pseudopotential energies (SALMON v2.3.0)

Abstraction: integrals over the grid are replaced by real numbers.
* `Ekin, Enl` : kinetic and nonlocal energies Σ f w ⟨ψ|T|ψ⟩, Σ f w ⟨ψ|V_nl|ψ⟩;
* `Uloc = ∫ V_psl ρ`, `UH = ∫ V_H ρ`, `Uxc = ∫ V_xc ρ` ;
* `EH` the Hartree energy (= UH/2 when V_H is the Hartree potential of ρ), `Exc`, `EII`.
The eigenvalue sum is `Σ f ε = Ekin + Enl + Uloc + UH + Uxc` because
`V_local = Vpsl + Vh + Vxc` (`src/common/hamiltonian.f90:1213`) and
`esp = ⟨ψ|H|ψ⟩` (`src/common/total_energy.f90:575,676`).
-/

noncomputable section
namespace SalmonAudit.GS

/-- Band-energy sum Σ f ε. -/
def bandSum (Ekin Enl Uloc UH Uxc : ℝ) : ℝ := Ekin + Enl + Uloc + UH + Uxc

/-- Kohn–Sham total energy (correct). -/
def EtotCorrect (Ekin Enl Uloc EH Exc EII : ℝ) : ℝ := Ekin + Enl + Uloc + EH + Exc + EII

/-- Code, isolated, `method_poisson='cg'/'dirichlet'` (`src/common/total_energy.f90:94-111,261-263`):
`Σfε − ½∫V_Hρ − ∫V_xcρ + E_xc + E_ion-ion`. -/
def EtotIsoCG (Ekin Enl Uloc UH Uxc Exc EII : ℝ) : ℝ :=
  bandSum Ekin Enl Uloc UH Uxc - (1/2) * UH - Uxc + Exc + EII

/-- GS-09 (OK). -/
theorem EtotIsoCG_ok (Ekin Enl Uloc EH Uxc Exc EII : ℝ) :
    EtotIsoCG Ekin Enl Uloc (2 * EH) Uxc Exc EII = EtotCorrect Ekin Enl Uloc EH Exc EII := by
  unfold EtotIsoCG bandSum EtotCorrect; ring

/-- Code, isolated, `method_poisson='ft'` (`src/common/total_energy.f90:226-239,261-267`):
`Σfε − ∫V_xcρ + E_xc + E_ion-ion + E_sum(1) + E_sum(2) + E_sum(3)` where
`E_sum(1)` = G-space Hartree energy `EHg` and `E_sum(2)+E_sum(3)` = G-space local e–ion energy `Elocg`. -/
def EtotIsoFT (Ekin Enl Uloc UH Uxc Exc EII EHg Elocg : ℝ) : ℝ :=
  bandSum Ekin Enl Uloc UH Uxc - Uxc + Exc + EII + EHg + Elocg

/-- GS-10 (BUG). Even in the most favourable case (G-space terms exactly equal to the real-space
ones: `EHg = EH`, `Elocg = Uloc`, `UH = 2EH`) the 'ft' isolated total energy double counts
(the actual size of the error additionally depends on GS-11/GS-33, which make EHg, Elocg wrong):
`E_code = E_correct + 2 E_H + E_loc`. -/
theorem EtotIsoFT_doublecount (Ekin Enl Uloc EH Uxc Exc EII : ℝ) :
    EtotIsoFT Ekin Enl Uloc (2 * EH) Uxc Exc EII EH Uloc
      = EtotCorrect Ekin Enl Uloc EH Exc EII + 2 * EH + Uloc := by
  unfold EtotIsoFT bandSum EtotCorrect; ring

/-- GS-10 (BUG) concrete witness: E_H = 1, E_loc = −3 (any system with 2E_H + E_loc ≠ 0). -/
theorem EtotIsoFT_ne :
    EtotIsoFT 0 0 (-3) 2 0 0 0 1 (-3) ≠ EtotCorrect 0 0 (-3) 1 0 0 := by
  unfold EtotIsoFT bandSum EtotCorrect; norm_num

/-- Code, valence electron–ion term in the isolated 'ft' path, per G
(`src/common/total_energy.f90:211`): `sysvol*coef*(−ρe·conj ρi)` with `sysvol = 8·det_a`
(`:178`), `ρe = FFT/(8N)` normalised to the doubled box (`src/poisson/poisson_isolated.f90:1147`)
but `ρi = S/det_A` normalised to the original box (`src/atom/pp/prep_pp.f90:769-771`).
`x = Re(ρe·conj S)`. -/
def EeiIsoFTCode (d c x : ℝ) : ℝ := (8 * d) * c * (-(x / d))

/-- Correct: both densities normalised to the doubled box `8d`. -/
def EeiIsoFTCorrect (d c x : ℝ) : ℝ := (8 * d) * c * (-(x / (8 * d)))

/-- GS-11 (BUG). Normalisation defect alone: the valence e–ion G-space term is 8× the value it would
have with consistently normalised densities (independent of GS-33, which corrupts ρe itself). -/
theorem EeiIsoFT_factor8 (d c x : ℝ) (hd : d ≠ 0) :
    EeiIsoFTCode d c x = 8 * EeiIsoFTCorrect d c x := by
  unfold EeiIsoFTCode EeiIsoFTCorrect; field_simp

/-- Code, ordinary isolated FT (`yn_ffte='n'`, no FFTW), `src/poisson/poisson_isolated.f90:1160-1166`:
`zrhoG_ele = ff4z` is read **before** `ff4z` is overwritten, so it holds the previous call's
`coef·ρ_prev(G)` (= V_H^{prev}(G)), not ρ(G). Hartree term (`src/common/total_energy.f90:207`):
`8V·coef·|zrhoG_ele|²/2`. Here `c = coef = 4π/G²`, `p = ρ(G)` (self-consistent: ρ_prev = ρ). -/
def EHisoFTCode (V c : ℝ) (p : ℂ) : ℝ := 8 * V * c * Complex.normSq (c * p) / 2

/-- Correct: `8V·coef·|ρ(G)|²/2`. -/
def EHisoFTCorrect (V c : ℝ) (p : ℂ) : ℝ := 8 * V * c * Complex.normSq p / 2

/-- GS-33 (BUG): even at self-consistency the code Hartree (and, via the same array,
the valence/core e–ion) G-term carries an extra factor `coef²` (resp. `coef`). -/
theorem EHisoFT_stale (V c : ℝ) (p : ℂ) :
    EHisoFTCode V c p = c ^ 2 * EHisoFTCorrect V c p := by
  unfold EHisoFTCode EHisoFTCorrect
  rw [Complex.normSq_mul, Complex.normSq_ofReal]; ring

/-- GS-33 witness: coef = 2, ρ(G) = 1, V = 1. -/
theorem EHisoFT_stale_ne : EHisoFTCode 1 2 1 ≠ EHisoFTCorrect 1 2 1 := by
  unfold EHisoFTCode EHisoFTCorrect; norm_num

/-- Code, periodic (`src/common/total_energy.f90:505`):
`E_kin + E_h + E_ion_loc + E_ion_nloc + E_xc + E_ion_ion`. GS-12 (OK). -/
theorem EtotPeriodic_ok (Ekin EH Eloc Enl Exc EII : ℝ) :
    Ekin + EH + Eloc + Enl + Exc + EII = EtotCorrect Ekin Enl Eloc EH Exc EII := by
  unfold EtotCorrect; ring

/-- Periodic Hartree energy per G (`src/common/total_energy.f90:445`):
`sysvol*coef*|ρe|²*0.5`, ρe(G) = (1/N)Σ_r ρ(r)e^{-iGr} (`src/poisson/poisson_periodic.f90:149,168`).
Correct `E_H = (V/2) Σ_{G≠0} 4π|ρ(G)|²/G²`. GS-13 (OK): identical expression (`p2 = |ρe|²`). -/
theorem EHperiodic_ok (V G2 p2 : ℝ) :
    V * (4 * Real.pi / G2) * (p2 * 0.5) = V / 2 * (4 * Real.pi * p2 / G2) := by
  ring

/-- Periodic local e–ion energy per G and atom (`src/common/total_energy.f90:448-455`):
`V·coef·(−ρe·conj ρi) + conj(ρe)·Ṽ_sr(G)·e^{−iG·R}` with `ρi = Z e^{−iGR}/V`;
correct `V·conj(ρe)·V_ion(G)`, `V_ion(G) = (Ṽ_sr(G) − 4πZ/G²) e^{−iGR}/V`
(`src/atom/pp/prep_pp.f90:879`). Real parts agree. GS-14 (OK). -/
theorem Eloc_periodic_ok (V c Z Vsr : ℝ) (rhoe ph : ℂ) (hV : V ≠ 0) :
    (V * c * (-(rhoe * star (Z * ph / V))) + star rhoe * Vsr * ph).re
      = (V * star rhoe * ((Vsr - c * Z) * ph / V)).re := by
  have hV' : (V : ℂ) ≠ 0 := by exact_mod_cast hV
  have key : V * star rhoe * ((Vsr - c * Z) * ph / V)
      = star rhoe * Vsr * ph - V * c * (star rhoe * (Z * ph / V)) := by
    field_simp
  rw [key]
  simp only [Complex.add_re, Complex.sub_re, Complex.mul_re, Complex.neg_re, Complex.neg_im,
    Complex.mul_im, Complex.ofReal_re, Complex.ofReal_im, Complex.star_def, Complex.conj_re,
    Complex.conj_im, Complex.div_re, Complex.div_im, Complex.ofReal_re, Complex.ofReal_im,
    Complex.conj_ofReal, map_mul, map_div₀]
  ring

/-- Nonlocal energy extraction (`src/common/total_energy.f90:612-626,744-753`):
`E_nl = Σ f w ⟨ψ|(H − (T + V_local))|ψ⟩`; with `H = T + V_local + V_nl` this is ⟨ψ|V_nl|ψ⟩.
GS-15 (OK). -/
theorem Enl_extract_ok {ι : Type*} (s : Finset ι) (f t v n : ι → ℝ) :
    ∑ i ∈ s, f i * ((t i + v i + n i) - (t i + v i)) = ∑ i ∈ s, f i * n i := by
  refine Finset.sum_congr rfl ?_; intro i _; ring

end SalmonAudit.GS
end
