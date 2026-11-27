import Mathlib

/-!
# XC audit: libxc interface (`src/xc/salmon_xc.f90:1557-1781`) and the GGA potential
assembly in `exchange_correlation` (`salmon_xc.f90:330-442`)

Conventions: libxc unpolarized takes total `ρ`, `σ = |∇ρ|²`, returns `zk = ε` (per particle),
`vrho = ∂(ρε)/∂ρ`, `vsigma = ∂(ρε)/∂σ`.  SALMON stores `rdedd = -2 vsigma ∇ρ` and adds
`∇·rdedd` (`calc_gradient_field` of each component, summed over `idir`) to `Vxc`.
-/

namespace SalmonAudit.XC

/-! ## XC-12  unpolarized GGA term: `rdedd_i = -2 vsigma ∂_iρ = -∂f/∂(∂_iρ)` (`salmon_xc.f90:1644`) -/

/-- For `f = F(σ)`, `σ = g₁²+g₂²+g₃²`: `∂f/∂g₁ = 2 F'(σ) g₁`; so `-2 vsigma g₁` is `-∂f/∂g₁`,
and `v = ∂f/∂ρ - ∇·(∂f/∂∇ρ) = vrho + ∇·rdedd` (the sign used at `salmon_xc.f90:397`). -/
theorem gga_rdedd_ok {F : ℝ → ℝ} {vsig g1 g2 g3 : ℝ}
    (hF : HasDerivAt F vsig (g1 ^ 2 + g2 ^ 2 + g3 ^ 2)) :
    HasDerivAt (fun x => F (x ^ 2 + g2 ^ 2 + g3 ^ 2)) (-(-2 * vsig * g1)) g1 := by
  have hs : HasDerivAt (fun x => x ^ 2 + g2 ^ 2 + g3 ^ 2) (2 * g1) g1 := by
    simpa using ((hasDerivAt_pow 2 g1).add_const (g2 ^ 2)).add_const (g3 ^ 2)
  exact (hF.comp g1 hs).congr_deriv (by ring)

/-! ## XC-13  spin-polarized GGA cross terms (`salmon_xc.f90:1693-1697, 1765-1773`) -/

/-- libxc ordering `σ = (∇ρ↑·∇ρ↑, ∇ρ↑·∇ρ↓, ∇ρ↓·∇ρ↓)`.  Varying the x-component `a` of ∇ρ↑
(with ∇ρ↓ x-component `b`, other components frozen in `cuu`, `cud`):
`∂f/∂a = 2 v₁ a + v₂ b`, which is minus the code's `rdedd_s(:,1,idir)`. -/
theorem gga_spin_rdedd_ok {F : ℝ × ℝ × ℝ → ℝ} {L : ℝ × ℝ × ℝ →L[ℝ] ℝ} {a b cuu cud sdd v1 v2 v3 : ℝ}
    (hF : HasFDerivAt F L (a ^ 2 + cuu, a * b + cud, sdd))
    (hL : ∀ x y z, L (x, y, z) = v1 * x + v2 * y + v3 * z) :
    HasDerivAt (fun t => F (t ^ 2 + cuu, t * b + cud, sdd)) (-(-2 * v1 * a - v2 * b)) a := by
  have h1 : HasDerivAt (fun t => t ^ 2 + cuu) (2 * a) a := by
    simpa using (hasDerivAt_pow 2 a).add_const cuu
  have h2 : HasDerivAt (fun t => t * b + cud) b a := by
    simpa using ((hasDerivAt_id a).mul_const b).add_const cud
  have h3 : HasDerivAt (fun _ : ℝ => sdd) 0 a := hasDerivAt_const a sdd
  have hp := h1.prodMk (h2.prodMk h3)
  have := hF.comp_hasDerivAt a hp
  refine this.congr_deriv ?_
  rw [hL]; ring

/-! ## XC-14  BUG: libxc unpolarized + NLCC energy density (`salmon_xc.f90:1580-1583, 1629-1631`)

`rho_1d = ρ + ρ_c` is passed to libxc, but `eexc = zk * rho` multiplies by the valence `ρ`.
Correct NLCC: `e_xc = (ρ+ρ_c) ε_xc(ρ+ρ_c)` (as the builtin PZ/PW/r2SCAN paths do). -/

def libxcEexcCode (zk ρv _ρc : ℝ) : ℝ := zk * ρv
def libxcEexcCorrect (zk ρv ρc : ℝ) : ℝ := zk * (ρv + ρc)

theorem libxc_eexc_nlcc_bug (zk ρv ρc : ℝ) :
    libxcEexcCode zk ρv ρc = libxcEexcCorrect zk ρv ρc - zk * ρc := by
  unfold libxcEexcCode libxcEexcCorrect; ring

theorem libxc_eexc_nlcc_counterexample : libxcEexcCode (-1) 1 1 ≠ libxcEexcCorrect (-1) 1 1 := by
  unfold libxcEexcCode libxcEexcCorrect; norm_num

/-! ## XC-15  BUG: NLCC omitted from the density gradient for GGA / meta-GGA
(`salmon_xc.f90:234-242` gradient of valence `rho_s(1)` only; `1586-1588` σ; `1484` r2SCAN)

The functional is evaluated at `n = ρ_v + ρ_c` but with `σ = |∇ρ_v|²` instead of `|∇ρ_v + ∇ρ_c|²`. -/

def sigmaCode (gv _gc : ℝ) : ℝ := gv ^ 2
def sigmaCorrect (gv gc : ℝ) : ℝ := (gv + gc) ^ 2

theorem nlcc_sigma_counterexample : sigmaCode 1 1 ≠ sigmaCorrect 1 1 := by
  unfold sigmaCode sigmaCorrect; norm_num

/-! ## XC-16  BUG: libxc meta-GGA potential drops `vtau` and `vlapl` (`salmon_xc.f90:1618-1650`)

Only `vrho` (and `vsigma` through `rdedd`) reach `Vxc`; `lvxc_tmp_1d` (`∂f/∂∇²ρ`) and
`tvxc_tmp_1d` (`∂f/∂τ`) are discarded, and no τ operator is installed (payload is set only by
builtin r2SCAN).  Model: one-parameter variation `ρ+λδρ`, `τ+λδτ` of `f(ρ,τ) = ρ·τ`;
the KS potential must reproduce `dE/dλ = vrho δρ + vtau δτ`, the code gives `vrho δρ`. -/

theorem mgga_first_variation (ρ τ dρ dτ : ℝ) :
    HasDerivAt (fun l => (ρ + l * dρ) * (τ + l * dτ)) (τ * dρ + ρ * dτ) 0 := by
  have h1 : HasDerivAt (fun l => ρ + l * dρ) dρ 0 := by
    simpa using ((hasDerivAt_id (0:ℝ)).mul_const dρ).const_add ρ
  have h2 : HasDerivAt (fun l => τ + l * dτ) dτ 0 := by
    simpa using ((hasDerivAt_id (0:ℝ)).mul_const dτ).const_add τ
  exact (h1.mul h2).congr_deriv (by simp; ring)

/-- code first-order energy change from its potential: only `vrho = τ`. -/
def mggaVarCode (_ρ τ dρ _dτ : ℝ) : ℝ := τ * dρ
def mggaVarCorrect (ρ τ dρ dτ : ℝ) : ℝ := τ * dρ + ρ * dτ

theorem mgga_vtau_dropped_counterexample : mggaVarCode 1 1 0 1 ≠ mggaVarCorrect 1 1 0 1 := by
  unfold mggaVarCode mggaVarCorrect; norm_num

/-! ## XC-17  BUG: τ fed to libxc meta-GGA in RT periodic omits the vector potential
(`salmon_xc.f90:626-633` uses `system%vec_k`; the Hamiltonian uses `vec_k + vec_Ac`,
`hamiltonian.f90:157`).  Model: one orbital `u = c` (constant), so `∇u = 0`:
`τ = ½|(∇ + i q) u|² = ½ q² c²`; code `q = k`, correct (velocity gauge) `q = k + A`. -/

noncomputable def tauCode (k _A c : ℝ) : ℝ := k ^ 2 * c ^ 2 / 2
noncomputable def tauCorrect (k A c : ℝ) : ℝ := (k + A) ^ 2 * c ^ 2 / 2

theorem rt_tau_missing_A_counterexample : tauCode 0 1 1 ≠ tauCorrect 0 1 1 := by
  unfold tauCode tauCorrect; norm_num

end SalmonAudit.XC
