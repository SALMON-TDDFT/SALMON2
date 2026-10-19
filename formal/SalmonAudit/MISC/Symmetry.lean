import Mathlib

/-!
# k-point / density / vector / force symmetrization (src/symmetry)

Conventions in SALMON (`src/symmetry/symmetry.f90:107-118`): `SymMatA` = point operation
in lattice (direct, reduced) coordinates read from `sym.dat`; `R_cart = A S_A A⁻¹`;
`SymMatB = B⁻¹ R_cart B` = the same operation in reciprocal-lattice coordinates,
which equals `(S_A⁻¹)ᵀ` for an orthogonal `R_cart`.

Concrete test system: fcc primitive cell `a1=(0,1,1)/2, a2=(1,0,1)/2, a3=(1,1,0)/2`
(Si-type cell) and the 4 operations of Oh that fix `a3` (what `yn_symmetry='yyn'`
keeps, `symmetry.f90:61-79`): E, mirror x↔y, mirror z→-z, C2 about [110].
Matrices below were obtained from `A⁻¹ R A` and `B⁻¹ R B` and are checked here
to satisfy `S_Bᵀ S_A = 1` and to be the lattice images of the Cartesian ops.
-/

namespace SalmonAudit.MISC.Symmetry

open Matrix

/-- fcc primitive cell, columns = lattice vectors (times 2, to stay integral). -/
def A2 : Matrix (Fin 3) (Fin 3) ℚ := !![0, 1, 1; 1, 0, 1; 1, 1, 0]

/-- Cartesian operations. -/
def R1 : Matrix (Fin 3) (Fin 3) ℚ := 1
def R2 : Matrix (Fin 3) (Fin 3) ℚ := !![0, 1, 0; 1, 0, 0; 0, 0, 1]
def R3 : Matrix (Fin 3) (Fin 3) ℚ := !![1, 0, 0; 0, 1, 0; 0, 0, -1]
def R4 : Matrix (Fin 3) (Fin 3) ℚ := !![0, 1, 0; 1, 0, 0; 0, 0, -1]

/-- `SymMatA` (lattice coordinates). -/
def SA1 : Matrix (Fin 3) (Fin 3) ℚ := 1
def SA2 : Matrix (Fin 3) (Fin 3) ℚ := !![0, 1, 0; 1, 0, 0; 0, 0, 1]
def SA3 : Matrix (Fin 3) (Fin 3) ℚ := !![0, -1, 0; -1, 0, 0; 1, 1, 1]
def SA4 : Matrix (Fin 3) (Fin 3) ℚ := !![-1, 0, 0; 0, -1, 0; 1, 1, 1]

/-- `SymMatB` (reciprocal-lattice coordinates). -/
def SB1 : Matrix (Fin 3) (Fin 3) ℚ := 1
def SB2 : Matrix (Fin 3) (Fin 3) ℚ := !![0, 1, 0; 1, 0, 0; 0, 0, 1]
def SB3 : Matrix (Fin 3) (Fin 3) ℚ := !![0, -1, 1; -1, 0, 1; 0, 0, 1]
def SB4 : Matrix (Fin 3) (Fin 3) ℚ := !![-1, 0, 1; 0, -1, 1; 0, 0, 1]

/-- The lattice matrices are the images of the Cartesian ops: `R A = A S_A`. -/
theorem SA_are_lattice_images :
    R1 * A2 = A2 * SA1 ∧ R2 * A2 = A2 * SA2 ∧ R3 * A2 = A2 * SA3 ∧ R4 * A2 = A2 * SA4 := by
  refine ⟨?_, ?_, ?_, ?_⟩ <;>
  · ext i j; fin_cases i <;> fin_cases j <;>
      simp [A2, R1, R2, R3, R4, SA1, SA2, SA3, SA4, Matrix.mul_apply, Fin.sum_univ_three]

/-- `S_B = (S_A⁻¹)ᵀ`, i.e. `S_Bᵀ S_A = 1`. -/
theorem SB_is_inverse_transpose :
    SB1.transpose * SA1 = 1 ∧ SB2.transpose * SA2 = 1 ∧ SB3.transpose * SA3 = 1 ∧
    SB4.transpose * SA4 = 1 := by
  refine ⟨?_, ?_, ?_, ?_⟩ <;>
  · ext i j; fin_cases i <;> fin_cases j <;>
      simp [SA1, SA2, SA3, SA4, SB1, SB2, SB3, SB4, Matrix.mul_apply, Fin.sum_univ_three,
        Matrix.one_apply]

/-- Current along `a3` direction, Cartesian `(1,1,0)`, in reduced coordinates
`A⁻¹ v = (0,0,2)`. -/
def vred : Fin 3 → ℚ := ![0, 0, 2]

/-- Code `sym_vector_xyz` (`src/symmetry/symmetry_vector.f90:24-32`): reduced vector
`(1/N) Σ S_B · (A⁻¹ v)`, then mapped back with `Amat`. -/
def codeSymRed : Fin 3 → ℚ :=
  (1 / 4 : ℚ) • (SB1 *ᵥ vred + SB2 *ᵥ vred + SB3 *ᵥ vred + SB4 *ᵥ vred)

/-- Correct: a real-space (polar) vector in lattice coordinates transforms with `S_A`. -/
def correctSymRed : Fin 3 → ℚ :=
  (1 / 4 : ℚ) • (SA1 *ᵥ vred + SA2 *ᵥ vred + SA3 *ᵥ vred + SA4 *ᵥ vred)

/-- MISC-14 (BUG): the invariant current is left unchanged by the correct average, but the
code returns reduced `(1,1,2)`, i.e. Cartesian `A·(1,1,2)/2 = (3/2,3/2,1)` instead of `(1,1,0)`. -/
theorem symvec_bug : correctSymRed = vred ∧ codeSymRed = ![1, 1, 2] := by
  constructor <;>
  · ext i; fin_cases i <;>
      simp [codeSymRed, correctSymRed, vred, SA1, SA2, SA3, SA4, SB1, SB2, SB3, SB4,
        Matrix.vecHead, Matrix.vecTail, Pi.smul_apply, Function.comp] <;> norm_num

theorem symvec_cartesian_bug :
    (1 / 2 : ℚ) • (A2 *ᵥ codeSymRed) = ![3 / 2, 3 / 2, 1] ∧
    (1 / 2 : ℚ) • (A2 *ᵥ correctSymRed) = ![1, 1, 0] := by
  have h := symvec_bug
  rw [h.1, h.2]
  constructor <;>
  · ext i; fin_cases i <;> simp [A2, vred, Matrix.mulVec, dotProduct, Fin.sum_univ_three] <;>
      norm_num

/-- Code `sym_vector_force_xyz` (`symmetry_vector.f90:51,63,67,76`):
`Faaa = Ainv·F`, `RotF = S_A·Faaa`, `F ← (1/N)Σ RotF` — never multiplied back by `Amat`.
Model: identity-only `sym.dat` (N = 1), simple cubic cell `A = a·1`. -/
def codeForce (a : ℚ) (F : Fin 3 → ℚ) : Fin 3 → ℚ := (1 : Matrix (Fin 3) (Fin 3) ℚ) *ᵥ (a⁻¹ • F)
/-- Correct: `A (1/N Σ S_A) A⁻¹ F` (Cartesian force). -/
def correctForce (a : ℚ) (F : Fin 3 → ℚ) : Fin 3 → ℚ :=
  a • ((1 : Matrix (Fin 3) (Fin 3) ℚ) *ᵥ (a⁻¹ • F))

/-- MISC-15 (BUG): forces come out in reduced (lattice) coordinates; with a = 10 bohr the
force is divided by 10. -/
theorem force_bug (a : ℚ) (ha : a ≠ 0) (F : Fin 3 → ℚ) : correctForce a F = F ∧
    codeForce a F = a⁻¹ • F := by
  constructor
  · simp [correctForce, smul_smul, ha]
  · simp [codeForce]

theorem force_counterexample : codeForce 10 ![1, 0, 0] ≠ correctForce 10 ![1, 0, 0] := by
  rw [(force_bug 10 (by norm_num) _).1, (force_bug 10 (by norm_num) _).2]
  intro h
  have := congrFun h 0
  norm_num at this

/-- MISC-16 (OK): `sym_rho` maps grid point `t` to `r_k = Σ_m t_m S_B(m,k) = (S_Bᵀ t)_k`
(`symmetry_rho.f90:201-203`); since `S_Bᵀ = S_A⁻¹`, this is `ρ({S|τ}⁻¹ x)`, the correct
pull-back for the group average. -/
theorem symrho_ok (SA SB : Matrix (Fin 3) (Fin 3) ℚ) (h : SB.transpose * SA = 1) :
    SB.transpose = SA⁻¹ := (Matrix.inv_eq_left_inv h).symm

/-- MISC-17 (OK): k-points in reduced reciprocal coordinates rotate with `S_B`
(`symmetry_kvector.f90:57`); check: `S_B k` pairs with `S_A r` invariantly,
`(S_B k)·(S_A r) = k·r`. -/
theorem kvec_ok (SA SB : Matrix (Fin 3) (Fin 3) ℚ) (h : SB.transpose * SA = 1)
    (k r : Fin 3 → ℚ) : (SB *ᵥ k) ⬝ᵥ (SA *ᵥ r) = k ⬝ᵥ r := by
  rw [dotProduct_comm, Matrix.dotProduct_mulVec, ← Matrix.mulVec_transpose,
    Matrix.mulVec_mulVec, h, Matrix.one_mulVec, dotProduct_comm]

end SalmonAudit.MISC.Symmetry
