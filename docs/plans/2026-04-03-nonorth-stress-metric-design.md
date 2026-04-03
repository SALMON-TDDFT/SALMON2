# Nonorthogonal Stress Metric Design

**Goal:** Fix SALMON's non-orthogonal primitive-cell stress so that the analytic `output_level='C'` tensor matches finite-difference derivatives and can be compared consistently against the orthogonal conventional-cell route.

**Context:** Builtin `r2scan` now runs on the primitive non-orthogonal Si cell, but the remaining mismatch is not specific to `r2scan`. Finite-difference checks on primitive-cell `PZ` show that the raw non-orthogonal stress tensor is not the Cartesian Cauchy stress. The mismatch is broad: `Kinetic`, `Hartree`, `Local`, `Nonlocal`, `XC`, and `Total` fail both shear and isotropic finite-difference checks, while `Ewald` already matches. This indicates that the current non-Ewald stress sectors are missing non-orthogonal metric contributions that are present in the Hamiltonian and energy evaluation.

**Observed evidence:**
- Primitive shear finite difference with symmetric primitive `4x4x4` k-grid gives `dE/dδ/V = 5.54 GPa`, while raw `Total xy = 1.74 GPa`.
- Primitive isotropic scaling finite difference gives `P_total = -66.66 GPa`, while analytic `P_total = -10.31 GPa`.
- In the same isotropic scan, `Ewald` matches within `0.2 GPa`, while all other sectors miss by tens to hundreds of GPa.
- The mismatch therefore cannot be explained by k-point mapping or a single output-side tensor transform.

**Approach:**
- Keep the orthogonal stress implementation unchanged.
- Treat the current primitive/non-orthogonal mismatch as a missing-terms problem inside each non-Ewald stress sector.
- Start with `Kinetic`, because its non-orthogonal Hamiltonian dependence on `rmatrix_B` and `coef_F` is explicit and finite-difference discrepancies are large and clean.
- Use `PZ` primitive finite-difference runs as the control harness before touching `r2scan`.

**Why Kinetic first:**
- `calc_stress_kin` currently accumulates `-Re[w_a^* w_b]` from `calc_gradient_psi(..., system%rmatrix_B)` but does not explicitly differentiate the non-orthogonal metric coefficients.
- The non-orthogonal Hamiltonian, by contrast, depends on both `B` and `F = B B^T`.
- If `Kinetic` is corrected first and begins to match finite differences, the missing-term pattern for the other sectors becomes much easier to localize.

**Chosen scope:**
- `PZ` primitive-cell control first
- `GS/stress only`
- `nspin=1 only`
- `nproc_rgrid=1 only`
- `4 MPI(k) x 1 OMP`
- `diamond Si` 2-atom primitive cell
- finite-difference validation for both:
  - isotropic scaling
  - Cartesian shear

**Non-goals for this milestone:**
- No change to orthogonal-cell formulas.
- No attempt to solve all non-orthogonal sectors in one patch.
- No attempt to reinterpret the current raw tensor by an output-only basis transform and declare success.

**Candidate approaches considered:**

1. **Output-side tensor transform only**
   - Reinterpret primitive raw stress in another basis and compare after transformation.
   - Pros: lowest code churn.
   - Cons: ruled out by runtime evidence, because `Ewald` already matches FD while the other sectors do not.

2. **Sector-by-sector nonorthogonal metric correction**
   - Add the missing non-orthogonal metric/strain terms to each analytic stress sector, starting from `Kinetic`.
   - Pros: consistent with the observed `Ewald` vs non-Ewald split, preserves orthogonal behavior.
   - Cons: requires careful derivation and repeated FD validation.

3. **Full stress re-derivation and refactor**
   - Rewrite the entire non-orthogonal stress path around a unified strain-derivative formalism.
   - Pros: cleanest long-term architecture.
   - Cons: too large for the current debugging loop.

**Recommendation:** Approach 2.

**Architecture details:**
- Add a non-orthogonal-only correction path to `calc_stress_kin`.
- Preserve the current orthogonal branch exactly.
- Validate `Kinetic` against primitive finite differences before propagating the same strategy to:
  - `Hartree`
  - `XC`
  - `Local`
  - `Nonlocal`
- Continue to use `Ewald` as the reference sector that already behaves correctly.

**Validation strategy:**
- Reuse the primitive finite-difference setup already established in `build-mpi-gcc15/runtime-checks`.
- For each sector under repair, compare:
  - analytic `xx/yy/zz/xy/yz/xz`
  - finite-difference derivatives for the matching strain path
- Keep orthogonal regression checks in the loop so the conventional-cell route does not move.

**Success criteria:**
- Primitive `PZ` non-orthogonal `Kinetic` matches finite differences within about `0.1 GPa`.
- The corrected `Total` stress moves materially toward FD without regressing orthogonal-cell results.
- Existing orthogonal `PZ` and `r2scan` stress tests remain green.
- The resulting correction pattern is reusable for later `r2scan` primitive stress cleanup.

**Main risks:**
- The exact non-orthogonal conjugate variable may be a strain on the normalized lattice basis rather than a naive Cartesian strain, so the derivation must stay aligned with the Hamiltonian path actually implemented.
- `Hartree`, `Local`, and `Nonlocal` may each require distinct missing terms, so a `Kinetic` fix may not immediately repair `Total`.
- Finite-difference noise can hide sign mistakes unless the same control inputs and convergence settings are reused across comparisons.
