# Primitive-Cell r2SCAN Nonorthogonal Design

**Goal:** Enable builtin `r2scan` ground-state and stress calculations for the 2-atom primitive Si cell on a non-orthogonal lattice, and compare the `output_level='C'` stress tensor decomposition against the existing 8-atom conventional-cell route.

**Context:** SALMON already accepts non-orthogonal lattice vectors through `al_vec1:3`, initializes `system%primitive_a`, and has a non-orthogonal stencil path in the main Hamiltonian. However, builtin `r2scan` currently enables an additional `tau`-dependent XC operator payload, and that operator is guarded to orthogonal stencils only. As a result, non-orthogonal primitive-cell `r2scan` stops before the stress path can be evaluated, even though `PZ` and `TBmBJ` can use the non-orthogonal GS path.

**Approach:**
- Keep the existing orthogonal `r2scan` implementation unchanged.
- Add a non-orthogonal branch for the XC `tau` operator in `add_xc_tau_operator`.
- Match its geometry to the existing non-orthogonal Laplacian stencil, including the `F(1:6)` metric coefficients and mixed-derivative terms.
- Reuse the existing `calc_gradient_psi(..., system%rmatrix_B)` and `stress_xc_tau` machinery so the GS and stress paths are discretized consistently.

**Why this is the first milestone:**
- The current blocker is not the stress formula itself, but the GS Hamiltonian path.
- `r2scan` stress cannot be meaningfully validated on the primitive cell until the non-orthogonal `tau` operator is available in the SCF loop.
- The stress implementation already works with full `3x3` tensors and uses non-orthogonal gradients through `rmatrix_B`, so opening the GS path is the highest-leverage step.

**Chosen scope:**
- `GS only`
- `nspin=1 only`
- `OpenACC off`
- `nproc_rgrid=1 only`
- `diamond Si` 2-atom primitive cell
- stress validation at `output_level='C'`

**Non-goals for this milestone:**
- No attempt to generalize non-orthogonal `tau` support to RT or MS modes.
- No attempt to support `nproc_rgrid > 1` for non-orthogonal lattices.
- No attempt to rework `TBmBJ` or libxc meta-GGA paths in the same change.

**Candidate approaches considered:**

1. **Targeted non-orthogonal `tau` operator branch**
   - Extend `add_xc_tau_operator` with a non-orthogonal path modeled after `zstencil_nonorthogonal`.
   - Pros: minimal blast radius, directly addresses the blocker, keeps stress and GS geometry aligned.
   - Cons: duplicates some stencil algebra now expressed separately for orthogonal and non-orthogonal cases.

2. **Refactor to a unified variable-coefficient stencil operator**
   - Introduce a common operator backend for orthogonal and non-orthogonal `vtau`.
   - Pros: cleaner long-term architecture.
   - Cons: too much surface area for the current milestone and harder to validate incrementally.

3. **Approximate primitive-cell `r2scan` without the `tau` operator**
   - Force a fallback route that uses `tau` only inside the XC kernel.
   - Pros: fastest path to a runnable case.
   - Cons: changes the functional definition and is not acceptable for comparison work.

**Recommendation:** Approach 1.

**Architecture details:**
- The orthogonal `tau` operator currently evaluates a variable-coefficient form of `-1/2 div(vtau grad psi)` along Cartesian axes.
- The non-orthogonal extension should use the same decomposition already present in `zstencil_nonorthogonal`:
  - direct-axis second-derivative contributions weighted by `F(1:3)`
  - mixed `yz/zx/xy` contributions weighted by `F(4:6)`
- `vtau` should be sampled with the same shadow-value assumptions already used in the orthogonal path.
- The mixed terms should be introduced through the same two-stage work-array pattern as `zstencil_nonorthogonal`, but with local coefficient averaging so the operator remains the variable-coefficient analogue of the constant-coefficient stencil.

**Validation strategy:**
- First prove that non-orthogonal primitive-cell `r2scan GS` converges.
- Then run primitive-cell `r2scan stress` at `output_level='C'`.
- Compare the primitive-cell stress tensor to the conventional-cell tensor after expressing them in the same Cartesian frame.
- Monitor every `C`-level sector, not only total pressure:
  - `Kinetic`
  - `Hartree`
  - `XC`
  - `Local`
  - `Nonlocal`
  - `Ewald`
  - `Total`

**Success criteria:**
- Primitive-cell `r2scan GS` runs to completion on a non-orthogonal lattice.
- Primitive-cell `r2scan stress` runs to completion at `output_level='C'`.
- For primitive vs conventional Si with the same functional, the transformed stress tensor sectors match within about `0.1 GPa` for the first milestone.
- Existing orthogonal-cell `PZ` and `r2scan` source tests remain green.

**Main risks:**
- The mixed-derivative variable-coefficient discretization may not be self-adjoint if implemented carelessly.
- A GS implementation that converges can still produce sector-level stress mismatches if the `tau` operator geometry and the `stress_xc_tau` geometry are inconsistent.
- Non-orthogonal primitive inputs may expose pre-existing assumptions in test tooling or sample preparation rather than in the physics code itself.
