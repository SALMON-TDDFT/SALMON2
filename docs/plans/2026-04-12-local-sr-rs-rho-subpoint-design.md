# Local-SR RS Rho-Subpoint Probe Design

**Goal:** Extend the diagnostic-only local-SR real-space subdivision probe so it can compare kernel-only subdivision against kernel-plus-density subdivision, then use the same probe family to study convergence for FHI, PSP, and VPS representative cases.

**Scope:** This design only changes diagnostic paths for local-SR RS stress analysis. It must not affect assembled stress, forces, energies, or non-stress code paths.

## Background

The current diagnostic probe supports `n x n x n` subcell sampling for the local-SR kernel while keeping `rho` fixed at the parent cell center. That probe already showed larger sensitivity than the earlier `C1` switch issue, especially for `si_pz_fhi_lda`, and `2x2x2 -> 4x4x4` is not yet converged.

What is still missing is separation of:

- kernel sampling sensitivity
- density sampling sensitivity

without perturbing the current production stress path.

## Recommended Approach

Add one new diagnostic-only input that controls how `rho` is sampled inside the existing RS subdivision probe:

- `center`
  - Keep the current behavior
  - Use the parent grid-point density for all subpoints
- `trilinear`
  - Build a temporary scalar density field with one-cell overlap
  - Evaluate density at each subpoint by trilinear interpolation in Cartesian grid coordinates

The probe stress tensor remains a diagnostic-only quantity stored in `stress_loc_sr_rs_subdiv_probe`. Each run chooses one rho mode, so comparisons are made across runs rather than by adding multiple probe tensors to one run.

This keeps the write scope small and avoids touching the production current RS stress assembly.

## Why This Approach

### Option A: Recommended

Keep a single probe tensor and switch rho sampling mode by input.

Pros:

- smallest change to existing diagnostic flow
- no impact on current stress path
- easy to verify with `n=1`, where `center` and `trilinear` should match
- campaign logic stays simple: compare runs by `(case, n, rho_mode)`

Cons:

- comparing two rho modes requires separate runs

### Option B: Dual probe outputs in one run

Accumulate both center-rho and trilinear-rho probes simultaneously.

Pros:

- direct in-run comparison

Cons:

- adds more diagnostic state, more output plumbing, and more code churn
- larger risk of confusing the existing output interpretation

### Option C: Modify the production RS integrator now

Pros:

- directly tests a candidate fix

Cons:

- too early for the current stage
- blurs diagnosis and physics change

## Data Flow

1. Build a temporary total-density scalar box on `mg%is_array:mg%ie_array`.
2. If `info%if_divide_rspace`, update its overlap with `srg_scalar`.
3. During RS probe evaluation:
   - for `center`, reuse the parent `rho_center`
   - for `trilinear`, interpolate `rho` at each subpoint from the overlap-aware scalar box
4. Continue to evaluate the local-SR kernel from the stress spline backend at each subpoint.
5. Write the probe mode and subdivision count into the sampled dump header and `variables.log`.

## Boundary Handling

- Interpolation is only diagnostic-only.
- The density box uses one-cell overlap so `ix+1`, `iy+1`, `iz+1` accesses are safe when the current grid point lies at the local edge.
- `n=1` must collapse to the current cell center exactly, so `center` and `trilinear` should agree to roundoff.

## Verification Strategy

### Source-level

- Add grep checks for the new rho probe input and validation.
- Add checks that `stress.f90` contains the new mode handling and trilinear helper.

### Numerical

1. FHI `n=1`, `rho=center` vs `rho=trilinear`
   - should match within roundoff
2. FHI `n=2,4`
   - compare kernel-only vs kernel-plus-density subdivision
3. Extend to representative `PSP` and `VPS`
   - at minimum `n=1,2,4` for the new rho mode

## Expected Outputs

- Existing sampled dump file, with header lines indicating:
  - `probe_subdiv_n`
  - `probe_rho_mode`
- Campaign-level compare tables keyed by:
  - case
  - `n`
  - `rho_mode`
- Focus metrics:
  - `P_loc_sr_rs_subdiv_probe - P_loc_sr_rs`
  - convergence from `n=1 -> 2 -> 4`
  - difference between `center` and `trilinear`

## Non-Goals

- No changes to the production `P_loc_sr_rs`
- No changes to `grad`
- No changes to total stress assembly
- No history/origlin analysis
