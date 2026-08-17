# Wannier90 Symmetry-Adapted Coordinate Optimization Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Replace iterative dense symmetry projection in Wannier90 localization with direct optimization in the complete symmetry-allowed multiplicity space, while preserving exact DMN covariance and bounded memory.

**Architecture:** A one-time analyzer reconstructs the finite-group representation from DMN generators, identifies isotypic multiplicities, and builds contraction/expansion maps for the commutant Lie algebra. Wannier90 then stores conjugate-gradient state in small multiplicity blocks and expands one accepted update per cycle. Delivery is staged: diagnostics first, analytic-versus-legacy comparison second, production activation third, legacy removal last.

**Tech Stack:** Fortran 2008, bundled Wannier90 patching through CMake, LAPACK/BLAS, MPI, Python regression drivers, standalone Wannier90 replay bundles.

---

## Task 1: Freeze the Si64 Replay and Receipt Contract

**Files:**
- Modify: `tests/dg/replay_dg_wannier90_bundle.py`
- Modify: `tests/dg/test_replay_dg_wannier90_bundle.py`
- Modify: `tests/dg/check_sawf_dmn_format.py`
- Reference: `cmakefiles/Builder/patches/apply_wannier90_generator_symmetry.cmake`

**Step 1: Write the failing receipt test**

Require the replay receipt to expose, without changing the optimizer result:

- group order and generator count;
- retained/Wannier dimensions;
- symmetry projector iteration count;
- gradient norm before projection;
- gradient norm after projection;
- covariance residual;
- total spread and localization cycle count;
- peak accounted symmetry workspace.

The test must reject a receipt missing any field and must keep exact numerical fields separate from human-readable log text.

**Step 2: Run the test to verify it fails**

Run: `python3 tests/dg/test_replay_dg_wannier90_bundle.py`

Expected: FAIL because the new diagnostic fields are absent.

**Step 3: Add the minimal parser/export path**

Extend the replay driver and bundled patch diagnostics only. Do not change projection tolerances, iteration limits, search directions, or Wannier updates.

**Step 4: Verify the focused tests**

Run:

```bash
python3 tests/dg/test_replay_dg_wannier90_bundle.py
python3 tests/dg/check_sawf_dmn_format.py
git diff --check
```

Expected: PASS.

**Step 5: Capture the baseline**

Run the existing Si64 replay once with the current legacy projector. Store the machine-readable receipt outside the source tree or under the existing ignored replay-output directory. Record the exact input fingerprint and executable fingerprint.

**Step 6: Commit**

```bash
git add tests/dg/replay_dg_wannier90_bundle.py tests/dg/test_replay_dg_wannier90_bundle.py tests/dg/check_sawf_dmn_format.py cmakefiles/Builder/patches/apply_wannier90_generator_symmetry.cmake
git commit -m "test: freeze Wannier90 symmetry replay receipts"
```

## Task 2: Add a Small Exact Commutant Diagnostic

**Files:**
- Modify: `cmakefiles/Builder/patches/apply_wannier90_generator_symmetry.cmake`
- Modify: `tests/dg/check_sawf_dmn_format.py`
- Create: `tests/dg/test_w90_symmetry_allowed_coordinates.py`

**Step 1: Write exact representation fixtures**

Add Z2, Z3, Z2 x Z2, and S3 fixtures covering:

- multiplicity one with no nontrivial internal rotation;
- two copies of one irrep with known `U(2)` freedom;
- a purely forbidden anti-Hermitian gradient;
- a purely allowed anti-Hermitian gradient;
- mixed allowed/forbidden components;
- malformed product, nonunitary generator, and nonassociative catalog rejection.

For these small fixtures only, build the linear constraints
`X D(g) - D(g) X = 0` and compute their null space. Treat this as an oracle, not the production algorithm.

**Step 2: Run the new test and confirm RED**

Run: `python3 tests/dg/test_w90_symmetry_allowed_coordinates.py`

Expected: FAIL because the diagnostic API does not exist.

**Step 3: Implement the test-only diagnostic**

Add a bounded helper that returns:

- allowed real Lie-algebra dimension;
- orthonormal coordinate basis for small N;
- reconstruction residual;
- allowed and forbidden gradient norms;
- deterministic fingerprint.

Reject dimensions above the explicit fixture limit before allocating the constraint matrix.

**Step 4: Verify**

Run:

```bash
python3 tests/dg/test_w90_symmetry_allowed_coordinates.py
python3 tests/dg/check_sawf_dmn_format.py
git diff --check
```

Expected: PASS, including MPI/rank-independent receipts where the fixture invokes MPI.

**Step 5: Commit**

```bash
git add cmakefiles/Builder/patches/apply_wannier90_generator_symmetry.cmake tests/dg/check_sawf_dmn_format.py tests/dg/test_w90_symmetry_allowed_coordinates.py
git commit -m "test: define exact Wannier symmetry coordinate oracle"
```

## Task 3: Implement Scalable Isotypic Analysis

**Files:**
- Modify: `cmakefiles/Builder/patches/apply_wannier90_generator_symmetry.cmake`
- Modify: `tests/dg/test_w90_symmetry_allowed_coordinates.py`
- Modify: `tests/dg/check_sawf_dmn_format.py`

**Step 1: Add failing scalable-analysis tests**

Require the production analyzer to reproduce the oracle multiplicities and allowed dimensions for all small fixtures. Add adverse tests for rank-disagreeing metadata, generator-word mismatch, representation-product mismatch, allocation failure, and checked extent overflow.

Add a memory assertion that persistent storage contains generators, component bases, multiplicity blocks, and at most one streamed full retained-space matrix; it must not contain a dense `N x N x |G|` tensor.

**Step 2: Run the tests and confirm RED**

Run: `python3 tests/dg/test_w90_symmetry_allowed_coordinates.py`

Expected: FAIL because only the small null-space oracle exists.

**Step 3: Implement deterministic group traversal**

Validate the associative finite-group catalog and deterministic generator words. Stream each represented group element from generators, verifying unitarity and products without retaining every full matrix.

Use checked wide extents, `STAT=` allocation handling, and collective MPI agreement before shape-dependent branches.

**Step 4: Implement isotypic and multiplicity coordinates**

Construct isotypic projectors from streamed characters, determine ranks with tolerance-separated spectral gates, and resolve equivalent copies into a tensor-product basis. Publish:

- irrep dimensions and multiplicities;
- total allowed real Lie-algebra dimension;
- reconstruction and generator covariance defects;
- persistent/transient accounted bytes;
- deterministic catalog/decomposition fingerprint.

Fail closed when a cluster boundary is ambiguous. Do not choose a LAPACK-dependent basis inside an unresolved multiplicity block; preserve that block as the optimization freedom.

**Step 5: Compare with the exact oracle**

Run:

```bash
python3 tests/dg/test_w90_symmetry_allowed_coordinates.py
python3 tests/dg/check_sawf_dmn_format.py
git diff --check
```

Expected: PASS with exact dimensions and tolerance-level reconstruction agreement.

**Step 6: Commit**

```bash
git add cmakefiles/Builder/patches/apply_wannier90_generator_symmetry.cmake tests/dg/check_sawf_dmn_format.py tests/dg/test_w90_symmetry_allowed_coordinates.py
git commit -m "feat: analyze Wannier symmetry multiplicity coordinates"
```

## Task 4: Add Analytic Gradient Contraction and Expansion

**Files:**
- Modify: `cmakefiles/Builder/patches/apply_wannier90_generator_symmetry.cmake`
- Modify: `tests/dg/test_w90_symmetry_allowed_coordinates.py`
- Modify: `tests/dg/replay_dg_wannier90_bundle.py`
- Modify: `tests/dg/test_replay_dg_wannier90_bundle.py`

**Step 1: Write contraction/expansion REDs**

For every exact fixture, require:

- allowed gradients survive contraction then expansion;
- forbidden gradients contract to zero;
- mixed gradients equal the exact oracle projection;
- expanded directions are anti-Hermitian and commute with all generators;
- contraction and expansion are adjoint under the Frobenius inner product;
- results are invariant under unitary changes of irrep-copy basis.

**Step 2: Implement block contraction**

For an isotypic block arranged as multiplicity x irrep, contract the dense gradient over the irrep index to an anti-Hermitian multiplicity block. Do not divide inside the innermost loop; precompute reciprocal irrep dimensions. Keep contiguous indices innermost and use BLAS where it avoids temporary transposes.

**Step 3: Implement one-shot expansion**

Expand a multiplicity search direction as `X_alpha tensor I` and transform it once to retained coordinates. Measure covariance after expansion; do not iteratively repair it.

**Step 4: Add compare-only mode**

During one localization cycle, compute both:

- the analytic allowed-space gradient;
- the legacy converged Reynolds projection.

Return their norm difference, allowed/forbidden fractions, timing, and workspace. Keep the legacy result as the applied update in this phase.

**Step 5: Verify fixtures and replay tooling**

Run:

```bash
python3 tests/dg/test_w90_symmetry_allowed_coordinates.py
python3 tests/dg/test_replay_dg_wannier90_bundle.py
python3 tests/dg/check_sawf_dmn_format.py
git diff --check
```

Expected: PASS.

**Step 6: Run Si64 diagnostic and stop for review**

Run the frozen Si64 replay in compare-only mode. Report:

- multiplicity table and allowed dimension;
- pre-projection, analytic-allowed, and forbidden gradient norms;
- analytic-versus-legacy difference;
- time and memory for both paths;
- covariance residual.

Stop here. Do not activate the new optimizer until the result demonstrates that useful symmetry-allowed freedom exists and the analytic projection agrees with the exact/legacy definitions within the specified tolerance.

**Step 7: Commit the diagnostic phase**

```bash
git add cmakefiles/Builder/patches/apply_wannier90_generator_symmetry.cmake tests/dg/test_w90_symmetry_allowed_coordinates.py tests/dg/replay_dg_wannier90_bundle.py tests/dg/test_replay_dg_wannier90_bundle.py
git commit -m "feat: compare analytic Wannier symmetry gradients"
```

## Task 5: Optimize Directly in Multiplicity Blocks

**Prerequisite:** Explicit approval after the Task 4 Si64 review checkpoint.

**Files:**
- Modify: `cmakefiles/Builder/patches/apply_wannier90_generator_symmetry.cmake`
- Modify: `tests/dg/test_w90_symmetry_allowed_coordinates.py`
- Modify: `tests/dg/replay_dg_wannier90_bundle.py`

**Step 1: Add optimizer REDs**

Require a known multiplicity-two fixture to lower a controlled localization objective while preserving covariance and unitarity. Require multiplicity-one/no-freedom fixtures to terminate cleanly without fake progress.

**Step 2: Store optimizer state in allowed coordinates**

Move conjugate-gradient history, preconditioning, direction updates, and line-search direction norms into the multiplicity blocks. Dense retained-space matrices may exist only for the current accepted update and existing Wannier90 overlap evaluation.

**Step 3: Apply accepted updates once**

Exponentiate each small anti-Hermitian multiplicity block, expand the block unitary, and apply it. Measure covariance and unitarity after application. A failed measurement is an error, not a request for another projection loop.

**Step 4: Verify exact fixtures**

Run:

```bash
python3 tests/dg/test_w90_symmetry_allowed_coordinates.py
python3 tests/dg/check_sawf_dmn_format.py
git diff --check
```

Expected: PASS.

**Step 5: Run the Si64 replay**

Compare with the frozen baseline:

- final spread must decrease when the allowed gradient is nonzero;
- covariance must stay within tolerance every cycle;
- no dense symmetry-projection iteration loop may execute;
- per-cycle symmetry time must be bounded by retained contractions plus small blocks;
- accounted memory must remain within the one-stream-workspace contract.

**Step 6: Commit**

```bash
git add cmakefiles/Builder/patches/apply_wannier90_generator_symmetry.cmake tests/dg/test_w90_symmetry_allowed_coordinates.py tests/dg/replay_dg_wannier90_bundle.py
git commit -m "feat: optimize Wannier functions in symmetry coordinates"
```

## Task 6: Remove the Legacy Iterative Projector and Integrate Production

**Files:**
- Modify: `cmakefiles/Builder/patches/apply_wannier90_generator_symmetry.cmake`
- Modify: `tests/dg/check_sawf_dmn_format.py`
- Modify: `tests/dg/check_dg_overlapping_wannier_route.py`
- Modify: `tests/dg/run_si64_overlapping_wannier_response_hhg.py`
- Modify: `docs/plans/2026-08-18-wannier90-symmetry-adapted-coordinate-optimization-design.md`

**Step 1: Add a route RED**

Reject source containing the legacy 100/1000-cycle gradient symmetrization loop or duplicated dense projector state. Require setup-time decomposition followed by block-coordinate localization.

**Step 2: Remove superseded code**

Delete the iterative Reynolds repair path, iteration-cap tuning, and diagnostics specific to convergence of that loop. Retain one-shot covariance validation.

**Step 3: Run focused verification**

Run:

```bash
python3 tests/dg/test_w90_symmetry_allowed_coordinates.py
python3 tests/dg/test_replay_dg_wannier90_bundle.py
python3 tests/dg/check_sawf_dmn_format.py
python3 tests/dg/check_dg_overlapping_wannier_route.py
python3 tests/dg/run_dg_overlapping_wannier_w90_mpi.py
git diff --check
```

Expected: PASS. Run MPI fixtures at 1/2/4/8 ranks using the existing runner configuration; do not increase rank count or OMP threads during debugging.

**Step 4: Run production Si64 through Wannier completion**

Use the same validated input and MPI rank count as the baseline. Capture RSS per rank, cycle timing, multiplicity receipts, spread, covariance, and the post-Wannier transition. Do not proceed to a larger or longer RT workload until Wannier completion and its immediate post-processing are stable.

**Step 5: Update the design status**

Append measured Si64 evidence and mark which performance and failure contracts were demonstrated. Document any intentionally deferred optimization without weakening correctness criteria.

**Step 6: Final verification**

Run the focused suite again from a clean build, inspect `git diff --check`, and confirm that only intended files are staged. Do not add generated `.wout` files or unrelated dirty worktree files.

**Step 7: Commit**

```bash
git add cmakefiles/Builder/patches/apply_wannier90_generator_symmetry.cmake tests/dg/check_sawf_dmn_format.py tests/dg/check_dg_overlapping_wannier_route.py tests/dg/run_si64_overlapping_wannier_response_hhg.py docs/plans/2026-08-18-wannier90-symmetry-adapted-coordinate-optimization-design.md
git commit -m "refactor: remove iterative Wannier symmetry projection"
```

## Non-Negotiable Constraints

- Do not raise iteration limits as a correctness fix.
- Do not weaken symmetry tolerance to obtain convergence.
- Do not retain all full group matrices or an `N x N x |G|` tensor.
- Do not introduce one Wannier90 invocation per character sector.
- Do not infer freedom from Wannier centers; derive it from the validated representation.
- Do not activate production updates before the Task 4 Si64 review checkpoint.
- Keep MPI rank count unchanged during comparisons and set OMP explicitly in every run receipt.
- Treat allocation, extent, metadata-agreement, and LAPACK failures collectively.
