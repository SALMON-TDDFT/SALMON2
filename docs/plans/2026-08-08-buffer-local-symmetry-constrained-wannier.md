# Buffer-Local Symmetry-Constrained Wannier Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add a buffer-local, symmetry-constrained periodic-spread optimizer that turns the accepted projected overlapping basis into a buffer-converged localized Wannier basis without storing full-system KS wavefunctions or a dense global gauge.

**Architecture:** A new pure localization module evaluates periodic localization matrices and applies accepted two-row unitary rotations on a sparse overlapping-pair graph.  Production builds pair and symmetry-pair orbits from the existing buffered tails and exact full-system maps, performs collective monotone sweeps before `S/H/X/V` publication, and records convergence in V3 evidence.  Wannier90 remains an optional small-fixture oracle.

**Tech Stack:** Fortran 2008, MPI, complex Hermitian algebra, existing periodic phase links, Python fixture runners, CMake, spglib, optional Wannier90 oracle.

---

### Task 1: Freeze the localization contract and periodic spread

**Files:**
- Create: `src/gs/dc/dg_overlapping_wannier_localization.f90`
- Create: `tests/dg/test_dg_overlapping_wannier_localization_mpi.f90`
- Create: `tests/dg/run_dg_overlapping_wannier_localization_mpi.py`
- Modify: `src/gs/dc/CMakeLists.txt`
- Modify: `tests/dg/check_dg_overlapping_wannier_route.py`

**Step 1: Write the RED test**

Define the wished-for API `evaluate_dg_periodic_localization(values, weights, phases, norm, moment, spread, ok, message)`.  Test a delta-localized pair, a delocalized pair, phase-origin invariance, nonfinite rejection, and dimension rejection.  Assert the localized gauge has smaller spread.

**Step 2: Verify RED**

Run:

```bash
python3 tests/dg/run_dg_overlapping_wannier_localization_mpi.py
```

Expected: compile failure because `dg_overlapping_wannier_localization` does not exist.

**Step 3: Implement the minimum evaluator**

Evaluate

```text
N_n = sum_p weight_p |W_np|^2
z_an = sum_p weight_p |W_np|^2 phase_ap
Omega = sum_n,a (1 - |z_an/N_n|^2)
```

with strict finite, positive-weight, positive-norm, and shape gates.  Do not add optimization yet.

**Step 4: Focused verification**

Run the fixture on 1, 2, 4, and 8 ranks and the route checker.  Expected: PASS with identical spread on every rank count.

**Step 5: Specification and quality review**

Confirm the functional matches the approved design, is bounded and origin independent, and does not claim finite-buffer MLWF status.  Resolve all Critical/Important findings, run `git diff --check`, and rerun focused verification.

**Step 6: Commit**

```bash
git add src/gs/dc/dg_overlapping_wannier_localization.f90 src/gs/dc/CMakeLists.txt \
  tests/dg/test_dg_overlapping_wannier_localization_mpi.f90 \
  tests/dg/run_dg_overlapping_wannier_localization_mpi.py \
  tests/dg/check_dg_overlapping_wannier_route.py
git commit -m "feat(dg): evaluate buffer-local periodic Wannier spread"
```

### Task 2: Add monotone two-Wannier localization rotations

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_localization.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_localization_mpi.f90`

**Step 1: Write the RED test**

Add `optimize_dg_wannier_pair(values, gradients, weights, phases, first, second, spread_tolerance, rotation, before, after, gradient, accepted, ok, message)`.  Use two deliberately mixed localized functions.  Require a unitary `2x2` rotation, `after < before`, unchanged pair density, rotated gradients, and no change for an already stationary pair.

**Step 2: Verify RED**

Run the localization MPI fixture.  Expected: compile failure for the missing optimizer.

**Step 3: Implement the minimum optimizer**

Parameterize a complex Givens/Jacobi rotation, compute the analytic pair-gradient of the periodic spread, and use a bounded backtracking line search.  Accept only finite rotations satisfying `after <= before + spread_tolerance`.  Rotate values and all three gradient components together.

**Step 4: Focused verification and reviews**

Run 1/2/4/8 ranks.  Review unitarity, phase conventions, stationary handling, and line-search termination.  Resolve Critical/Important findings and rerun.

**Step 5: Commit**

```bash
git add src/gs/dc/dg_overlapping_wannier_localization.f90 \
  tests/dg/test_dg_overlapping_wannier_localization_mpi.f90
git commit -m "feat(dg): minimize periodic spread with local pair rotations"
```

### Task 3: Build sparse pairs and dense-representation symmetry generators

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_localization.f90`
- Modify: `src/gs/dc/dg_overlapping_wannier_symmetry.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_localization_mpi.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_fragment_symmetry_mpi.f90`

**Step 1: Write the RED tests**

Require `build_dg_overlapping_pair_graph` to retain only pairs sharing buffer support above tolerance.  Require `build_dg_symmetry_constrained_generator` to group-average a sparse anti-Hermitian pair seed as `sum_g D(g) K D(g)^dagger / |G|`, preserve a nontrivial dense-mixing representation, reject nonunitary/nonclosed input, and reduce to the original seed for identity-only symmetry.

**Step 2: Verify RED**

Run both MPI fixtures.  Expected: compile failure for the missing graph/orbit APIs.

**Step 3: Implement graph and orbit construction**

Use local support products plus `MPI_Allreduce`; store only `(first, second)` integer pairs.  Build each symmetry-compatible generator using the exact dense Wannier representation from the instantaneous full-system symmetry.  Validate anti-Hermiticity and commutators with every retained `D(g)`.  Restrict exponentiation to the generator's connected support block.  Never restore a parent operation absent from the instantaneous catalog.

**Step 4: Focused verification and reviews**

Run both fixtures on 1/2/4/8 ranks.  Review pair canonicalization, dense orbital mixing, anti-Hermiticity, commutator residuals, deterministic ordering, collective consistency, and linear-memory behavior.  Resolve Critical/Important findings and rerun.

**Step 5: Commit**

```bash
git add src/gs/dc/dg_overlapping_wannier_localization.f90 \
  src/gs/dc/dg_overlapping_wannier_symmetry.f90 \
  tests/dg/test_dg_overlapping_wannier_localization_mpi.f90 \
  tests/dg/test_dg_overlapping_wannier_fragment_symmetry_mpi.f90
git commit -m "feat(dg): constrain local generators by dense symmetry"
```

### Task 4: Implement symmetry-tied monotone localization sweeps

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_localization.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_localization_mpi.f90`

**Step 1: Write the RED test**

Add `localize_dg_overlapping_wannier_basis`.  Test monotone spread history, convergence below the pair-gradient tolerance, metric orthogonality, pair-density conservation, exact symmetry covariance, identity-only operation, and deterministic equality on 1/2/4/8 ranks.  Add a forced-stall case that must reject publication.

**Step 2: Verify RED**

Run the fixture.  Expected: compile failure for the missing sweep API.

**Step 3: Implement the minimum sweep driver**

For each sparse pair seed, build `K_sym`, exponentiate its connected block, and apply the complete symmetry-compatible update transactionally.  Roll it back if total spread rises or symmetry/orthogonality exceeds tolerance.  Stop only on the gradient tolerance or maximum iterations.

**Step 4: Focused verification and reviews**

Run the fixture under bounds checking and floating-point traps on 1/2/4/8 ranks.  Review transactionality, collective failure paths, reproducibility, and memory scaling.  Resolve Critical/Important findings and rerun.

**Step 5: Commit**

```bash
git add src/gs/dc/dg_overlapping_wannier_localization.f90 \
  tests/dg/test_dg_overlapping_wannier_localization_mpi.f90
git commit -m "feat(dg): sweep symmetry-constrained Wannier localization"
```

### Task 5: Integrate localization into the accepted GS route and V3 evidence

**Files:**
- Modify: `src/gs/main_dft.f90`
- Modify: `src/gs/dc/dg_overlapping_wannier_checkpoint.f90`
- Modify: `src/io/salmon_global.f90`
- Modify: `src/io/inputoutput.f90`
- Modify: `tests/dg/check_dg_overlapping_wannier_route.py`
- Modify: `tests/dg/check_dg_fragment_symmetry_production.py`
- Modify: `tests/dg/test_dg_overlapping_wannier_checkpoint_mpi.f90`

**Step 1: Write RED production contracts**

Require localization after global tail materialization and before metric/SCF/operator publication.  Require V3 evidence for initial spread, final spread, maximum pair-gradient, iteration count, convergence, and the exact symmetry-orbit fingerprint.  Require failure on non-convergence.

**Step 2: Verify RED**

Run the route, symmetry-production, and checkpoint fixtures.  Expected: failure because production and V3 lack localization evidence.

**Step 3: Implement production wiring**

Add narrowly scoped OW localization tolerances and maximum iterations to the retained `dc` namelist.  Run the optimizer on `ow_box_values/ow_box_gradients`, refresh centers from periodic moments, recompute the basis fingerprint, and gate checkpoint publication.  Preserve conventional SALMON, normal DC LCFO, EigenExa, and optional Wannier90 paths.

**Step 4: Focused verification and reviews**

Run all three fixtures plus obsolete-route and input contracts.  Review restart compatibility, fingerprint invalidation, zero-iteration rejection, and collective error handling.  Resolve Critical/Important findings and rerun.

**Step 5: Commit**

```bash
git add src/gs/main_dft.f90 src/gs/dc/dg_overlapping_wannier_checkpoint.f90 \
  src/io/salmon_global.f90 src/io/inputoutput.f90 tests/dg
git commit -m "feat(dg): publish symmetry-localized OW ground states"
```

### Task 6: Clean-first overlay build and Si64 GS verification

**Files:**
- Modify as required by failures only
- Record evidence under the existing Si64 focused-test output convention

**Step 1: Create the clean-first overlay**

Archive committed `HEAD`, overlay only the reviewed worktree diff, and configure Release with MPI, ScaLAPACK, EigenExa, spglib ON, and Wannier90 OFF.

**Step 2: Build and run focused GS**

Run the ideal Si64 buffer-supported OW GS on 8 ranks.  Require monotone spread, localization convergence, exact group closure, accepted SCF, and V3 checkpoint publication.  Compare peak memory with the pre-localization baseline and reject any dense full-system wavefunction/gauge allocation.

**Step 3: Specification and quality reviews**

Review the numerical evidence and source diff.  Resolve every Critical/Important item.  Rebuild clean-first and rerun the focused GS after any source change.

**Step 4: Commit**

Commit only reviewed fixes and evidence-contract changes with a focused message.

### Task 7: Polarization-derived linear response and HHG validation

**Files:**
- Modify as required: `tests/dg/run_si64_overlapping_wannier_response_hhg.py`
- Modify as required: `tests/dg/check_si64_overlapping_wannier_response_hhg.py`

**Step 1: Run RT from the accepted V3 checkpoint**

Use generalized-eigenvalue Exp coefficient RT.  Run field-off, weak impulse/linear response, and a long laser pulse with the largest verified stable `dt`.  Keep polarization as the primary observable and current as secondary evidence.

**Step 2: Produce the HHG figure**

Compute the spectrum from polarization, use a semi-log plot, and label harmonic orders.  Do not introduce lattice displacement into the ideal inversion-symmetry gate.

**Step 3: Physical gates**

Require field-off stationarity, cubic-axis agreement, small transverse polarization, odd-harmonic peaks, and suppression of even orders.  Explicitly classify whether each even order is a peak or a dip relative to neighboring bins.  Treat results as qualitative because the system is small.

**Step 4: Reviews and final clean-first verification**

Perform specification and code-quality reviews, resolve all Critical/Important findings, rerun the full focused suite, and run `git diff --check`.

**Step 5: Commit and push**

Commit the verified RT/HHG changes, push `codex/wpw-s-orthogonal-complement` to `origin`, and push the same commit to the configured upstream branch.  Report exact commit IDs, test commands, physical diagnostics, and the absolute path to the semi-log HHG figure.
