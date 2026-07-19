# WPW Density-Carrying Seed Implementation Plan

> **For Codex:** REQUIRED SUB-SKILL: Use `executing-plans` to implement this plan task-by-task.

**Goal:** Build a qualified W+P fixed-H checkpoint from the converged DC density-carrying fragment ensemble without running the Flux-inclusive dense LCFO eigensolver.

**Architecture:** Construct W and P overlaps for occupied fragment orbitals, solve the distributed metric equations `S C=B`, rank-qualify the occupied block, and remove it from deterministic extra states. Solve the frozen H0 problem at zero interface strength, continue the complete interface operator to one, and publish only after frozen-state and representation checks pass.

**Tech Stack:** Fortran 2008, MPI, BLAS/LAPACK, SALMON DC-LCFO/WPW bounded operators, Python source-contract tests.

**Repository constraints:** Work in the existing dirty worktree. Do not reset, checkout, recreate Tasks 0–10, commit, push, or open a PR. New runs must use dedicated directories and must not overwrite sample results.

---

### Task 0: Re-establish the protected workspace

**Files:**
- Read: `docs/plans/2026-07-19-wpw-density-carrying-seed-design.md`
- Read: `docs/plans/2026-07-19-wpw-density-carrying-seed-session-handoff.md`

1. Run `pwd`, `git branch --show-current`, `git rev-parse HEAD`, and
   `git status --short`.
2. Require branch `codex/singlescale-vortex-observables` and HEAD
   `125bdb30d1cd1e49ab09835888797c015f75aa6d`.
3. Confirm that the worktree remains intentionally dirty. Do not clean it.
4. Run `git diff --check`.
5. Run `cmake --build build-mpi-eigenexa -j4`.
6. Run the focused tests listed in the handoff document and save their exact
   results before editing.

### Task 1: Repair and qualify the frozen-H snapshot

**Files:**
- Modify: `src/gs/dc/lcfo_flux.f90`
- Test: `tests/dg/check_dg_wpw_fixed_h_relaxation.py`
- Test: `tests/dg/run_dg_wpw_production_operator_mpi.py`

1. Extend the source contract so frozen validation requires WW values and IDs,
   WP/PP volume/nonlocal/face values, bounded H0/interface caches, layout IDs,
   and provenance.
2. Run the contract and confirm RED.
3. Audit the partially added `wpw_frozen_*` arrays from the previous session.
   Do not assume that grouped deallocation or deep assignment is transactional.
4. Add `stat=` allocation paths and collective failure synchronization. Ensure
   every root enters the same production-communicator collectives.
5. Make validation shape-safe: test allocation and shape before comparing
   array values.
6. Add an MPI mutation oracle for one H0 array, one interface array, and one
   transport ID. Each mutation must be detected collectively without deadlock.
7. Verify rollback leaves the last valid operator and callbacks intact.

### Task 2: Replace the false projection oracle with `S C=B`

**Files:**
- Modify: `src/gs/dc/dg_wpw_matrix_free_scf.f90`
- Modify: `tests/dg/test_dg_wpw_matrix_free_scf_mpi.f90`
- Test: `tests/dg/test_dg_wpw_matrix_free_scf.py`

1. Replace the current oracle that normalizes coefficient guesses.
2. Build a two-rank fixture with nonzero W and P rows and a known Hermitian
   positive-definite W+P metric containing WP coupling.
3. Choose known coefficients `C_ref`, form `B=S C_ref`, and require the new
   projection routine to recover the same S-metric occupied projector.
4. Rotate `C_ref` by an occupied unitary matrix and require projector
   invariance.
5. Add a rank-deficient RHS case and require fail-closed behavior.
6. Add a nonfinite RHS on one rank and require collective failure.
7. Use separate input/output buffers in callback tests; never alias an
   `intent(in)` and `intent(out)` actual argument.
8. Run the test and confirm RED.
9. Implement a bounded distributed metric solve for `S C=B`. Report the true
   relative equation residual. A fixed iteration cap must return failure, not a
   partially converged result.
10. Apply rank-revealing S-orthonormalization only after the solve.
11. Run the MPI oracle and expect PASS.

### Task 3: Build both W and P overlap blocks from the DC ensemble

**Files:**
- Modify: `src/gs/dc/lcfo_flux.f90`
- Test: `tests/dg/check_dg_wpw_fixed_h_relaxation.py`
- Test: add or extend an MPI fixture under `tests/dg/`

1. Rename routines and diagnostics so they say `density_carrying_fragment_seed`,
   not `global_lcfo_occupied_projection`.
2. Add a source contract requiring both `b_W` and `b_P` accumulation.
3. Confirm RED because the current code explicitly zeros the P block.
4. Define the source occupation map from the converged DC fragment ensemble.
   Check occupation counts and total charge collectively; do not assume that
   every fragment has the same count.
5. Accumulate `W^dagger phi` and `P^dagger phi` on fragment-owned core points
   with the production quadrature weights.
6. Reduce fragment contributions to fragment roots, then solve `S C=B` across
   the production-root communicator.
7. Record metric residual, captured norm, effective rank, S-orthogonality, and
   charge. Reject nonfinite or inconsistent results.
8. State explicitly in output and manifest that the source is the direct sum
   of DC fragment occupied orbitals, not `coef_wf` from Flux diagonalization.
9. Add a small two-fragment oracle whose P overlap is required for the correct
   answer.

### Task 4: Complete deterministic extra states without polluting occupancy

**Files:**
- Modify: `src/gs/dc/dg_wpw_matrix_free_scf.f90`
- Modify: `src/gs/dc/lcfo_flux.f90`
- Test: `tests/dg/test_dg_wpw_matrix_free_scf_mpi.f90`

1. Add a failing test with an extra seed containing a deliberate occupied
   component.
2. Require `Q_occ^dagger S Q_extra=0` after completion.
3. Require occupied projector invariance before and after extra completion.
4. Rank-qualify only the extra block after occupied projection.
5. Do not run a final full-block normalization that can rotate occupied and
   extra columns together unless the oracle proves occupied-projector
   invariance.
6. Reject insufficient extra rank; do not replace occupied columns randomly.

### Task 5: Implement zero-interface fixed-H relaxation

**Files:**
- Modify: `src/gs/dc/lcfo_flux.f90`
- Modify: `src/gs/dc/dg_wpw_matrix_free_scf.f90`
- Test: `tests/dg/test_dg_wpw_matrix_free_scf_mpi.f90`

1. Set the bounded operator to `lambda_interface=0` transactionally.
2. Validate the frozen state before and after every solver stage.
3. Start block CG/LOBPCG from the qualified density-carrying seed.
4. Keep density, Hartree, XC, local potential, and all H0 blocks fixed.
5. Forbid every call to `wpw_potential_step` in this route.
6. Use occupied trace, maximum generalized residual, S-orthogonality, and
   occupied-projector change as acceptance metrics.
7. Reconstruct density only for diagnostics and establish the zero-interface
   representation baseline.
8. Fail closed at the iteration cap.

### Task 6: Continue the complete interface operator to one

**Files:**
- Modify: `src/common/dg_wpw_bounded_operator.f90`
- Modify: `src/gs/dc/lcfo_flux.f90`
- Test: `tests/dg/test_dg_wpw_production_operator_mpi.f90`

1. Preserve the existing λ=0, 0.5, 1 dense oracle.
2. Add rejected-step and rollback cases.
3. Snapshot accepted coefficients, lambda, solver merit, and frozen invariant
   state before every trial.
4. On rejection, restore the last accepted state and reduce the lambda step.
5. Scale only `ww_face_self`, `ww_cross_value`, and `wp_h_face`.
6. Prove H0 caches and transport IDs are bitwise unchanged across all trials.
7. Require convergence at `lambda_interface=1`.

### Task 7: Qualify and publish the checkpoint

**Files:**
- Modify: `src/common/dg_wpw_checkpoint.f90`
- Modify: `src/gs/dc/lcfo_flux.f90`
- Modify: `tests/dg/check_dg_wpw_checkpoint_publication.py`

1. Extend checkpoint metadata with fixed-H mode, density-carrying seed
   provenance, metric residual, captured norm, projection rank, projection
   charge, final interface lambda, tolerance profile, and frozen fingerprints.
2. Confirm the publication contract is RED.
3. Publish only after lambda one, converged fixed-H residuals, bounded density
   diagnostics, and exact frozen-state validation.
4. Ensure failure leaves no manifest that could be mistaken for a valid
   checkpoint.
5. Run checkpoint read-back and identity validation.

### Task 8: Review and preflight before Si64

1. Run all focused tests listed in the handoff document.
2. Run `cmake --build build-mpi-eigenexa -j4` and `git diff --check`.
3. Perform a findings-first review focused on collective ordering, metric
   projection mathematics, frozen-state coverage, rollback, and provenance.
4. Fix every P0/P1 finding and rerun the focused suite.
5. Run only a short low-cost preflight. Confirm MPI ranks, OpenMP threads,
   estimated memory, output path, and free capacity first.
6. Do not start the full Si64 production checkpoint until the preflight proves
   that projection, zero-interface solve, and at least one continuation trial
   finish with finite diagnostics.

### Task 9: Resume the production sequence

After the fixed-H checkpoint passes review:

1. Run the Si64 2x2x2 DG-DC GS/checkpoint calculation in a dedicated run
   directory.
2. Stop if GS convergence, charge, energy, provenance, or checkpoint identity
   fails.
3. Run checkpoint-backed field-off RT and review finite observables, midpoint
   convergence, S-norm/charge drift, and field-off Pz/Jz drift.
4. Only after field-off acceptance, run the 20-step Sin2 laser smoke from the
   same checkpoint.
5. Do not run or interpret long-time HHG.
