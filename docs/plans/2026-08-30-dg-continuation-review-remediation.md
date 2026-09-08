# DG Continuation Review Remediation Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Correct and simplify the production WF+PW DG continuation before resuming the complete GS checkpoint and RT connection.

**Architecture:** Keep one explicit contained continuation driver and its fixed variational payload.  Batch immutable setup data once, reuse the DC distributed density/potential path, run expensive physical acceptance only for candidate stages, and represent the reconstructed DG residual as one volume channel plus the three actual SIPG face-functional channels.  Preserve every legacy route and remove no shared facility until all callers are proven absent.

**Tech Stack:** Fortran 2008, MPI, SALMON DC/WF+PW/SIPG infrastructure, ScaLAPACK, standalone Python MPI runners.

---

Work only in `/Users/otobetoshihito/SALMON-dev/SALMON2_RTDG/.worktrees/wpw-s-orthogonal-complement`.
Do not create another worktree.  Preserve all pre-existing dirty changes.  Use
`git add -p` for every file that was dirty before a task.  Before each commit,
run `git diff --cached --check` and inspect the complete `git diff --cached`.
Do not use a timeout.

### Task 1: Correct nonlocal strong action and collective failure

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_nonlocal.f90`
- Modify: `src/gs/main_dft.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_nonlocal_mpi.f90`
- Modify: `tests/dg/run_dg_overlapping_wannier_physical_matrices_mpi.py`

**Step 1: Write RED tests**

Add a nonunit-volume fixture that independently evaluates projected
`hvol*rinv_uvu` matrix entries and pointwise `rinv_uvu` action.  Inject missing
projector support on one rank and require communicator-wide rejection without
entering a later collective.

**Step 2: Run RED**

Run `python3 tests/dg/run_dg_overlapping_wannier_physical_matrices_mpi.py`.
Expected: wrong strong-action scale and asymmetric-failure case fail.

**Step 3: Implement minimally**

Return separately named `matrix_strength` and `action_strength`.  Replace all
rank-local early returns before projector collectives with `local_bad` followed
by `MPI_Allreduce(MPI_MAX)` and one collective return.

**Step 4: Run GREEN**

Run the runner at 1, 2, 4, and 8 ranks.  Run `cmake --build build-hybrid-commit -j2`.

**Step 5: Commit**

Commit only these hunks as `fix(dg): correct distributed nonlocal strong action`.

### Task 2: Batch immutable interior materialization

**Files:**
- Modify: `src/gs/dc/dg_hybrid_production_face_traces.f90`
- Modify: `tests/dg/test_dg_hybrid_production_face_traces_mpi.f90`
- Modify: `tests/dg/run_dg_hybrid_production_face_traces_mpi.py`
- Modify: `src/gs/main_dft.f90`

**Step 1: Write RED tests**

Instrument the fixture receipt with request and response counts.  Multiple
basis columns on one fragment must send each destination's point IDs once and
return value, three gradients, and kinetic action in the same response.  Keep
irregular fragments and periodic points.  Require no collective face gather.

**Step 2: Run RED**

Run `python3 tests/dg/run_dg_hybrid_production_face_traces_mpi.py`.
Expected: request-count assertion fails for the basis-by-basis implementation.

**Step 3: Implement minimally**

Replace the two materializers with one grouped routine.  Group effective basis
columns by `(owner, fragment)`, exchange point IDs once per participating peer,
and pack five complex values per point and column.  Delete the separate kinetic
communication loop and update the production call site.

**Step 4: Run GREEN**

Run the face runner at 1, 2, 4, and 8 ranks and the SIPG runner at 1, 2, and 4.

**Step 5: Commit**

Commit as `refactor(dg): batch fixed basis materialization`.

### Task 3: Remove inner-loop fixed work and gate expensive acceptance

**Files:**
- Modify: `src/gs/dc/dg_hybrid_broken_volume.f90`
- Modify: `tests/dg/test_dg_hybrid_broken_volume_mpi.f90`
- Modify: `src/gs/main_dft.f90`
- Modify: `tests/dg/check_dg_hybrid_continuation_route.py`
- Modify: `tests/dg/run_dg_hybrid_continuation_scf_mpi.py`

**Step 1: Write RED tests**

Require a local-potential-only row assembler to match the local block of the
existing combined assembler without evaluating gradients.  Count symmetry,
Hermiticity, and reconstructed-action calls: they must be zero for failed cheap
iterations and exactly once for an accepted candidate.  Force convergence on
the former final SCF iteration and require a separate final refresh to run.

**Step 2: Run RED**

Run the broken-volume, continuation-SCF, and route runners.  Expected: missing
local-only API, excessive acceptance calls, and refresh-boundary failure.

**Step 3: Implement minimally**

Add the local-only projection.  Allocate density and action workspaces outside
the loop.  Use `controller%controls%density_damping`; remove the driver-local
constant.  Compute cheap gates first, invoke expensive gates only for a
candidate, and move lambda-one refresh to one explicit post-loop solve/check.

**Step 4: Run GREEN**

Run all three runners and the full build.

**Step 5: Commit**

Commit as `refactor(dg): simplify continuation candidate checks`.

### Task 4: Measure complete SIPG reconstructed-action residuals

**Files:**
- Modify: `src/gs/dc/dg_hybrid_real_space_residual.f90`
- Modify: `src/gs/dc/dg_hybrid_production_face_traces.f90`
- Modify: `src/gs/main_dft.f90`
- Modify: `tests/dg/test_dg_hybrid_real_space_residual_mpi.f90`
- Modify: `tests/dg/run_dg_hybrid_real_space_residual_mpi.py`
- Modify: `tests/dg/test_dg_hybrid_production_face_traces_mpi.f90`

**Step 1: Write RED tests**

Construct a basis with zero coefficient residual and zero trace-change residual
but a nonzero missing SIPG boundary action.  Fail consistency,
adjoint-consistency, and penalty channels separately.  Require independent
normalized residuals and prove that `R_T=0` cannot make them pass.

**Step 2: Run RED**

Run both real-space and face runners.  Expected: missing face-action residual API.

**Step 3: Implement minimally**

Apply the three existing frozen SIPG face components to reconstructed occupied
traces.  Compare their boundary functionals with the generalized eigen-equation
boundary contribution without inventing a volume lifting.  Keep volume and
three face residuals separate and require all four in candidate acceptance.

**Step 4: Run GREEN**

Run both runners at 1, 2, and 4 ranks, face materialization also at 8 ranks,
then run the SIPG operator regression.

**Step 5: Commit**

Commit as `fix(dg): gate complete reconstructed SIPG action`.

### Task 5: Support zero gaps and fractional occupations

**Files:**
- Modify: `src/gs/main_dft.f90`
- Modify: `src/gs/dc/dg_hybrid_continuation_controller.f90`
- Modify: `tests/dg/test_dg_hybrid_continuation_controller_mpi.f90`
- Modify: `tests/dg/test_dg_hybrid_continuation_scf_mpi.f90`
- Modify: `tests/dg/run_dg_hybrid_continuation_scf_mpi.py`

**Step 1: Write RED tests**

Add an exactly degenerate zero-gap fixture and a fractionally occupied crossing.
Require acceptance when the occupation kernel and `S`-projector are continuous.
Require a shrinking meaningful gap to reduce the next lambda step without
becoming an acceptance gate.  Permit a retained basis with no extra empty state
when the configured occupations need none.

**Step 2: Run RED**

Run controller and continuation runners.  Expected: hard gap and `nstate+1`
requirements reject valid cases.

**Step 3: Implement minimally**

Derive the solved state count from configured occupations/smearing and available
basis size.  Set `occupation_ok` from finite, electron-count-consistent
occupations and subspace continuity.  Remove every positive-gap acceptance
condition.  Compute `gap_shrinking` only when a meaningful occupied/unoccupied
separation exists.

**Step 4: Run GREEN**

Run both runners at all supported decompositions and generalized eigensystem regression.

**Step 5: Commit**

Commit as `fix(dg): continue through occupied subspace crossings`.

### Task 6: Reuse distributed DC Hartree and simplify projector distribution

**Files:**
- Modify: `src/gs/dc/dcdft.f90`
- Modify: `src/gs/dc/dg_overlapping_wannier_nonlocal.f90`
- Modify: `src/gs/main_dft.f90`
- Modify: `tests/dg/check_dg_hybrid_divided_dc_controls.py`
- Modify: `tests/dg/test_dg_overlapping_wannier_nonlocal_mpi.f90`

**Step 1: Write RED tests**

Require the continuation potential update to call the established DC
distributed-density/Hartree path and forbid a full-grid density allocation per
rank.  Add repeated projector keys in nonsorted rank order and require
deterministic canonical ownership without all-rank complete-overlap replication
or quadratic key search.

**Step 2: Run RED**

Run the DC controls contract and physical nonlocal runner.  Expected: replicated
density and Allgatherv/quadratic-key assertions fail.

**Step 3: Implement minimally**

Expose the smallest existing DC density-to-Hartree entry point needed by the
continuation.  Reuse its layouts and halo exchange.  Sort projector keys,
reduce to canonical owners, and send complete overlaps only to row owners and
ranks whose local projector support needs the strong action.

**Step 4: Run GREEN**

Run both focused tests, existing DC route contracts, and full build.

**Step 5: Commit**

Commit as `refactor(dg): reuse distributed DC potential and projectors`.

### Task 7: Remove duplicate checkpoint metric representation

**Files:**
- Modify: `src/common/dg_hybrid_sparse_operators.f90`
- Modify: `src/rt/dg/rt_dg_hybrid_checkpoint.f90`
- Modify: `tests/dg/test_rt_dg_hybrid_checkpoint_mpi.f90`
- Modify: `tests/dg/run_rt_dg_hybrid_checkpoint_mpi.py`

**Step 1: Write RED tests**

Use a sparse metric graph and a different operator-union graph.  Require exact
round-trip of the authoritative metric and Hamiltonian graph without an
operator-side metric payload.  Corrupt either graph and require rejection.

**Step 2: Run RED**

Run `python3 tests/dg/run_rt_dg_hybrid_checkpoint_mpi.py`.
Expected: the format still requires and serializes `operators%metric_values`.

**Step 3: Implement minimally**

Remove operator-side metric serialization from the new complete format and its
hash.  Keep legacy version readers byte-compatible.  Store and hash the metric
and operator graphs independently.  Replace the current uncommitted fixture
that labels a truncated metric with the full metric fingerprint.

**Step 4: Run GREEN and protected regressions**

Run checkpoint at 1, 2, 4, and 8 ranks, occupied-only checkpoint regression,
hybrid metric solver, continuation route contracts, all remediation runners,
and `cmake --build build-hybrid-commit -j2`.

**Step 5: Commit**

Use `git add -p` for the pre-dirty checkpoint files.  Commit as
`refactor(rt): store one authoritative DG metric`.

### Task 8: Review checkpoint and resume the GS-to-RT plan

**Files:**
- Modify only review-driven files if a test first demonstrates a defect.

**Step 1: Request code review**

Use `superpowers:requesting-code-review` on all Task 1--7 commits.  Resolve
Critical and Important findings with TDD and task-scoped commits.

**Step 2: Run fresh remediation verification**

Run every focused runner named above and all protected route contracts without
a timeout.  Preserve user dirty files and generated evidence.

**Step 3: Resume checkpoint Task 8**

Continue at Task 8 of
`docs/plans/2026-08-27-wpw-dg-continuation-ground-state.md`.  Do not claim
overall completion and do not invoke branch finishing until its Tasks 8--12,
including eight-rank Si64 GS-to-RT zero-field acceptance, are complete.
