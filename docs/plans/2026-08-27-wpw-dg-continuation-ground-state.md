# WF+PW DG Continuation Ground-State Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Starting from the exact converged DC density, converge the occupied WF+PW subspace, density, and complete SIPG interface observables through an adaptive lambda continuation to the fully self-consistent lambda-one DG ground state, publish one complete atomic checkpoint, and start stationary zero-field hybrid RT from that exact payload.

**Architecture:** Freeze one symmetry-closed WF+PW catalog for a continuation attempt, preassemble the coefficient-independent complete SIPG interface blocks, and update only the density-dependent volume operator and projector-derived state.  A transactional controller jointly converges density, occupied `S`-projector, and interface observables at each uniform lambda, with adaptive step rejection and complete rollback.  The accepted lambda-one operator and state are serialized together and consumed directly by a new, isolated hybrid branch in `main_tddft.f90`.

**Tech Stack:** Fortran 2008, MPI, SALMON DC/Wannier90/WF+PW infrastructure, SIPG weak form, ScaLAPACK or EigenExa generalized eigensolver, BLAS/LAPACK, standalone Python MPI runners.

---

Use the fixed worktree
`/Users/otobetoshihito/SALMON-dev/SALMON2_RTDG/.worktrees/wpw-s-orthogonal-complement`.
Do not create another worktree.  Preserve every pre-existing uncommitted change.
The runners in `tests/dg` are invoked directly; do not register them with
CTest and do not modify `tests/CMakeLists.txt`.

For every task, first add a failing test and record its RED result, then add
only the minimum production implementation and record the GREEN result.  Stage
only the files named by that task.  If a named file was dirty before the task,
use `git add -p` and select only task-specific hunks.  Before every commit run
both `git diff --cached --check` and `git diff --cached`, and remove any
unrelated hunk from the index.  Never discard, restore, or overwrite user
changes.

The old DC+LCFO/Wannier90, overlapping-Wannier, ordinary GS, and ordinary RT
branches are protected paths.  Every new production call must be inside the
explicit hybrid WF+PW DG-continuation flag.

### Task 1: Specify the DC-density seed and fixed catalog

**Files:**
- Create: `tests/dg/test_dg_hybrid_continuation_state_mpi.f90`
- Create: `tests/dg/run_dg_hybrid_continuation_state_mpi.py`
- Create: `src/gs/dc/dg_hybrid_continuation_state.f90`
- Modify: `src/gs/dc/CMakeLists.txt`

**Step 1: Write the failing seed test**

Construct a distributed DC density, a symmetry-closed basis catalog, and a
different density reconstructible from a trial coefficient matrix.  Require
initialization to copy the DC density exactly, preserve its fingerprint, set
lambda to zero, and reject duplicate core ownership, zero fingerprints, or a
catalog mutation after initialization.

Require

```fortran
seed_residual = norm2(seed_density-dc_density)/max(1d0,norm2(dc_density))
```

to be zero within the configured seed tolerance.  Explicitly fail an
implementation that initializes from the trial coefficient density.

**Step 2: Run RED**

Run: `python3 tests/dg/run_dg_hybrid_continuation_state_mpi.py`

Expected: compile failure because `dg_hybrid_continuation_state` is absent.

**Step 3: Implement the minimal state type**

Add types for the immutable catalog fingerprint set, mixed density, occupied
projector, interface observables, lambda controller state, epochs, residual
receipt, and an accepted snapshot.  Initialization takes `dc_density` as an
explicit required argument and never calls a density reconstruction callback.

**Step 4: Run GREEN**

Run: `python3 tests/dg/run_dg_hybrid_continuation_state_mpi.py`

Expected: PASS at 1, 2, and 4 ranks with identical global fingerprints.

**Step 5: Commit**

Stage the three new files and only the module-list hunk from
`src/gs/dc/CMakeLists.txt`; inspect the cached diff and commit:

`feat(dg): seed continuation from converged DC density`

### Task 2: Assemble complete WF+PW SIPG interface blocks

**Files:**
- Create: `src/gs/dc/dg_hybrid_sipg_operator.f90`
- Create: `tests/dg/test_dg_hybrid_sipg_operator_mpi.f90`
- Create: `tests/dg/run_dg_hybrid_sipg_operator_mpi.py`
- Modify: `src/gs/dc/CMakeLists.txt`

**Step 1: Write the failing two-fragment test**

Use analytic linear basis functions on two adjacent fragments.  Compare all
four minus/minus, minus/plus, plus/minus, and plus/plus blocks against a dense
reference containing consistency, adjoint-consistency, and penalty terms.
Require nonzero off-diagonal coupling, Hermiticity, reciprocal-face
cancellation, exactly one canonical owner, and uniform scaling of every block
by one scalar lambda.  Include a periodic face and reject face-local lambdas.

**Step 2: Run RED**

Run: `python3 tests/dg/run_dg_hybrid_sipg_operator_mpi.py`

Expected: compile failure because the projected SIPG assembler is absent.

**Step 3: Implement the minimum assembler**

Reuse `evaluate_dg_nodal_sipg_face`.  Exchange fixed basis values and outward
normal derivatives, assemble both row directions for each canonical face, and
store the result in row-owned sparse form.  Keep consistency, adjoint, and
penalty diagnostic norms separately.  Do not use current eigenvector
coefficients while constructing the operator.

**Step 4: Run GREEN**

Run:

`python3 tests/dg/run_dg_hybrid_sipg_operator_mpi.py`

Expected: PASS at 1, 2, and 4 ranks; the fixture directly exercises the
existing nodal SIPG evaluator as well as the new projected assembler.

**Step 5: Commit**

Commit only the task files as `feat(dg): assemble complete hybrid SIPG blocks`.

### Task 3: Add occupied-projector and interface residual algebra

**Files:**
- Create: `src/gs/dc/dg_hybrid_continuation_residuals.f90`
- Create: `tests/dg/test_dg_hybrid_continuation_residuals_mpi.f90`
- Create: `tests/dg/run_dg_hybrid_continuation_residuals_mpi.py`
- Modify: `src/gs/dc/CMakeLists.txt`

**Step 1: Write failing gauge-invariance tests**

For a positive-definite complex metric, form an occupied cluster and apply
independent column phases and a dense unitary rotation inside a degenerate
cluster.  Require invariant occupied `S`-projector and projector-change
residual.  Construct the distinct occupation-weighted density matrix
`Gamma=C f C^dagger` and use it for density, electron number, energy, and
occupied face density matrices.  Require invariance under rotations within an
equally occupied degenerate cluster and reject symmetry-incompatible unequal
occupations.  Demonstrate that raw coefficient differences are nonzero and
must not be used.  Compare `R_H`, `R_rho`, `R_T`, and `R_S` with dense
references and reject nonfinite or rank-deficient inputs.

**Step 2: Run RED**

Run: `python3 tests/dg/run_dg_hybrid_continuation_residuals_mpi.py`

Expected: compile failure because the residual module is absent.

**Step 3: Implement minimum projector/trace operations**

Implement distributed `S`-metric occupied-projector comparison, the distinct
occupation-weighted density matrix, optional Procrustes alignment for
deterministic output, gauge-invariant face density matrices, and independently
normalized residual channels.  Do not implement raw eigenvector mixing.

**Step 4: Run GREEN**

Run the new runner at 1, 2, and 4 ranks; expect PASS and rank-invariant norms.

**Step 5: Commit**

Commit only the task files as `feat(dg): measure projector and interface fixed points`.

### Task 4: Implement transactional adaptive lambda control

**Files:**
- Create: `src/gs/dc/dg_hybrid_continuation_controller.f90`
- Create: `tests/dg/test_dg_hybrid_continuation_controller_mpi.f90`
- Create: `tests/dg/run_dg_hybrid_continuation_controller_mpi.py`
- Modify: `src/gs/dc/CMakeLists.txt`

**Step 1: Write failing controller fixtures**

Test monotone inexact tolerances, acceptance only when every residual channel
passes, bounded step growth after easy stages, and rejection on residual
growth, gap collapse, crossing, projector discontinuity, or symmetry failure.
Mutate density, potential, projector, trace, occupation, eigenvalue, mixing
history, and derived-cache epochs before rejection; require bitwise restoration
of the accepted snapshot and step reduction.  Require the next trial to apply
one identical lambda to all faces.

**Step 2: Run RED**

Run: `python3 tests/dg/run_dg_hybrid_continuation_controller_mpi.py`

Expected: compile failure because the controller is absent.

**Step 3: Implement the minimum state machine**

Implement propose, accept, reject/restore, tolerance scheduling, rollback
limits, and collective decision routines.  Density and interface damping are
separate controls with `alpha_trace <= alpha_density` by default.  Keep their
histories and residuals separate.

**Step 4: Run GREEN**

Run the new runner at 1, 2, and 4 ranks; expect identical decisions and PASS.

**Step 5: Commit**

Commit only the task files as `feat(dg): control adaptive DG continuation`.

### Task 5: Build the coupled fixed-point driver

**Files:**
- Create: `src/gs/dc/dg_hybrid_continuation_scf.f90`
- Create: `tests/dg/test_dg_hybrid_continuation_scf_mpi.f90`
- Create: `tests/dg/run_dg_hybrid_continuation_scf_mpi.py`
- Modify: `src/gs/dc/CMakeLists.txt`

**Step 1: Write the failing synthetic SCF test**

Use a small density-dependent generalized eigenproblem with a known lambda-one
fixed point.  Start from a supplied DC density and occupied projector.  Require
the callback order `volume -> full H/S solve -> projector -> density/trace ->
residuals -> separate mixing`, no lambda advance while any gate fails, at
least one forced rollback, and convergence to the dense reference.  Rotate the
occupied eigenvectors randomly at every solve to prove gauge stability.

**Step 2: Run RED**

Run: `python3 tests/dg/run_dg_hybrid_continuation_scf_mpi.py`

Expected: compile failure because the coupled driver is absent.

**Step 3: Implement the callback-driven driver**

Compose the state, residual, and controller modules.  Preserve the exact DC
density before the first volume build.  Mix density and trace separately; use
the occupied projector for tracking.  At lambda one, perform an unmixed full
refresh and require every final gate to pass again.

**Step 4: Run GREEN and existing solver regression**

Run:

```text
python3 tests/dg/run_dg_hybrid_continuation_scf_mpi.py
python3 tests/dg/run_dg_hybrid_generalized_eigensystem_mpi.py
python3 tests/dg/run_dg_hybrid_scf_mpi.py
```

Expected: all PASS.

**Step 5: Commit**

Commit only the task files as `feat(dg): converge coupled DG fixed points`.

### Task 6: Add symmetry covariance and real-space DG residual gates

**Files:**
- Create: `src/gs/dc/dg_hybrid_continuation_acceptance.f90`
- Create: `tests/dg/test_dg_hybrid_continuation_acceptance_mpi.f90`
- Create: `tests/dg/run_dg_hybrid_continuation_acceptance_mpi.py`
- Modify: `src/gs/dc/CMakeLists.txt`

**Step 1: Write failing acceptance tests**

Construct symmetry-related faces, occupied clusters, and retained basis
blocks.  Require covariance of `H_volume`, `H_interface`, `S`, and the complete
occupied projector for every operation, while allowing individual degenerate
eigenvectors to mix.  Require symmetry-compatible occupations in each
degenerate occupied block.  Measure retained-space leakage
`||(1-Q_ret)D(g)Q_ret||`; show that a face-local lambda, a basis missing one
symmetry partner, and a cutoff splitting a multiplet or reciprocal star fail.
Also show that a smaller but symmetry-complete excitation space passes the
closure test, while being reported separately as not proving observable-level
excitation-cutoff convergence.  Add a truncated basis whose coefficient
residual is zero but whose reconstructed real-space complete-DG residual is
large; require rejection.

**Step 2: Run RED**

Run: `python3 tests/dg/run_dg_hybrid_continuation_acceptance_mpi.py`

Expected: compile failure because the acceptance oracle is absent.

**Step 3: Implement minimum acceptance oracles**

Evaluate normalized covariance defects of the zero-field operators, retained
space, occupied projector, and occupation density matrix using the verified
basis representation.  Never require an individual eigenvector to be
invariant.  Add a callback that reconstructs selected occupied and near-gap
states and applies volume plus complete SIPG action in real space.  Aggregate
maxima collectively and fail closed on omitted operations, split symmetry
blocks, or faces.

**Step 4: Run GREEN**

Run the new runner and existing fragment-symmetry runner; expect PASS.

**Step 5: Commit**

Commit only the task files as `feat(dg): gate symmetry and real-space DG residuals`.

### Task 7: Replace the one-shot divided branch with continuation

**Files:**
- Modify: `src/gs/main_dft.f90`
- Modify: `src/gs/dc/dcdft.f90`
- Modify: `src/io/salmon_global.f90`
- Modify: `src/io/inputoutput.f90`
- Modify: `tests/dg/check_dg_hybrid_divided_dc_controls.py`
- Create: `tests/dg/check_dg_hybrid_continuation_route.py`

**Step 1: Write failing source-contract tests**

Require the explicit hybrid branch to pass the converged DC total density and
its fingerprint into continuation initialization, construct the complete SIPG
operator, run the coupled continuation, and publish only after the lambda-one
refresh.  Forbid `solve_dg_hybrid_generalized_once_and_publish` and occupied-only
checkpoint publication in this branch.  Require all legacy branches to retain
their former calls and flag conditions.

**Step 2: Run RED**

Run:

```text
python3 tests/dg/check_dg_hybrid_continuation_route.py
python3 tests/dg/check_dg_hybrid_divided_dc_controls.py
```

Expected: the new contract FAILS on the one-shot branch.

**Step 3: Wire the production callbacks minimally**

Reuse the accepted WF+PW basis, DC potential and density infrastructure, and
distributed generalized solver.  Assemble complete cross-fragment SIPG rows.
Keep the catalog frozen for the attempt.  Do not alter DC+LCFO/Wannier90 or
overlapping-Wannier branches.

**Step 4: Run GREEN and route regressions**

Run the two contracts plus
`check_dg_hybrid_self_consistent_route.py` and the existing
overlapping-Wannier route check.  Expected: all PASS.

**Step 5: Commit**

These files are already dirty.  Use `git add -p` for every modified file,
inspect both cached checks, and commit only new continuation hunks as
`feat(dg): run WF+PW DG continuation from DC density`.

### Task 8: Publish one complete atomic GS checkpoint

**Files:**
- Modify: `src/rt/dg/rt_dg_hybrid_checkpoint.f90`
- Modify: `tests/dg/test_rt_dg_hybrid_checkpoint_mpi.f90`
- Modify: `tests/dg/run_rt_dg_hybrid_checkpoint_mpi.py`
- Modify: `src/gs/main_dft.f90`

**Step 1: Extend the checkpoint test and observe RED**

Require one file to round-trip the exact sparse `H_DG(0)`, `S_DG`, basis
ownership/catalog, cutoff/selection metadata, coefficients, occupations,
eigenvalues, density, interface observables, DC seed fingerprint,
continuation receipt, and all fingerprints.  Flip one byte independently in
each payload class and require rejection.  Require an interrupted write to
leave the previous accepted file intact.

Run: `python3 tests/dg/run_rt_dg_hybrid_checkpoint_mpi.py`

Expected: FAIL because the existing format lacks the full state/catalog.

**Step 2: Extend the versioned atomic format minimally**

Hash metadata and every serialized payload value together.  Write to a unique
temporary file, close and verify it, then atomically rename.  The reader must
return the serialized matrices rather than regenerate them.  Keep old occupied
checkpoint routines available for protected legacy callers.

**Step 3: Wire publication after the final refresh only**

Replace the occupied-only writer only inside the new hybrid continuation
branch.  A failed final gate must not write or replace a checkpoint.

**Step 4: Run GREEN**

Run the checkpoint runner at 1, 2, and 4 ranks and the route contract; expect
PASS.

**Step 5: Commit**

Use partial staging for dirty `main_dft.f90`; commit only task hunks as
`feat(dg): checkpoint complete DG ground state`.

### Task 9: Connect the exact checkpoint payload to hybrid RT

**Files:**
- Modify: `src/rt/main_tddft.f90`
- Create: `src/rt/dg/rt_dg_hybrid_initialization.f90`
- Modify: `src/rt/CMakeLists.txt`
- Create: `tests/dg/test_rt_dg_hybrid_initialization_mpi.f90`
- Create: `tests/dg/run_rt_dg_hybrid_initialization_mpi.py`

**Step 1: Write the failing GS-to-RT identity test**

Write a complete synthetic GS checkpoint and initialize RT from it.  Require
bitwise-identical sparse matrix payloads and matching complete-payload
fingerprint.  At RT startup re-evaluate generalized residual, `S`
orthogonality, electron number, Hermiticity, zero-field basis/operator
covariance, and covariance of the complete occupied projector or occupation
density matrix.  Do not test individual eigenvector symmetry.  Reject a test
that reconstructs numerically equal matrices under a different payload
identity.

**Step 2: Run RED**

Run: `python3 tests/dg/run_rt_dg_hybrid_initialization_mpi.py`

Expected: compile or source-contract failure because `main_tddft.f90` has no
hybrid checkpoint branch.

**Step 3: Implement isolated hybrid initialization**

Add an explicit default-off hybrid RT branch before ordinary initialization.
Read the complete checkpoint once, validate its exact payload, redistribute
the stored basis/state according to its catalog, and initialize the existing
hybrid metric solver and propagator with the stored `H_DG(0)` and `S_DG`.

**Step 4: Run GREEN and legacy RT regression**

Run the new runner, `run_rt_dg_hybrid_checkpoint_mpi.py`,
`run_rt_dg_hybrid_metric_solver_mpi.py`, and
`run_rt_dg_hybrid_length_gauge_mpi.py`; expect all PASS.

**Step 5: Commit**

Commit only the task files as `feat(rt): load exact hybrid DG ground state`.

### Task 10: Add zero-field RT stationarity acceptance

**Files:**
- Create: `src/rt/dg/rt_dg_hybrid_stationarity.f90`
- Modify: `src/rt/CMakeLists.txt`
- Create: `tests/dg/test_rt_dg_hybrid_stationarity_mpi.f90`
- Create: `tests/dg/run_rt_dg_hybrid_stationarity_mpi.py`
- Modify: `src/rt/main_tddft.f90`

**Step 1: Write failing stationary and phase-rotation tests**

Propagate a known generalized eigenstate with zero external field.  Require
bounded drift in density, total energy, occupied `S`-projector, electron
number, and DG Hamiltonian residual.  Apply arbitrary occupied
phases and a degenerate-space rotation between samples; require the projector
test to pass while a deliberately changed occupied subspace fails.

**Step 2: Run RED**

Run: `python3 tests/dg/run_rt_dg_hybrid_stationarity_mpi.py`

Expected: compile failure because the stationarity evaluator is absent.

**Step 3: Implement the minimum evaluator and sampling hook**

Reuse the projector/residual algebra from GS.  Record initial invariants from
the checkpoint payload and compare them at configured RT samples.  Do not use
raw coefficient differences as an acceptance measure.  Do not add a separate
RT symmetry-drift gate: the initial zero-field operator and occupied-space
symmetry were already accepted, and stationary projector/density checks cover
the zero-field case.  Driven-state symmetry is outside this plan.

**Step 4: Run GREEN**

Run the new runner at 1, 2, and 4 ranks and the existing hybrid RT runners;
expect PASS.

**Step 5: Commit**

Commit only task files as `test(rt): verify zero-field hybrid stationarity`.

### Task 11: Add the Si64 end-to-end runner and protected-route regression

**Files:**
- Create: `tests/dg/data/si64_overlapping_wannier_rt/input_hybrid_dg_continuation.in`
- Create: `tests/dg/data/si64_overlapping_wannier_rt/input_hybrid_dg_zero_field_rt.in`
- Create: `tests/dg/run_dg_hybrid_si64_continuation_rt.py`
- Modify: `tests/dg/check_dg_hybrid_continuation_route.py`

**Step 1: Write RED acceptance parsing first**

Make the runner require lambda-one convergence receipts for every GS residual,
the DC seed-density identity, exact GS-to-RT payload identity, zero-field GS
operator/occupied-space symmetry at handoff, and zero-field stationarity
receipts for density, total energy, occupied projector, electron number, and
DG residual.  Do not require a separate RT-state symmetry receipt.  Feed it the preserved
old one-shot log and confirm rejection because it has no complete DG
continuation or RT payload evidence.

**Step 2: Add only the new hybrid input files**

Do not modify the already-dirty `input_hybrid_scf.in`.  Set the explicit new
GS and RT flags, zero external field for RT, and production final tolerances.
The runner invokes Python/SALMON directly and is not registered with CTest.

**Step 3: Run focused source/input checks GREEN**

Run the runner's parser-only mode and the route contract.  Expected: PASS for
synthetic complete receipts and rejection of incomplete receipts.

**Step 4: Run protected-route regressions**

Run the existing DC+LCFO/Wannier90 and overlapping-Wannier source/fixture
runners specified by their current Python entry points.  Expected: PASS with
no changed legacy branch behavior.

**Step 5: Commit**

Commit only the four task files as `test(dg): add Si64 DG continuation RT acceptance`.

### Task 12: Fresh verification, review, and branch completion

**Files:**
- No production changes expected.
- Preserve outputs under a new `verification-si64-dg-continuation-<date>/`
  directory without staging generated data.

**Step 1: Invoke verification-before-completion**

Read and follow `superpowers:verification-before-completion`.  Run every new
focused runner and every protected-route regression freshly; do not rely on a
previous result.

**Step 2: Run the final Si64 GS-to-RT calculation**

From a clean output directory run the new end-to-end Python runner with:

```text
MPI ranks: 8
OMP_NUM_THREADS=1
timeout: none
```

Do not use `timeout`, `gtimeout`, scheduler wall-time termination, or a Python
subprocess timeout.  Preserve the complete log whether the run passes, fails,
or is manually interrupted.

**Step 3: Verify the final evidence**

Require simultaneous evidence for:

- exact DC density as the lambda-zero seed;
- accepted adaptive stages ending at lambda one;
- final refreshed `R_H`, `R_rho`, `R_T`, `R_S`, electron count, symmetry, and
  real-space DG residual;
- complete checkpoint payload and matching GS/RT fingerprint;
- symmetry closure of the retained WF+PW space and covariance of the complete
  occupied ground-state projector, without individual-state symmetry;
- zero-field density, energy, projector, electron, and Hamiltonian
  stationarity, without a separate RT-state symmetry gate.

Any missing receipt means the implementation is not complete.

**Step 4: Request code review**

Use `superpowers:requesting-code-review`, resolve all critical or important
findings with TDD and task-scoped commits, then rerun affected verification.

**Step 5: Complete the branch workflow**

Only after all fresh checks pass, invoke
`superpowers:finishing-a-development-branch` and present its integration
options.  Do not claim completion before this step.
