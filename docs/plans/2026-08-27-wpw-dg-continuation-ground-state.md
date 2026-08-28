# WF+PW DG Continuation Ground-State Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Starting from the exact converged DC density, converge the occupied WF+PW subspace, density, and complete SIPG interface observables through an adaptive lambda continuation to the fully self-consistent lambda-one DG ground state, publish one complete atomic checkpoint, and start stationary zero-field hybrid RT from that exact payload.

**Architecture:** Freeze one WF+PW catalog closed under the actual full-system symmetry group and one metric for a continuation attempt; when no nontrivial symmetry exists, use the identity group rather than a separate or disabled path.  One concrete solver owns the mutable fixed-point state, the last accepted state, the MPI collective schedule, rollback, and acceptance.  It preassembles the coefficient-independent complete SIPG interface blocks, updates the density-dependent volume operator, refreshes occupied projector/density/interface traces after every solve, and advances one uniform adaptive lambda to a fully refreshed lambda-one fixed point.  No generic callback controller or production adapter is introduced.  The complete basis/operator/state payload is consumed by an isolated RT branch that updates Hartree/XC once at the start of every explicit time step.

**Tech Stack:** Fortran 2008, MPI, SALMON DC/Wannier90/WF+PW infrastructure, SIPG weak form, ScaLAPACK or EigenExa generalized eigensolver, BLAS/LAPACK, standalone Python MPI runners.

---

Use the fixed worktree
`/Users/otobetoshihito/SALMON-dev/SALMON2_RTDG/.worktrees/wpw-s-orthogonal-complement`.
Do not create another worktree.  Preserve every pre-existing uncommitted change.
Invoke the Python runners in `tests/dg` directly.

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

Use the explicit new flags `yn_dg_hybrid_continuation_scf` and
`yn_rt_dg_hybrid_continuation`; do not reinterpret an existing legacy flag.
The accepted file name is `hybrid_dg_ground_state.chk`.

The planned public interfaces are fixed before implementation:

```fortran
call initialize_dg_hybrid_continuation(comm, dc_density, dc_density_fingerprint, &
  catalog, dc_occupied_coefficients, occupations, supported_scope, state, ok, message)
call assemble_dg_hybrid_sipg_operator(comm, catalog, penalty_factor, &
  interface_operator, diagnostics, ok, message)
call run_dg_hybrid_continuation_scf(comm, catalog, state, controls, &
  accepted_state, ok, message)
call write_rt_dg_hybrid_ground_state_checkpoint(comm, path, payload, &
  payload_fingerprint, ok, message)
call initialize_rt_dg_hybrid_from_checkpoint(comm, path, rt_state, ok, message)
call update_rt_dg_hybrid_density_operator(comm, rt_state, density, &
  hamiltonian, ok, message)
```

The checkpoint format receives a new version and magic; the existing hybrid
and occupied-only readers remain unchanged for their current callers.

### Task 1: Specify the DC-density seed and fixed catalog

**Files:**
- Create: `tests/dg/test_dg_hybrid_continuation_state_mpi.f90`
- Create: `tests/dg/run_dg_hybrid_continuation_state_mpi.py`
- Create: `src/gs/dc/dg_hybrid_continuation_state.f90`
- Modify: `src/common/dg_hybrid_production_pw_basis.f90`
- Modify: `src/common/dg_hybrid_windowed_pw_basis.f90`
- Modify: `src/common/dg_hybrid_windowed_pw_types.f90`
- Modify: `src/common/dg_hybrid_reciprocal_catalog.f90`
- Modify: `tests/dg/test_dg_hybrid_production_pw_basis_mpi.f90`
- Modify: `src/gs/dc/CMakeLists.txt`
- Modify: `src/common/CMakeLists.txt`

**Step 1: Write the failing seed test**

Construct a distributed DC density, a symmetry-closed basis catalog, and a
different density reconstructible from a trial coefficient matrix.  Require
initialization to copy the DC density exactly, preserve its fingerprint, set
lambda to zero, and reject duplicate core ownership, zero fingerprints, or a
catalog mutation after initialization.

Give the initial selection a cutoff that splits a reciprocal star and a WF
symmetry multiplet.  Require deterministic expansion to the complete orbit
before catalog freezing, followed by recomputed ownership, distribution, and
fingerprints.  Require failure when the supplied action maps cannot produce a
finite closed selection.

Add an identity-only catalog with inequivalent fragments, unequal local basis
dimensions, irregular face geometry, and no nonidentity operation.  Require
normal initialization, no added orbit member, and execution of the same
closure path rather than a symmetry-check bypass.  Record zero nonidentity
operations, successful authoritative symmetry analysis and its provenance,
and `identity_only=.true.` in its receipt.  Separately provide an unexecuted
analysis, failed analysis, absent provenance, empty operation list, and
malformed group; require collective rejection rather than identity fallback.

Add a collective closure routine taking the requested WF block IDs and their
group action, plus requested PW packet IDs and their packet action.  It returns
the sorted effective IDs, a parent/reason entry for every added member, and a
fingerprint covering both requested and effective selections.  Pass these
actions from the already accepted basis symmetry representation; do not infer
WF multiplets from eigenvalue proximity.

Normalize a successfully completed authoritative analysis with no nontrivial
operation to one explicit identity action before closure.  Do not require
equivalent fragments, equal fragment basis dimensions, uniform face areas, or
a regular neighbor graph for that identity-only case.  For a known nontrivial
physical group, require fragment geometry, topology, basis actions, and face
orbits to be covariant; reject a fragmentation that loses an operation and do
not relabel it as identity-only.

Require collective failure before seed initialization unless
`theory=='dft'`, `system%Nspin==1`, `yn_spinorbit=='n'`,
`.not.PLUS_U_ON`, `yn_hse=='n'`, `yn_fix_func=='n'`,
`yn_jm=='n'`, the boundary is periodic, and all active
`xc_func%xctype` entries are the built-in PZ, PZM, PW, or PBE constants.
Hash every selector and `xctype` entry into the supported-scope receipt.

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
projector, freshly derived interface observables, lambda controller state,
epochs, residual receipt, an immutable DC seed snapshot, and a separate
accepted snapshot.  Initialization takes `dc_density` as an explicit required
argument and a collective supported-scope receipt, never calls a density
reconstruction helper, and marks lambda zero
as unaccepted until its full fixed-point gates pass.  Reuse the existing
reciprocal-star and basis-action maps to close the retained selection before
freezing the catalog.  Store requested cutoff/selection separately from the
effective retained blocks, added orbit members, and their closure action map.
Do not change the basis during continuation.

**Step 4: Run GREEN**

Run: `python3 tests/dg/run_dg_hybrid_continuation_state_mpi.py`

Expected: PASS at 1, 2, and 4 ranks with identical global fingerprints.

**Step 5: Commit**

Stage only the Task 1 files.  Any already-dirty basis file must use
`git add -p` for its symmetry-closure hunks; likewise stage only the new
module-list hunks from both CMake files.  Inspect the cached diff and commit:

`feat(dg): seed continuation from converged DC density`

### Task 2: Assemble complete WF+PW SIPG interface blocks

**Files:**
- Create: `src/gs/dc/dg_hybrid_sipg_operator.f90`
- Create: `tests/dg/test_dg_hybrid_sipg_operator_mpi.f90`
- Create: `tests/dg/run_dg_hybrid_sipg_operator_mpi.py`
- Modify: `src/common/dg_hybrid_sparse_metric.f90`
- Modify: `src/common/dg_hybrid_sparse_operators.f90`
- Modify: `tests/dg/test_dg_hybrid_sparse_operators_mpi.f90`
- Modify: `src/gs/dc/CMakeLists.txt`

**Step 1: Write the failing two-fragment test**

Use analytic complex linear basis functions on two adjacent fragments.  Fix
one canonical normal from minus to plus, define jumps and averages with that
same normal on both traces, and compare all
four minus/minus, minus/plus, plus/minus, and plus/plus blocks against a dense
reference containing consistency, adjoint-consistency, and penalty terms.
The reference is derived independently from SALMON's
`-0.5d0*nabla^2`: require an outer factor `0.5d0` on consistency,
adjoint-consistency, and penalty contributions exactly once, with
`penalty_factor` interpreted as the dimensionless eta inside the bracket.
Require nonzero off-diagonal coupling, Hermiticity, reciprocal-face
cancellation, exactly one canonical owner, and uniform scaling of every block
by one scalar lambda.  Include a periodic face and reject face-local lambdas.
Include an interface Hamiltonian entry for which the metric entry is zero and
require a metric CSR graph independent of the operator-union CSR graph.  All
Hamiltonian components and position operators share a coupling-envelope graph
constructed from basis support, nonlocal support, and face topology, using
explicit component zeros rather than separate exchange graphs.  Change a
generic component value on an initially zero allowed edge and require it to
use the unchanged graph and exchange schedule; the physical Hartree/XC
density-perturbation test belongs to Task 9.  Add a crossing nonlocal-projector
fixture and require exactly-once volume accounting, separate from SIPG kinetic
faces.

**Step 2: Run RED**

Run: `python3 tests/dg/run_dg_hybrid_sipg_operator_mpi.py`

Expected: compile failure because the projected SIPG assembler is absent.

**Step 3: Implement the minimum assembler**

Reuse `evaluate_dg_nodal_sipg_face`.  Convert its two outward-normal traces to
the documented canonical-normal convention, apply complex conjugation to the
test-function trace and derivative, assemble both row directions for each
canonical face, and
store the result in row-owned sparse form.  Keep consistency, adjoint, and
penalty diagnostic norms separately.  Do not use current eigenvector
coefficients while constructing the operator.  Give the metric and operator
coupling envelope independent row offsets, columns, structure fingerprints, and
communication schedules.  Do not require the metric graph to contain every
SIPG Hamiltonian edge and do not create a separate graph per component.  The
nodal evaluator returns the unscaled bracket action; multiply both its value
and normal-action outputs by `0.5d0` exactly once when assembling the SALMON
kinetic operator.  Keep raw nodal diagnostic energies labeled bracket units;
multiply consistency, adjoint, and penalty diagnostics by `0.5d0` in the
physical projected receipt and checkpoint.  Assert both raw and physical
penalty values in the fixture.

**Step 4: Run GREEN**

Run:

`python3 tests/dg/run_dg_hybrid_sipg_operator_mpi.py`

Then run: `python3 tests/dg/run_dg_hybrid_sparse_operators_mpi.py`

Expected: PASS at 1, 2, and 4 ranks; the fixture directly exercises the
existing nodal SIPG evaluator as well as the new projected assembler.

**Step 5: Commit**

Commit only the task files as `feat(dg): assemble complete hybrid SIPG blocks`.

### Task 3: Add occupied-projector and interface residual algebra

**Files:**
- Create: `src/common/dg_hybrid_continuation_residuals.f90`
- Create: `tests/dg/test_dg_hybrid_continuation_residuals_mpi.f90`
- Create: `tests/dg/run_dg_hybrid_continuation_residuals_mpi.py`
- Modify: `src/common/CMakeLists.txt`

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

Use a nonorthogonal dense representation satisfying `D^dagger S D=S` and
require `D Q D^-1=Q` for the mixed-index occupied projector.  Separately test
`D^dagger H D=H` and the selected coefficient convention for the occupation
kernel.  For `C -> D C`, require explicitly
`Gamma -> D Gamma D^dagger`; reject use of the projector similarity rule for
`Gamma`.  This pins forward versus pullback transformations.

**Step 2: Run RED**

Run: `python3 tests/dg/run_dg_hybrid_continuation_residuals_mpi.py`

Expected: compile failure because the residual module is absent.

**Step 3: Implement minimum projector/trace operations**

Implement shared GS/RT distributed `S`-metric occupied-projector comparison, the distinct
occupation-weighted density matrix, optional Procrustes alignment for
deterministic output, gauge-invariant face density matrices, and independently
normalized residual channels.  Do not implement raw eigenvector mixing.

**Step 4: Run GREEN**

Run the new runner at 1, 2, and 4 ranks; expect PASS and rank-invariant norms.

**Step 5: Commit**

Commit only the task files as `feat(dg): measure projector and interface fixed points`.

### Task 4: Implement the internal adaptive lambda state machine

**Files:**
- Create: `src/gs/dc/dg_hybrid_continuation_controller.f90`
- Create: `tests/dg/test_dg_hybrid_continuation_controller_mpi.f90`
- Create: `tests/dg/run_dg_hybrid_continuation_controller_mpi.py`
- Modify: `src/gs/dc/CMakeLists.txt`

**Step 1: Write failing controller fixtures**

Test monotone inexact tolerances, acceptance only when every residual channel
passes, bounded step growth after easy stages, and rejection on residual
growth, crossing with failed cluster-aware occupation, projector
discontinuity, or symmetry failure.  A shrinking gap alone reduces the next
step but does not reject an otherwise acceptable stage.
Mutate density, potential, projector, derived trace cache, occupation, eigenvalue, mixing
history, and derived-cache epochs before rejection; require bitwise restoration
of the accepted snapshot and step reduction.  Require the next trial to apply
one identical lambda to all faces.

Pin the defaults and exact decisions: initial/minimum/maximum step
`0.125/0.015625/0.5`, growth/shrink `1.5/0.5`, residual-growth limit `4`,
density damping `0.5`, minimum projector overlap `0.9`, and eight rollbacks.
For every residual use the documented linear-in-lambda intermediate tolerance
clamped below by its final tolerance.  Treat a small gap by cluster-aware
occupation and projector overlap; do not add a separate gap-cutoff rejection.
After the first inner iteration, define each residual-growth ratio within the
same trial against `max(previous_inner_channel, numerical_floor)` and use
the maximum channel ratio for the decision.  Require growth above the limit
for two consecutive inner updates before rejection.  Use the previous
accepted stage only to propose the initial lambda step.  An easy stage is one
that converges within half of the configured iteration limit without rollback
and passes the projector-overlap gate; only then may the step grow.  Read the
iteration limit and intermediate/final channel tolerances from `controls` and
pin these transitions in the fixture rather than adding another controller.

**Step 2: Run RED**

Run: `python3 tests/dg/run_dg_hybrid_continuation_controller_mpi.py`

Expected: compile failure because the controller is absent.

**Step 3: Implement the minimum state machine**

Implement propose, accept, reject/restore, tolerance scheduling, rollback
limits, and collective decision routines.  Density damping is the only
nonlinear mixing control.  Interface traces are invalidated on every state
change and rebuilt from the current occupation density matrix; keep `R_T` as
an independent residual history without an `alpha_trace` control.  Document
this as the explicit global-SIPG-matrix replacement of trace damping.

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

Use a small density-dependent generalized eigenproblem with known lambda-zero
and lambda-one fixed points.  Start from a supplied DC density deliberately
different from the lambda-zero fixed point and an occupied projector.  Require
the very first Hamiltonian build to see the exact DC density, but forbid a
lambda increase until lambda zero itself passes every acceptance gate.  Require
the operation order `volume -> full H/S solve -> projector -> density/trace ->
residuals -> density mixing`, no lambda advance while any gate fails, at
least one forced rollback, and convergence to the dense reference.  Rotate the
occupied eigenvectors randomly at every solve to prove gauge stability.

**Step 2: Run RED**

Run: `python3 tests/dg/run_dg_hybrid_continuation_scf_mpi.py`

Expected: compile failure because the coupled driver is absent.

**Step 3: Implement the concrete fixed-point driver**

Compose the state, residual, and controller modules.  Preserve the exact DC
density before the first volume build.  Mix density only, use the occupied
projector for tracking, and fully refresh interface traces from `Gamma` after
each solve.  Numerical test kernels may be passed privately by the fixture,
but they are not a production backend API.  At lambda one, perform an unmixed
full refresh and require every final gate to pass again.

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
- Create: `src/common/dg_hybrid_continuation_acceptance.f90`
- Create: `tests/dg/test_dg_hybrid_continuation_acceptance_mpi.f90`
- Create: `tests/dg/run_dg_hybrid_continuation_acceptance_mpi.py`
- Modify: `src/gs/dc/dg_hybrid_continuation_scf.f90`
- Modify: `tests/dg/test_dg_hybrid_continuation_scf_mpi.f90`
- Modify: `src/common/CMakeLists.txt`

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
residual is zero but whose reconstructed-grid residual under the actual
discrete DG action is large; require rejection.

Add an identity-only, non-equivalent-fragment fixture.  Require the complete
candidate-acceptance path to run and pass its identity covariance checks;
reject implementations that treat zero nonidentity operations as either an
error or permission to omit the symmetry check.  Keep the nontrivial-group
fixtures to prove that identity normalization does not weaken real symmetry
enforcement.

Run the coupled driver with a fixture that passes coefficient-space gates but
fail first the symmetry oracle and then the reconstructed-grid oracle.  Require
the candidate stage to remain unaccepted in both cases; this test fails if the
production driver can bypass either oracle.

**Step 2: Run RED**

Run:

```text
python3 tests/dg/run_dg_hybrid_continuation_acceptance_mpi.py
python3 tests/dg/run_dg_hybrid_continuation_scf_mpi.py
```

Expected: the acceptance runner fails because the oracle is absent, and the
SCF runner fails because production candidate acceptance can bypass it.

**Step 3: Implement minimum acceptance oracles**

Evaluate normalized covariance defects of the zero-field operators, retained
space, occupied projector, and occupation density matrix using the verified
basis representation.  Never require an individual eigenvector to be
invariant.  For an otherwise acceptable stage, use the concrete helper that
reconstructs every occupied state, applies the actual discrete volume and
complete SIPG action, lifts it to the production grid, and evaluates the
existing real-space quadrature norm.  Report the three face contributions
separately with the existing face weights; do not invent an additional
combined DG norm.  Repeat this expensive check after the final lambda-one
refresh.  Aggregate maxima collectively and fail closed on omitted
operations, split symmetry blocks, or faces.
The operation list must contain at least the identity.  Zero nonidentity
operations is valid and is reported distinctly from a missing or malformed
operation list.  Acceptance also requires successful authoritative symmetry
analysis and its provenance; missing required data or analysis failure cannot be
represented as an identity-only result.

Wire this oracle directly into `dg_hybrid_continuation_scf`.  Invoke it only
after the inexpensive inner gates pass and again after the final lambda-one
refresh.  Missing payload or a failed oracle rejects the stage collectively.

**Step 4: Run GREEN**

Run:

```text
python3 tests/dg/run_dg_hybrid_continuation_acceptance_mpi.py
python3 tests/dg/run_dg_hybrid_continuation_scf_mpi.py
python3 tests/dg/check_dg_fragment_symmetry_production.py
```

Expected: all PASS.

**Step 5: Commit**

Commit only the task files as `feat(dg): gate symmetry and real-space DG residuals`.

### Task 7a: Materialize immutable production face traces

**Files:**
- Modify: `src/gs/dc/dg_hybrid_fragment_basis.f90`
- Create: `src/gs/dc/dg_hybrid_production_face_traces.f90`
- Create: `tests/dg/test_dg_hybrid_production_face_traces_mpi.f90`
- Create: `tests/dg/run_dg_hybrid_production_face_traces_mpi.py`
- Modify: `src/gs/dc/CMakeLists.txt`

**Step 1: Write the failing production face test**

Construct two irregular neighboring fragments and one periodic image face.
Require one canonical owner, one minus-to-plus normal, complete value and
normal-derivative traces for every retained basis function, and a collective
fingerprint covering topology, geometry, basis IDs, quadrature, values, and
derivatives.  Feed the payload to `assemble_dg_hybrid_sipg_face` and require a
nonzero Hermitian cross-fragment block.  Reject duplicate ownership, missing
neighbors, inconsistent periodic shifts, incomplete point correspondence,
nonfinite traces, and an effective selection that is not closed under the
accepted action.

**Step 2: Run RED**

Run: `python3 tests/dg/run_dg_hybrid_production_face_traces_mpi.py`

Expected: compile failure because `dg_hybrid_production_face_traces` is absent.

**Step 3: Implement the minimum immutable payload**

Add an optional face-trace component to the fragment basis and build it only
for the explicit continuation route.  Reconstruct gradients with the existing
SALMON stencil, convert both sides to one canonical normal, validate exactly-once
face ownership collectively, and freeze the payload and its fingerprint.
Legacy callers remain valid without the optional component.

**Step 4: Run GREEN**

Run:

```text
python3 tests/dg/run_dg_hybrid_production_face_traces_mpi.py
python3 tests/dg/run_dg_hybrid_sipg_operator_mpi.py
```

Expected: PASS at 1, 2, and 4 ranks.

**Step 5: Commit**

Stage only these Task 7a files and commit as
`feat(dg): materialize production SIPG face traces`.

### Task 7a.5: Expose the production symmetry and selection boundary

**Files:**
- Create: `docs/plans/2026-08-29-dg-production-selection-boundary-design.md`
- Modify: `src/common/dg_hybrid_windowed_pw_types.f90`
- Modify: `src/common/dg_hybrid_production_pw_basis.f90`
- Modify: `tests/dg/test_dg_hybrid_production_pw_basis_mpi.f90`
- Modify: `tests/dg/run_dg_hybrid_production_pw_basis_mpi.py`

**Step 1: Write failing MPI tests**

Require an authoritative production symmetry receipt and exported fragment,
spatial, reciprocal, and PW-packet actions.  Test a nontrivial group and an
explicit identity-only group.  Require a known nonidentity physical operation
that does not map whole fragments to fail instead of being silently removed.
Request one member of a multi-packet orbit, close it with
`close_dg_hybrid_selection`, and require catalog freezing to accept the closed
effective IDs and reject the original non-closed IDs.

**Step 2: Run RED**

Run `python3 tests/dg/run_dg_hybrid_production_pw_basis_mpi.py`.

Expected: compilation fails because the phased production analysis and
effective-selection API is absent.

**Step 3: Implement the minimal phased API**

Add a receipt/action type to `dg_hybrid_windowed_pw_types`.  Separate
fragment-covariance analysis and packet-universe construction from catalog
freezing.  Stable packet IDs enumerate `(fragment, reciprocal star)` in the
existing order.  Validate that effective IDs are unique, in the universe, and
closed under the returned packet action.  Recompute catalog ownership and
selection fingerprints from effective IDs.  Preserve the existing convenience
entry point by selecting the complete universe.

**Step 4: Run GREEN and regressions**

Run:

```text
python3 tests/dg/run_dg_hybrid_production_pw_basis_mpi.py
python3 tests/dg/check_dg_hybrid_windowed_pw_route.py
python3 tests/dg/run_dg_hybrid_production_face_traces_mpi.py
```

Expected: all PASS.

**Step 5: Commit**

Use `git add -p` for every already-dirty file.  Run
`git diff --cached --check` and inspect `git diff --cached`.  Commit only this
task as `feat(dg): expose production symmetry selection boundary`.

### Task 7a.6: Build and connect one concrete production continuation solver

**Files:**
- Create: `docs/plans/2026-08-29-dg-concrete-continuation-solver-design.md`
- Modify: `src/gs/dc/dg_hybrid_continuation_scf.f90`
- Modify: `tests/dg/test_dg_hybrid_continuation_scf_mpi.f90`
- Modify: `tests/dg/run_dg_hybrid_continuation_scf_mpi.py`
- Modify: `src/gs/main_dft.f90`
- Modify: `src/gs/dc/dcdft.f90`
- Modify: `src/io/salmon_global.f90`
- Modify: `src/io/inputoutput.f90`
- Modify: `tests/dg/check_dg_hybrid_divided_dc_controls.py`
- Create: `tests/dg/check_dg_hybrid_continuation_route.py`
- Delete before commit if still untracked: `src/gs/dc/dg_hybrid_production_continuation_adapter.f90`
- Delete before commit if still untracked: `tests/dg/test_dg_hybrid_production_continuation_adapter_mpi.f90`
- Delete before commit if still untracked: `tests/dg/run_dg_hybrid_production_continuation_adapter_mpi.py`

**Step 1: Preserve the rejected experiment and write the full-path failing test**

Keep the current test output and review findings in the existing verification
records.  Do not commit the experimental adapter.  Extend the continuation
fixture so it exercises the same state-transition routine used by the
contained production driver, rather than manually calling phases.  Use a
nonorthogonal two-fragment problem with nonzero cross-fragment SIPG blocks.
Require the physical basis-space projector
`C_occ C_occ^dagger S`, uniform lambda on every canonical face, complete
rollback, and a fully refreshed lambda-one state.

**Step 2: Run RED**

Run `python3 tests/dg/run_dg_hybrid_continuation_scf_mpi.py`.

Expected: FAIL because the existing callback path does not own a complete
production state transition and currently forms the wrong projector.

Add a source-contract RED test requiring the explicit continuation branch to
pass the exact converged DC density and the already assembled fixed `S` and
complete `H_interface` payload.  Forbid a runtime catalog builder, adapter,
callback table, one-shot final
LCFO, and occupied-only checkpoint publication in that branch.

**Step 3: Add rank-local failure, stale-state, and route RED cases**

Inject a volume-kernel failure on one rank and require communicator-wide
failure without deadlock.  Make Hermiticity, electron number, symmetry, and
each numeric residual fail independently.  Rotate the occupied eigenvectors
by phases and a degenerate-space unitary and require projector invariance.
Reject any final state whose density, trace, operator, or epoch predates the
last solve.

Run both the MPI fixture and `check_dg_hybrid_continuation_route.py`.  The MPI
fixture must fail on the missing concrete catalog/solver contract and the
route test must fail because production is not connected.

**Step 4: Implement the minimum concrete solver**

Complete symmetry closure, basis construction, metric assembly, and SIPG
assembly before entering SCF.  Keep only the completed fixed arrays and their
fingerprints, one current mutable state, and one deep copy of the last accepted
state.  Implement one contained production continuation driver where the
existing SALMON state is already available.  Do not add a runtime catalog
builder, `class(*)`, abstract backend, or procedure table.  Call the existing
assembly, distributed eigensolver, density, trace, and residual routines in
one fixed order.  Convert every rank-local
failure to collective consensus before entering another collective.  Mix only
density.  Evaluate acceptance as a pure operation on the current fully
refreshed state.  Roll back the whole state atomically.

Inside the explicit default-off continuation branch, build the supported-scope
receipt, close WF blocks and PW packets, materialize only the effective
selection, recompute ownership and fingerprints, finish the fixed basis and
matrix payload, and start the contained driver from the exact converged
`dc%rho_tot`.
Do not alter protected routes and do not publish a checkpoint in this task.

Selection closure and face assembly retain their existing focused tests.  The
production task only connects their completed outputs; it does not add another
selection or materialization layer.

**Step 5: Run focused RED then GREEN at all decompositions**

Run the continuation fixture at 1, 2, 4, and 8 ranks.  Before implementation
record failure for the new cases; after the minimal implementation require all
cases to pass.

**Step 6: Run focused regressions**

Run:

```text
python3 tests/dg/run_dg_hybrid_continuation_scf_mpi.py
python3 tests/dg/run_dg_hybrid_continuation_acceptance_mpi.py
python3 tests/dg/run_dg_hybrid_sipg_operator_mpi.py
python3 tests/dg/run_dg_hybrid_production_face_traces_mpi.py
python3 tests/dg/check_dg_hybrid_divided_dc_controls.py
python3 tests/dg/check_dg_hybrid_continuation_route.py
python3 tests/dg/check_dg_hybrid_self_consistent_route.py
python3 tests/dg/check_dg_overlapping_wannier_route.py
```

Expected: all PASS.

**Step 7: Commit**

Use `git add -p` for dirty files.  Run `git diff --cached --check` and inspect
`git diff --cached`.  Commit only this task as
`feat(dg): consolidate concrete production continuation solver`.

### Task 8: Publish one complete atomic GS checkpoint

**Files:**
- Modify: `src/rt/dg/rt_dg_hybrid_checkpoint.f90`
- Modify: `tests/dg/test_rt_dg_hybrid_checkpoint_mpi.f90`
- Modify: `tests/dg/run_rt_dg_hybrid_checkpoint_mpi.py`
- Modify: `src/gs/main_dft.f90`

**Step 1: Extend the checkpoint test and observe RED**

Require one file to round-trip the exact sparse `H_DG(0)`, `S_DG`, the metric
CSR graph and operator-union CSR graph, actual distributed basis values, grid IDs and
weights, partition data, face values and normals, nonlocal distribution,
ownership/catalog, requested cutoff/selection, effective symmetry-closed
selection, added orbit members and closure action/reason metadata, coefficients, occupations,
eigenvalues, density, interface observables, separately identified fixed and
initial Hartree/XC Hamiltonian components, DC seed fingerprint,
continuation receipt, exchange-correlation functional, pseudopotential and
energy-decomposition provenance, and all fingerprints.  Corrupt one representative value
from metadata, basis, matrix, and state payloads and require rejection.  Require
an interrupted write to leave the previous accepted file intact.
Include authoritative-analysis completion and provenance, the actual-group
operation count, nonidentity-operation count, and identity-only normalization
flag in the hashed catalog metadata.

Run: `python3 tests/dg/run_rt_dg_hybrid_checkpoint_mpi.py`

Expected: FAIL because the existing format lacks the full state/catalog.

**Step 2: Extend the versioned atomic format minimally**

Hash metadata and every serialized payload value together.  Write to a unique
temporary file, close and verify it, then atomically rename.  The reader must
return the serialized matrices and basis payload rather than regenerate them.
Remove the existing requirement that metric and Hamiltonian share identical
row degrees and column IDs.  Hash the metric graph, operator-union graph, and
all basis data, including requested-versus-effective selection provenance.
Keep old occupied checkpoint routines available for protected
legacy callers.

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
- Modify: `src/rt/dg/rt_dg_hybrid_length_gauge.f90`
- Modify: `src/rt/dg/rt_dg_hybrid_sparse_exchange.f90`
- Create: `src/rt/dg/rt_dg_hybrid_initialization.f90`
- Create: `src/rt/dg/rt_dg_hybrid_density_update.f90`
- Modify: `src/rt/CMakeLists.txt`
- Create: `tests/dg/test_rt_dg_hybrid_initialization_mpi.f90`
- Create: `tests/dg/run_rt_dg_hybrid_initialization_mpi.py`
- Modify: `tests/dg/test_rt_dg_hybrid_length_gauge_mpi.f90`

**Step 1: Write the failing GS-to-RT identity test**

Write a complete synthetic GS checkpoint and initialize RT from it.  Require
bitwise identity of the canonical serialized global payload and its
complete-payload fingerprint.  After MPI redistribution, compare rows by
global row and column IDs and require exact values, without requiring local CSR
arrays from different rank counts to be bitwise identical.  At RT startup
re-evaluate generalized residual, `S`
orthogonality, electron number, Hermiticity, zero-field basis/operator
covariance, and covariance of the complete occupied projector or occupation
density matrix.  Do not test individual eigenvector symmetry.  Reject a test
that reconstructs numerically equal matrices under a different payload
identity.

In the same fixture, rebuild Hartree/XC from the stored initial density and
require the resulting `H_DG(0)` to match the stored payload.  Perturb the
density after initialization and require the RT Hamiltonian to change; this
must fail for an implementation that freezes `H_DG(0)` throughout RT.
Require RT initialization to reject a checkpoint whose supported-scope
receipt is absent or requests spin, spin-orbit, DFT+U, HSE/exact exchange, or
another state-dependent Hamiltonian channel.  Require local RT selectors to
match the hashed receipt, with `theory` restricted to `tddft_response` or
`tddft_pulse`, `yn_fix_func=='n'`, `yn_jm=='n'`, periodic
boundaries, and the same allowed built-in `xc_func%xctype` entries.

**Step 2: Run RED**

Run: `python3 tests/dg/run_rt_dg_hybrid_initialization_mpi.py`

Expected: compile or source-contract failure because `main_tddft.f90` has no
hybrid checkpoint branch.

**Step 3: Implement isolated hybrid initialization**

Add an explicit default-off hybrid RT branch before ordinary initialization.
Read the complete checkpoint once, validate its exact payload, redistribute
the stored basis/state according to its catalog, and initialize the existing
hybrid metric solver and propagator with the stored `H_DG(0)` and `S_DG`.
Separate stored time-independent kinetic/SIPG/ionic/nonlocal components from
Hartree/XC.  Before time step zero, run the normal density-dependent potential
update and require reconstruction of the stored complete Hamiltonian.  During
RT, rebuild Hartree/XC from `rho(t)` once at the beginning of each explicit
SALMON time step before calling the existing one-step length-gauge propagator.
Do not add an inner midpoint or predictor-corrector SCF.  Keep one stable
operator-union structure fingerprint and a separate value fingerprint for
each update; key sparse-exchange schedules only by the structure fingerprint
so a value update does not rebuild communication metadata.  Reject any value
outside the precomputed coupling envelope.  Exercise an initially zero edge
that becomes nonzero after a density perturbation.

**Step 4: Run GREEN and legacy RT regression**

Run:

```text
python3 tests/dg/run_rt_dg_hybrid_initialization_mpi.py
python3 tests/dg/run_rt_dg_hybrid_checkpoint_mpi.py
python3 tests/dg/run_rt_dg_hybrid_length_gauge_mpi.py
```

Expected: all PASS.  The full RT regression set is deferred to Task 11 and
final verification.

**Step 5: Commit**

Commit only the task files as `feat(rt): load exact hybrid DG ground state`.

### Task 10: Add zero-field RT stationarity acceptance

**Files:**
- Create: `src/rt/dg/rt_dg_hybrid_stationarity.f90`
- Create: `src/common/dg_hybrid_total_energy.f90`
- Modify: `src/common/CMakeLists.txt`
- Modify: `src/rt/CMakeLists.txt`
- Create: `tests/dg/test_rt_dg_hybrid_stationarity_mpi.f90`
- Create: `tests/dg/run_rt_dg_hybrid_stationarity_mpi.py`
- Modify: `src/rt/main_tddft.f90`

**Step 1: Write failing stationary and phase-rotation tests**

Propagate a known generalized eigenstate with zero external field.  Require
bounded drift in density, total energy, occupied `S`-projector, electron
number, and DG Hamiltonian residual.  Apply arbitrary occupied
phases and a degenerate-space rotation between samples; require the projector
test to pass while a deliberately changed occupied subspace fails.  Instrument
the Hartree/XC update callback and require exactly one call at the beginning of
each explicit RT step; a frozen-Hamiltonian propagator and an accidental inner
SCF loop must both fail the fixture.

Compare the hybrid energy evaluator against SALMON's existing DFT energy
decomposition for a no-interface fixture, then add an analytic two-fragment
fixture and require exactly one complete SIPG face-energy contribution.
Perturb Hartree and XC independently to prove that the evaluator applies the
existing double-counting corrections and is not `Tr(Gamma H)`.
Obtain `E_ion_nloc` from the existing nonlocal projector action while
discarding the ordinary-grid kinetic value.  Set the correctly normalized DG
kinetic field, call the unchanged `calc_Total_Energy_periodic`, and verify
every final `s_dft_energy` component against an independently summed
periodic fixture.  Require isolated boundaries to fail the hybrid scope gate.

**Step 2: Run RED**

Run: `python3 tests/dg/run_rt_dg_hybrid_stationarity_mpi.py`

Expected: compile failure because the stationarity evaluator is absent.

**Step 3: Implement the minimum evaluator and sampling hook**

Leave `src/common/total_energy.f90` and all protected callers unchanged.
Reuse the common projector/residual algebra.  Record initial invariants from
the checkpoint payload and compare them at configured RT samples.  Do not use
raw coefficient differences as an acceptance measure.  Reconstruct the
production-grid density and potentials, obtain the nonlocal component from
the existing projector action, set broken-volume plus correctly normalized
complete SIPG kinetic energy from `Gamma_occ`, and call the existing
`calc_Total_Energy_periodic` to compute and sum Hartree, XC, local-ionic,
and ion--ion terms.  Verify functional and pseudopotential provenance.  Do
not add a separate
RT symmetry-drift gate: the initial zero-field operator and occupied-space
symmetry were already accepted, and stationary projector/density checks cover
the zero-field case.  Driven-state symmetry is outside this plan.

**Step 4: Run GREEN**

Run the new runner at 1, 2, and 4 ranks and then:

```text
python3 tests/dg/run_rt_dg_hybrid_metric_solver_mpi.py
```

Expected: PASS.  Checkpoint and length-gauge regressions already ran in Task 9;
all protected RT runners run together in Task 11 and final verification.

**Step 5: Commit**

Commit only task files as `feat(rt): verify zero-field hybrid stationarity`.

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
The runner invokes Python/SALMON directly.

**Step 3: Run focused source/input checks GREEN**

Run the runner's parser-only mode and the route contract.  Expected: PASS for
synthetic complete receipts and rejection of incomplete receipts.

**Step 4: Run protected-route regressions**

Run the existing protected entry points and all hybrid RT runners available
after Tasks 9 and 10:

```text
python3 tests/dg/check_dg_hybrid_self_consistent_route.py
python3 tests/dg/check_dg_overlapping_wannier_route.py
python3 tests/dg/check_dg_fragment_symmetry_production.py
python3 tests/dg/test_replay_dg_wannier90_bundle.py
python3 tests/dg/run_rt_dg_hybrid_checkpoint_mpi.py
python3 tests/dg/run_rt_dg_hybrid_metric_solver_mpi.py
python3 tests/dg/run_rt_dg_hybrid_length_gauge_mpi.py
python3 tests/dg/run_rt_dg_hybrid_initialization_mpi.py
python3 tests/dg/run_rt_dg_hybrid_stationarity_mpi.py
```

Expected: PASS with no changed legacy branch behavior.

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
- separate acceptance of the self-consistent lambda-zero fixed point before
  any positive lambda stage;
- accepted adaptive stages ending at lambda one;
- final refreshed `R_H`, `R_rho`, `R_T`, `R_S`, electron count, symmetry, and
  real-space DG residual;
- complete checkpoint payload and matching GS/RT fingerprint;
- closure of the retained WF+PW space and covariance of the complete occupied
  ground-state projector under the actual group, accepting the identity-only
  group without bypassing the common checks and without individual-state
  symmetry;
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
