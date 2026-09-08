# WF+PW Variational DG Continuation Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Replace the invalid ordinary-`hpsi` final projection with a self-consistent WF+PW DG continuation built from broken-volume, complete nonlocal, and SIPG variational matrices.

**Architecture:** Freeze one retained WF+PW basis and assemble `S`, broken kinetic, exactly-once nonlocal, and complete SIPG rows before continuation. During SCF update only the density-dependent local-potential rows, mix density only, and advance one adaptive uniform lambda after every fixed point converges. Keep the numerical state and lambda decisions in small modules, but put the one physical SCF loop in the explicit production branch where SALMON's Hartree, XC, density, and basis objects already exist. Do not add a callback fixture or adapter layer.

**Tech Stack:** Fortran 2008, MPI, SALMON density and potential infrastructure, ScaLAPACK, SIPG weak form, standalone Python MPI runners.

---

Work only in
`/Users/otobetoshihito/SALMON-dev/SALMON2_RTDG/.worktrees/wpw-s-orthogonal-complement`.
Do not create another worktree. Preserve every existing uncommitted change.
For a pre-dirty file use `git add -p` and stage only task-specific hunks.
Before every commit run `git diff --cached --check` and inspect the complete
`git diff --cached`. Run Python test runners directly; do not add CTest
registration. Never use a time cutoff.

The protected ordinary GS/RT, DC+LCFO, Wannier90, and overlapping-Wannier
routes must remain unchanged. Every production call added below is guarded by
`yn_dg_hybrid_continuation_scf == 'y'`.

### Task 1: Assemble the fixed broken-volume kinetic and local matrices

**Files:**
- Create: `src/gs/dc/dg_hybrid_broken_volume.f90`
- Create: `tests/dg/test_dg_hybrid_broken_volume_mpi.f90`
- Create: `tests/dg/run_dg_hybrid_broken_volume_mpi.py`
- Modify: `src/gs/dc/CMakeLists.txt`

**Step 1: Write the failing analytic test**

Use two fragments with unequal interior point counts and unequal basis counts.
Provide complex basis values and analytic gradients on exactly-once owned
interior points. Require

```text
T_ij = 0.5 * sum_K sum_r w_r conj(grad phi_i) dot grad phi_j
V_ij =       sum_K sum_r w_r conj(phi_i) v_local(r) phi_j
```

for every owned row. Require Hermiticity, decomposition independence at 1, 2,
and 4 ranks, and zero kinetic and local cross-fragment blocks even when basis
tails are present in a neighboring buffer. Reject duplicate/missing interior
point ownership, duplicate/missing basis-row ownership, inconsistent fragment
IDs, nonfinite values, and gradients whose shapes do not match the basis.

**Step 2: Run RED**

Run: `python3 tests/dg/run_dg_hybrid_broken_volume_mpi.py`

Expected: compile failure because `dg_hybrid_broken_volume` is absent.

**Step 3: Implement the minimum row-owned assembler**

Add one routine:

```fortran
call assemble_dg_hybrid_broken_volume_rows(comm, global_basis_count, &
  row_ids, basis_fragment, interior_ids, interior_fragment, weights, &
  basis_values, basis_gradients, local_potential, kinetic_rows, &
  local_rows, diagnostics, ok, message)
```

Accumulate a matrix element only when both basis IDs belong to the owned
interior fragment. Use conjugation on the test function. Validate ownership
collectively before entering row reductions. Do not call `hpsi`, sample face
traces, or add a surface term.

**Step 4: Run GREEN**

Run the new runner and require PASS at 1, 2, and 4 ranks.

**Step 5: Commit**

Stage only the four task files and commit:

`feat(dg): assemble hybrid broken-volume rows`

### Task 2: Assemble crossing nonlocal projectors exactly once

**Files:**
- Modify: `src/gs/dc/dg_hybrid_broken_volume.f90`
- Modify: `tests/dg/test_dg_hybrid_broken_volume_mpi.f90`
- Modify: `tests/dg/run_dg_hybrid_broken_volume_mpi.py`
- Modify: `src/gs/dc/CMakeLists.txt` only if a new dependency is required

**Step 1: Extend the test and observe RED**

Add one projector contained in a fragment and one projector whose support
crosses the interface. Duplicate the crossing projector metadata on both
neighbor ranks, but assign one canonical projector owner. Require the dense
reference

```text
V_NL(i,j) = sum_a,mu,nu <phi_i|beta_a,mu> D_a(mu,nu)
                           <beta_a,nu|phi_j>
```

including nonzero cross-fragment blocks, independent of rank decomposition.
Require exactly-once accounting and reject zero, duplicate, or disagreeing
canonical ownership. Verify that changing lambda does not change `V_NL`.

**Step 2: Run RED**

Run the broken-volume runner. Expected: FAIL because the returned fixed-volume
operator lacks the crossing nonlocal contribution.

**Step 3: Add the minimum nonlocal input**

Extend the module with a separate routine that consumes canonical projector
IDs, owners, strengths/coupling blocks, and distributed basis-projector
overlaps. Reuse existing nonlocal row algebra where its contract matches.
Do not truncate a projector at a fragment boundary and do not place it in the
SIPG operator.

**Step 4: Run GREEN and regression**

Run:

```text
python3 tests/dg/run_dg_hybrid_broken_volume_mpi.py
python3 tests/dg/run_dg_overlapping_wannier_nonlocal_mpi.py
```

Expected: PASS.

**Step 5: Commit**

Commit only task hunks as:

`feat(dg): retain complete nonlocal volume coupling`

### Task 3: Freeze and compose the variational DG payload

**Files:**
- Create: `src/gs/dc/dg_hybrid_variational_payload.f90`
- Create: `tests/dg/test_dg_hybrid_variational_payload_mpi.f90`
- Create: `tests/dg/run_dg_hybrid_variational_payload_mpi.py`
- Modify: `src/gs/dc/CMakeLists.txt`

**Step 1: Write the failing composition test**

Construct row-owned `S`, broken kinetic, nonlocal, local, and SIPG matrices.
Require

```text
H(lambda,rho) = T_broken + V_NL + V_local(rho) + lambda * H_SIPG
```

with one scalar lambda. At lambda zero require no SIPG contribution; at lambda
one require the complete SIPG contribution. Require volume kinetic cross
blocks to remain zero and the kinetic cross blocks of `H` to equal the SIPG
cross blocks. Require fixed-matrix fingerprints and values to remain unchanged
when local rows are replaced. Reject face-local lambda arrays, inconsistent
row layouts, non-Hermitian components, and a fixed fingerprint mutation.

**Step 2: Run RED**

Run: `python3 tests/dg/run_dg_hybrid_variational_payload_mpi.py`

Expected: compile failure because the payload module is absent.

**Step 3: Implement three small types and one composition routine**

Add immutable `s_dg_hybrid_fixed_payload`, mutable
`s_dg_hybrid_variational_iterate`, and accepted snapshot types. The fixed type
owns `S`, `T_broken`, `V_NL`, and `H_SIPG`; the iterate owns `V_local` and
complete `H`. Compose rows without communication after collective contract
validation. Do not add callbacks or a runtime catalog.

**Step 4: Run GREEN**

Run the new runner at 1, 2, and 4 ranks, followed by the SIPG and face runners.

**Step 5: Commit**

Commit only task files as:

`feat(dg): freeze variational continuation payload`

### Task 4: Implement the concrete DC-style density-to-potential kernel

**Files:**
- Create: `src/gs/dc/dg_hybrid_variational_potential.f90`
- Create: `tests/dg/test_dg_hybrid_variational_potential_mpi.f90`
- Create: `tests/dg/run_dg_hybrid_variational_potential_mpi.py`
- Modify: `src/gs/dc/CMakeLists.txt`
- Modify: `src/gs/dc/dcdft.f90` only to expose an existing DC redistribution
  helper if it cannot be called directly

**Step 1: Write the failing decomposed-potential test**

Use irregular fragments with unequal core sizes. Require exactly-once
fragment core density assembly into the total grid, a deterministic reference
Hartree operation on that total density, redistribution of the Hartree field
to fragment buffers, and independent fragment-local XC evaluation. Verify
that the combined local field equals `V_H + V_xc + V_ion_local` at every
fragment point. Include a semilocal fixture whose gradient stencil uses the
fragment halo. Reject duplicate or missing core ownership and incomplete
semilocal halos.

**Step 2: Run RED**

Run: `python3 tests/dg/run_dg_hybrid_variational_potential_mpi.py`

Expected: compile failure because the concrete potential kernel is absent.

**Step 3: Implement the minimum DC-style kernel**

Reuse the mappings and communication order of `calc_rho_total_dcdft` and
`calc_vlocal_fragment_dcdft`. Assemble the total density, call the existing
total-system Hartree FFT, return Hartree values to fragments, evaluate the
supported local/semilocal XC on fragment buffers, and combine the fixed local
ionic field. Expose a concrete routine over SALMON grid, Poisson, XC, and DC
state types; do not use procedure arguments or a callback table.

**Step 4: Run GREEN and DC regressions**

Run the new runner at 1, 2, and 4 ranks, then the divided-DC control and
protected DC route checks.

**Step 5: Commit**

Commit only task files as:

`feat(dg): update hybrid potential from divided density`

### Task 5: Remove the callback SCF layer and test the concrete loop contract

**Files:**
- Delete: `src/gs/dc/dg_hybrid_continuation_scf.f90`
- Delete: `tests/dg/test_dg_hybrid_continuation_scf_mpi.f90`
- Modify: `tests/dg/run_dg_hybrid_continuation_scf_mpi.py`
- Modify: `src/gs/dc/CMakeLists.txt`

**Step 1: Replace the fixture-only RED contract**

Inspect the contained production loop and require the first local build to see
the exact supplied DC density. Require gradual updates

```text
rho_next = rho_in + alpha_rho * (rho_out - rho_in)
```

with neither raw coefficients nor interface traces mixed. Require operation
order `local -> compose -> solve -> Gamma/Q -> density/trace -> residuals ->
density mix`. Require lambda zero to converge before positive lambda, uniform
lambda, a forced complete rollback, gauge-invariant projector tracking, and
one fully refreshed lambda-one final state.

Require operation order `Hartree -> fragment XC -> local projection -> compose
-> solve -> Gamma/Q -> density/trace -> residuals -> density mix`. Add a source
assertion forbidding abstract callback interfaces, `class(*)`, procedure-pointer
tables, and a fixture-only public entry point.

**Step 2: Run RED**

Run: `python3 tests/dg/run_dg_hybrid_continuation_scf_mpi.py`

Expected: FAIL because the production branch does not yet contain the loop and
the callback fixture still exists.

**Step 3: Remove the obsolete abstraction**

Delete the callback-only SCF module and its synthetic callback test. Keep the
existing controller, residual, acceptance, eigensystem, density, potential,
and payload modules as independently tested numerical kernels. The following
task supplies their only physical orchestration loop.

**Step 4: Run GREEN at every focused decomposition**

Run the source-contract runner, then the residual, controller, acceptance,
potential, density, and generalized-eigensystem runners.

**Step 5: Commit**

Commit only the three task files as:

`refactor(dg): remove callback continuation fixture`

### Task 6: Implement the production WF+PW continuation loop

**Files:**
- Modify: `src/gs/main_dft.f90`
- Modify: `tests/dg/check_dg_hybrid_continuation_route.py`
- Modify: `tests/dg/check_dg_hybrid_divided_dc_controls.py` only if the
  explicit flag contract changes

**Step 1: Strengthen the route RED test**

Require the explicit continuation branch to:

- seed from `dc%rho_tot`;
- materialize the fixed retained basis and interior gradients;
- assemble broken kinetic, complete nonlocal, metric, and SIPG rows;
- contain the single explicit lambda/SCF loop;
- exclude ordinary `hpsi` projection, the one-shot generalized solve, and the
  occupied-only checkpoint from the continuation driver.

Require those exclusions only inside the explicit continuation branch so the
protected legacy path remains byte-for-byte unchanged.

**Step 2: Run RED**

Run the route test. Expected: FAIL because the contained production driver and
broken-volume connection are absent.

**Step 3: Add the contained production driver**

Materialize global basis values and gradients only on each rank's owned
fragment interior. Assemble all fixed matrices before SCF. Seed the loop with
the exact converged DC density. During every inner iteration assemble that
density for the existing total-system Hartree FFT, redistribute only the
Hartree field, evaluate local/semi-local XC on fragment buffers, and project
the combined local field. Compose, solve, rebuild Gamma/Q, density and traces,
evaluate residuals, then mix density only. Use the compact controller for
uniform adaptive lambda and complete rollback. In the explicit
continuation branch, skip `solve_dg_hybrid_generalized_once_and_publish` and
`write_rt_dg_hybrid_occupied_checkpoint`. Do not change their legacy callers.

**Step 4: Run GREEN and protected regressions**

Run:

```text
python3 tests/dg/check_dg_hybrid_continuation_route.py
python3 tests/dg/check_dg_hybrid_divided_dc_controls.py
python3 tests/dg/run_dg_hybrid_broken_volume_mpi.py
python3 tests/dg/run_dg_hybrid_variational_payload_mpi.py
python3 tests/dg/run_dg_hybrid_continuation_scf_mpi.py
python3 tests/dg/run_dg_hybrid_production_face_traces_mpi.py
python3 tests/dg/run_dg_hybrid_sipg_operator_mpi.py
python3 tests/dg/check_dg_hybrid_self_consistent_route.py
python3 tests/dg/check_dg_overlapping_wannier_route.py
```

Expected: PASS, followed by `cmake --build build-hybrid-commit -j2` PASS.

**Step 5: Commit**

Use `git add -p` for `src/gs/main_dft.f90` and any other pre-dirty file. Commit
only task hunks as:

`feat(dg): connect variational production continuation`

### Task 7: Resume the checkpoint, RT, and Si64 acceptance tasks

**Files:**
- Modify the Task 8--12 files already listed in
  `docs/plans/2026-08-27-wpw-dg-continuation-ground-state.md`.

**Step 1: Reconcile the remaining plan**

Before Task 8, update references from a generic volume operator to the exact
payload components `T_broken`, `V_NL`, `V_local`, and `H_SIPG`. The checkpoint
must store the complete lambda-one Hamiltonian and its component fingerprints.

**Step 2: Execute Tasks 8--11 with TDD**

Keep their GS-to-RT identical-payload and zero-field stationarity acceptance
criteria. Do not add an RT driven-symmetry gate.

**Step 3: Run fresh final verification**

Invoke `verification-before-completion`, run every focused and protected
runner freshly, then run Si64 with eight MPI ranks, `OMP_NUM_THREADS=1`, and no
timeout. Preserve every log.

**Step 4: Review and branch completion**

Invoke `requesting-code-review`, resolve important findings with TDD, rerun
affected verification, and only then invoke `finishing-a-development-branch`.
