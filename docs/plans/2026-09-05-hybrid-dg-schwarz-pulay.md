# Hybrid DG Schwarz + Pulay Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Include the complete nearest-neighbor DG coupling from the first divided-Hybrid density epoch, retain one MPI rank per fragment and bounded CG work, use persistent Hybrid-specific Pulay mixing without invalidating reusable ordinary-DC seeds, and perform exactly one terminal LCFO solve.

**Architecture:** Each rank owns the coefficient rows for one fragment but shares a dynamically sized state-column inventory.  An immutable schedule derived from the production basis directory and SIPG face catalog exchanges only adjacent-fragment coefficient blocks.  A transactional block-Jacobi/Schwarz update applies the full DG H/S operator, globally S-orthonormalizes the common columns, assigns 300 K occupations with one chemical potential, and publishes density only after collective validation.  SALMON's existing mixer implementations are selected through a post-seed Hybrid control, and the already-existing terminal LCFO path remains the sole full diagonalization.

**Tech Stack:** Fortran 2008, MPI, BLAS/LAPACK, SALMON density mixing, CMake, Python contract runners, bounds/FPE-checked MPI tests.

---

## Execution constraints

- Work only in `/Users/otobetoshihito/SALMON-dev/SALMON2_RTDG/.worktrees/wpw-s-orthogonal-complement`.
- Do not create another worktree and do not clean, overwrite, or discard existing dirty files or verification logs.
- Keep production at one MPI rank per fragment.  A different MPI size or rank-fragment map must reject DC-seed reuse.
- Never hard-code the Si64 state count (including 384); derive the common column count from electrons, occupations, degeneracy, guard policy, and available candidates.
- Keep `dg_hybrid_fragment_cg_steps` as a hard per-density-epoch cap (default 3), with early exit allowed.
- Do not add a fallback to self-block CG, complete coefficient all-gather, repeated complete diagonalization, or simple mixing.
- Commit only files intentionally changed by each task.  Preserve every unrelated dirty change.

### Task 1: Add the Hybrid-only mixing selector without changing the DC-seed fingerprint

**Files:**

- Modify: `src/io/salmon_global.f90`
- Modify: `src/io/inputoutput.f90`
- Modify: `src/gs/main_dft.f90`
- Modify: `tests/dg/check_dg_hybrid_divided_dc_controls.py`
- Modify: `tests/dg/check_dg_dc_seed_route.py`
- Modify: `tests/dg/check_dg_hybrid_localization_first_inputs.py`
- Modify: `tests/dg/data/si64_overlapping_wannier_rt/input_hybrid_divided_lcfo.in`

**Step 1: Write failing input and fingerprint contracts**

Extend the Python checks to require:

- `character(16) :: dg_hybrid_divided_mixing` in global input state;
- default `pulay`, lowercase normalization, broadcast, variables-log echo, and validation of only `inherit`, `simple`, `pulay`, `broyden`;
- the divided callback resolves `inherit` to `method_mixing` and otherwise uses the Hybrid selector;
- the ordinary DC convergence fingerprint continues to hash `method_mixing` and `mixing%mixrate`, but does not hash `dg_hybrid_divided_mixing`;
- the Si64 divided fixture requests `dg_hybrid_divided_mixing='pulay'` while keeping its ordinary seed-compatible `method_mixing` unchanged.

**Step 2: Run the contracts and confirm failure**

Run:

```bash
python3 tests/dg/check_dg_hybrid_divided_dc_controls.py
python3 tests/dg/check_dg_dc_seed_route.py
python3 tests/dg/check_dg_hybrid_localization_first_inputs.py
```

Expected: failure because the new input is absent.

**Step 3: Implement the input and selector**

Add the namelist field, default, normalization, broadcast, echo, and collective validation.  Resolve the selected method once on entering divided Hybrid SCF and log:

```text
[DG-HYBRID-MIXING] method=pulay mixrate=... seed_fingerprint_unchanged=T
```

Keep the ordinary DC convergence-fingerprint code byte-for-byte independent of the new selector.

**Step 4: Re-run the contracts**

Expected: all three pass.

**Step 5: Commit the task**

```bash
git add src/io/salmon_global.f90 src/io/inputoutput.f90 src/gs/main_dft.f90 tests/dg/check_dg_hybrid_divided_dc_controls.py tests/dg/check_dg_dc_seed_route.py tests/dg/check_dg_hybrid_localization_first_inputs.py tests/dg/data/si64_overlapping_wannier_rt/input_hybrid_divided_lcfo.in
git commit -m "feat(dg): select divided Hybrid density mixing"
```

### Task 2: Publish a dynamic common state-column inventory

**Files:**

- Create: `src/gs/dc/dg_hybrid_schwarz_state.f90`
- Modify: `src/gs/dc/CMakeLists.txt`
- Create: `tests/dg/test_dg_hybrid_schwarz_state_mpi.f90`
- Create: `tests/dg/run_dg_hybrid_schwarz_state_mpi.py`

**Step 1: Write the failing MPI test**

Cover 2, 4, and 8 ranks with unequal local basis sizes.  Require one rank per fragment; common ordered column IDs; initial count derived from global electron target, 300 K occupation tail, degeneracy completion, and guard; transactional extension on every rank; stable nonzero fingerprint; and collective rejection of duplicate IDs, stale generation, insufficient candidates, rank-fragment mismatch, or a rank-local publication failure.

**Step 2: Confirm the test fails to compile**

```bash
python3 tests/dg/run_dg_hybrid_schwarz_state_mpi.py
```

Expected: missing module/API.

**Step 3: Implement the state contract**

Define inventory and accepted-state types containing generation, common IDs, local row coefficients `C_f(nb_f,ntrial)`, owner/mapping fingerprints, and occupation provenance.  Build and extend them collectively without publishing partially valid allocations.  Do not pad with material-specific constants or bind columns by fragment-local eigenstate ordinal.

**Step 4: Run the MPI test**

Expected: pass on 2, 4, and 8 ranks under bounds/FPE checks.

**Step 5: Commit**

```bash
git add src/gs/dc/dg_hybrid_schwarz_state.f90 src/gs/dc/CMakeLists.txt tests/dg/test_dg_hybrid_schwarz_state_mpi.f90 tests/dg/run_dg_hybrid_schwarz_state_mpi.py
git commit -m "feat(dg): publish common Schwarz state inventory"
```

### Task 3: Build an immutable nearest-neighbor DG schedule

**Files:**

- Create: `src/gs/dc/dg_hybrid_schwarz_operator.f90`
- Modify: `src/gs/dc/CMakeLists.txt`
- Create: `tests/dg/test_dg_hybrid_schwarz_operator_mpi.f90`
- Create: `tests/dg/run_dg_hybrid_schwarz_operator_mpi.py`
- Modify: `tests/dg/run_dg_hybrid_production_face_traces_mpi.py`

**Step 1: Write schedule failure tests**

Construct small line/ring fixtures from the same row-owner directory and interface rows used by production.  Assert exact neighbor lists, remote basis-column requests, deterministic ordering, and fingerprints.  Add collective failures for missing reciprocal faces, duplicate requests, a non-neighbor destination, stale basis generation, stale face fingerprint, and changed rank-fragment mapping.

**Step 2: Confirm failure**

```bash
python3 tests/dg/run_dg_hybrid_schwarz_operator_mpi.py
```

Expected: missing schedule builder.

**Step 3: Implement schedule construction**

Create the schedule once per basis generation from `basis_owner`, `basis_fragment`, `basis_local_slot`, the fixed-payload directory fingerprint, and nonzero cross-fragment SIPG interface entries.  Store only adjacent peers and the precise send/receive row slots.  Authenticate all immutable provenance collectively.

**Step 4: Verify schedule and production face tests**

```bash
python3 tests/dg/run_dg_hybrid_schwarz_operator_mpi.py
python3 tests/dg/run_dg_hybrid_production_face_traces_mpi.py
```

Expected: both pass on 2, 4, and 8 ranks.

**Step 5: Commit**

```bash
git add src/gs/dc/dg_hybrid_schwarz_operator.f90 src/gs/dc/CMakeLists.txt tests/dg/test_dg_hybrid_schwarz_operator_mpi.f90 tests/dg/run_dg_hybrid_schwarz_operator_mpi.py tests/dg/run_dg_hybrid_production_face_traces_mpi.py
git commit -m "feat(dg): build immutable Schwarz neighbor schedule"
```

### Task 4: Apply full DG H and S through neighbor-only exchange

**Files:**

- Modify: `src/gs/dc/dg_hybrid_schwarz_operator.f90`
- Modify: `tests/dg/test_dg_hybrid_schwarz_operator_mpi.f90`
- Modify: `tests/dg/run_dg_hybrid_schwarz_operator_mpi.py`
- Modify: `tests/dg/check_dg_hybrid_divided_dc_controls.py`

**Step 1: Add a failing numerical oracle**

For unequal fragment basis sizes and multiple common columns, compare distributed H/S actions with explicit dense application assembled from `metric_rows`, `kinetic_rows`, `nonlocal_rows`, `interface_rows`, and iteration-local potential rows.  Run 2, 4, and 8 ranks and include complex coefficients.  Instrument peer traffic and assert that no non-neighbor and no complete coefficient all-gather occurs.

**Step 2: Confirm numerical failure**

The existing self-block action must differ whenever an off-diagonal face block is nonzero.

**Step 3: Implement the distributed application**

At each call, freeze the accepted step-start coefficients; apply local interior/nonlocal/potential/SIPG-self rows; exchange only requested neighbor coefficient blocks with deterministic tags; add off-diagonal face terms; and return only locally owned output rows.  Apply S through the same ownership/provenance contract.  Collectively reject missing, duplicate, stale, nonfinite, or out-of-schedule messages.

**Step 4: Run operator and static contracts**

```bash
python3 tests/dg/run_dg_hybrid_schwarz_operator_mpi.py
python3 tests/dg/check_dg_hybrid_divided_dc_controls.py
```

Expected: dense-oracle equality within the test tolerance on 2, 4, and 8 ranks; communication contract passes.

**Step 5: Commit**

```bash
git add src/gs/dc/dg_hybrid_schwarz_operator.f90 tests/dg/test_dg_hybrid_schwarz_operator_mpi.f90 tests/dg/run_dg_hybrid_schwarz_operator_mpi.py tests/dg/check_dg_hybrid_divided_dc_controls.py
git commit -m "feat(dg): apply full operator by neighbor exchange"
```

### Task 5: Add transactional bounded Schwarz-CG updates

**Files:**

- Create: `src/gs/dc/dg_hybrid_schwarz_solver.f90`
- Modify: `src/gs/dc/CMakeLists.txt`
- Create: `tests/dg/test_dg_hybrid_schwarz_solver_mpi.f90`
- Create: `tests/dg/run_dg_hybrid_schwarz_solver_mpi.py`

**Step 1: Write failing solver tests**

Use a small positive-definite generalized problem distributed by rows.  Require one to three preconditioned CG steps per density epoch, early exit, neighbor coefficients fixed at the start of each block-Jacobi step, global `C^dagger S C=I`, monotonically accepted Rayleigh/residual criteria, and whole-step rollback after an injected rank-local H, S, preconditioner, or orthogonalization failure.

**Step 2: Confirm missing solver failure**

```bash
python3 tests/dg/run_dg_hybrid_schwarz_solver_mpi.py
```

**Step 3: Implement the solver**

Reuse the existing fragment preconditioner and bounded-CG acceptance ideas, but operate on distributed coefficient rows and use communicator-wide reductions only for small state-space matrices, norms, and energies.  Stage every new coefficient block until collective validation.  Record attempted/accepted step counts and rollback reason; never invoke a complete eigensolver.

**Step 4: Run the test**

Expected: pass on 2, 4, and 8 ranks with bounds/FPE checks.

**Step 5: Commit**

```bash
git add src/gs/dc/dg_hybrid_schwarz_solver.f90 src/gs/dc/CMakeLists.txt tests/dg/test_dg_hybrid_schwarz_solver_mpi.f90 tests/dg/run_dg_hybrid_schwarz_solver_mpi.py
git commit -m "feat(dg): add bounded Schwarz CG updates"
```

### Task 6: Use one 300 K occupation epoch and dynamic common extension

**Files:**

- Modify: `src/gs/dc/dg_hybrid_schwarz_state.f90`
- Modify: `src/gs/dc/dg_hybrid_schwarz_solver.f90`
- Modify: `tests/dg/test_dg_hybrid_schwarz_state_mpi.f90`
- Modify: `tests/dg/test_dg_hybrid_schwarz_solver_mpi.f90`

**Step 1: Add failing occupation tests**

Require global Rayleigh energies, one 300 K chemical potential, occupations summing to the original global electron target, degeneracy-safe inclusion at the boundary, and simultaneous common-column extension when the upper occupation tail is unresolved.  Confirm that early SCF target drift is reported but not by itself fatal, matching ordinary DC convergence behavior.  Exhausted candidate capacity must fail collectively.

**Step 2: Confirm failure**

Run both new MPI runners; expected failure because occupations are still fragment-local.

**Step 3: Implement global occupation/extension**

Compute global Rayleigh energies from distributed H/S products, call the established finite-temperature occupation policy at 300 K, and publish one chemical potential and occupation vector.  Extend every rank's common inventory transactionally from its ordered DC/PW candidate catalog when required.  Reset the accepted solver search state after extension while preserving the last accepted density.

**Step 4: Re-run tests**

Expected: both state and solver runners pass on 2, 4, and 8 ranks.

**Step 5: Commit**

```bash
git add src/gs/dc/dg_hybrid_schwarz_state.f90 src/gs/dc/dg_hybrid_schwarz_solver.f90 tests/dg/test_dg_hybrid_schwarz_state_mpi.f90 tests/dg/test_dg_hybrid_schwarz_solver_mpi.f90
git commit -m "feat(dg): assign global thermal Schwarz occupations"
```

### Task 7: Make Pulay history persistent and explicitly resettable

**Files:**

- Create: `src/gs/dc/dg_hybrid_divided_mixing.f90`
- Modify: `src/gs/dc/CMakeLists.txt`
- Modify: `src/gs/main_dft.f90`
- Create: `tests/dg/test_dg_hybrid_divided_mixing_mpi.f90`
- Create: `tests/dg/run_dg_hybrid_divided_mixing_mpi.py`
- Modify: `tests/dg/check_dg_hybrid_divided_dc_controls.py`

**Step 1: Write failing history tests**

Require the selected method and `mixrate` to persist across density epochs; Pulay history length to grow on an unchanged basis/inventory; no reset for an ordinary accepted epoch; and a diagnosed reset only for basis generation change, common-inventory extension, or collective rollback.  Confirm `inherit`, simple, Pulay, and Broyden dispatch and reject unsupported values without fallback.

**Step 2: Confirm failure**

```bash
python3 tests/dg/run_dg_hybrid_divided_mixing_mpi.py
python3 tests/dg/check_dg_hybrid_divided_dc_controls.py
```

**Step 3: Implement the adapter**

Move Hybrid mixing lifecycle and reset bookkeeping out of the large main callback into a small module.  Reuse `copy_density`, `simple_mixing`, `pulay`, and `wrapper_broyden`; keep the same `s_mixing` instance across epochs.  Log method, mix rate, history length, and any reset reason.  Expose callbacks to gather/scatter the physical core density without changing ordinary DC mixing configuration or seed fingerprint.

**Step 4: Re-run tests**

Expected: history and static contracts pass.

**Step 5: Commit**

```bash
git add src/gs/dc/dg_hybrid_divided_mixing.f90 src/gs/dc/CMakeLists.txt src/gs/main_dft.f90 tests/dg/test_dg_hybrid_divided_mixing_mpi.f90 tests/dg/run_dg_hybrid_divided_mixing_mpi.py tests/dg/check_dg_hybrid_divided_dc_controls.py
git commit -m "feat(dg): persist divided Hybrid Pulay history"
```

### Task 8: Wire Schwarz DG into the production divided SCF

**Files:**

- Modify: `src/gs/main_dft.f90`
- Modify: `src/gs/dc/dg_hybrid_divided_scf.f90`
- Modify: `src/gs/dc/dg_hybrid_fragment_thermal.f90`
- Modify: `tests/dg/test_dg_hybrid_divided_scf_mpi.f90`
- Modify: `tests/dg/run_dg_hybrid_divided_scf_mpi.py`
- Modify: `tests/dg/check_dg_hybrid_divided_dc_controls.py`

**Step 1: Write failing route and integration tests**

Require production setup to create the common inventory and schedule once per basis generation; every divided-SCF H/S callback to use the Schwarz operator; `fragment_count==MPI size`; CG steps never to exceed the user cap; logs to report neighbor exchanges, accepted CG steps, extensions, density convergence, electron defect, and mixing history.  Assert the divided route no longer calls `extract_dg_hybrid_fragment_self_block` for its iterative solve.

**Step 2: Confirm current route fails**

```bash
python3 tests/dg/check_dg_hybrid_divided_dc_controls.py
python3 tests/dg/run_dg_hybrid_divided_scf_mpi.py
```

Expected: self-block-only route detected.

**Step 3: Wire the new components**

Initialize distributed state from admitted WF+PW candidates, construct the immutable neighbor schedule after freezing the payload, and replace `solve_dg_hybrid_bounded_fragments`/local thermal callbacks with the bounded Schwarz solver.  Reconstruct core density from the globally occupied common columns, then invoke persistent Hybrid mixing.  Preserve the accepted density and coefficients across a failed collective trial.

**Step 4: Run focused integration tests**

```bash
python3 tests/dg/run_dg_hybrid_schwarz_state_mpi.py
python3 tests/dg/run_dg_hybrid_schwarz_operator_mpi.py
python3 tests/dg/run_dg_hybrid_schwarz_solver_mpi.py
python3 tests/dg/run_dg_hybrid_divided_mixing_mpi.py
python3 tests/dg/run_dg_hybrid_divided_scf_mpi.py
python3 tests/dg/check_dg_hybrid_divided_dc_controls.py
```

Expected: all pass on their configured 2/4/8-rank matrices.

**Step 5: Commit**

```bash
git add src/gs/main_dft.f90 src/gs/dc/dg_hybrid_divided_scf.f90 src/gs/dc/dg_hybrid_fragment_thermal.f90 tests/dg/test_dg_hybrid_divided_scf_mpi.f90 tests/dg/run_dg_hybrid_divided_scf_mpi.py tests/dg/check_dg_hybrid_divided_dc_controls.py
git commit -m "feat(dg): run divided SCF with Schwarz coupling"
```

### Task 9: Preserve the one-terminal-LCFO contract

**Files:**

- Modify: `tests/dg/check_dg_hybrid_divided_lcfo_route.py`
- Modify: `tests/dg/check_dg_hybrid_divided_dc_controls.py`
- Modify: `tests/dg/run_dg_hybrid_si64_divided_lcfo.py`
- Modify: `src/gs/main_dft.f90` only if the tests expose a route defect

**Step 1: Strengthen the failing/static contract**

Assert that no complete generalized eigensolve occurs inside the divided density loop; after convergence the total potential is refreshed once, final complete H/S rows are assembled once, and `solve_dg_hybrid_generalized_once_and_publish` is called exactly once.  Assert no density mixing/update follows LCFO unless the existing explicit short-global-refinement mode is selected.

**Step 2: Run the route checks**

```bash
python3 tests/dg/check_dg_hybrid_divided_lcfo_route.py
python3 tests/dg/check_dg_hybrid_divided_dc_controls.py
```

Expected: pass if wiring preserved the current route; otherwise fail at the exact misplaced call.

**Step 3: Make the minimum correction if required**

Keep the LCFO result authoritative for symmetry, occupations, checkpoint state, and real-time initialization.  Do not add a post-LCFO density update to the default route.

**Step 4: Re-run route checks**

Expected: both pass.

**Step 5: Commit**

```bash
git add tests/dg/check_dg_hybrid_divided_lcfo_route.py tests/dg/check_dg_hybrid_divided_dc_controls.py tests/dg/run_dg_hybrid_si64_divided_lcfo.py src/gs/main_dft.f90
git commit -m "test(dg): enforce one terminal Hybrid LCFO solve"
```

Before committing, omit `src/gs/main_dft.f90` from `git add` if it required no Task 9 change.

### Task 10: Verify Si64 simple versus Pulay with one reused seed

**Files:**

- Modify: `tests/dg/run_dg_hybrid_si64_divided_lcfo.py`
- Create only new timestamped output directories under `/tmp`; never replace earlier logs

**Step 1: Build the release executable**

```bash
cmake --build build-hybrid-release -j2
```

Expected: successful build.

**Step 2: Validate the reusable seed before both runs**

Use `/tmp/si64-task8-bounded-smoke-20260905/dc-seed` only if its publication ID is `7047888166118007469`, MPI size is 8, and mapping fingerprint is `254086644876463474`.  Have the runner fail before computation if either Hybrid variant changes these values or executes ordinary DC SCF.

**Step 3: Run a simple-mixing baseline**

Create a new timestamped `/tmp/si64-task8-schwarz-simple-*` directory, select `dg_hybrid_divided_mixing='simple'`, keep the seed-side `method_mixing` and mix rate unchanged, and run eight ranks.  Preserve the full log even on failure.

**Step 4: Run Pulay from the identical seed**

Create a different timestamped `/tmp/si64-task8-schwarz-pulay-*` directory, change only the Hybrid selector to `pulay`, and run eight ranks.  Preserve the full log.

**Step 5: Compare evidence**

Require both logs to show the same seed publication/mapping fingerprints, full neighbor DG from epoch 1, CG counts in `[1,3]`, no repeated full diagonalization, no NaN/rollback/fallback, and exactly one terminal LCFO if converged.  Compare density-convergence histories and electron defects rather than judging from norm preservation alone.  If a run does not converge, report its measured trend and retain the directory; do not loosen tolerances silently.

**Step 6: Commit runner changes**

```bash
git add tests/dg/run_dg_hybrid_si64_divided_lcfo.py
git commit -m "test(dg): compare Schwarz simple and Pulay mixing"
```

### Task 11: Remove the obsolete self-block production route and run final verification

**Files:**

- Modify: `src/gs/main_dft.f90`
- Modify: `src/gs/dc/dg_hybrid_divided_operator.f90`
- Modify: `tests/dg/test_dg_hybrid_divided_operator_mpi.f90`
- Modify: `tests/dg/run_dg_hybrid_divided_operator_mpi.py`
- Modify: relevant route checks that explicitly inventory production entry points

**Step 1: Prove the old route is unused**

```bash
rg -n "extract_dg_hybrid_fragment_self_block|solve_dg_hybrid_bounded_fragments" src tests/dg
```

Expected: references only in obsolete implementation/tests, never in the production divided route.  If another supported diagnostic genuinely uses the extractor, retain the extractor but delete only the obsolete divided-SCF callback.

**Step 2: Remove only dead production code and update focused tests**

Delete the superseded self-block divided-SCF callback and unused saved state.  Do not remove the full-cell diagnostic oracle or any route still exercised by a supported mode.

**Step 3: Run focused and regression checks**

```bash
python3 tests/dg/run_dg_hybrid_divided_operator_mpi.py
python3 tests/dg/run_dg_hybrid_production_face_traces_mpi.py
python3 tests/dg/run_dg_hybrid_schwarz_state_mpi.py
python3 tests/dg/run_dg_hybrid_schwarz_operator_mpi.py
python3 tests/dg/run_dg_hybrid_schwarz_solver_mpi.py
python3 tests/dg/run_dg_hybrid_divided_mixing_mpi.py
python3 tests/dg/run_dg_hybrid_divided_scf_mpi.py
python3 tests/dg/check_dg_hybrid_divided_lcfo_route.py
python3 tests/dg/check_dg_hybrid_divided_dc_controls.py
python3 tests/dg/check_dg_dc_seed_route.py
python3 tests/dg/check_dg_hybrid_localization_first_inputs.py
cmake --build build-hybrid-release -j2
git diff --check
```

Expected: every test/build/check passes.

**Step 4: Review the complete change**

Inspect `git diff --stat`, the scoped diffs, and all new logs.  Verify explicitly that there is no fixed state count, no coefficient all-gather, no rank-sharing of one fragment, no repeated LCFO solve, no post-LCFO density update by default, no seed-fingerprint dependency on the Hybrid mixer, and no fallback path.

**Step 5: Commit cleanup**

```bash
git add src/gs/main_dft.f90 src/gs/dc/dg_hybrid_divided_operator.f90 tests/dg/test_dg_hybrid_divided_operator_mpi.f90 tests/dg/run_dg_hybrid_divided_operator_mpi.py tests/dg/check_dg_hybrid_divided_lcfo_route.py tests/dg/check_dg_hybrid_divided_dc_controls.py
git commit -m "refactor(dg): remove obsolete divided self-block route"
```

Only add files that actually changed in Task 11.
