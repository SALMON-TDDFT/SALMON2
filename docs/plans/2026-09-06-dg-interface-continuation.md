# Fixed-Density DG Interface Continuation Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Replace the nonconvergent production density-mixing loop with a diagnostic fixed-DC-density continuation that scales the complete SIPG interface operator from zero to one before one terminal LCFO solve.

**Architecture:** Introduce a small collective continuation-state module that owns the lambda schedule and rollback rules.  Pass lambda explicitly into the Schwarz Hamiltonian action and apply the same scaling to the local preconditioner; the production driver holds the restored DC density and potential fixed throughout the schedule.  At lambda one, assemble the unscaled full operator and diagonalize exactly once without a post-LCFO density update.

**Tech Stack:** Fortran 2008, MPI, SALMON DG/Schwarz modules, Python source-contract tests, CMake/CTest.

**Workspace constraint:** Work only in `/Users/otobetoshihito/SALMON-dev/SALMON2_RTDG/.worktrees/wpw-s-orthogonal-complement`. Do not create another worktree. Preserve all existing dirty changes and verification logs; stage only task-owned hunks.

---

### Task 1: Collective fixed-step continuation state

**Files:**
- Create: `src/gs/dc/dg_hybrid_interface_continuation.f90`
- Modify: `src/gs/dc/CMakeLists.txt`
- Create: `tests/dg/test_dg_hybrid_interface_continuation_mpi.f90`
- Create: `tests/dg/run_dg_hybrid_interface_continuation_mpi.py`

**Step 1: Write the failing MPI test**

Cover rates `0.2` and `0.3`.  Assert that the initial point is exactly zero, accepted points are monotone, the terminal point is exactly one, and all ranks have identical state fingerprints.  Add rejection cases for rate `<=0`, rate `>1`, rank-disagreeing controls, and basis/mapping fingerprint disagreement.  Assert that a rejected step leaves lambda and the accepted-step counter unchanged.

**Step 2: Run the test to verify RED**

Run:

```bash
python3 tests/dg/run_dg_hybrid_interface_continuation_mpi.py
```

Expected: FAIL because the continuation module is absent.

**Step 3: Implement the minimal state machine**

Define `s_dg_hybrid_interface_continuation` with `valid`, `basis_generation`, `mapping_fingerprint`, `rate`, `lambda`, `step_index`, `accepted_steps`, `finished`, and `fingerprint`.  Export:

```fortran
initialize_dg_hybrid_interface_continuation(comm,basis_generation,mapping_fingerprint,rate,state,ok,message)
accept_dg_hybrid_interface_point(comm,basis_generation,mapping_fingerprint,local_accept,state,ok,message)
```

Initialization sets `lambda=0`.  A successful acceptance at lambda below one advances with `min(1d0,lambda+rate)`.  Acceptance at lambda one marks `finished`.  A rejected point returns a collective diagnostic failure without advancing state.  Validate every control collectively and fingerprint all persistent fields.

**Step 4: Run GREEN verification**

Run the MPI runner at `-np 1,2,4,8`; expect PASS at every size.  Run `git diff --check`.

**Step 5: Commit and checkpoint**

Commit only the four Task 1 files:

```bash
git commit -m "feat(dg): add collective interface continuation state"
```

Stop for the first review checkpoint.

---

### Task 2: Scale the complete SIPG operator consistently

**Files:**
- Modify: `src/gs/dc/dg_hybrid_schwarz_operator.f90`
- Modify: `tests/dg/test_dg_hybrid_schwarz_operator_mpi.f90`
- Modify: `tests/dg/run_dg_hybrid_schwarz_operator_mpi.py`

**Step 1: Extend the failing operator tests**

Pass an explicit `interface_scale` to the Hamiltonian action.  For the same input vectors, verify:

```text
H(0)   = kinetic + nonlocal + local-potential
H(0.5) = H(0) + 0.5 * interface
H(1)   = existing full Hamiltonian action
```

Check that the metric action is unchanged, `H(lambda)` remains Hermitian for zero, half, and one, and invalid or rank-disagreeing lambda is rejected collectively.

**Step 2: Run the test to verify RED**

Run `python3 tests/dg/run_dg_hybrid_schwarz_operator_mpi.py`; expect the new interface-scale contract to fail.

**Step 3: Implement one-point scaling**

Add `interface_scale` to `apply_dg_hybrid_schwarz_hamiltonian`.  Validate `0 <= interface_scale <= 1` and collective rank agreement.  Form the action as:

```fortran
kinetic_action + nonlocal_action + local_potential_action + interface_scale*interface_action
```

Do not scale the metric rows or any volume contribution.  Scale the complete already-assembled SIPG row block, never individual jump/penalty components separately.

**Step 4: Update all call sites and run GREEN verification**

Update tests and production callers to pass `1d0` until Task 3 owns the live lambda.  Run the Schwarz operator, solver, production-face, and divided-SCF MPI runners at their existing rank matrices.  Expect PASS.

**Step 5: Commit**

```bash
git commit -m "feat(dg): scale complete Schwarz interface action"
```

---

### Task 3: Drive Schwarz updates at fixed DC density

**Files:**
- Modify: `src/gs/main_dft.f90`
- Modify: `tests/dg/check_dg_hybrid_fragment_wannier_route.py`
- Create: `tests/dg/check_dg_hybrid_fixed_density_continuation_route.py`

**Step 1: Write the failing source-contract test**

Require the production divided entry to initialize the continuation state, update the potential from the restored DC density exactly once before the continuation loop, pass the current lambda to both the Schwarz Hamiltonian and preconditioner assembly, cap each point by `dg_hybrid_fragment_cg_steps`, and carry the accepted coefficients forward.

Reject production calls to density assembly, density mixing, Pulay/Broyden, or `run_dg_hybrid_divided_scf` inside this fixed-density route.  Require the existing exact MPI-rank/fragment and seed compatibility guards.

**Step 2: Run the tests to verify RED**

Run:

```bash
python3 tests/dg/check_dg_hybrid_fixed_density_continuation_route.py
python3 tests/dg/check_dg_hybrid_fragment_wannier_route.py
```

Expected: the new fixed-density contract fails while the existing route test remains green.

**Step 3: Implement the production continuation loop**

In `run_dg_hybrid_divided_ground_state_for_main`:

1. Keep `initial_density` from the ordinary-DC seed immutable.
2. Call the distributed potential update once with that density.
3. Initialize lambda using the configured `mixing%mixrate` value; do not initialize or mutate density-mixing history.
4. At each lambda, run one capped Schwarz epoch and common 300 K occupation assignment.
5. On acceptance, advance lambda and retain the updated coefficients.
6. On rollback/nonfinite diagnostics, retain the last accepted coefficients, report the last accepted lambda collectively, and stop without pretending convergence.

Store the live lambda in an explicit saved production-state value used only by the operator callbacks.  Never infer it from iteration count.

**Step 4: Make the preconditioner match the operator**

Assemble each local preconditioner block with:

```fortran
kinetic_rows + nonlocal_rows + local_potential_rows + lambda*interface_rows
```

Pass the same lambda to `apply_dg_hybrid_schwarz_hamiltonian`.  Leave `bounded_fixed_payload` immutable.

**Step 5: Run GREEN verification**

Run both source-contract tests, all Schwarz MPI runners, divided-operator runners, seed-route checks, and `cmake --build build-hybrid-release -j2`.  Expect PASS and a successful link.

**Step 6: Commit**

Stage only Task 3 hunks from the already-dirty files and commit:

```bash
git commit -m "feat(dg): continue interface at fixed DC density"
```

---

### Task 4: Continuation diagnostics and one terminal LCFO

**Files:**
- Modify: `src/gs/main_dft.f90`
- Modify: `tests/dg/check_dg_hybrid_divided_lcfo_route.py`
- Modify: `tests/dg/check_dg_hybrid_fixed_density_continuation_route.py`

**Step 1: Write failing diagnostic and terminal-state assertions**

Require one record per lambda containing lambda, accepted CG steps, residual, orthogonality defect, electron defect, Rayleigh/energy trace, interface-action norm, and accepted/rollback status.  Assert that terminal LCFO is unreachable before lambda one, is called exactly once at lambda one, uses the full unscaled interface rows, and is not followed by a density or potential update.

**Step 2: Run RED verification**

Run both source-contract checks; expect missing diagnostics and lambda-one guards.

**Step 3: Implement diagnostics and terminal guard**

Compute diagnostics collectively without gathering coefficient matrices.  Include lambda and continuation fingerprint in the log line.  Require `state%finished` and exact `lambda==1d0` before final row composition and `solve_dg_hybrid_generalized_once_and_publish`.  Label the published result as fixed-density/non-self-consistent in output.

**Step 4: Run GREEN verification**

Run route checks, generalized-eigensystem tests, Schwarz MPI tests, `git diff --check`, and the release build.  Expect PASS.

**Step 5: Commit**

```bash
git commit -m "test(dg): certify fixed-density terminal LCFO"
```

---

### Task 5: Remove the obsolete production density-mixing branch

**Files:**
- Modify: `src/gs/main_dft.f90`
- Modify: `tests/dg/check_dg_hybrid_fixed_density_continuation_route.py`
- Retain unless proven unused elsewhere: `src/gs/dc/dg_hybrid_divided_scf.f90`
- Retain unless proven unused elsewhere: `src/gs/dc/dg_hybrid_divided_mixing.f90`

**Step 1: Add a failing obsolete-route test**

Reject old production callbacks and saved state that only supported divided density mixing: potential-update-per-epoch, core-density assembly for feedback, production Pulay/Broyden history, and density-mixing rollback state.  Do not remove shared modules or tests if another supported route still imports them.

**Step 2: Run RED verification**

Run the fixed-density source-contract test; expect it to identify the remaining obsolete production symbols.

**Step 3: Remove only proven-dead production wiring**

Delete unused imports, saved arrays, callbacks, and messages from `main_dft.f90`.  Preserve diagnostic helpers and continuation code used by the supported RT/continuation routes.  Use repository-wide symbol searches before deleting either shared module.

**Step 4: Run GREEN and build verification**

Run all DG route checks, all affected MPI runners at `1,2,4,8` where supported, the release build, and `git diff --check`.  Expect PASS.

**Step 5: Commit**

```bash
git commit -m "refactor(dg): remove obsolete density-mixed production route"
```

---

### Task 6: Si64 fixed-seed diagnostic sweep

**Files:**
- Create: a timestamped directory under `/tmp` for each run
- Modify only if needed for reusable automation: `tests/dg/run_dg_hybrid_si64_schwarz.py`
- Update: `docs/plans/2026-09-06-dg-interface-continuation-design.md` with measured results

**Step 1: Reuse the exact compatible DC seed**

Use the Task 10 Si64 seed only with the same eight MPI ranks and exact rank-fragment mapping.  Confirm the publication ID, ownership-map fingerprint, and ordinary-DC skip before accepting the run.

**Step 2: Run the fixed schedule**

Run lambda `0,0.2,0.4,0.6,0.8,1.0`, three CG steps per point, 300 K occupations, and no density update.  Save the complete log and a compact CSV/table of diagnostics.

**Step 3: Interpret without loosening tolerances**

Report whether lambda zero reduces the residual and identify the first interval where residual or energy behavior degrades.  If lambda zero remains near the previous `~2e3` residual, record that the SIPG interface is not the primary cause and return to the Schwarz metric/preconditioner investigation.  Do not tune tolerances to force a pass.

**Step 4: Run the full regression suite**

Run the affected Python/MPI checks, release build, and `git diff --check`.  Confirm exactly one terminal LCFO and no post-LCFO density update in the run log.

**Step 5: Document and commit**

Append commands, seed identity, lambda table, interpretation, and log directories to the design/results document.  Stage only the report or runner changes and commit:

```bash
git commit -m "test(dg): diagnose Si64 interface continuation"
```
