# Hybrid Terminal LCFO Refinement Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Make the divided Hybrid WF+PW ground-state route perform one complete LCFO generalized eigensolve in the normal case and at most three additional Hartree+XC/local-potential refinements, while removing the whole-system eigensolve from every local continuation step.

**Architecture:** Keep the converged DC density immutable during the one-fragment-per-rank Schwarz/interface phase and apply the full DG/SIPG coupling from its first local step. Move all complete-system diagonalization, 300 K occupations, density reconstruction, and final spectral certification into a terminal refinement driver. The driver reuses the existing history-aware density mixer, changes only the local-potential matrix, exits after convergence, and publishes the last finite valid state with a named warning after three unsuccessful additional solves. Preserve the byte-level distributed-v5 format; store the new audit data in a separately authenticated receipt bound to the v5 publication fingerprint.

**Tech Stack:** Fortran 2008, MPI, ScaLAPACK `PZHEEVD`, SALMON DC/LCFO modules, Python source-contract tests, CMake/CTest, SHA-256 checkpoint authentication.

---

Implementation uses the existing ordinary clone
`/tmp/salmon-v230-dg-integration-20260907/repository` on branch
`codex/wpw-s-orthogonal-complement-v2.3.0`. Do not create another worktree and
do not modify the original dirty worktree.

## Task 1: Add a testable terminal-refinement state machine

**Files:**

- Create: `src/gs/dc/dg_hybrid_terminal_refinement.f90`
- Modify: `src/gs/dc/CMakeLists.txt`
- Create: `tests/dg/test_dg_hybrid_terminal_refinement_mpi.f90`
- Create: `tests/dg/run_dg_hybrid_terminal_refinement_mpi.py`

**Step 1: Write the failing MPI policy test**

Follow `tests/dg/run_dg_hybrid_divided_mixing_mpi.py`. The Fortran test must
cover, identically on 1, 2, 4, and 8 ranks:

- a first sample satisfying density and energy tolerances exits with one total solve;
- sequences converge after total solve counts 2, 3, and 4;
- a fourth failed sample marks `exhausted=.true.`, `converged=.false.`, and `publish_last_valid=.true.`;
- no request for a fifth solve is possible;
- nonfinite metrics and rank-disagreeing controls/samples are rejected collectively;
- the receipt reports total/additional solve counts, final changes, and a nonzero fingerprint.

Run:

```bash
python3 tests/dg/run_dg_hybrid_terminal_refinement_mpi.py
```

Expected: FAIL because the module does not exist.

**Step 2: Implement the minimal policy module**

Add control and receipt types plus collective initialize/observe routines.
Count the mandatory terminal solve as solve 1 and hard-cap additional solves
at three. A finite fourth failed sample selects warning publication, not an
error. Keep eigensolver and density-mixer code out of this module.

**Step 3: Verify and commit**

```bash
python3 tests/dg/run_dg_hybrid_terminal_refinement_mpi.py
git diff --check
git add src/gs/dc/dg_hybrid_terminal_refinement.f90 src/gs/dc/CMakeLists.txt tests/dg/test_dg_hybrid_terminal_refinement_mpi.f90 tests/dg/run_dg_hybrid_terminal_refinement_mpi.py
git commit -m "feat(dg): add terminal LCFO refinement policy"
```

## Task 2: Lock down the local/global solve boundary

**Files:**

- Modify: `tests/dg/check_dg_hybrid_divided_lcfo_route.py`
- Create: `tests/dg/check_dg_hybrid_terminal_refinement_route.py`
- Modify: `src/gs/main_dft.f90`

**Step 1: Write failing source-contract tests**

Assert that the local Schwarz/interface section contains no
`solve_dg_hybrid_generalized_*` call, never replaces `initial_density` with a
reconstructed LCFO density, reaches full DG/SIPG coupling, and uses the bounded
one-fragment-per-rank local solver. Assert that the terminal section contains
one syntactic LCFO solve call inside one loop governed by the new policy.
Also assert that `run_dg_hybrid_concrete_continuation` can no longer be selected
as a production route, because it diagonalizes inside `stage_pass`.

Run both route tests. Expected: the new test FAILS because the legacy route is
still selectable and the adaptive loop is absent.

**Step 2: Refactor route selection**

Make both supported Hybrid GS selectors enter
`run_dg_hybrid_divided_ground_state_for_main`, preserving input compatibility.
Remove the production call to `run_dg_hybrid_concrete_continuation`; retain its
body until Task 7 proves no remaining references. Add named local and terminal
source ranges so tests verify the correct scope.

**Step 3: Verify and commit**

```bash
python3 tests/dg/check_dg_hybrid_divided_lcfo_route.py
python3 tests/dg/check_dg_hybrid_terminal_refinement_route.py
git diff --check
git add src/gs/main_dft.f90 tests/dg/check_dg_hybrid_divided_lcfo_route.py tests/dg/check_dg_hybrid_terminal_refinement_route.py
git commit -m "refactor(dg): isolate terminal LCFO from local continuation"
```

## Task 3: Implement adaptive terminal LCFO refinement

**Files:**

- Modify: `src/gs/main_dft.f90`
- Modify: `src/gs/dc/dg_hybrid_ground_state_types.f90`
- Modify: `tests/dg/test_dg_hybrid_generalized_eigensystem_mpi.f90`
- Modify: `tests/dg/run_dg_hybrid_generalized_eigensystem_mpi.py`
- Modify: `tests/dg/check_dg_hybrid_terminal_refinement_route.py`

**Step 1: Add failing executable cases**

Extend the generalized-eigensystem test with callback-counting density/energy
sequences that produce total solve counts 1, 2, 3, and 4. Check 300 K
occupations and electron conservation for every final state. Check that an
exhausted case retains the coefficients, occupations, eigenvalues, and
operator fingerprint from solve 4.

**Step 2: Build the terminal loop**

In `run_dg_hybrid_divided_ground_state_for_main`:

1. Initialize density history from the immutable DC density.
2. Project the local potential and compose `H=T+Vnl+Vsipg+Vlocal`.
3. Call the existing ScaLAPACK solve and 300 K occupation path.
4. Reconstruct distributed density and compute weighted relative density and total-energy changes.
5. Feed the sample to the new policy.
6. If requested, call `mix_dg_overlapping_wannier_density_history`, update Hartree+XC/local potential only, and repeat.
7. Exit early on convergence or retain solve 4 on finite exhaustion.

Use `dg_dc_gs_final_density_tolerance` and the established GS final-energy
tolerance; do not add a new input. Store actual total/additional solve counts
and convergence/exhaustion flags in ground-state audit fields.

**Step 3: Verify and commit**

```bash
python3 tests/dg/run_dg_hybrid_terminal_refinement_mpi.py
python3 tests/dg/run_dg_hybrid_generalized_eigensystem_mpi.py
python3 tests/dg/check_dg_hybrid_terminal_refinement_route.py
git diff --check
git add src/gs/main_dft.f90 src/gs/dc/dg_hybrid_ground_state_types.f90 tests/dg/test_dg_hybrid_generalized_eigensystem_mpi.f90 tests/dg/run_dg_hybrid_generalized_eigensystem_mpi.py tests/dg/check_dg_hybrid_terminal_refinement_route.py
git commit -m "feat(dg): refine terminal LCFO at most three times"
```

## Task 4: Prove only the local potential changes

**Files:**

- Modify: `src/gs/dc/dg_hybrid_terminal_refinement.f90`
- Create: `tests/dg/test_dg_hybrid_terminal_operator_guard_mpi.f90`
- Create: `tests/dg/run_dg_hybrid_terminal_operator_guard_mpi.py`
- Modify: `src/gs/main_dft.f90`

**Step 1: Write the failing immutable-component test**

Accept changes to local rows and the recomposed Hamiltonian only. Reject a
one-bit mutation to metric, kinetic, nonlocal, SIPG, basis generation, row
ownership, fixed-payload fingerprint, or immutable DC seed density, including
a mutation on just one rank.

**Step 2: Implement and wire the epoch guard**

Record collective fingerprints at terminal entry and validate them before and
after every additional solve. Guard failures are fatal; the nonfatal exhaustion
policy applies only to finite density/energy nonconvergence.

**Step 3: Verify and commit**

```bash
python3 tests/dg/run_dg_hybrid_terminal_operator_guard_mpi.py
python3 tests/dg/check_dg_hybrid_terminal_refinement_route.py
git diff --check
git add src/gs/dc/dg_hybrid_terminal_refinement.f90 src/gs/main_dft.f90 tests/dg/test_dg_hybrid_terminal_operator_guard_mpi.f90 tests/dg/run_dg_hybrid_terminal_operator_guard_mpi.py
git commit -m "fix(dg): guard fixed operators during LCFO refinement"
```

## Task 5: Add an authenticated refinement receipt without changing v5

**Files:**

- Create: `src/rt/dg/rt_dg_hybrid_refinement_receipt.f90`
- Modify: `src/rt/dg/CMakeLists.txt`
- Modify: `src/gs/main_dft.f90`
- Create: `tests/dg/test_rt_dg_hybrid_refinement_receipt_mpi.f90`
- Create: `tests/dg/run_rt_dg_hybrid_refinement_receipt_mpi.py`
- Modify: `tests/dg/check_rt_dg_hybrid_distributed_v5.py`

**Step 1: Write failing round-trip and compatibility tests**

Test a version-1 receipt containing the v5 publication fingerprint, total and
additional counts, density/energy changes, convergence/exhaustion flags,
early-exit reason, and SHA-256 digest. Collectively reject truncation, metric
mutation, wrong v5 binding, and rank disagreement. Separately assert that v5
still uses exactly eight acceptance receipts and existing fixtures remain
readable.

**Step 2: Implement atomic authenticated publication**

Write the receipt after successful v5 publication using a temporary file plus
rename and full SHA-256. Bind it to v5 publication and MPI rank/fragment
identity. On exhaustion rank zero prints exactly once:

```text
[DG-HYBRID-REFINEMENT-WARNING] maximum additional LCFO solves exhausted; publishing last finite valid state
```

Do not change `s_rt_dg_hybrid_v5_shard` or its reader. RT may report the
companion receipt but must accept an unconverged flag and legacy v5 files with
no companion.

**Step 3: Verify and commit**

```bash
python3 tests/dg/run_rt_dg_hybrid_refinement_receipt_mpi.py
python3 tests/dg/run_rt_dg_hybrid_checkpoint_v5_mpi.py
python3 tests/dg/check_rt_dg_hybrid_distributed_v5.py
git diff --check
git add src/rt/dg/rt_dg_hybrid_refinement_receipt.f90 src/rt/dg/CMakeLists.txt src/gs/main_dft.f90 tests/dg/test_rt_dg_hybrid_refinement_receipt_mpi.f90 tests/dg/run_rt_dg_hybrid_refinement_receipt_mpi.py tests/dg/check_rt_dg_hybrid_distributed_v5.py
git commit -m "feat(dg): authenticate terminal refinement receipt"
```

## Task 6: Build and run focused production checks

**Files:**

- Modify: `tests/dg/run_dg_fragment_wf_production_smoke.py`
- Modify: `tests/dg/run_dg_hybrid_si64_divided_lcfo.py`
- Modify: `tests/dg/check_si64_overlapping_wannier_gate.py`

**Step 1: Add log and receipt assertions**

Require the good-seed H4 smoke to report exactly one global solve. Require all
production runs to report a count in `[1,4]`, 300 K occupation, and a matching
receipt. Keep the Si64 analyzer read-only and accept both legacy v5 without the
companion and new v5 with it.

**Step 2: Build and run focused checks**

```bash
cmake -S . -B /tmp/salmon-v230-terminal-refinement-build -DCMAKE_BUILD_TYPE=Release
cmake --build /tmp/salmon-v230-terminal-refinement-build -j 4
python3 tests/dg/run_salmon_v230_task6_production_compile.py
```

Then run the documented H4 production smoke and Si64 read-only analyzer. Reuse
a checksum-compatible Wannier90 cache instead of regenerating it.

**Step 3: Commit**

```bash
git add tests/dg/run_dg_fragment_wf_production_smoke.py tests/dg/run_dg_hybrid_si64_divided_lcfo.py tests/dg/check_si64_overlapping_wannier_gate.py
git commit -m "test(dg): cover adaptive terminal LCFO production"
```

## Task 7: Remove the superseded whole-system continuation route

**Files:**

- Modify: `src/gs/main_dft.f90`
- Modify: `src/gs/dc/CMakeLists.txt`
- Modify/Delete: continuation-only modules and tests proven unreferenced
- Modify: `tests/dg/check_obsolete_dg_routes_removed.py`
- Modify: `tests/dg/check_dg_hybrid_terminal_refinement_route.py`

**Step 1: Add the failing obsolete-route audit**

Forbid `run_dg_hybrid_concrete_continuation` and any controller used only by
it after a complete reference audit. Require all supported selectors to reach
the divided local-plus-terminal route. Do not delete shared residual, spectral,
acceptance, or checkpoint code.

**Step 2: Delete only proven-dead code**

Remove the old stage loop that diagonalizes at every lambda/SCF iteration and
unreferenced helpers/tests. Preserve input compatibility by mapping the old
selector to the new route with a concise alias diagnostic.

**Step 3: Verify and commit**

```bash
python3 tests/dg/check_obsolete_dg_routes_removed.py
python3 tests/dg/check_dg_hybrid_terminal_refinement_route.py
python3 tests/dg/check_dg_overlapping_wannier_route.py
git diff --check
git add -A src/gs/main_dft.f90 src/gs/dc tests/dg
git commit -m "refactor(dg): remove repeated global continuation solves"
```

## Task 8: Run full regression and certify the branch

**Files:**

- Modify: `docs/plans/2026-09-10-hybrid-terminal-lcfo-refinement-design.md`
- Create: `docs/plans/2026-09-10-hybrid-terminal-lcfo-refinement-results.md`

**Step 1: Run focused MPI gates and the complete inventory**

Run all new tests plus the existing generalized eigensystem, divided
mixing/operator, Schwarz solver/state, v5 checkpoint/initialization, production
smoke, and system-identity tests. Then enumerate the same 30 certified DG gates
and add the three new gates. Record command, exit status, log path, and SHA-256.
Expected: `PASS_FRESH_ALL_33`; if the inventory count changes, document the
exact addition/removal rather than forcing the number.

**Step 2: Run physical regressions**

- H4 current-binary GS handoff and zero-field Exp RT startup;
- fresh Si8 checkpoint miss, hit, and recovery with identical rank mapping;
- preserved Si64 read-only analysis without unnecessary Wannier90 work;
- conventional bulk-Si GS/RT to prove the non-DG route is unchanged.

Record solve/refinement counts, convergence flag, density/energy changes,
occupation/electron defect, fixed fingerprints, v5 digest, and companion digest.

**Step 3: Review and commit certification**

Run `git diff --check`, inspect status, and perform independent specification
and quality reviews. Resolve Critical and Important findings and rerun affected
gates. Commit the evidence summary:

```bash
git add docs/plans/2026-09-10-hybrid-terminal-lcfo-refinement-design.md docs/plans/2026-09-10-hybrid-terminal-lcfo-refinement-results.md
git commit -m "docs(dg): certify adaptive terminal LCFO refinement"
```

Do not update the original repository branch reference until the user accepts
the final checkpoint. Never modify or clean the original dirty worktree.
