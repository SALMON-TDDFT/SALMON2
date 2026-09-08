# Hybrid DC Symmetry Handoff Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Allow the eight-rank Si64 Hybrid continuation route to defer physical symmetry recovery to the final LCFO solve while retaining exact DC/DG decomposition identity and authoritative Wannier provenance.

**Architecture:** Add an explicit LCFO-deferred production-selection mode.  It uses identity bookkeeping only for fragment-local PW packet ownership, records and hashes the already-certified Wannier basis fingerprint, and never interprets that bookkeeping action as the physical symmetry group.  The existing final LCFO actual-group occupied-projector acceptance remains authoritative.

**Tech Stack:** Fortran 2008, MPI, Python test runners, CMake SALMON build.

---

### Task 1: Specify LCFO-deferred production selection

**Files:**
- Modify: `src/common/dg_hybrid_windowed_pw_types.f90`
- Modify: `src/common/dg_hybrid_production_pw_basis.f90`
- Modify: `tests/dg/test_dg_hybrid_production_pw_basis_mpi.f90`
- Modify: `tests/dg/run_dg_hybrid_production_pw_basis_mpi.py`

**Step 1: Write the failing MPI test**

Add a test fixture whose physical row action contains
`mixed_row_action(:,3)=[1,3,2,4]`, which cuts both two-point fragments.  Call a
new `analyze_dg_hybrid_lcfo_selection` entry point with a nonzero certified
Wannier fingerprint.  Require:

```fortran
call require(ok,'LCFO-deferred production selection rejected split fragment action: '//trim(message))
call require(selection%lcfo_symmetry_deferred,'LCFO deferral provenance is missing')
call require(selection%wannier_symmetry_fingerprint==701_int64,'Wannier provenance was not retained')
call require(selection%operation_count==1.and.selection%identity_only,&
  'fragment-local bookkeeping must use one explicit identity action')
call require(all(selection%requested_packet_ids==selection%packet_ids),&
  'LCFO-deferred preparation must retain the complete PW packet catalog')
```

Call the same entry point with a zero Wannier fingerprint and require rejection
containing `Wannier symmetry provenance`.

**Step 2: Run RED**

Run:

```text
python3 tests/dg/run_dg_hybrid_production_pw_basis_mpi.py
```

Expected: compilation fails because the new fields and entry point do not yet
exist.

**Step 3: Implement the minimal selection mode**

Extend `s_dg_hybrid_production_selection` with:

```fortran
logical :: lcfo_symmetry_deferred=.false.
integer(8) :: wannier_symmetry_fingerprint=0_8
```

Add `analyze_dg_hybrid_lcfo_selection` beside the existing strict analyzer.
It must:

1. collectively reject a zero or rank-inconsistent Wannier fingerprint;
2. construct one distributed identity row action from `core_ids` and one
   identity reciprocal rotation;
3. call the existing analyzer using only that operational identity;
4. retain every packet by copying `packet_ids` to `requested_packet_ids`;
5. set `lcfo_symmetry_deferred=.true.` and record the Wannier fingerprint;
6. recompute `analysis_fingerprint` after those fields and the complete packet
   selection are set.

Hash both new fields in `production_selection_fingerprint`.  Make
`freeze_dg_hybrid_production_selection` reject a deferred receipt with a zero
Wannier fingerprint.  Do not weaken or modify the existing strict analyzer;
its known-nonidentity split-fragment test must continue to fail.

**Step 4: Run GREEN**

Run:

```text
python3 tests/dg/run_dg_hybrid_production_pw_basis_mpi.py
```

Expected: PASS on 1, 2, 4, and 8 ranks, including both the strict rejection and
LCFO-deferred acceptance fixtures.

**Step 5: Commit task-scoped hunks**

Use `git add -p` because all four files were dirty before this task.  Stage only
the LCFO-deferred fields, analyzer, fingerprint validation, and associated
tests.  Commit as:

```text
feat(dg): defer Hybrid fragment symmetry to LCFO
```

### Task 2: Wire the Hybrid continuation handoff

**Files:**
- Modify: `src/gs/main_dft.f90`
- Modify: `tests/dg/check_dg_hybrid_continuation_route.py`

**Step 1: Write RED route assertions**

Require the Hybrid preparation callback to call
`analyze_dg_hybrid_lcfo_selection` with `basis_fingerprint_arg`.  Require it to
keep the existing DC-derived fragment origins, extents, and core ownership as
the DG inputs.  Reject source patterns that call the strict physical-group
analyzer from `prepare_dg_hybrid_divided_production_basis`.

Also require a diagnostic receipt containing:

```text
[HYBRID-LCFO-SYMMETRY-HANDOFF]
```

with nonzero Wannier and production fingerprints.

**Step 2: Run RED**

Run:

```text
python3 tests/dg/check_dg_hybrid_continuation_route.py
```

Expected: FAIL because the callback still invokes strict fragment covariance
analysis and emits no LCFO handoff receipt.

**Step 3: Implement the minimal wiring**

In `prepare_dg_hybrid_divided_production_basis`, replace only the production
selection call with `analyze_dg_hybrid_lcfo_selection`, passing the existing
`basis_fingerprint_arg` as authoritative Wannier provenance.  Keep
`core_fragment_ids` derived from the DC fragment origins and extents; do not
introduce a second DG partition.

After successful preparation, rank zero must emit:

```fortran
'[HYBRID-LCFO-SYMMETRY-HANDOFF] wannier_fingerprint=',basis_fingerprint_arg,&
  ' production_fingerprint=',pw_fingerprint_arg
```

Do not remove any Wannier construction checks or any final LCFO occupied-space
acceptance checks.

**Step 4: Run GREEN and focused MPI tests**

Run:

```text
python3 tests/dg/check_dg_hybrid_continuation_route.py
python3 tests/dg/run_dg_hybrid_production_pw_basis_mpi.py
python3 tests/dg/run_dg_hybrid_continuation_acceptance_mpi.py
```

Expected: PASS.

**Step 5: Build**

Run:

```text
cmake --build build-hybrid-commit -j2
```

Expected: successful SALMON build.

**Step 6: Commit**

Stage only `main_dft.f90` and the route-contract hunks.  Commit as:

```text
fix(dg): hand Hybrid symmetry recovery to LCFO
```

### Task 3: Fresh protected and Si64 verification

**Files:**
- Modify only if a new failing test proves a defect.
- Preserve generated evidence under `verification-si64-dg-continuation-20260830/` without staging it.

**Step 1: Run focused and protected verification freshly**

Run every focused runner and protected route listed in Task 12 of
`docs/plans/2026-08-27-wpw-dg-continuation-ground-state.md`, including:

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

Expected: PASS.  The strict production analyzer must still reject a known
nonidentity operation that cuts fragments; only the explicitly deferred LCFO
handoff may bypass that intermediate condition.

**Step 2: Run the final Si64 calculation**

From a new clean output directory run:

```text
python3 tests/dg/run_dg_hybrid_si64_continuation_rt.py \
  build-hybrid-commit/salmon \
  verification-si64-dg-continuation-20260830/si64-run \
  --ranks 8
```

Use `OMP_NUM_THREADS=1` through the runner and no timeout.  Preserve complete
GS and RT logs on pass, failure, or manual interruption.

**Step 3: Verify final evidence**

Require the new LCFO symmetry-handoff receipt, all existing Wannier affine
proofs, lambda-zero and lambda-one continuation receipts, final actual-group
occupied-projector covariance, matching GS/RT payload fingerprints, and all
zero-field RT stationarity receipts.  Missing evidence is failure.

**Step 4: Request code review**

Use `superpowers:requesting-code-review` on the two implementation commits.
Resolve Critical and Important findings with TDD and task-scoped commits, then
rerun affected verification.

**Step 5: Finish only after fresh success**

Invoke `superpowers:verification-before-completion`, then
`superpowers:finishing-a-development-branch`.  Present integration options and
do not merge automatically.
