# Direct Retained-Frame Wannier Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Replace production spectral-basin trial construction with the already orthonormal retained frame and its directly matching Wannier symmetry representation.

**Architecture:** Keep the adapted occupied plus occupied-orthogonal s/p retained frame unchanged, expose a small row-owned direct-frame preparation primitive, and feed its identity coefficients and band representation into the single Wannier90 call. Remove basin construction, eigensolves, propagation, and basin-induced DMN actions from `main_dft.f90`, while retaining reusable basin primitives and focused historical tests.

**Tech Stack:** Fortran 2008, MPI, BLAS/LAPACK, EigenExa, Wannier90 SAWF/DMN adapter, Python source-contract and MPI fixture runners.

---

### Task 1: Add the production-route RED

**Files:**
- Modify: `tests/dg/check_dg_overlapping_wannier_route.py:220-255`
- Modify: `tests/dg/check_dg_overlapping_wannier_route.py:820-845`
- Modify: `tests/dg/check_dg_overlapping_wannier_route.py:1620-1645`

**Step 1: Write the failing source-contract assertions**

Replace the existing complement-basin production assertions with assertions
that the production ground-state subroutine contains none of these calls:

```python
for obsolete_call in (
    "build_dg_periodic_spectral_basins",
    "prepare_dg_spectral_basin_operators",
    "project_dg_prepared_spectral_basin_operator",
    "diagonalize_dg_spectral_basin_operator",
    "select_dg_spectral_basin_channel_ranks",
    "propagate_dg_spectral_basin_orbit_channels",
    "build_dg_spectral_channel_generator_actions",
    "compose_dg_occupied_complement_trial_rows",
):
    assert not re.search(rf"call\\s+{obsolete_call}\\b", ow_ground_state_body)
```

Require one call to the new direct-frame primitive, one Wannier90 call, and
the existing post-Wannier route order. Require the DMN target action to be
copied from the retained band action rather than built from basin labels.

**Step 2: Run the route test and verify RED**

Run:

```bash
python3 tests/dg/check_dg_overlapping_wannier_route.py
```

Expected: FAIL because the production subroutine still calls spectral-basin
construction and does not call the direct-frame primitive.

**Step 3: Commit the RED**

```bash
git add tests/dg/check_dg_overlapping_wannier_route.py
git commit -m "test(dg): require direct retained Wannier frame"
```

### Task 2: Add a row-owned direct-frame primitive with focused REDs

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_construction.f90:105-125`
- Modify: `src/gs/dc/dg_overlapping_wannier_construction.f90:122-258`
- Modify: `tests/dg/test_dg_overlapping_wannier_construction_mpi.f90`
- Modify: `tests/dg/run_dg_overlapping_wannier_construction_mpi.py`

**Step 1: Write the failing MPI fixture**

Add a focused fixture for a six-state retained frame split unevenly across MPI
ranks. The fixture must call:

```fortran
call prepare_dg_direct_retained_wannier_frame(comm,row_ids,6,band_action_rows,&
  retained_fingerprint,operation_fingerprint,1d-12,trial_rows,wannier_action_rows,&
  frame_defect,action_defect,frame_fingerprint,workspace_bytes,ok,message)
```

Assert:

- `trial_rows(p,row_ids(p)) == 1` and every other entry is zero;
- its distributed Gram is identity;
- `wannier_action_rows == band_action_rows` locally;
- a nontrivial 2x2 rotation within a repeated representation block still
  intertwines the copied band/target actions;
- duplicate, missing, and out-of-range row IDs reject collectively;
- rank-disagreeing dimensions/fingerprints reject collectively; and
- the output fingerprint is identical for MPI 1/2/4/8.

**Step 2: Run the fixture and verify RED**

Run:

```bash
python3 tests/dg/run_dg_overlapping_wannier_construction_mpi.py
```

Expected: compile failure because
`prepare_dg_direct_retained_wannier_frame` does not exist.

**Step 3: Implement the minimal primitive**

Add a public routine that:

1. collectively agrees `global_row_count`, operation count, tolerance, and
   both nonzero provenance fingerprints;
2. validates exact global ownership with occurrence counts, including
   same-rank duplicates;
3. checks default-integer and int64 extents before allocation;
4. allocates `trial_rows(nlocal,nstate)` and
   `wannier_action_rows(nlocal,nstate,noperation)` with `stat=` and collective
   allocation consensus;
5. fills row-owned identity coefficients;
6. copies the already row-owned band actions into the target actions;
7. measures the distributed identity Gram and copied-action difference;
8. produces a decomposition-independent streamed fingerprint bound to row
   IDs, retained-frame provenance, operation provenance, and copied action;
9. reports a checked conservative owned-workspace receipt; and
10. cleans partial outputs on any collective failure.

Do not perform an eigensolve, all-gather a dense frame, or allocate basin
metadata.

**Step 4: Run focused MPI tests**

Run:

```bash
python3 tests/dg/run_dg_overlapping_wannier_construction_mpi.py
```

Expected: PASS on MPI 1/2/4/8, including one common direct-frame fingerprint.

**Step 5: Commit**

```bash
git add src/gs/dc/dg_overlapping_wannier_construction.f90 \
  tests/dg/test_dg_overlapping_wannier_construction_mpi.f90 \
  tests/dg/run_dg_overlapping_wannier_construction_mpi.py
git commit -m "feat(dg): prepare direct retained Wannier frame"
```

### Task 3: Replace the production basin route

**Files:**
- Modify: `src/gs/main_dft.f90:60-80`
- Modify: `src/gs/main_dft.f90:590-730`
- Modify: `src/gs/main_dft.f90:1308-1788`

**Step 1: Run the route RED again**

Run:

```bash
python3 tests/dg/check_dg_overlapping_wannier_route.py
```

Expected: FAIL on the old basin calls.

**Step 2: Remove production basin state and calls**

Delete production declarations, allocations, loops, diagnostics, and cleanup
for basin labels, spectra, selected ranks, orbit representatives, prepared
operators, and basin-induced action construction. Remove now-unused basin
procedures from the `use` list in `main_dft.f90`; do not delete their reusable
module implementations in this task.

**Step 3: Prepare the direct trial frame**

After the retained symmetry action has been assembled and validated, call the
new primitive with the retained row IDs and row-owned band action. Use the
returned identity coefficients as `spectral_trial_rows` (rename production
variables to `direct_trial_rows` where this makes ownership clear) and the
copied actions as the DMN target actions.

Bind the direct-frame fingerprint into `w90_input_fingerprint` and the final
checkpoint provenance. Keep the single `run_dg_w90_gamma_library` call and all
post-Wannier sector/Gamma/cocycle processing unchanged.

**Step 4: Avoid redundant spatial materialization**

Where `materialize_dg_row_owned_sector_on_spatial_grid` would multiply the
retained spatial frame by identity, pass/reuse `global_closed_core` directly
as the Wannier anchor view if its storage orientation matches. If the adapter
requires the transposed layout, retain only the existing unavoidable adapter
buffer and document its receipt; do not construct a second coefficient-space
identity product.

**Step 5: Run route checks**

Run:

```bash
python3 tests/dg/check_dg_overlapping_wannier_route.py
python3 tests/dg/check_obsolete_dg_routes_removed.py
```

Expected: PASS; production has zero basin calls and exactly one Wannier90 call.

**Step 6: Compile the production overlay**

Run the repository's existing overlapping-Wannier overlay/build command used
by the focused runners. Expected: successful compilation with no new unused
or undefined production symbols.

**Step 7: Commit**

```bash
git add src/gs/main_dft.f90 tests/dg/check_dg_overlapping_wannier_route.py
git commit -m "fix(dg): bypass redundant spectral basin trials"
```

### Task 4: Regression and numerical verification

**Files:**
- Modify only if a genuine regression is exposed:
  `src/gs/dc/dg_overlapping_wannier_construction.f90`
- Do not modify the user-owned Si64 input/checker files.

**Step 1: Run focused suites**

Run:

```bash
python3 tests/dg/run_dg_overlapping_wannier_construction_mpi.py
python3 tests/dg/run_dg_overlapping_wannier_eigenexa_mpi.py
python3 tests/dg/run_dg_overlapping_wannier_w90_mpi.py
python3 tests/dg/check_dg_overlapping_wannier_route.py
python3 tests/dg/check_obsolete_dg_routes_removed.py
git diff --check
```

Expected: PASS. Remove generated `w90_converged_fixture.wout` and
`w90_exhausted_fixture.wout` from the worktree after the W90 fixture; never
stage them.

**Step 2: Run Si64 through Wannier90 completion**

Use the existing Si64 MPI-rank count and input without editing user-owned
files. Record:

- maximum RSS per rank before and after Wannier90;
- direct-frame Gram/action defects;
- Wannier90 iteration, convergence status, and final spread;
- post-Wannier unitarity, Gamma, and cocycle defects; and
- whether execution passes the former basin Gram failure point.

Expected: no spectral-basin diagnostics, entry into Wannier90, and either
successful Wannier90 completion or a bounded explicit convergence failure.
An unbounded run is not acceptable.

**Step 3: Apply the design decision gate**

If Wannier90 localizes with finite acceptable spread and all symmetry receipts
pass, retain `D_wann = D_band`. If it is demonstrably overconstrained, stop;
capture the evidence and write a separate design for an algebraic regular
target representation. Do not restore basin selection or weaken symmetry
checks in this task.

**Step 4: Commit any evidence-only test adjustment**

Commit only repository-owned test or source changes. Leave user-owned Si64
files and generated W90 logs untouched.

### Task 5: Final review and handoff

**Files:**
- Review: `src/gs/main_dft.f90`
- Review: `src/gs/dc/dg_overlapping_wannier_construction.f90`
- Review: `tests/dg/check_dg_overlapping_wannier_route.py`
- Review: focused MPI fixtures and runners

**Step 1: Inspect the final diff**

Confirm the diff contains no user-owned Si64 edits and no generated W90 logs.
Confirm no unrelated dirty-file content is staged.

**Step 2: Request code review**

Use `superpowers:requesting-code-review` and ask reviewers specifically about:

- whether `D_wann = D_band` uses the correct action orientation;
- whether the direct-frame provenance reaches checkpoint output;
- MPI collective safety and exact row ownership;
- whether any basin work remains in production; and
- whether identity materialization duplicates an `N x M` buffer.

**Step 3: Address findings and rerun verification**

Use `superpowers:receiving-code-review` for actionable findings. Rerun every
focused check and `git diff --check` after corrections.

**Step 4: Commit final corrections**

Create narrowly scoped commits and report Si64 progress, spread, symmetry
receipts, and measured memory rather than claiming success from source checks
alone.

