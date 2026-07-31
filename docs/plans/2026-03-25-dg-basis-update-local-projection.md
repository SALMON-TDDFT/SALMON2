# DG Basis Update Local Projection Implementation Plan

> **Historical/removed:** This document describes an obsolete experimental DG route
> removed on 2026-07-31. It is retained only as an implementation record and is
> not executable guidance.

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Replace basis-update-time global dense mixed re-diagonalization with local overlap-based projection so DG RT basis updates can run without keeping global dense mixed operators resident.

**Architecture:** Treat basis update as a state-transfer problem from old basis to new basis. Use fragment-fragment and fragment-PW overlap blocks (`FF/FP/PP`) to compute the new propagated coefficients. Dense mixed operators become optional fallback/debug state, not normal runtime dependencies.

**Tech Stack:** Fortran, MPI, OpenMP, LAPACK, SALMON DG-RT mixed fragment/PW runtime

---

### Task 1: Inventory and isolate the normal basis-update dense path

**Files:**
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/docs/plans/2026-03-25-dg-basis-update-local-projection.md`
- Inspect: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_fragment_basis_update.f90`
- Inspect: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_plane_wave.f90`
- Inspect: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_basis_projection.f90`

**Step 1: Write the failing source audit**

Record the routines that still require normal-update access to:

- `H_mat_mixed`
- `S_mat_mixed_prop`
- global mixed dense diagonalization during basis update

**Step 2: Run the audit to verify it fails**

Run:

```bash
rg -n "H_mat_mixed|S_mat_mixed_prop|diagonalize_mixed_basis|diagonalize_full_system_dg" /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg
```

Expected: normal basis-update paths are still reported

**Step 3: Record the replacement map**

Document which routines will be:

- removed from the normal path
- kept as fallback only
- replaced by local overlap projection

**Step 4: Re-run the audit**

Confirm this task remains intentionally red until Tasks 3-5 land.

**Step 5: Commit**

```bash
git add docs/plans/2026-03-25-dg-basis-update-local-projection.md
git commit -m "docs: inventory dg basis-update dense dependencies"
```

### Task 2: Add mixed overlap transfer helpers from FF/FP/PP

**Files:**
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_fragment_ops.f90`
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_fragment_types.f90`

**Step 1: Write the failing test/check**

Add a source-level check or temporary assertion plan for a helper contract that can build/apply mixed overlap transfer using only:

- fragment overlap blocks
- `S_mat_frag_pw`
- PW identity block

**Step 2: Run check to verify it fails**

Confirm no such helper exists yet.

**Step 3: Write minimal implementation**

Add helpers that can:

- assemble the overlap action for old->new mixed transfer
- apply that action to a batch of coefficient columns
- optionally build only the local dense scratch required for a small solve

Do not introduce new persistent dense mixed arrays.

**Step 4: Run build**

Run:

```bash
make -C /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/build -j4
```

Expected: build succeeds

**Step 5: Commit**

```bash
git add src/rt/dg/rt_dg_fragment_ops.f90 src/rt/dg/rt_dg_fragment_types.f90
git commit -m "feat: add local mixed overlap transfer helpers"
```

### Task 3: Replace fragment-only update transfer with local projection

**Files:**
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_fragment_basis_update.f90`
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_basis_projection.f90`
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_fragment_ops.f90`

**Step 1: Write the failing test/check**

Identify the normal fragment-only update path that still calls full-system re-diagonalization for coefficient transfer.

**Step 2: Run check to verify it fails**

Confirm the current path still reaches:

- `diagonalize_full_system_dg`
- `diagonalize_and_update_basis`

for the normal fragment-only update transfer.

**Step 3: Write minimal implementation**

Replace that handoff with:

- compute new-vs-old fragment overlap
- solve the local projection transfer
- update propagated coefficients without rebuilding a new global eigenbasis

Keep any dense work strictly local to the projection solve.

**Step 4: Run build**

Run:

```bash
make -C /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/build -j4
```

Expected: build succeeds

**Step 5: Commit**

```bash
git add src/rt/dg/rt_dg_fragment_basis_update.f90 src/rt/dg/rt_dg_basis_projection.f90 src/rt/dg/rt_dg_fragment_ops.f90
git commit -m "refactor: use local projection for dg basis transfer"
```

### Task 4: Replace mixed basis update transfer with FF/FP/PP projection

**Files:**
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_fragment_basis_update.f90`
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_plane_wave.f90`
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_fragment_ops.f90`

**Step 1: Write the failing test/check**

Identify the mixed update path that still depends on dense mixed diagonalization to transfer propagated coefficients.

**Step 2: Run check to verify it fails**

Confirm the normal mixed update path still reaches `diagonalize_mixed_basis`.

**Step 3: Write minimal implementation**

Replace the mixed update transfer with overlap-based projection using:

- new fragment/new PW basis
- old fragment/old PW coefficients
- `FF/FP/PP` overlap data

Keep the old dense mixed diagonalization only behind an explicit fallback guard if needed.

**Step 4: Run build**

Run:

```bash
make -C /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/build -j4
```

Expected: build succeeds

**Step 5: Commit**

```bash
git add src/rt/dg/rt_dg_fragment_basis_update.f90 src/rt/dg/rt_dg_plane_wave.f90 src/rt/dg/rt_dg_fragment_ops.f90
git commit -m "refactor: project mixed dg basis updates without dense evp"
```

### Task 5: Remove normal-path dependence on persistent mixed dense operators

**Files:**
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_plane_wave.f90`
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_integrator_derivative.f90`
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_integrator_unitarity.f90`
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_fragment_basis_update.f90`

**Step 1: Write the failing test/check**

Add a check list for any remaining normal RT dependency on:

- persistent `H_mat_mixed`
- persistent `S_mat_mixed_prop`

outside explicit fallback/debug paths.

**Step 2: Run check to verify it fails**

Confirm the current branch still reports those normal-path references.

**Step 3: Write minimal implementation**

Make normal RT use:

- `H_mat_frag_pw`
- `H_mat_pw_diag`
- `S_mat_frag_pw`
- fragment block operators

directly, without requiring persistent mixed dense state.

**Step 4: Run build**

Run:

```bash
make -C /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/build -j4
```

Expected: build succeeds

**Step 5: Commit**

```bash
git add src/rt/dg/rt_dg_plane_wave.f90 src/rt/dg/rt_dg_integrator_derivative.f90 src/rt/dg/rt_dg_integrator_unitarity.f90 src/rt/dg/rt_dg_fragment_basis_update.f90
git commit -m "refactor: drop normal-path dense mixed dg state"
```

### Task 6: Validate memory-focused behavior and record remaining fallback debt

**Files:**
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/docs/plans/2026-03-25-dg-basis-update-local-projection.md`
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/docs/plans/2026-03-25-dg-basis-update-local-projection-design.md`

**Step 1: Write the validation checklist**

Record the required runtime checks:

- basis update trigger in fragment-only mode
- basis update trigger in mixed mode
- `nproc_frag > 1`
- memory before/after update trigger

**Step 2: Run available local verification**

Run local build verification and note any runtime gaps that require Fugaku.

**Step 3: Update the docs**

Document:

- which dense mixed dependencies were removed from the normal path
- which fallback paths remain
- what still requires explicit runtime validation on target hardware

**Step 4: Commit**

```bash
git add docs/plans/2026-03-25-dg-basis-update-local-projection.md docs/plans/2026-03-25-dg-basis-update-local-projection-design.md
git commit -m "docs: record dg basis-update local projection validation"
```
