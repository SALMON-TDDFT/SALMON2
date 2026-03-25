# DG Fragment Dense Elimination Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Eliminate remaining persistent dense fragment-only operator storage from DG-RT by making fragment block matrices authoritative while preserving the current mixed-basis dense path for a later phase.

**Architecture:** Phase 1 migrates fragment-only `S/H/momentum` consumers from dense arrays to block-based application or temporary local dense scratch. Dense overlap broadcast is removed from the steady-state path. Mixed fragment-PW dense operators remain for a separate follow-up phase.

**Tech Stack:** Fortran, MPI, OpenMP, CMake, SALMON DG-RT runtime modules

---

### Task 1: Inventory dense fragment-only operator consumers

**Files:**
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/docs/plans/2026-03-25-dg-fragment-dense-elimination.md`
- Inspect: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_integrator_derivative.f90`
- Inspect: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_integrator_unitarity.f90`
- Inspect: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_basis_projection.f90`
- Inspect: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_fragment_basis_update.f90`
- Inspect: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_fragment_basis_update_soi.f90`

**Step 1: Write the failing check**

Create a lightweight source audit script that fails if the above files still contain direct dense fragment-only accesses that are in scope for phase 1.

**Step 2: Run check to verify it fails**

Run a grep/perl-based audit and confirm the current branch still reports dense `S_mat`, `S_mat_prop`, `H_mat`, or `momentum_mat` usage in fragment-only code paths.

**Step 3: Record the initial inventory**

List each remaining dense consumer and mark whether it can move to:

- block apply helper
- temporary local dense scratch
- deferred mixed-basis path

**Step 4: Re-run the check**

Confirm the audit script still fails until each targeted use site is migrated.

**Initial inventory (2026-03-25)**

- Audit script: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/.worktrees/codex/dg-fragment-dense-elimination/tools/check_dg_fragment_dense_elimination.py`
- `rt_dg_integrator_unitarity.f90`
  - Dense fragment-only `S_mat_prop[_c]` and `S_mat[_c]` `matmul` calls remain.
  - Migration target: block apply helper for `n_pw == 0`; keep `S_mat_mixed_prop` for deferred mixed-basis use.
- `rt_dg_integrator_derivative.f90`
  - Dense fragment-only overlap copies into `S_eval` remain.
  - Migration target: temporary local dense scratch assembled from block overlap data for the eigensolver path.
  - Dense fragment-only `H_mat[_c]` and `momentum_mat[_c]` reads also remain.
  - Migration target: later phase; likely block apply or temporary local dense scratch depending on the call site.
- `rt_dg_basis_projection.f90`
  - Dense `H_mat(:, :, ispin)` copy remains.
  - Migration target: temporary local dense scratch for LAPACK input.
- `rt_dg_fragment_basis_update.f90`
  - Dense fragment-only `S_mat_prop[_c]` copies, `S_mat` stability checks, `matmul(S_mat, x)`, and `H_mat` construction/copies remain.
  - Migration target: block apply helper for operator action; temporary local dense scratch where eigensolver input still requires dense storage.
- `rt_dg_fragment_basis_update_soi.f90`
  - Dense fragment-only `S_mat_prop[_c]`, `S_mat[_c]` stability checks, and `H_mat` dense staging remain.
  - Migration target: block apply helper plus temporary local dense scratch; keep mixed-basis dense paths deferred.
- Task 1 remains intentionally red after this inventory update, which is the expected state until Tasks 4 and 5 migrate the remaining basis-update and projection consumers.

**Step 5: Commit**

```bash
git add docs/plans/2026-03-25-dg-fragment-dense-elimination.md
git commit -m "docs: inventory dg dense fragment-only consumers"
```

### Task 2: Remove steady-state dense overlap broadcast

**Files:**
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_fragment_ops.f90`
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_fragment_hamiltonian.f90`
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_fragment_hamiltonian_soi.f90`
- Test: source audit script from Task 1

**Step 1: Write the failing test/check**

Extend the audit to fail when `ensure_overlap_prop_available` still broadcasts persistent dense overlap arrays in the steady-state path.

**Step 2: Run check to verify it fails**

Run the audit and confirm it flags the current dense overlap broadcast.

**Step 3: Write minimal implementation**

Refactor overlap propagation so that:

- block overlap state is authoritative
- dense overlap copies are rebuilt only as temporary local scratch when unavoidable
- subgroup broadcast no longer depends on persistent dense overlap arrays

**Step 4: Run checks and build**

Run:

```bash
make -C /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/build -j4
```

Expected: build succeeds

**Step 5: Commit**

```bash
git add src/rt/dg/rt_dg_fragment_ops.f90 src/rt/dg/rt_dg_fragment_hamiltonian.f90 src/rt/dg/rt_dg_fragment_hamiltonian_soi.f90
git commit -m "refactor: remove dense overlap steady-state broadcast"
```

### Task 3: Migrate unitarity and derivative paths off dense fragment-only overlap

**Files:**
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_integrator_unitarity.f90`
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_integrator_derivative.f90`
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_fragment_ops.f90`
- Test: source audit script from Task 1

**Step 1: Write the failing test/check**

Add checks for direct fragment-only `matmul(dg_frag%S_mat...)` and equivalent dense overlap use in these files.

**Step 2: Run check to verify it fails**

Run the audit and confirm both files are still flagged.

**Step 3: Write minimal implementation**

Introduce or extend helper routines so these paths use:

- block apply for fragment-only overlap action
- dense mixed-basis arrays only when `n_pw > 0` requires the deferred mixed path

**Step 4: Run verification**

Run the audit and build:

```bash
make -C /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/build -j4
```

Expected: audit passes for these files, build succeeds

**Step 5: Commit**

```bash
git add src/rt/dg/rt_dg_integrator_unitarity.f90 src/rt/dg/rt_dg_integrator_derivative.f90 src/rt/dg/rt_dg_fragment_ops.f90
git commit -m "refactor: use block overlap operators in dg integrators"
```

### Task 4: Migrate basis projection and basis update fragment-only dense accesses

**Files:**
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_basis_projection.f90`
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_fragment_basis_update.f90`
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_fragment_basis_update_soi.f90`
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_fragment_ops.f90`
- Test: source audit script from Task 1

**Step 1: Write the failing test/check**

Add checks for direct dense fragment-only `H_mat`, `S_mat`, and `momentum_mat` accesses that are in phase 1 scope.

**Step 2: Run check to verify it fails**

Run the audit and confirm these files are still flagged.

**Step 3: Write minimal implementation**

Replace dense fragment-only accesses with:

- block apply helpers
- temporary dense work arrays assembled locally inside the one routine that needs LAPACK input

Do not remove mixed-basis dense use in code paths where `n_pw > 0`.

**Step 4: Run verification**

Run the audit and build:

```bash
make -C /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/build -j4
```

Expected: audit passes for phase 1 targets, build succeeds

**Step 5: Commit**

```bash
git add src/rt/dg/rt_dg_basis_projection.f90 src/rt/dg/rt_dg_fragment_basis_update.f90 src/rt/dg/rt_dg_fragment_basis_update_soi.f90 src/rt/dg/rt_dg_fragment_ops.f90
git commit -m "refactor: remove dense fragment-only basis operator dependencies"
```

### Task 5: Remove persistent dense fragment-only allocation where no longer needed

**Files:**
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_fragment.f90`
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_fragment_soi.f90`
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_fragment_hamiltonian.f90`
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_fragment_hamiltonian_soi.f90`
- Test: source audit script from Task 1

**Step 1: Write the failing test/check**

Add a check that fails if steady-state fragment-only paths still allocate persistent dense `H_mat`, `S_mat`, `S_mat_prop`, or `momentum_mat` unnecessarily.

**Step 2: Run check to verify it fails**

Confirm current allocation sites are still reported.

**Step 3: Write minimal implementation**

Reduce allocation policy so these dense arrays are created only when:

- a temporary dense rebuild is explicitly requested
- a deferred mixed-basis or dense-only path still requires them

**Step 4: Run verification**

Run:

```bash
make -C /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/build -j4
```

Expected: build succeeds and phase 1 audit passes

**Step 5: Commit**

```bash
git add src/rt/dg/rt_dg_fragment.f90 src/rt/dg/rt_dg_fragment_soi.f90 src/rt/dg/rt_dg_fragment_hamiltonian.f90 src/rt/dg/rt_dg_fragment_hamiltonian_soi.f90
git commit -m "refactor: drop persistent dense fragment-only operator storage"
```

### Task 6: Validate runtime behavior and document remaining mixed-basis debt

**Files:**
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/docs/plans/2026-03-25-dg-fragment-dense-elimination.md`
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/docs/plans/2026-03-25-dg-fragment-dense-elimination-design.md`

**Step 1: Write the failing validation checklist**

Record the required runtime checks:

- non-SOI small DG-RT case
- SOI small DG-RT case
- `nproc_frag > 1`
- memory comparison before and after

**Step 2: Run available local verification**

Run the build and any available lightweight checks. If runtime MPI validation is not available locally, explicitly record that gap.

**Step 3: Write minimal documentation update**

Add a short section summarizing:

- what fragment-only dense allocations were removed
- what mixed-basis dense allocations remain
- what runtime validation is still required on production MPI systems

**Step 4: Re-run verification**

Re-run the build and confirm docs are accurate relative to current code.

**Step 5: Commit**

```bash
git add docs/plans/2026-03-25-dg-fragment-dense-elimination.md docs/plans/2026-03-25-dg-fragment-dense-elimination-design.md
git commit -m "docs: record dg dense elimination phase 1 validation"
```
