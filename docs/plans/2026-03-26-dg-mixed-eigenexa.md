# DG Mixed Basis EigenExa Implementation Plan

> **Historical/removed:** This document describes an obsolete experimental DG route
> removed on 2026-07-31. It is retained only as an implementation record and is
> not executable guidance.

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Propagate `n_pw > 0` DG-RT cases in an exact orthonormal mixed fragment+PW basis so the RT hot path no longer carries raw `FP/PF` overlap and Hamiltonian coupling every step.

**Architecture:** Keep raw fragment/PW quantities only as rebuild views. Introduce a canonical propagation basis `X_mix` with coefficients `coef_mix`, regenerate raw coefficients only for density and potential rebuild, and transform rebuilt raw operators back into mixed space before derivative evaluation. Basis updates must rebuild the mixed basis and reproject coefficients into the new orthonormal subspace.

**Tech Stack:** Fortran 2008, existing SALMON DG RT modules, EigenExa `eigen_sx`, LAPACK/BLAS, current mixed fragment/PW assembly routines.

---

### Task 1: Finalize mixed-basis state in `s_dg_fragment_rt`

**Files:**
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_fragment_types.f90`
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_fragment.f90`
- Test: `make -C /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/build -j4`

**Step 1: Write the failing test**

Define the expected state by inspection:
- `mixed_basis_ready`
- `mixed_basis_dim`
- `mixed_transform`
- `coef_mix`
must exist, be invalidated on teardown, and be sized by retained mixed dimension per spin.

**Step 2: Run test to verify it fails**

Run:
```bash
make -C /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/build -j4
```
Expected: either missing symbols for `coef_mix` or compile success before the new data model is added.

**Step 3: Write minimal implementation**

Add to `s_dg_fragment_rt`:
- `complex(8), allocatable :: coef_mix(:,:,:)`
- retained mixed dimension metadata
- cleanup/invalidation in finalization paths

**Step 4: Run test to verify it passes**

Run:
```bash
make -C /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/build -j4
```
Expected: PASS

**Step 5: Commit**

```bash
git -C /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2 add src/rt/dg/rt_dg_fragment_types.f90 src/rt/dg/rt_dg_fragment.f90
git -C /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2 commit -m "Store mixed-basis propagation state"
```

### Task 2: Make `diagonalize_mixed_basis` return a propagation basis

**Files:**
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_plane_wave.f90`
- Test: `make -C /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/build -j4`

**Step 1: Write the failing test**

Define the required postconditions:
- `mixed_transform` stores retained raw-to-mixed eigenvectors
- `coef_mix` is populated from the current raw coefficients
- `mixed_basis_ready=.true.` only after both transform and coefficients are valid

**Step 2: Run test to verify it fails**

Run:
```bash
make -C /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/build -j4
```
Expected: current code compiles but does not yet populate `coef_mix`.

**Step 3: Write minimal implementation**

In `diagonalize_mixed_basis(...)`:
- keep `mixed_transform`
- compute `coef_mix = T^H * c_raw`
- keep raw `coef` / `coef_pw` as rebuild views
- set retained dimensions per spin from overlap rank truncation

**Step 4: Run test to verify it passes**

Run:
```bash
make -C /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/build -j4
```
Expected: PASS

**Step 5: Commit**

```bash
git -C /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2 add src/rt/dg/rt_dg_plane_wave.f90
git -C /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2 commit -m "Populate mixed propagation coefficients"
```

### Task 3: Add raw/mixed coefficient projection helpers

**Files:**
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_fragment_ops.f90`
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_fragment_types.f90`
- Test: `make -C /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/build -j4`

**Step 1: Write the failing test**

Define two helper behaviors:
- expand mixed coefficients to raw `(coef, coef_pw)`
- compress raw `(coef, coef_pw)` back to `coef_mix`

**Step 2: Run test to verify it fails**

Run:
```bash
make -C /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/build -j4
```
Expected: helper routines not found yet.

**Step 3: Write minimal implementation**

Add routines such as:
- `expand_mixed_coef_to_raw(...)`
- `project_raw_coef_to_mixed(...)`

These should use `mixed_transform` and support all active spins/states.

**Step 4: Run test to verify it passes**

Run:
```bash
make -C /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/build -j4
```
Expected: PASS

**Step 5: Commit**

```bash
git -C /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2 add src/rt/dg/rt_dg_fragment_ops.f90 src/rt/dg/rt_dg_fragment_types.f90
git -C /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2 commit -m "Add mixed coefficient projection helpers"
```

### Task 4: Reproject after basis update

**Files:**
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_fragment_basis_update.f90`
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_plane_wave.f90`
- Test: `make -C /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/build -j4`

**Step 1: Write the failing test**

Define the required sequence:
- save old mixed coefficients
- rebuild raw fragment basis
- re-diagonalize mixed basis
- project old state into the new mixed basis

**Step 2: Run test to verify it fails**

Run:
```bash
make -C /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/build -j4
```
Expected: current code invalidates mixed basis but does not reproject.

**Step 3: Write minimal implementation**

Implement:
- old mixed-to-raw expansion before invalidation
- new raw-to-mixed projection after rebuild
- logging or guards when retained mixed dimension changes

**Step 4: Run test to verify it passes**

Run:
```bash
make -C /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/build -j4
```
Expected: PASS

**Step 5: Commit**

```bash
git -C /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2 add src/rt/dg/rt_dg_fragment_basis_update.f90 src/rt/dg/rt_dg_plane_wave.f90
git -C /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2 commit -m "Reproject mixed coefficients after basis updates"
```

### Task 5: Build mixed-space operators from raw rebuilds

**Files:**
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_integrator_stage_update.f90`
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_fragment_ops.f90`
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_fragment_types.f90`
- Test: `make -C /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/build -j4`

**Step 1: Write the failing test**

Define expected behavior:
- stage update rebuilds raw `H` and `M`
- mixed-ready path forms transformed operators `T^H H_raw T`, `T^H M_raw T`
- transformed operators are cached for derivative evaluation

**Step 2: Run test to verify it fails**

Run:
```bash
make -C /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/build -j4
```
Expected: no transformed mixed-space operator storage exists yet.

**Step 3: Write minimal implementation**

Add:
- mixed-space operator cache fields
- helper routines to form transformed operators
- stage-update hook to refresh them after each self-consistent rebuild

**Step 4: Run test to verify it passes**

Run:
```bash
make -C /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/build -j4
```
Expected: PASS

**Step 5: Commit**

```bash
git -C /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2 add src/rt/dg/rt_dg_integrator_stage_update.f90 src/rt/dg/rt_dg_fragment_ops.f90 src/rt/dg/rt_dg_fragment_types.f90
git -C /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2 commit -m "Cache transformed mixed-space operators"
```

### Task 6: Switch derivative evaluation to the orthonormal mixed basis

**Files:**
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_integrator_derivative.f90`
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_fragment_ops.f90`
- Test: `make -C /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/build -j4`

**Step 1: Write the failing test**

Define expected behavior:
- when `mixed_basis_ready` is true, derivative uses `coef_mix`
- overlap solve is bypassed because propagation overlap is identity
- raw `FP/PF` coupling path remains fallback only

**Step 2: Run test to verify it fails**

Run:
```bash
make -C /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/build -j4
```
Expected: current derivative still uses raw `coef/cof_pw` and overlap solves.

**Step 3: Write minimal implementation**

Implement mixed-ready branch in `calculate_time_derivative(...)`:
- read `coef_mix`
- apply transformed mixed-space `H` and `M`
- skip `solve_overlap_operator_batch*`
- write `dcoef_dt_mix`

Keep existing raw mixed branch unchanged as fallback.

**Step 4: Run test to verify it passes**

Run:
```bash
make -C /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/build -j4
```
Expected: PASS

**Step 5: Commit**

```bash
git -C /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2 add src/rt/dg/rt_dg_integrator_derivative.f90 src/rt/dg/rt_dg_fragment_ops.f90
git -C /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2 commit -m "Propagate mixed RT states in orthonormal basis"
```

### Task 7: Verify startup, basis update, and mixed RT behavior

**Files:**
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/docs/plans/2026-03-26-dg-mixed-eigenexa-design.md`
- Test: `make -C /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/build -j4`
- Test: project-specific small `n_pw > 0` smoke case

**Step 1: Write the failing test**

Define verification checklist:
- `T^H S T ≈ I`
- mixed basis dimension is stable or clearly logged after updates
- mixed RT path runs without raw overlap solve
- observables remain finite and electron count drift is acceptable

**Step 2: Run test to verify it fails**

Run:
```bash
make -C /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/build -j4
```
Expected: PASS for build, pending runtime verification.

**Step 3: Write minimal implementation**

Add minimal orthogonality and dimension diagnostics guarded behind a debug flag and update the design doc with final verification notes.

**Step 4: Run test to verify it passes**

Run:
```bash
make -C /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/build -j4
```
Expected: PASS

**Step 5: Commit**

```bash
git -C /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2 add docs/plans/2026-03-26-dg-mixed-eigenexa-design.md
git -C /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2 commit -m "Document mixed-basis RT verification"
```
