# DG F-P Orthonormalization Implementation Plan

> **Historical/removed:** This document describes an obsolete experimental DG route
> removed on 2026-07-31. It is retained only as an implementation record and is
> not executable guidance.

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Make the mixed fragment + plane-wave basis explicitly orthonormal between `F` and `P` at startup and after basis updates, so runtime propagation no longer depends on dense `FP/PF` overlap couplings.

**Architecture:** Introduce a reusable orthogonalization pipeline in the plane-wave module that transforms raw PW basis functions into a fragment-orthogonalized, internally orthonormal PW basis. Rebuild mixed overlap/Hamiltonian data in that transformed basis and update coefficient transfer paths to use the new representation consistently.

**Tech Stack:** Fortran 90, LAPACK (`dsyev`, `zheev`, linear solves), existing DG mixed-basis modules.

---

### Task 1: Add orthogonalized PW basis state to the DG fragment type

**Files:**
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_fragment_types.f90`
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_fragment.f90`

**Step 1: Write the failing implementation note**

Document in comments which arrays now represent raw PW basis versus orthogonalized PW basis, and remove the misleading `[UNUSED]` comment.

**Step 2: Add minimal fields**

Add persistent storage for orthogonalized real-space PW basis and, if needed, a compact transformation matrix that maps raw PW coefficients to orthogonalized PW coefficients.

**Step 3: Clean up deallocation paths**

Ensure all new arrays are freed in `cleanup_dg_fragment_rt`.

**Step 4: Build check**

Run:

```bash
make -C /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/build -j4
```

Expected: build succeeds.

### Task 2: Implement reusable F-P orthonormalization in the plane-wave module

**Files:**
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_plane_wave.f90`

**Step 1: Add a dedicated helper**

Create a helper that:

- computes `S_FP`,
- builds or loads `S_FF`,
- solves for fragment projection coefficients,
- subtracts fragment components from the raw PW basis,
- builds and orthonormalizes the residual PW basis,
- drops numerically dependent directions.

**Step 2: Reuse existing regularization policy**

Use the same eigenvalue-flooring/conditioning safeguards already used for overlap regularization when inverting or factorizing the fragment overlap.

**Step 3: Add diagnostics**

Print compact diagnostics for:

- number of kept PW modes,
- max `|S_FP|` after orthogonalization,
- deviation of `S_PP` from identity.

**Step 4: Build check**

Run the same `make` command and verify success.

### Task 3: Rebuild mixed startup operators in the orthogonalized basis

**Files:**
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_plane_wave.f90`
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_fragment.f90`

**Step 1: Update startup path**

Wire the startup path so `prepare_mixed_basis_startup` and `diagonalize_mixed_basis` operate on the orthogonalized PW basis, not the raw PW basis.

**Step 2: Ensure stored mixed operators are transformed**

Store `S_mat_frag_pw`, `H_mat_frag_pw`, and PW diagonal terms for the transformed basis.

**Step 3: Keep fallback behavior only where necessary**

Do not delete old code until the transformed path is verified, but make the orthogonalized path the default for mixed startup.

**Step 4: Build check**

Run the same `make` command and verify success.

### Task 4: Re-run orthogonalization after adaptive basis updates

**Files:**
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_fragment_basis_update_soi.f90`
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_plane_wave.f90`

**Step 1: Insert the rebuild point**

After fragment basis refresh in the mixed basis-update path, rerun PW orthonormalization and mixed operator rebuild before coefficient projection.

**Step 2: Preserve coefficient meaning**

If the orthogonalized PW basis changes dimension or ordering, update the coefficient transfer path so the old propagated PW coefficients are projected into the new PW basis consistently.

**Step 3: Add a basis-update diagnostic**

Print one compact message showing old/new PW count and post-update `max |S_FP|`.

**Step 4: Build check**

Run the same `make` command and verify success.

### Task 5: Simplify runtime overlap/momentum handling for orthogonalized mixed basis

**Files:**
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_fragment_ops.f90`
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_integrator_derivative.f90`

**Step 1: Gate the old `FP/PF` path**

When the orthogonalized mixed basis is active, skip dense `FP/PF` overlap contributions and treat any residual `S_mat_frag_pw` as diagnostics/fallback.

**Step 2: Align restricted overlap solve**

Make the local overlap solve assume orthogonalized `P`, so mixed runs do not need the ad hoc `S_FP rhs_P` correction on the main path.

**Step 3: Keep a safe fallback**

If orthogonalization is unavailable or diagnostics exceed tolerance, fall back to the current mixed path.

**Step 4: Build check**

Run the same `make` command and verify success.

### Task 6: Verify numerics and performance regressions

**Files:**
- Modify if needed: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_plane_wave.f90`

**Step 1: Verify compile**

Run:

```bash
make -C /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/build -j4
```

Expected: success.

**Step 2: Verify overlap diagnostics**

Run the existing mixed RT case and confirm logs show:

- post-orthogonalization `max |S_FP|` near tolerance,
- `before-overlap-solve` advances further or faster than before.

**Step 3: Verify no obvious regression**

Check that occupied norms and early-step energies remain reasonable compared with the current mixed-basis baseline.

**Step 4: Commit**

```bash
git add /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/docs/plans/2026-03-26-dg-fp-orthogonalization-design.md /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/docs/plans/2026-03-26-dg-fp-orthogonalization.md /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_fragment_types.f90 /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_fragment.f90 /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_plane_wave.f90 /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_fragment_basis_update_soi.f90 /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_fragment_ops.f90 /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_integrator_derivative.f90
git commit -m "feat: orthogonalize fragment and plane-wave bases"
```
