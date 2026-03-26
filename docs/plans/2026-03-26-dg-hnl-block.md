# DG Nonlocal H Block Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Replace the dense DG nonlocal pseudopotential fragment cache with direct block storage so RT `H0` and observables can apply the nonlocal operator through owner-local block lists.

**Architecture:** Introduce `H_nl_blocks` with the same fragment-pair topology as `H_mat_blocks`, build nonlocal PP contributions directly into those blocks inside the A-dependent refresh path, reduce them with the existing block reduction routine, and switch fragment-only RT and observables from dense nonlocal matmul to block apply. Dense nonlocal storage remains only as a fallback if a mixed/dense path still needs it.

**Tech Stack:** Fortran 90, SALMON DG RT block operator infrastructure, existing MPI block reductions.

---

### Task 1: Add nonlocal block storage to the DG fragment state

**Files:**
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_fragment_types.f90`
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_fragment.f90`
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_fragment_soi.f90`

**Step 1: Add nonlocal block fields**

Add to `s_dg_fragment_rt`:

- `type(matrix_block_info), allocatable :: H_nl_blocks(:)`
- `integer, allocatable :: H_nl_block_map(:,:)`
- `integer :: n_H_nl_blocks = 0`

If the existing `H_local_block_ids` can be reused because the topology is identical, do not add a second local list yet.

**Step 2: Extend finalization**

Deallocate the new nonlocal block arrays in both non-SOI and SOI finalizers.

**Step 3: Run a build-only sanity check**

Run:

```bash
make -C /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/build -j4
```

Expected: build succeeds.

### Task 2: Build nonlocal PP directly into blocks

**Files:**
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_fragment_ops.f90`

**Step 1: Add a direct block builder**

Create a new helper that mirrors the current nonlocal PP build loops but writes into:

- `dg_frag%H_nl_blocks(iblk)%val(io, jo, ispin)`

Use `find_matrix_block_runtime` or the stored block map to locate the block.

**Step 2: Add a strict topology assertion**

If a nonlocal contribution finds no valid block id, stop with a clear diagnostic. This protects against silently dropping couplings outside the assumed block topology.

**Step 3: Zero and reduce the nonlocal blocks**

Initialize `H_nl_blocks` before filling, then use `reduce_matrix_blocks(...)` on `H_nl_blocks` after local accumulation.

**Step 4: Keep dense nonlocal storage only as fallback if still needed**

If a dense fallback must remain temporarily, separate the block builder from the dense path so the fragment-only RT path no longer depends on `H_nl_cache`.

### Task 3: Integrate nonlocal block refresh into the existing A-dependent path

**Files:**
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_fragment_ops.f90`

**Step 1: Update `ensure_nonlocal_pp_matrix_A`**

Make this routine ensure that block storage exists and is refreshed whenever the nonlocal operator is refreshed.

Responsibilities:

- initialize `H_nl_blocks` and `H_nl_block_map` if needed
- rebuild `H_nl_blocks` when the current rebuild condition triggers
- avoid rebuilding when the cache policy says reuse is valid

**Step 2: Rebuild local block ids if topology changes**

If the nonlocal block topology is guaranteed to match `H_mat_blocks`, reuse the existing local block ids. Otherwise, add a nonlocal-specific local block list and rebuild it here.

### Task 4: Switch RT `H0` nonlocal apply from dense matmul to block apply

**Files:**
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_integrator_derivative.f90`
- Modify: `/Users/otobetoshito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_fragment_ops.f90`

**Step 1: Add block apply for `H_nl_blocks`**

Use `apply_matrix_blocks_batch` with the owner-local block list on `H_nl_blocks`.

**Step 2: Remove the dense nonlocal fragment matmul from the fragment-only block path**

Replace:

```fortran
matmul(dg_frag%H_nl_cache(1:n_frag, 1:n_frag, ispin), coef_all(1:n_frag, :))
```

with block apply accumulation.

**Step 3: Keep mixed/dense fallbacks safe**

If a dense path still exists for mixed basis or complex dense mode, guard it carefully and do not break existing fallback behavior.

### Task 5: Switch observables to the nonlocal block path

**Files:**
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_observables.f90`

**Step 1: Update fragment-only Hamiltonian application**

Where observables currently do:

```fortran
matmul(dg_frag%H_nl_cache(...), coef_frag_all(...))
```

replace it with block apply on `H_nl_blocks`.

**Step 2: Keep dense fallback only where required**

If observables still need dense nonlocal storage in a mixed/dense branch, keep that branch explicit.

### Task 6: Verify compile and call-path coverage

**Files:**
- Test: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/build`

**Step 1: Build**

Run:

```bash
make -C /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/build -j4
```

Expected: build succeeds.

**Step 2: Confirm block-path references**

Run:

```bash
rg -n "H_nl_blocks|H_nl_block_map|H_nl_cache|apply_matrix_blocks_batch" /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg
```

Expected:

- RT fragment-only path uses `H_nl_blocks`
- observables fragment-only path uses `H_nl_blocks`
- remaining `H_nl_cache` references are only fallback or diagnostics

**Step 3: Residual risk review**

Check carefully for:

- dropped nonlocal couplings because of block topology mismatch
- stale nonlocal blocks after `A(t)` changes
- dense fallback paths still allocating more memory than intended

**Step 4: Commit**

```bash
git add /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/docs/plans/2026-03-26-dg-hnl-block-design.md /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/docs/plans/2026-03-26-dg-hnl-block.md
git commit -m "docs: plan DG nonlocal H block storage"
```
