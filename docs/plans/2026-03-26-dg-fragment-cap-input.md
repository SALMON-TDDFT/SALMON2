# DG Fragment Cap Input Implementation Plan

> **Historical/removed:** This document describes an obsolete experimental DG route
> removed on 2026-07-31. It is retained only as an implementation record and is
> not executable guidance.

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** DG-Fragment RT の `n_mat` cap を namelist から制御できるようにし、固定 cap と占有軌道数倍率 cap を選べるようにする。

**Architecture:** `&dg_fragment` に cap 制御入力を追加し、`rt_dg_fragment.f90` と `rt_dg_fragment_soi.f90` の既存 cap 適用点を namelist + 環境変数 override 対応に置き換える。`occ_multiple` では fragment ごとの occupied-state 数を見積もって `index_basis` を局所的に間引き、その後 global basis index を再圧縮する。

**Tech Stack:** Fortran, SALMON namelist input, DG-Fragment RT

---

### Task 1: Add input variables to global state

**Files:**
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/io/salmon_global.f90`
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/io/inputoutput.f90`

**Step 1: Add new global variables**

Add:

```fortran
character(16) :: dg_nmat_cap_mode
integer       :: dg_nmat_cap_fixed
real(8)       :: dg_nmat_cap_multiple
```

**Step 2: Register them in `&dg_fragment`**

Add the three variables to the `namelist/dc/`-adjacent DG fragment input block in `inputoutput.f90`.

**Step 3: Set defaults**

Use:

```fortran
dg_nmat_cap_mode = 'none'
dg_nmat_cap_fixed = 0
dg_nmat_cap_multiple = 0.0d0
```

**Step 4: Build**

Run: `make -C /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/build -j4`

Expected: build succeeds

**Step 5: Commit**

```bash
git add src/io/salmon_global.f90 src/io/inputoutput.f90
git commit -m "feat: add DG fragment cap input controls"
```

### Task 2: Refactor fixed cap logic behind one helper

**Files:**
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_fragment.f90`
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_fragment_soi.f90`

**Step 1: Extract cap source selection**

Implement logic that chooses:

```fortran
if (SALMON_DG_NMAT_CAP set) then
  effective_mode = 'fixed'
  effective_fixed = env value
else
  effective_mode = dg_nmat_cap_mode
  effective_fixed = dg_nmat_cap_fixed
end if
```

**Step 2: Preserve current fixed-cap behavior**

Keep existing semantics for:

```fortran
dg_frag%n_mat = min(dg_frag%n_mat, effective_fixed)
```

and the later `index_basis` cleanup and dense re-compression.

**Step 3: Add startup logging**

Print selected mode and cap on root rank.

**Step 4: Build**

Run: `make -C /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/build -j4`

Expected: build succeeds

**Step 5: Commit**

```bash
git add src/rt/dg/rt_dg_fragment.f90 src/rt/dg/rt_dg_fragment_soi.f90
git commit -m "refactor: route DG cap selection through input controls"
```

### Task 3: Implement occupied-state counting per fragment

**Files:**
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_fragment.f90`
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_fragment_soi.f90`

**Step 1: Add helper to estimate `nocc_frag`**

For each occupied state and spin:

```fortran
w(ifrag) = sum(abs(coef(local basis of ifrag, io, ispin))**2)
```

Assign that state to the fragment with maximum `w(ifrag)` and accumulate one count there.

**Step 2: Define occupied-state range**

Use the existing total-state/occupation conventions in DG startup. Limit counting to occupied states only.

**Step 3: Log summary**

On root rank, print min/max/avg `nocc_frag`.

**Step 4: Build**

Run: `make -C /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/build -j4`

Expected: build succeeds

**Step 5: Commit**

```bash
git add src/rt/dg/rt_dg_fragment.f90 src/rt/dg/rt_dg_fragment_soi.f90
git commit -m "feat: estimate occupied DG states per fragment"
```

### Task 4: Apply `occ_multiple` cap to fragment basis indices

**Files:**
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_fragment.f90`
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_fragment_soi.f90`

**Step 1: Compute fragment-local cap**

For each fragment and spin:

```fortran
cap_frag = min(n_basis(ifrag,ispin), int(floor(dg_nmat_cap_multiple * nocc_frag(ifrag,ispin))))
cap_frag = max(1, cap_frag)
```

**Step 2: Zero out excess `index_basis` entries**

Keep only the first `cap_frag` fragment-local basis entries and set the rest to zero.

**Step 3: Reuse existing global re-compression**

After zeroing, run the existing path that recomputes `n_mat` and compresses surviving global basis indices into a dense contiguous range.

**Step 4: Build**

Run: `make -C /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/build -j4`

Expected: build succeeds

**Step 5: Commit**

```bash
git add src/rt/dg/rt_dg_fragment.f90 src/rt/dg/rt_dg_fragment_soi.f90
git commit -m "feat: support occupied-multiple DG cap mode"
```

### Task 5: Add validation logging and smoke checks

**Files:**
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_fragment.f90`
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_fragment_soi.f90`
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/docs/plans/2026-03-26-dg-fragment-cap-input-design.md`

**Step 1: Add root-rank summary output**

Print:

- selected mode
- fixed cap or multiple
- min/max fragment-local cap
- final `n_mat_max`

**Step 2: Run fixed-mode check**

Run one small case with:

- `dg_nmat_cap_mode='fixed'`
- `dg_nmat_cap_fixed=<known value>`

Expected: startup `n_mat_max` matches the old `SALMON_DG_NMAT_CAP` behavior

**Step 3: Run occupied-multiple check**

Run one small case with:

- `dg_nmat_cap_mode='occ_multiple'`
- `dg_nmat_cap_multiple=2.0d0`

Expected: final `n_mat_max` is smaller than uncapped run and startup succeeds

**Step 4: Final build**

Run: `make -C /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/build -j4`

Expected: build succeeds

**Step 5: Commit**

```bash
git add src/rt/dg/rt_dg_fragment.f90 src/rt/dg/rt_dg_fragment_soi.f90 docs/plans/2026-03-26-dg-fragment-cap-input-design.md
git commit -m "chore: validate DG cap input modes"
```
