# Nonorthogonal Stress Metric Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Repair the primitive non-orthogonal stress formulation, starting from the kinetic sector, so analytic stress matches finite-difference checks without changing the existing orthogonal implementation.

**Architecture:** Keep `stencil%if_orthogonal=.true.` behavior fixed. Add non-orthogonal-only corrections sector-by-sector, beginning with `calc_stress_kin`, and validate each correction against primitive-cell finite differences before moving to the next sector.

**Tech Stack:** Fortran, `ctest` source checks, local Open MPI runtime checks, SALMON primitive-cell finite-difference smoke runs.

---

### Task 1: Freeze the current nonorthogonal evidence in a source-facing design note

**Files:**
- Create: `docs/plans/2026-04-03-nonorth-stress-metric-design.md`

**Step 1: Record the runtime evidence**

Document:
- primitive shear FD mismatch
- primitive isotropic FD mismatch
- `Ewald` matching while the other sectors do not

**Step 2: Record the implementation hypothesis**

Write down that the current non-Ewald sectors likely miss non-orthogonal metric derivatives already present in the Hamiltonian path.

**Step 3: Commit**

Skip commit unless explicitly requested.

### Task 2: Add a failing source test that protects the orthogonal branch

**Files:**
- Modify: `tests/source/CMakeLists.txt`
- Create: `tests/source/check_nonorth_stress_keeps_orthogonal_branch.sh`

**Step 1: Write the failing test**

Check that:
- `calc_stress_kin` still contains an orthogonal-only path or branch guard
- any new nonorthogonal logic is gated behind `if (.not. stencil%if_orthogonal)` or equivalent
- no existing orthogonal formulas are deleted

**Step 2: Run the test to verify it fails**

Run: `ctest --output-on-failure -R check_nonorth_stress_keeps_orthogonal_branch`

Expected: FAIL until the guard expectations are encoded in source.

**Step 3: Add the minimal source-test registration**

Update `tests/source/CMakeLists.txt` so the new check is discoverable.

### Task 3: Add a runtime helper for primitive nonorthogonal FD comparison

**Files:**
- Create: `tests/tools/run_primitive_nonorth_stress_fd.sh`
- Create: `tests/source/check_primitive_nonorth_stress_fd_tool.sh`
- Modify: `tests/source/CMakeLists.txt`

**Step 1: Write the failing source test**

Check that the new tool advertises:
- primitive-cell input generation
- isotropic scaling mode
- Cartesian shear mode
- `4 MPI(k) x 1 OMP`
- `stress_fd_detail='C'`

**Step 2: Create the helper script**

The tool should:
- write primitive Si inputs
- vary either lattice scale or a selected shear parameter
- run SALMON with the local MPI binary
- collect `*_info.data` and `*_stress_energy.data`

**Step 3: Run the source checks**

Run: `ctest --output-on-failure -R 'check_nonorth_stress_keeps_orthogonal_branch|check_primitive_nonorth_stress_fd_tool'`

Expected: PASS.

### Task 4: Add a failing kinetic-sector comparison check

**Files:**
- Create: `tests/tools/check_primitive_nonorth_kinetic_fd.py`

**Step 1: Write the failing comparator**

Implement a small analysis tool that:
- reads primitive FD outputs
- extracts analytic `P_kin` and `xy` kinetic components
- computes matching finite differences
- returns nonzero if the discrepancy exceeds tolerance

**Step 2: Run it on current primitive FD data**

Run it against the existing primitive control directories.

Expected: FAIL with a large kinetic mismatch.

### Task 5: Implement nonorthogonal-only kinetic correction

**Files:**
- Modify: `src/common/stress.f90`

**Step 1: Preserve the orthogonal path**

Do not alter the current orthogonal `calc_stress_kin` algebra.

**Step 2: Add the nonorthogonal correction path**

Introduce a branch for `stencil%if_orthogonal=.false.` that accounts for the non-orthogonal metric dependence missing from the current `-Re[w_a^* w_b]` accumulation.

**Step 3: Keep the implementation narrow**

Do not touch `Hartree`, `XC`, `Local`, or `Nonlocal` yet.

**Step 4: Re-run the failing kinetic comparator**

Run the comparator from Task 4.

Expected: the kinetic discrepancy shrinks materially, ideally to within `0.1 GPa`.

### Task 6: Verify orthogonal regressions

**Files:**
- Runtime only

**Step 1: Run source regressions**

Run: `ctest --output-on-failure -R 'check_builtin_r2scan|check_stress_|check_nonorth_stress_keeps_orthogonal_branch|check_primitive_nonorth_stress_fd_tool'`

Expected: PASS.

**Step 2: Rebuild local binaries**

Run:
- `cmake --build build-gcc15-env -j4`
- `cmake --build build-mpi-gcc15 -j4`

Expected: both succeed.

**Step 3: Re-run orthogonal smoke**

Run the existing orthogonal `PZ` and `r2scan` stress smoke scripts.

Expected: no change in orthogonal-sector results.

### Task 7: Re-run primitive isotropic and shear FD

**Files:**
- Runtime only

**Step 1: Primitive isotropic control**

Run the primitive isotropic FD helper and collect:
- `P_kin`
- `P_total`
- sector table at `output_level='C'`

**Step 2: Primitive shear control**

Run the primitive shear FD helper and collect:
- `xy` kinetic and total components
- corresponding FD derivatives

**Step 3: Decide whether the next sector is `Hartree` or `Local`**

Choose the next sector by largest remaining controlled discrepancy.

### Task 8: Summarize the corrected and remaining mismatch

**Files:**
- Optionally update: `docs/plans/2026-04-03-nonorth-stress-metric-design.md`

**Step 1: Report**

Include:
- kinetic before/after FD mismatch
- total before/after FD mismatch
- whether orthogonal stress stayed fixed
- which sector is next

**Step 2: Commit**

Recommended:

```bash
git add docs/plans/2026-04-03-nonorth-stress-metric-design.md \
  docs/plans/2026-04-03-nonorth-stress-metric.md \
  tests/source/CMakeLists.txt \
  tests/source/check_nonorth_stress_keeps_orthogonal_branch.sh \
  tests/source/check_primitive_nonorth_stress_fd_tool.sh \
  tests/tools/run_primitive_nonorth_stress_fd.sh \
  tests/tools/check_primitive_nonorth_kinetic_fd.py \
  src/common/stress.f90
git commit -m "fix: start nonorthogonal stress metric correction"
```
