# r2SCAN vsigma Diagnostic Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add a diagnostic-only `vsigma` payload and stress decomposition path for builtin `r2scan`, so the current `rdedd` gradient stress can be compared against a no-division `vsigma` formulation on the same runtime outputs.

**Architecture:** Extend the builtin `r2scan` payload with one scalar field `vsigma`, where `vsigma = -vgrad / (2 |grad rho|)` when `|grad rho|` is finite and zero otherwise. Keep the production XC stress unchanged, but store an alternate tensor `stress_xc_grad_vsigma` and print it at `output_level='C'` next to the current `XC-grad`. This keeps the comparison local, reversible, and tied to the same SCF state and discretization.

**Tech Stack:** Fortran, CMake/ctest source tests, local Open MPI runtime checks.

---

### Task 1: Add a failing source test for the new diagnostic split

**Files:**
- Modify: `tests/source/CMakeLists.txt`
- Create: `tests/source/check_builtin_r2scan_stress_vsigma_diag.sh`

**Step 1: Write the failing test**

Check for:
- `vsigma` in `s_xc_operator_payload`
- `stress_xc_grad_vsigma` in `s_dft_system`
- `XC-grad-vsigma` output row in `write.f90`

**Step 2: Run test to verify it fails**

Run: `ctest --output-on-failure -R check_builtin_r2scan_stress_vsigma_diag`

Expected: FAIL because none of the new symbols exist yet.

**Step 3: Commit**

Skip commit unless explicitly requested.

### Task 2: Add diagnostic storage for `vsigma`

**Files:**
- Modify: `src/common/structures.f90`
- Modify: `src/common/stress.f90`

**Step 1: Add minimal fields**

Add:
- `vsigma_has_shadow_values`
- `type(s_scalar) :: vsigma`
- `real(8) :: stress_xc_grad_vsigma(3,3)`

**Step 2: Zero-initialize and symmetrize**

Reset and symmetrize `stress_xc_grad_vsigma` in `calc_stress`.

**Step 3: Run focused source tests**

Run: `ctest --output-on-failure -R 'check_builtin_r2scan_stress_detail_split|check_builtin_r2scan_stress_vsigma_diag'`

Expected: still FAIL until the write path is added.

### Task 3: Populate `vsigma` in builtin `r2scan`

**Files:**
- Modify: `src/xc/builtin_r2scan.f90`
- Modify: `src/xc/salmon_xc.f90`

**Step 1: Derive `vsigma` from existing `vgrad`**

Use:
- `grad = |grad rho|`
- `vsigma = -vgrad / (2 * grad)` for finite `grad`
- `0` otherwise

This matches the current `rdedd = -vgrad * grad(rho) / |grad rho|`.

**Step 2: Thread through payload**

Allocate/copy/update overlap for `payload%vsigma` just like the existing scalar payload path.

**Step 3: Run source tests**

Run: `ctest --output-on-failure -R 'check_builtin_r2scan|check_stress_'`

Expected: PASS.

### Task 4: Assemble and print alternate gradient stress

**Files:**
- Modify: `src/common/stress.f90`
- Modify: `src/io/write.f90`

**Step 1: Build alternate tensor**

Using existing `grho_local`, assemble:
- `stress_xc_grad_vsigma(a,b) += 2 * vsigma * grho_local(a) * grho_local(b)`

Keep the existing `stress_xc` production path unchanged.

**Step 2: Print in detail output**

For `r2scan` detail output, print:
- `XC-grad-rdedd`
- `XC-grad-vsigma`

Keep `XC-local`, `XC-grad`, and `XC-tau` available.

**Step 3: Run source tests**

Run: `ctest --output-on-failure -R 'check_builtin_r2scan|check_stress_'`

Expected: PASS.

### Task 5: Rebuild and compare `18/22/26`

**Files:**
- Runtime only

**Step 1: Rebuild**

Run:
- `cmake --build build-mpi-gcc15 -j4`
- `cmake --build build-gcc15-env -j4`

**Step 2: Re-run runtime checks**

Run:
- `OMP_NUM_THREADS=1 mpiexec -n 4 ../../salmon < inputfile > outputfile 2>&1`

for:
- `build-mpi-gcc15/runtime-checks/exercise4_r2scan_stress_18_k4mpi_C`
- `build-mpi-gcc15/runtime-checks/exercise4_r2scan_stress_22_k4mpi_C`
- `build-mpi-gcc15/runtime-checks/exercise4_r2scan_stress_26_k4mpi_C`

**Step 3: Compare norms**

Extract and compare `dev norm` for:
- `XC-grad-rdedd`
- `XC-grad-vsigma`

Expected:
- If `vsigma` is materially smaller, the `|grad rho|` division is implicated.
- If both are the same, the issue is deeper than representation choice.

### Task 6: Summarize the outcome

**Files:**
- None required

**Step 1: Report**

Include:
- which representation is smaller
- whether the improvement is systematic across `18/22/26`
- whether the next step should be `h` or a deeper formula/sign review

**Step 2: Commit**

Skip commit unless explicitly requested.
