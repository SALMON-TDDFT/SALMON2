# Primitive-Cell r2SCAN Nonorthogonal Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Enable builtin `r2scan` GS and stress for the 2-atom primitive Si cell on a non-orthogonal lattice, then verify `output_level='C'` stress tensor sectors against the existing 8-atom conventional-cell route.

**Architecture:** Extend the existing orthogonal `tau`-operator payload path with a non-orthogonal branch that mirrors the geometry of `zstencil_nonorthogonal`. Keep the current orthogonal path untouched, preserve the existing `r2scan` stress assembly, and validate the primitive-cell implementation by comparing transformed stress tensors sector-by-sector in Cartesian coordinates.

**Tech Stack:** Fortran, CMake/ctest source tests, local Open MPI runtime checks, SALMON samples and runtime smoke scripts.

---

### Task 1: Add a failing source test for the new non-orthogonal `tau` branch

**Files:**
- Modify: `tests/source/CMakeLists.txt`
- Create: `tests/source/check_builtin_r2scan_nonorth_tau_operator.sh`

**Step 1: Write the failing test**

Check for:
- removal of the unconditional non-orthogonal stop in `add_xc_tau_operator`
- presence of a non-orthogonal branch in the XC `tau` operator path
- continued presence of the orthogonal branch

**Step 2: Run test to verify it fails**

Run: `ctest --output-on-failure -R check_builtin_r2scan_nonorth_tau_operator`

Expected: FAIL because the orthogonal-only guard is still present.

**Step 3: Commit**

Skip commit unless explicitly requested.

### Task 2: Add a primitive-cell runtime driver input path

**Files:**
- Modify: `tests/tools/run_exercise4_xc_mpi_smoke.sh`
- Create: `tests/tools/run_primitive_si_r2scan_stress_smoke.sh`

**Step 1: Prepare the runtime target shape**

Add a dedicated primitive-cell smoke driver that:
- writes a 2-atom primitive Si input
- uses `al_vec1:3`
- keeps `nproc_rgrid=1,1,1`
- runs with `4 MPI(k) x 1 OMP`
- enables `yn_out_stress='y'`, `yn_out_stress_decomp='y'`, `stress_fd_detail='C'`

**Step 2: Keep it isolated**

Do not change the existing conventional-cell smoke script semantics.

**Step 3: Run the driver in dry form if possible**

Run: `bash tests/tools/run_primitive_si_r2scan_stress_smoke.sh --help`

Expected: the script is parseable and documents the intended runtime path.

### Task 3: Implement a non-orthogonal `tau` operator branch

**Files:**
- Modify: `src/common/hamiltonian.f90`

**Step 1: Preserve the current orthogonal code path**

Keep the existing orthogonal branch unchanged.

**Step 2: Add the non-orthogonal branch**

Implement a variable-coefficient form consistent with `zstencil_nonorthogonal`, using:
- direct-axis terms aligned with `F(1:3)`
- mixed `yz/zx/xy` terms aligned with `F(4:6)`
- the same `mg%idx/idy/idz` periodic indexing model
- local averaging of `vtau` on links used by each difference term

**Step 3: Keep current guards that still apply**

Retain:
- `OpenACC` stop
- `single-scale Maxwell-TDDFT` stop
- shadow-value bounds checks

Only relax the orthogonal-only stop.

**Step 4: Run the focused failing test**

Run: `ctest --output-on-failure -R check_builtin_r2scan_nonorth_tau_operator`

Expected: PASS.

### Task 4: Add a source test for primitive-cell smoke coverage

**Files:**
- Modify: `tests/source/CMakeLists.txt`
- Create: `tests/source/check_builtin_r2scan_primitive_smoke.sh`

**Step 1: Write the failing test**

Check that the primitive smoke driver:
- writes `al_vec1`, `al_vec2`, `al_vec3`
- uses `xc='r2scan'`
- requests stress output at detail level `C`
- keeps `nproc_rgrid(1:3)=1`

**Step 2: Run it to verify it fails**

Run: `ctest --output-on-failure -R check_builtin_r2scan_primitive_smoke`

Expected: FAIL until the new script exists.

**Step 3: Re-run after the script is added**

Run: `ctest --output-on-failure -R 'check_builtin_r2scan_nonorth_tau_operator|check_builtin_r2scan_primitive_smoke'`

Expected: PASS.

### Task 5: Verify source regressions stay green

**Files:**
- Runtime only

**Step 1: Run focused source regressions**

Run: `ctest --output-on-failure -R 'check_builtin_r2scan|check_stress_'`

Expected: PASS.

**Step 2: Rebuild the local binaries**

Run:
- `cmake --build build-gcc15-env -j4`
- `cmake --build build-mpi-gcc15 -j4`

Expected: both builds succeed.

### Task 6: Run primitive-cell GS and stress smoke

**Files:**
- Runtime only

**Step 1: Run primitive `r2scan` GS**

Run the primitive driver with:
- `4 MPI(k) x 1 OMP`
- non-orthogonal primitive Si
- `output_level='C'`

Expected: GS converges without hitting the orthogonal `tau` stop.

**Step 2: Run primitive `r2scan` stress**

Use the same driver with stress output enabled.

Expected: stress completes and writes the full sector decomposition.

**Step 3: Run primitive `PZ` stress as a control**

Expected: a non-orthogonal primitive baseline exists for comparison.

### Task 7: Compare primitive vs conventional stress tensors

**Files:**
- Runtime only

**Step 1: Extract sector tensors**

From `output_level='C'`, extract:
- `Kinetic`
- `Hartree`
- `XC`
- `Local`
- `Nonlocal`
- `Ewald`
- `Total`

for both:
- primitive `2-atom` Si
- conventional `8-atom` Si

**Step 2: Transform tensors into the same Cartesian basis**

Do not compare raw lattice-frame numbers if the primitive input is expressed in a different basis.

**Step 3: Check tolerances**

Target:
- each sector within about `0.1 GPa`
- if not, identify the first sector that diverges materially

### Task 8: Summarize and decide next scope

**Files:**
- Optionally update: `docs/plans/2026-04-03-primitive-r2scan-nonorthogonal-design.md`

**Step 1: Report**

Include:
- whether primitive `r2scan GS` now runs
- whether primitive `r2scan stress` now runs
- which stress sectors match conventional-cell results
- any residual mismatch and its likely origin

**Step 2: Commit**

Recommended:

```bash
git add src/common/hamiltonian.f90 tests/source/CMakeLists.txt \
  tests/source/check_builtin_r2scan_nonorth_tau_operator.sh \
  tests/source/check_builtin_r2scan_primitive_smoke.sh \
  tests/tools/run_primitive_si_r2scan_stress_smoke.sh \
  docs/plans/2026-04-03-primitive-r2scan-nonorthogonal-design.md \
  docs/plans/2026-04-03-primitive-r2scan-nonorthogonal.md
git commit -m "feat: add primitive-cell r2scan nonorthogonal tau support"
```
