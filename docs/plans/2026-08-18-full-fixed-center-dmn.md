# Full Fixed-Center DMN Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Make Wannier90's one-pass site-symmetry average an exact Reynolds projection by publishing the complete closed fixed-center group.

**Architecture:** Stream one representation for each canonical fixed-center operation from `fixed_center_symmetry_map`, append it to DMN, and validate the same complete operation list at publication. Preserve the existing row-owned representation assembly and one-operation memory bound.

**Tech Stack:** Fortran 2008, MPI, Wannier90 library DMN interface, Python/Fortran regression tests.

---

### Task 1: Add RED coverage for complete-group DMN publication

**Files:**
- Modify: `tests/dg/check_dg_overlapping_wannier_route.py`
- Modify: `tests/dg/check_sawf_dmn_format.py`

**Step 1: Write the failing route assertions**

Require the production DMN count and loop bound to use `size(fixed_center_operations)`, require representation assembly from `fixed_center_symmetry_map`, reject `global_affine_generators` within the DMN transaction, and reject `require_closed_group=.false.`.

**Step 2: Add the numerical Reynolds regression**

Use the exact two-dimensional representation `D=diag(1,exp(2*pi*i/3))`. Verify that averaging an off-diagonal matrix over `{I,D}` is not invariant and is not idempotent, while averaging over `{I,D,D^2}` is invariant to numerical tolerance.

**Step 3: Run tests to verify RED**

Run:

```bash
python3 tests/dg/check_dg_overlapping_wannier_route.py
python3 tests/dg/check_sawf_dmn_format.py
```

Expected: the route checker fails because production still publishes identity plus affine generators. The standalone Reynolds assertions pass and document the mathematical failure.

### Task 2: Stream the complete fixed-center group

**Files:**
- Modify: `src/gs/main_dft.f90:1550-1630`

**Step 1: Set the DMN operation count**

Pass `fixed_center_group_order` to `begin_sawf_dmn`.

**Step 2: Assemble the matching operation**

Loop `fixed_center_operation=1,fixed_center_group_order` and call `assemble_dg_distributed_basis_symmetry_overlap_rows` with `fixed_center_symmetry_map(:,fixed_center_operation:fixed_center_operation)`.

**Step 3: Preserve the exact identity convention**

Pass `fixed_center_operation==fixed_center_identity_operation` to `append_sawf_dmn_operation`.

**Step 4: Publish the matching closed group**

Pass all `fixed_center_operations` to `finish_sawf_dmn` and omit the optional false closure flag.

**Step 5: Run tests to verify GREEN**

Run:

```bash
python3 tests/dg/check_dg_overlapping_wannier_route.py
python3 tests/dg/check_sawf_dmn_format.py
```

Expected: PASS.

### Task 3: Focused MPI and build verification

**Files:**
- No source changes expected.

**Step 1: Run W90 MPI coverage**

```bash
python3 tests/dg/run_dg_overlapping_wannier_w90_mpi.py
```

Expected: PASS for MPI 1/2/4/8.

**Step 2: Run SAWF/DMN checks**

Run the focused SAWF DMN format, fragment symmetry, and route checks used by the worktree.

Expected: PASS.

**Step 3: Build production**

Rebuild the existing production build directory without changing MPI rank or OpenMP policy.

Expected: successful link.

**Step 4: Check patch hygiene**

```bash
git diff --check
git status --short
```

Expected: no whitespace errors; unrelated pre-existing changes remain untouched.

### Task 4: Si64 production confirmation

**Files:**
- Reuse: `/Users/otobetoshihito/SALMON-dev/verification/20260818-si64-reviewed-complex-gauge-mpi8-omp1/inputfile`

**Step 1: Start a fresh run directory**

Copy only the required input and pseudopotential files, use MPI 8 and OMP 1, and record peak RSS per rank.

**Step 2: Confirm DMN diagnostics**

Verify the DMN header reports 12 symmetry operations and the pre-Wannier closed-group transaction succeeds.

**Step 3: Confirm the post-Wannier gate**

Allow Wannier90 to finish. Verify the post-Wannier generator covariance gate passes. If it fails, preserve the measured defect and artifacts; do not relax tolerance or remove the check.

