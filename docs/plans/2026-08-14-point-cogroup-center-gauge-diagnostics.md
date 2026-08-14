# Point-Cogroup Center-Gauge Diagnostics Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Measure whether the first affine center-orbit failure is caused by a non-monomial retained point-cogroup gauge or by center data inconsistent with an otherwise center-permuting representation.

**Architecture:** Extend the existing center validator to return the first failing operation, then diagnose only that operation with the existing row-owned symmetry-overlap assembler.  Compare the retained representation weights against tolerance-compatible affine center targets and publish collective monomial, center-block leakage, unitarity, and workspace receipts without weakening the center gate.

**Tech Stack:** Fortran 2008, MPI, BLAS-backed distributed overlap assembly, Python MPI runner and route checker, CMake.

---

### Task 1: Add the center-gauge diagnostic RED fixtures

**Files:**
- Modify: `tests/dg/test_dg_overlapping_wannier_construction_mpi.f90`

**Step 1: Import the new public diagnostic API and declare fixture storage**

Add `diagnose_dg_point_center_gauge` to the module `use` list.  Add small row-owned basis/map/center fixtures, scalar defect receipts, an `integer(int64)` workspace receipt, and a `failed_operation` integer.

**Step 2: Require the center validator to return the first failing operation**

Extend the existing broken two-center call with `failed_operation=failed_operation` and assert `failed_operation==2`.

**Step 3: Add three diagnostic fixtures**

Use a four-orbital distributed orthonormal basis with globally unique spatial rows:

1. An exact center permutation.  Require monomial defect, center-block leakage, and unitarity defect below `1d-12`.
2. A unitary rotation mixing two different mapped-center blocks.  Require a small unitarity defect but center-block leakage greater than `1d-3`.
3. A unitary rotation confined to two orbitals sharing the same center.  Require nonzero elementwise monomial defect and center-block leakage below `1d-12`.

Also pass duplicate row ownership and rank-disagreeing operation metadata and require collective rejection.

**Step 4: Emit rank-independent receipt evidence**

On rank zero print a stable line containing all diagnostic receipts.  Extend `tests/dg/run_dg_overlapping_wannier_construction_mpi.py` to compare the evidence across MPI 1/2/4/8.

**Step 5: Run the RED**

Run:

```bash
python3 tests/dg/run_dg_overlapping_wannier_construction_mpi.py
```

Expected: compilation fails because `diagnose_dg_point_center_gauge` and the optional `failed_operation` output do not exist.

**Step 6: Commit the RED**

```bash
git add tests/dg/test_dg_overlapping_wannier_construction_mpi.f90 \
  tests/dg/run_dg_overlapping_wannier_construction_mpi.py
git commit -m "test(dg): cover point center gauge diagnostics"
```

### Task 2: Return the first failing center operation

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_construction.f90:2347-2419`
- Test: `tests/dg/test_dg_overlapping_wannier_construction_mpi.f90`

**Step 1: Add the optional output**

Add `integer,intent(out),optional::failed_operation`.  Initialize it to zero before contract validation and set it to `operation` immediately before returning the existing mismatch message.  Do not change matching or tolerance behavior.

**Step 2: Run the focused fixture**

Run:

```bash
python3 tests/dg/run_dg_overlapping_wannier_construction_mpi.py
```

Expected: the existing validator assertions pass; compilation still fails or the new diagnostic assertions fail because the diagnostic primitive is absent.

**Step 3: Commit the validator extension**

```bash
git add src/gs/dc/dg_overlapping_wannier_construction.f90 \
  tests/dg/test_dg_overlapping_wannier_construction_mpi.f90
git commit -m "diag(dg): return failing center operation"
```

### Task 3: Implement the row-owned center-gauge diagnostic

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_construction.f90`
- Test: `tests/dg/test_dg_overlapping_wannier_construction_mpi.f90`

**Step 1: Define the collective API**

Add a public subroutine with inputs:

```fortran
subroutine diagnose_dg_point_center_gauge(comm,local_basis,weights,point_map,&
    integer_rotation,fractional_translation,centers,tolerance,&
    monomial_defect,center_block_leakage,unitarity_defect,&
    workspace_peak_bytes,ok,message)
```

`local_basis` follows the existing `(nstate,nlocal_spatial)` convention and `point_map` is the one-operation row-owned spatial permutation.

**Step 2: Validate before any indexed or shape-dependent work**

Collectively agree `nstate`, local spatial extent, tolerance, rotation, and translation.  Require finite basis/weights/centers/translation, a valid rotation shape, and an exactly-once global spatial permutation.  Use checked `int64` products/additions for all allocation and byte receipts.  Use `allocate(stat=...)` followed by collective allocation consensus and a single cleanup path.

**Step 3: Assemble one retained representation**

Call `assemble_dg_distributed_basis_symmetry_overlap_rows` with `reshape(point_map,[nlocal,1])`.  Preserve its row IDs and one row-owned representation slice.  Measure `D^H D-I` collectively and reject a nonunitary result using the supplied tolerance.

**Step 4: Measure center compatibility and leakage**

For each source orbital, compute

```fortran
mapped_center = modulo(R * centers(:,source) + tau, 1d0)
```

A target is compatible when the maximum periodic component residual is at most `tolerance`.  Accumulate `abs(D(target,source))**2` outside compatible targets, reduce it with `MPI_SUM`, and publish the maximum source leakage.

**Step 5: Measure elementwise monomiality independently**

For each source column and target row, obtain the maximum `abs(D)**2`.  Publish

```text
max(1-min_column_max, 1-min_row_max)
```

This remains nonzero for allowed rotations inside repeated-center blocks, while center-block leakage remains zero.

**Step 6: Publish a checked collective workspace receipt**

Account for the overlap assembler receipt plus compatibility, column maxima, leakage, row-ID, and temporary arrays.  Reduce the peak with `MPI_MAX` so every rank returns the same value.

**Step 7: Run GREEN**

Run:

```bash
python3 tests/dg/run_dg_overlapping_wannier_construction_mpi.py
```

Expected: PASS on MPI 1/2/4/8 with identical diagnostic evidence.

**Step 8: Commit the primitive**

```bash
git add src/gs/dc/dg_overlapping_wannier_construction.f90 \
  tests/dg/test_dg_overlapping_wannier_construction_mpi.f90 \
  tests/dg/run_dg_overlapping_wannier_construction_mpi.py
git commit -m "diag(dg): measure point center gauge leakage"
```

### Task 4: Connect the diagnostic to the production failure path

**Files:**
- Modify: `src/gs/main_dft.f90:1691-1700`
- Modify: `tests/dg/check_dg_overlapping_wannier_route.py`

**Step 1: Add the route RED**

Require the production source to capture `failed_operation`, invoke the diagnostic with exactly that column of `global_symmetry_map`, `global_point_integer_rotations`, and `global_point_fractional_translations`, and preserve the original error stop.

**Step 2: Run the route RED**

Run:

```bash
python3 tests/dg/check_dg_overlapping_wannier_route.py
```

Expected: FAIL because the production diagnostic call is absent.

**Step 3: Add failure-path reporting**

When `verify_dg_wannier_center_affine_orbits` rejects, call the new primitive before the existing error stop.  On rank zero print:

```text
[OW-GS-DIAGNOSTIC] point_center_gauge operation=...
center_bottleneck_residual=... monomial_defect=...
center_block_leakage=... representation_unitarity_defect=...
workspace_peak_bytes=...
```

Keep the original center mismatch message authoritative.  If the diagnostic itself fails, print its message and still stop at the center gate.

**Step 4: Run the route checker**

Run:

```bash
python3 tests/dg/check_dg_overlapping_wannier_route.py
```

Expected: PASS.

**Step 5: Commit production integration**

```bash
git add src/gs/main_dft.f90 tests/dg/check_dg_overlapping_wannier_route.py
git commit -m "diag(dg): report production point center gauge leakage"
```

### Task 5: Verify the focused implementation

**Files:**
- Verify only; do not stage unrelated user files.

**Step 1: Run focused MPI tests**

```bash
python3 tests/dg/run_dg_overlapping_wannier_construction_mpi.py
```

Expected: PASS on MPI 1/2/4/8.

**Step 2: Run the route checker**

```bash
python3 tests/dg/check_dg_overlapping_wannier_route.py
```

Expected: PASS.

**Step 3: Build the release executable**

```bash
cmake --build /Users/otobetoshihito/SALMON-dev/verification/20260814-core-ownership-fix-build -j 1
```

Expected: PASS.

**Step 4: Check formatting and scope**

```bash
git diff --check
git status --short
```

Expected: no whitespace errors.  Preserve the four existing user-modified Si64 test/input files and the two untracked W90 fixture outputs without staging them.

**Step 5: Commit any verification-only corrections**

If corrections were required, commit only the diagnostic implementation and tests.  Otherwise do not create an empty commit.

### Task 6: Run Si64 only after focused verification

**Files:**
- Runtime output only under `/Users/otobetoshihito/SALMON-dev/verification/`.

**Step 1: Create a fresh runtime directory and record provenance**

Copy the previously validated Si64 input assets, record the new commit and binary hash, and retain MPI rank count 8.

**Step 2: Run to the existing center gate**

Run the release executable with input redirected from `inputfile`.  Monitor rank RSS and progress without changing MPI rank count.

**Step 3: Interpret the new receipt**

- Large center-block leakage with small unitarity defect confirms missing point-center gauge fixing.
- Small leakage with a large center residual indicates that measured centers or affine center action, rather than internal gauge mixing, is inconsistent.
- Do not change the center tolerance based on this run.

**Step 4: Stop and design the actual repair**

Use the measured contract failure to choose a separate point-gauge repair design.  Do not bundle that repair into this diagnostic change.
