# Distributed Total-Density Buffer Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Feed the stitched overlap-density assembly with the conserved distributed total-system density without replicating a global grid array.

**Architecture:** Introduce a checked MPI request/reply redistribution primitive for a row-owned real scalar field. Production maps the local `dc%rho_tot_s` slab to global physical IDs and requests its fragment buffer IDs directly.

**Tech Stack:** Fortran 2008, MPI Alltoall/Alltoallv, Python route checks, existing MPI fixture runners.

---

### Task 1: Add a failing redistribution test

**Files:**
- Modify: `tests/dg/test_dg_overlapping_wannier_construction_mpi.f90`
- Modify: `tests/dg/run_dg_overlapping_wannier_construction_mpi.py`

1. Add a fixture with cyclic row ownership and overlapping/reordered request IDs.
2. Compare returned local values with the analytic value indexed by global ID.
3. Add adverse ownership, ID, metadata, and finite-value cases.
4. Run the construction MPI test and confirm failure because the primitive is absent.

### Task 2: Implement the direct scalar redistribution

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_construction.f90`

1. Add the public primitive and collective contract checks.
2. Build owner lookup from exactly-once distributed IDs.
3. Exchange request IDs and original positions with `MPI_Alltoallv`.
4. Return values to requesters and restore original request order.
5. Add checked extent/workspace receipt and collective cleanup.
6. Run MPI 1/2/4/8 and confirm GREEN.

### Task 3: Connect production and enforce the route

**Files:**
- Modify: `src/gs/main_dft.f90`
- Modify: `tests/dg/check_dg_overlapping_wannier_route.py`

1. Make the route test require direct total-density materialization and reject fragment `rho_s` assignment.
2. Run the route test and confirm RED.
3. Construct local total-grid physical IDs from `dc%mg_tot` and call the primitive.
4. Remove the fragment-density sampling assignment.
5. Run the route test and focused MPI tests until GREEN.

### Task 4: Verify production

**Files:**
- No source changes unless a new independently reproduced defect is found.

1. Build the production executable.
2. Run Si64 with MPI 8 and OMP 1.
3. Verify `stitched_overlap_density electrons=256` within the existing tolerance.
4. Record the next completed gate or exact next failure.
5. Run `git diff --check` and review the scoped diff.
