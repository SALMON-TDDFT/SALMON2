# Wannier Partition Count Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add an input-level Wannier partition count so the DG fragment grid and the Wannier construction grid can be specified independently.

**Architecture:** Keep `num_fragment` as the DG/RT fragment decomposition.  Add `num_wannier_cluster` as the number of Wannier construction regions.  The existing `wannier_cluster_size` remains supported for compatibility, but when `num_wannier_cluster` is specified it is converted internally to an effective cluster size.

**Tech Stack:** Fortran namelist parsing in SALMON, DC-LCFO Wannier export code, existing MPI 2x2x2 diagnostic.

---

### Task 1: Add Input Variable

**Files:**
- Modify: `src/io/salmon_global.f90`
- Modify: `src/io/inputoutput.f90`

**Steps:**
1. Add `num_wannier_cluster(3)` to global DC input state.
2. Add it to the `&dc` namelist.
3. Default to `0,0,0` so old inputs remain unchanged.
4. Broadcast and log it with the other DC/Wannier parameters.

### Task 2: Derive Effective Cluster Size

**Files:**
- Modify: `src/io/inputoutput.f90`
- Modify: `src/gs/dc/lcfo_flux.f90`

**Steps:**
1. If `num_wannier_cluster` is positive, validate that it divides `num_fragment`.
2. Convert it to `wannier_cluster_size = num_fragment / num_wannier_cluster`.
3. Keep the old `wannier_cluster_size` validation path for legacy inputs.
4. Update logs to print both the input partition count and the effective cluster size.

### Task 3: 2x2x2 Diagnostic Input

**Files:**
- Runtime input under `/tmp/dg_wannier_trace_diamond/inputfile_gs`

**Steps:**
1. Set `num_fragment = 2,2,2`.
2. Set `num_wannier_cluster = 1,1,1`.
3. Keep `OMP_NUM_THREADS=2` and `mpirun -np 8`.
4. Run with `SALMON_WANNIER90_COMMAND=reuse` when `c64_dc.chk` already exists.

### Task 4: Verification

**Commands:**
- `cmake --build build-mpi-eigenexa-wannier-lib -j 4`
- `git diff --check`
- `mpirun -np 8 ...`

**Expected:**
- The log reports `num_wannier_cluster=1 1 1`.
- The effective `cluster_size` is `2 2 2`.
- Wannier cluster manifest has one cluster covering fragments `1:2 1:2 1:2`.
