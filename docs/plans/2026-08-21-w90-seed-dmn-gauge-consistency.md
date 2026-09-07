# Wannier90 Seed–DMN Gauge Consistency Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Make the localized SALMON seed A matrix and Wannier90 DMN symmetry representation use one identical gauge.

**Architecture:** Add one bounded distributed A-only assembly entry point before Wannier90 setup, then replace the raw overlap by its full-rank unitary polar factor Q. Reuse coordinator-owned Q in DMN construction, post-Wannier covariance validation, and the existing M/A assembler so the spatial overlap and gauge choice are each evaluated only once.

**Tech Stack:** Fortran 2008, MPI, Wannier90 library adapter, Python route checks.

---

### Task 1: Add failing gauge-consistency contracts

**Files:**
- Modify: `tests/dg/test_dg_overlapping_wannier_w90_mpi.f90`
- Modify: `tests/dg/check_dg_overlapping_wannier_route.py`

1. Add an MPI oracle for the polar gauge of `A=<retained|trial>` and precomputed-Q reuse.
2. Require production DMN and post-Wannier validation to use `Q^H D Q` and the same Q payload.
3. Run the focused and route tests and confirm failure for the missing behavior.

### Task 2: Implement bounded A assembly and reuse

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_w90.f90`

1. Add a collectively safe A-only assembler using the existing tiled overlap convention and a coordinator-side full-rank polar decomposition.
2. Add optional validated precomputed-A reuse to the M/A assembler.
3. Run the focused MPI test on 1/2/4/8 ranks.

### Task 3: Wire one gauge through production

**Files:**
- Modify: `src/gs/main_dft.f90`

1. Compute Q after retaining the localized seed anchors and before DMN publication.
2. Publish `Q^H D Q`, retained D, and Q in every DMN operation.
3. Reuse Q during M assembly and post-Wannier covariance validation.
4. Remove the identity-DMN fallback and obsolete spectral-A variables from this path.
5. Run route, build, focused MPI, and diff checks.

### Task 4: Run Si64 physical validation

**Files:**
- Use: `tests/dg/data/si64_overlapping_wannier_rt/input_hybrid_scf.in`

1. Run eight MPI ranks with `OMP_NUM_THREADS=1` and correct geometry units.
2. Record Wannier90 iterations, spread/gradient trajectory, peak RSS, post-Wannier symmetry receipts, and Hybrid-SCF entry.
3. Stop and diagnose rather than increasing iteration count if the trajectory still oscillates.
