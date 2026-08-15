# Point-Closed Periodic-Center Gauge Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Make the reference-sector periodic-center gauge explicitly invariant under the compact point-cogroup representation.

**Architecture:** Assemble row-owned point overlaps into small `m x m x npoint` matrices, close the six periodic-position Hermitian matrices under those representations, and jointly diagonalize the closed set with objective and block-monomial gates.

**Tech Stack:** Fortran 2008, MPI, BLAS/LAPACK, Python route contracts, focused MPI fixtures.

---

### Task 1: Add point-closed canonicalizer REDs

1. Extend the W90 MPI fixture with a two-center point-swap representation.
2. Require point-rotation covariance and block-monomial closure.
3. Add nonunitary and rank-disagreeing point-representation rejects.
4. Run MPI fixture and observe the missing API failure.

### Task 2: Close and validate the joint objective

1. Add point representations and provenance to the canonicalizer API.
2. Validate payload agreement and unitarity collectively.
3. Build `D_p^H H_q D_p` and reuse the bounded Jacobi kernel.
4. Gate normalized final objective and final center-block leakage.
5. Update checked workspace and fingerprint receipts.
6. Run W90 MPI 1/2/4/8.

### Task 3: Connect compact production point overlaps

1. Add a route RED requiring point-overlap assembly before canonicalization.
2. Assemble the reference-sector point matrices with the existing streamed overlap routine.
3. Gather only the compact row-owned representation and pass it to the canonicalizer.
4. Release overlap workspace before character iteration.
5. Run route checks and Release build with `-j 1`.

### Task 4: Verify Si64

1. Run a fresh MPI-8 Si64 calculation with one thread per rank.
2. Compare joint objective, alignment time, center leakage, and operation-2 residual.
3. Keep the original center gate authoritative and do not relax tolerance.

