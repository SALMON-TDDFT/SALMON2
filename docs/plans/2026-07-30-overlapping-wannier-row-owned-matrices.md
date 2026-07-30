# Overlapping-Wannier Row-Owned Matrices Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Remove persistent full \(N_W^2\) overlap and Hamiltonian arrays from every MPI rank while preserving the accepted overlapping-Wannier algebra and route isolation.

**Architecture:** Reuse the solver's `row_ids`, `hrows`, and `srows` distribution. Unique-core ranks form partial destination rows and reduce them to row owners; exact metric spectral diagnostics use only a root-temporary full matrix, while checkpoints persist row shards.

**Tech Stack:** Fortran 2008, MPI collectives/reductions, LAPACK on rank zero for the metric gate, Python source contracts, focused MPI fixtures, CMake.

## Global Constraints

- Work only in `/Users/otobetoshihito/SALMON-dev/SALMON2_RTDG/.worktrees/wpw-s-orthogonal-complement`.
- Keep the parent worktree read-only and do not import its uncommitted implementation.
- This work remains inside Task 9; Task 10 stays blocked until the complete Si64 ground-state matrix gate passes.
- Do not enable direct real-space DG, LCFO, EigenExa, WPW, normal checkpoint publication, or conventional RT.
- Preserve normal DC LCFO plus EigenExa behavior when the route is disabled.
- Every task requires recorded RED, focused verification, specification review, code-quality review, Critical/Important resolution, `git diff --check`, and a clean-first parent-prerequisite overlay build.

---

### Task 1: Assemble weak operators directly into owned rows

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_operators.f90`
- Modify: `src/gs/main_dft.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_operators_mpi.f90`
- Create: `tests/dg/check_dg_overlapping_wannier_row_storage.py`

**Interfaces:**
- Consumes: `row_ids(:)` containing a unique partition of global Wannier rows and local unique-core `values`, `gradients`, and potentials.
- Produces: `assemble_dg_overlapping_wannier_weak_operator_rows(...,row_ids,...,kinetic_rows,potential_rows,...)`, with both outputs shaped `(size(row_ids),nwann)`.

- [x] **Step 1: Write and run the RED source and MPI contracts**

Require the new row interface, reject `allocate(...nwann,nwann...)` and
full-matrix `MPI_Allreduce` inside production weak assembly, and compare its
rows against the existing direct reference on 1, 2, 4, and 8 ranks.

- [x] **Step 2: Implement destination-row reduction**

Compute local-core partial rows in bounded row batches. Use the collective
row ownership map to reduce each batch to its owner and return only owned
rows. Detect duplicate/missing row owners before payload collectives.

- [x] **Step 3: Switch `ow_build_hamiltonian` to owned rows**

Add kinetic, local, and nonlocal contributions without allocating
`global_h`; compute global Hermiticity and finite-value diagnostics from
distributed rows.

- [x] **Step 4: Run focused verification and reviews**

Run the source contract, operator fixture on 1/2/4/8 ranks, SCF fixture on
1/2/4/8 ranks, route checker, `git diff --check`, specification review, and
code-quality review. Resolve every Critical/Important finding.

### Task 2: Return the metric as owned rows

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_metric.f90`
- Modify: `src/gs/main_dft.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_metric_mpi.f90`
- Modify: `tests/dg/check_dg_overlapping_wannier_row_storage.py`

**Interfaces:**
- Consumes: the same exact `row_ids(:)` ownership used by the coefficient solver.
- Produces: `srows(size(row_ids),nwann)`, broadcast spectrum and retained transformation, with any full metric confined to rank zero and released before return.

- [x] **Step 1: Extend and run the RED contracts**

Require non-root ranks never allocate the full metric and verify identical
spectrum, rank rejection, condition number, and owned rows on 1/2/4/8 ranks.

- [x] **Step 2: Implement root-temporary spectral validation**

Reduce row blocks to owners, gather owned rows only to rank zero, solve the
Hermitian metric there, broadcast diagnostics and retained transformation,
then deallocate the root full matrix before returning.

- [x] **Step 3: Remove persistent `ow_s`**

Store `ow_srows` only, pass it directly to fingerprinting and SCF, and
compute metric diagnostics collectively from owned rows.

- [x] **Step 4: Run focused verification and reviews**

Run metric, solver, SCF, and checkpoint fixtures on 1/2/4/8 ranks, the
source contract, route checker, `git diff --check`, both reviews, and resolve
every Critical/Important finding.

### Task 3: Persist and restore row-sharded overlap

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_checkpoint.f90`
- Modify: `src/gs/main_dft.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_checkpoint_mpi.f90`
- Modify: `tests/dg/check_dg_overlapping_wannier_row_storage.py`

**Interfaces:**
- Consumes: exact global `row_ids(:)` and `srows(:,:)`.
- Produces: one shard per rank containing row IDs and overlap rows; restart returns only the caller's owned rows after collective uniqueness validation.

- [x] **Step 1: Extend and run the RED checkpoint contract**

Reject a checkpoint implementation that broadcasts or allocates a full
overlap on every rank. Exercise uneven and zero-row owners and reject
duplicate, missing, reordered-without-ID, and corrupt shards.

- [x] **Step 2: Implement row-sharded checkpoint I/O**

Write row IDs and row payload per rank, include them in provenance hashes,
and validate the global row partition before accepting restart.

- [x] **Step 3: Run focused verification and reviews**

Run checkpoint, SCF, solver, metric, and operator fixtures on 1/2/4/8 ranks,
restart reuse, source/route contracts, `git diff --check`, both reviews, and
resolve every Critical/Important finding.

### Task 4: Re-run the Task 9 smoke and memory gate

**Files:**
- Modify: `tests/dg/check_si64_overlapping_wannier_gate.py`
- Modify: `tests/dg/run_si64_overlapping_wannier_gate.py`
- Modify: `docs/plans/2026-07-27-si64-overlapping-wannier-results.md`

**Interfaces:**
- Consumes: row-owned metric, Hamiltonian, solver, and checkpoint path.
- Produces: raw evidence that the accepted smoke retains its numerical gates and no rank persistently allocates a full \(N_W^2\) operator.

- [x] **Step 1: Write and run the RED memory evidence check**

Require row counts, local operator bytes, peak RSS, restart identity, and
absence of replicated full-matrix markers.

- [x] **Step 2: Run the accepted Si64 smoke and restart**

Use 8 MPI ranks and 1 OpenMP thread. Verify unchanged route isolation,
metric/SCF/charge gates, immutable manifest, and reduced per-rank matrix
storage.

- [x] **Step 3: Complete reviews and overlay build**

Perform specification and code-quality reviews, resolve every
Critical/Important finding, run all focused Task 9 checks and a clean-first
parent-prerequisite overlay build. Do not claim the full Task 9 matrix gate
or begin Task 10.
