# Batched Symmetry-Constrained Wannier Localization Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Replace edge-wise localization transactions with one exact-symmetry-commuting batched descent transaction per sweep.

**Architecture:** Collect every sparse-edge gradient into a single group-projected anti-Hermitian generator.  Use one bounded exponential line search per sweep, with collective spread acceptance and identical transforms on every fragment rank.

**Tech Stack:** Fortran 2008, MPI, LAPACK/BLAS, Python MPI fixture runner, CMake.

---

### Task 1: Gate MPI production preprocessing and spread-evaluation complexity

**Files:**
- Modify: `tests/dg/test_dg_overlapping_wannier_localization_mpi.f90`
- Modify: `tests/dg/run_dg_overlapping_wannier_localization_mpi.py`
- Modify: `src/gs/dc/dg_overlapping_wannier_localization.f90`

**Step 1: Write RED tests**

Compile through a generated `config.h` rather than `-DUSE_MPI`.  Add rank-owned graph data whose union exists only after `MPI_Allreduce`.  Add an optional `spread_evaluations` result to the sweep call and require
`spread_evaluations <= 1 + 40 * maximum_iterations`, independent of edge count.

**Step 2: Verify RED**

Run `python3 tests/dg/run_dg_overlapping_wannier_localization_mpi.py`.
Expected: the missing optional argument first fails compilation; before the preprocessing fix, the distributed graph assertion fails on 2+ ranks.

**Step 3: Implement the minimal preprocessing and counter contract**

Include `config.h` in the module.  Increment the counter only around collective periodic-spread evaluations and publish it through the optional result.

**Step 4: Verify focused GREEN**

Run the MPI fixture on 1/2/4/8 ranks and `python3 tests/dg/check_dg_overlapping_wannier_route.py`.

### Task 2: Batch the symmetry-constrained descent

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_localization.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_localization_mpi.f90`

**Step 1: Confirm the evaluation-count RED remains attributable to edge-wise line search**

Run the focused fixture and retain the expected counter failure.

**Step 2: Accumulate the batched generator**

For every pair in the local retained block, compute both collective gradient components and assemble
the full anti-Hermitian negative gradient.  Reynolds-project that matrix through the exact group
once per sweep.  The bounded support graph remains diagnostic only because an arbitrary sparse
graph is not invariant under a dense representation.  Reject nonfinite, non-anti-Hermitian, or
noncommuting projected gradients.

**Step 3: Apply one transactional line search per sweep**

Normalize the accumulated generator, exponentiate once per trial, update the full values and gradients, evaluate collective spread, and either publish the transform or restore the backup.  Reject a nonstationary sweep when all 40 trials fail.

**Step 4: Verify GREEN and invariants**

Run MPI 1/2/4/8.  Require convergence, strict spread reduction, unitary transform, exact group commutation, identical rank result, and bounded spread-evaluation count.

**Step 5: Review**

Perform specification and code-quality reviews.  Resolve all Critical and Important findings, especially descent sign, normalization, rollback completeness, finite checks, and false convergence.

### Task 3: Genuine Si64 and clean-first verification

**Files:**
- Modify only if a verified defect is found.

**Step 1: Run genuine Si64 focused verification**

Use the tracked Si64 OW input and 8 MPI ranks.  Require one graph diagnostic, decreasing spread, localization convergence, exact closure, and substantially fewer spread evaluations than edge-wise localization.

**Step 2: Run retained-route focused suites**

Run point-group, construction, symmetry, localization, checkpoint, route, and response/HHG contract fixtures on their specified ranks.

**Step 3: Perform clean-first parent-prerequisite overlay build**

Create a new clean source/build root from HEAD.  Configure MPI, ScaLAPACK, EigenExa, and spglib; build `eigenexa-project-build -j1` first; then build all with `-j4`.  Wannier90 remains off for this accepted internal retained route.

**Step 4: Review and commit**

Repeat specification and code-quality reviews, resolve all Critical/Important findings, run `git diff --check`, then commit only the localization/MPI changes.
