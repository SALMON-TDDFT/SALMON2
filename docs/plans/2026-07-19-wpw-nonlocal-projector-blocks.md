# WPW Nonlocal Projector Blocks Implementation Plan

> **For Codex:** Execute with test-driven development in the existing protected worktree.

**Goal:** Add complete bounded projector-mediated WW/WP/PP nonlocal pseudopotential blocks to the production WPW operator.

**Architecture:** Basis owners create sparse `(basis ID, projector ID, overlap)` records using SALMON projector conventions, exchange them only with bounded support peers, and row owners form outer-product blocks.  Fixed nonlocal blocks remain separate from SCF-updated local potential blocks.

**Tech Stack:** Fortran 2008, MPI, SALMON pseudopotential grids, Python test runners, CMake.

---

### Task 1: RED dense projector outer-product oracle

**Files:**
- Create: `src/common/dg_wpw_nonlocal_projector.f90`
- Create: `tests/dg/test_dg_wpw_nonlocal_projector.f90`
- Create: `tests/dg/run_dg_wpw_nonlocal_projector.py`

1. Define sparse overlap records with stable basis/projector IDs.
2. Test signed projector weights and complex overlaps against a direct dense oracle.
3. Reject duplicate, missing, and nonfinite records.
4. Run before implementation and confirm the missing module/API failure.

### Task 2: Implement bounded outer-product assembly

1. Validate sorted unique overlap records and finite weights.
2. Intersect projector IDs for requested row/column basis pairs.
3. Assemble WW/WP/PP values without changing requested sparse-entry order.
4. Run the dense oracle until GREEN.

### Task 3: RED/GREEN bounded overlap exchange

**Files:**
- Modify: `src/common/dg_wpw_nonlocal_projector.f90`
- Create: `tests/dg/test_dg_wpw_nonlocal_projector_mpi.f90`
- Create: `tests/dg/run_dg_wpw_nonlocal_projector_mpi.py`

1. Route owned overlap records to only declared support peers.
2. Require the two-rank assembled blocks to match the global dense oracle.
3. Reject duplicate ownership and truncated payloads collectively.

### Task 4: Integrate GS projector overlaps

**Files:**
- Modify: `src/gs/dc/lcfo_flux.f90`
- Modify: `src/rt/dg/rt_dg_wpw_production_builder.f90`
- Modify: `src/common/dg_wpw_production_context.f90`
- Modify: relevant CMake lists.

1. Build owned W/P overlaps from `dc%ppg_tot` with the existing `hvol`, projector values, and weights.
2. Exchange sparse records across the bounded support peers.
3. Build authoritative fixed WW/WP/PP nonlocal blocks.
4. Compare owned-owned WW against the existing LCFO nonlocal block.
5. Add fixed nonlocal terms during initial build and every potential rebuild.

### Task 5: Verification and Si64 checkpoint gate

1. Remove temporary LOBPCG success traces.
2. Run projector oracles, bounded H/S, matrix-free SCF, routing, checkpoint, and provenance tests.
3. Run the full MPI/EigenExa build and `git diff --check`.
4. Run a fresh 8-rank full-basis Si64 preflight.
5. Do not launch production GS unless Ritz values/residuals are finite and physically consistent.

No commit, push, or pull request is allowed.
