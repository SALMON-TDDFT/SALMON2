# WPW Fragment-Root Ownership Remediation Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Replace the provisional all-rank arithmetic P-column ownership with a bounded fragment-root layout that preserves every window-overlap volume and canonical-face contribution.

**Architecture:** Fragment roots form a production communicator ordered by fragment ID and own all G columns for their fragment. Spatial fragments integrate uniquely owned core points using a bounded window-support halo, then route sparse candidates to P-row owners through fixed neighbor schedules; scanners remain communication-free and callbacks use the same bounded support schedule.

**Tech Stack:** Fortran 2008, MPI, SALMON D&C fragment communicators, linked MPI fixtures, Python source contracts, CMake.

---

### Task R1: Define fragment-root column ownership

**Files:**
- Modify: `src/rt/dg/rt_dg_wpw_column_layout.f90`
- Modify: `src/common/dg_wpw_production_context.f90`
- Create: `tests/dg/test_dg_wpw_fragment_root_layout_mpi.f90`
- Create: `tests/dg/run_dg_wpw_fragment_root_layout_mpi.py`
- Create: `tests/dg/check_dg_wpw_fragment_root_layout.py`

1. Write a four-rank RED fixture representing two fragments with two ranks each. Require only `id_frag==0` ranks to publish a context, production rank `K-1`, and owned IDs `(K-1)*n_G+1:K*n_G` even when total-rank ordering is not the ownership rule.
2. Add RED cases for duplicate/missing fragment roots, wrong production rank, overflow, invalid stable IDs, and accidental global owner arrays.
3. Run the linked fixture and source contract; require failure from the current arithmetic layout.
4. Add an explicit ownership-kind field and O(1) fragment-root owner/pair functions. Keep the old arithmetic initializer reachable only from existing mathematical fixtures.
5. Add a two-phase context initializer: all-rank preflight/split participation followed by root-only validation. Nonroots receive an inactive context and cannot bind callbacks.
6. Run the new fixture, existing sparse-builder/layout tests, full build, and `git diff --check`; require GREEN.

### Task R2: Build bounded window-support schedules

**Files:**
- Create: `src/common/dg_wpw_fragment_support.f90`
- Modify: `src/common/CMakeLists.txt`
- Create: `tests/dg/test_dg_wpw_fragment_support.f90`
- Create: `tests/dg/run_dg_wpw_fragment_support.py`
- Modify: `tests/dg/check_dg_wpw_production_consumer.py`

1. Write RED fixtures for face, edge, and corner overlaps; periodic wrapping; two-fragment duplicate images; one `q_K` per fragment in normalization; and buffer extent beyond the supported stencil.
2. Require construction from structured fragment coordinates/local decomposition inputs without `do ifrag=1,n_frag` and without MPI.
3. Implement bounded support records containing fragment ID, periodic image shift, overlap box, and owner rank. Keep canonical-face records separate.
4. Sort deterministically, retain distinct image regions, and expose a deduplicated fragment list to the normalized-window evaluator.
5. Run support/window/G-mode fixtures, source contracts, full build, and `git diff --check`; require GREEN.

### Task R3: Prepare the volume-support W halo

**Files:**
- Create: `src/rt/dg/rt_dg_wpw_volume_halo_provider.f90`
- Modify: `src/rt/CMakeLists.txt`
- Create: `tests/dg/test_dg_wpw_volume_halo_mpi.f90`
- Create: `tests/dg/run_dg_wpw_volume_halo_mpi.py`

1. Write a RED MPI fixture in which a remote-owned W overlaps a local core at a corner. Require exact value/gradient delivery and no remote reconstruction.
2. Add RED cases for stale epoch, bad checksum, truncated count, nonfinite payload, unexpected sender/image, invalid outbound header, and a peer-hang guard.
3. Implement a fixed neighbor schedule: validate bounded headers, post receives before sends, complete every posted request, then collectively validate status. Do not return early during an active exchange.
4. Store only overlap-box W values/gradients and stable IDs. Bind a persistent read-only provider whose epoch must match quadrature.
5. Run the new MPI fixture, trace-provider fixtures, full build, and `git diff --check`; require GREEN.

### Task R4: Integrate uniquely owned core volume points

**Files:**
- Create: `src/gs/dc/dg_wpw_rank_local_quadrature.f90`
- Modify: `src/gs/dc/CMakeLists.txt`
- Modify: `src/gs/dc/lcfo_flux.f90`
- Create: `tests/dg/test_dg_wpw_rank_local_quadrature_mpi.f90`
- Create: `tests/dg/run_dg_wpw_rank_local_quadrature_mpi.py`

1. Write a RED linked fixture with nontrivial normalized windows, local and remote W support, local potential, nonidentity WP/PP overlap, and a dense mathematical oracle.
2. Require each core point to be visited exactly once, remote W to come only from the current volume-halo epoch, and P values to use `chi_K exp(iG.r)/sqrt(Omega_cell)`.
3. Add fail-closed RED cases for missing support, duplicate core ownership, nonfinite input, invalid column IDs, and incomplete fragment reduction.
4. Implement local point assembly through `evaluate_windowed_kg_point`, `assemble_wpw_volume_point`, and immediate deterministic packing. Reduce fragment-rank partial buffers to the fragment root only.
5. Connect `lcfo_flux` after basis gradients are available and before they are deallocated. Do not scan all fragments or allocate global WP/PP.
6. Run the new fixture, Task 5 operator fixture, WPW contracts, full build, and `git diff --check`; require GREEN.

### Task R5: Route staged candidates to P-row owners

**Files:**
- Create: `src/rt/dg/rt_dg_wpw_candidate_halo.f90`
- Modify: `src/rt/CMakeLists.txt`
- Modify: `src/rt/dg/rt_dg_wpw_production_builder.f90`
- Modify: `src/rt/dg/rt_dg_wpw_sparse_builder.f90`
- Create: `tests/dg/test_dg_wpw_candidate_halo_mpi.f90`
- Create: `tests/dg/run_dg_wpw_candidate_halo_mpi.py`

1. Write a RED MPI fixture where a spatial fragment produces WP and PP rows owned by its face, edge, and corner neighbors. Compare the routed owned blocks with a dense oracle.
2. Add RED cases for unexpected sender, unsupported ID, missing expected message, duplicate image record, duplicate final pair, bad epoch, invalid payload, and peer-hang prevention.
3. Implement fixed-schedule candidate delivery using O(1) fragment-root ownership. Aggregate distinct image-region contributions by final stable-ID pair in deterministic order.
4. Permit non-owned rows only in unpublished staging buffers. Keep sparse-builder rejection of non-owned final candidates.
5. Route complete canonical-face blocks after scanning; never communicate from the scanner or publish a partial face.
6. Run candidate, operator, canonical-face, two-fragment-periodic, full build, and `git diff --check`; require GREEN.

### Task R6: Bind distributed H/S callbacks

**Files:**
- Create: `src/rt/dg/rt_dg_wpw_coefficient_halo.f90`
- Modify: `src/common/dg_wpw_production_context.f90`
- Modify: `src/rt/dg/rt_dg_wpw_production_builder.f90`
- Modify: `tests/dg/test_dg_wpw_production_operator_mpi.f90`
- Create: `tests/dg/test_dg_wpw_fragment_root_apply_mpi.f90`

1. Write a RED MPI fixture requiring H/S action on distributed owned W/P slices with neighbor support exchange and reverse accumulation to owners. Compare against the dense oracle.
2. Add RED cases for stale layout/epoch, missing coefficient slice, callback after invalidation, wrong fingerprint, nonfinite coefficient, and incomplete reverse exchange.
3. Implement the bounded coefficient schedule from final sparse block support. Exchange only required W/P slices and reverse-reduce output through the same schedule.
4. Bind callbacks only after communicator, quadrature, trace, candidate, operator, and coefficient schedules share one valid epoch/fingerprint.
5. Run apply/operator/matrix-free SCF fixtures, all WPW contracts, full build, and `git diff --check`; require GREEN.

### Task R7: Complete Task 5B/5C production reachability

**Files:**
- Modify: `src/gs/dc/lcfo_flux.f90`
- Modify: `src/gs/main_dft.f90`
- Modify: `src/gs/dc/dg_wpw_matrix_free_scf.f90`
- Modify: `tests/dg/check_dg_wpw_production_consumer.py`
- Modify: `tests/dg/test_dg_wpw_matrix_free_scf.py`

1. Extend the existing intentional RED contract to require the fragment-root initializer, window/volume halo, candidate halo, sparse builder, callback binder, and matrix-free solver in the explicit production route.
2. Connect the complete pipeline in `lcfo_flux`; nonroots participate only in required fragment collectives and never expose a context.
3. Route `main_dft` to bounded matrix-free DG-DC only after the root pipeline is valid. Keep the dense fixed-basis solver unreachable.
4. Add fragment-count scaling checks proving per-root metadata/storage depends only on owned and bounded support blocks.
5. Run all Task 5 focused tests, every `check_dg_wpw_*.py` contract, MPI fixtures, full MPI/EigenExa build, and `git diff --check`.
6. Review the entire ownership remediation findings-first. Add a RED regression before every valid fix and repeat until no Critical or Important finding remains.

### Task R8: Propagate the ownership fingerprint to Tasks 6-10

**Files:**
- Modify later Task 6 checkpoint schema and Task 7 RT handoff files when those tasks are executed.
- Modify: `docs/plans/2026-07-15-wpw-production-task5-10-remediation.md` only if the user later authorizes changing that intentionally untracked approved plan; otherwise treat this plan as its ordered Task 5B supplement.

1. Add the ownership kind/version, owned fragment, owned/support ID hashes, geometry/window-support hash, and neighbor-schedule hash to the shared checkpoint fingerprint.
2. Require RT to recreate the same fragment-root layout without redistribution or metric reselection.
3. Rerun checkpoint identity, field-off stationarity, and Si64 provenance gates under the fragment-root layout.
4. Do not commit, push, or create a pull request.
