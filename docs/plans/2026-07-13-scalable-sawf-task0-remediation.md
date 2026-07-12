# Scalable SAWF Task 0 Remediation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Complete Task 0 with production-connected, same-supercell SAWF template generation for bulk, defect, interface, surface/vacuum, and amorphous local environments.

**Architecture:** Canonically fingerprint the exact supercell and every core+buffer environment, form orbits only from exact fingerprints and the actual supercell group, and independently generate every inequivalent environment. Align independently generated buffered bases with metric-whitened U(N) Procrustes before constructing WW/WP/DG blocks; cache only the final validated basis inside the same supercell restart lineage.

**Tech Stack:** Fortran 2008, MPI, LAPACK, existing LCFO/Wannier90 import, Python 3 focused drivers, CMake.

## Global Constraints

- Reuse is restricted to one exact supercell and its restart lineage; cross-supercell import is forbidden.
- Actual-supercell symmetry is authoritative; parent symmetry is diagnostic provenance only.
- Local-environment equality is exact canonical equality, not tolerance-based structural similarity.
- Rank-local failures print fragment, orbit, operation, rank, singular values, and residual before any collective reduction.
- Gauge alignment uses common buffered-volume Gram and cross-overlap matrices, not face-only overlap.
- Basis/state coefficients are rotated before WW/WP/DG face construction.
- Scientific controls remain in `&dc`, never environment variables.

---

### Task 0A: Canonical supercell and local-environment fingerprints

**Files:**
- Modify: `src/gs/dc/lcfo_wannier_sawf_templates.f90`
- Modify: `src/gs/dc/lcfo_flux.f90`
- Create: `tests/dg/test_sawf_environment_fingerprint.py`

**Interfaces:**
- Produces: `build_sawf_supercell_fingerprint(...)`, `build_sawf_local_environment_fingerprint(...)`, and `classify_sawf_environment_orbits(...)`.

- [ ] Write a compiled Fortran test with bulk-equivalent, defect, two-material interface, surface/vacuum, and amorphous fixtures. Require only exact same-supercell bulk environments to share an orbit.
- [ ] Run `python3 tests/dg/test_sawf_environment_fingerprint.py`; expect failure because the fingerprint APIs do not exist.
- [ ] Implement canonical fixed-width serialization of lattice/PBC, all species/coordinates, pseudopotential content hashes, grid, buffer, band window, shells, XC and schema; implement core+buffer local serialization including vacuum occupancy.
- [ ] Replace `sawf_defect_intersects=.false.` in `generate_sawf_dmn` with classification output and actual-group orbit traversal.
- [ ] Run the compiled test and `python3 tests/dg/test_sawf_local_materialization.py`; expect PASS.
- [ ] Commit with `git commit -m "feat: classify exact SAWF local environments"`.

### Task 0B: Production representative and independent basis materialization

**Files:**
- Modify: `src/gs/dc/lcfo_flux.f90`
- Modify: `src/gs/dc/lcfo_wannier_sawf_templates.f90`
- Create: `tests/dg/test_sawf_production_materialization.py`

**Interfaces:**
- Consumes: Task 0A orbit/representative/regeneration arrays.
- Produces: buffered complex basis/state arrays for every local environment and explicit representative/operation provenance.

- [ ] Write a compiled/runtime test requiring a production call from `generate_sawf_dmn` through representative generation, actual grid action and `D_wann`, while inequivalent environments call the independent generator.
- [ ] Run the test; expect failure because current materialization is test-only/dead code.
- [ ] Add fragment-root distribution of representative buffered bases and actual symmetry point maps; materialize equivalent environments and regenerate every inequivalent environment.
- [ ] Validate direct-orbital covariance for a noncommuting operation pair and reject incomplete atom/grid maps before collective communication.
- [ ] Run materialization, band-representation and MPI fragment-map tests; expect PASS.
- [ ] Commit with `git commit -m "feat: materialize production SAWF environments"`.

### Task 0C: Buffered-volume gauge connection before operators

**Files:**
- Modify: `src/gs/dc/lcfo_flux.f90`
- Modify: `src/gs/dc/lcfo_wannier_sawf_templates.f90`
- Modify: `tests/dg/test_sawf_gauge_application.py`
- Create: `tests/dg/test_sawf_production_gauge_order.py`

**Interfaces:**
- Consumes: Task 0B buffered complex bases.
- Produces: aligned basis/state coefficients, per-neighbor U(N), rank and projector residuals.

- [ ] Add a production-order test requiring buffered Gram/cross-overlap, whitening and gauge rotation before every WW/WP/face builder call; require a non-unit face Gram not to control acceptance.
- [ ] Run the test; expect failure because production calls are absent.
- [ ] Form common-buffer quadrature matrices, apply relative metric cutoff, compute U(N) polar alignment and validate aligned occupied-projector/subspace residuals.
- [ ] Rotate basis/state coefficients first; rebuild derivatives, nonlocal data and face traces from the rotated basis. Keep post-hoc covariance transforms diagnostic-only.
- [ ] Emit rank-local diagnostics and stop before reduction on rank loss or residual failure.
- [ ] Run gauge tests and MPI build; expect PASS.
- [ ] Commit with `git commit -m "feat: connect buffered SAWF gauges before operators"`.

### Task 0D: Three-buffer operator convergence

**Files:**
- Modify: `src/gs/dc/lcfo_flux.f90`
- Modify: `src/gs/dc/lcfo_wannier_sawf_templates.f90`
- Rewrite: `tests/dg/test_sawf_buffer_convergence.py`

**Interfaces:**
- Consumes: aligned production bases and namelist `wannier_sawf_buffer_steps`.
- Produces: measured largest-two-buffer residual record and admission result.

- [ ] Replace token/NumPy checks with a compiled driver that constructs three buffers and independently evaluates centers, projector, overlap, WW, WP, and all six face blocks.
- [ ] Run the test; expect failure because the production driver/gate is absent.
- [ ] Execute all three configured buffers, compare the largest two for every component, record component maxima and reject if any exceeds `wannier_sawf_buffer_tolerance`.
- [ ] Run buffer and focused SAWF tests; expect PASS.
- [ ] Commit with `git commit -m "feat: gate SAWF operator buffer convergence"`.

### Task 0E: Safe same-supercell checkpoint publication and reload

**Files:**
- Modify: `src/gs/dc/lcfo_wannier_sawf_templates.f90`
- Modify: `src/gs/dc/lcfo_flux.f90`
- Rewrite: `tests/dg/test_sawf_template_checkpoint.py`

**Interfaces:**
- Consumes: Task 0A fingerprints and Task 0C/0D validated payload.
- Produces: atomic same-supercell cache with checksum and mandatory post-load revalidation.

- [ ] Add corrupt, truncated, foreign-supercell, nonfinite, excessive-dimension and interrupted-publication tests.
- [ ] Run the test; expect failures against native `status='replace'` checkpoint handling.
- [ ] Write versioned fixed-width payload to a temporary file, include bounded dimensions and checksum, flush/close, then atomically publish; forbid foreign supercell IDs.
- [ ] On reload validate checksum/finiteness/dimensions, then rerun group, rank, gauge and operator gates before reuse.
- [ ] Run checkpoint tests and build; expect PASS.
- [ ] Commit with `git commit -m "feat: publish validated same-supercell SAWF cache"`.

### Task 0F: Actual monolithic-versus-hierarchical equivalence and review gate

**Files:**
- Rewrite: `tests/dg/test_sawf_global_local_equivalence.py`
- Modify: `src/gs/dc/lcfo_flux.f90`
- Modify: `docs/plans/2026-07-12-scalable-sawf-contract.md`

**Interfaces:**
- Consumes: complete production hierarchical basis/operator.
- Produces: Task-0 admission report for occupied projector, S, fixed H and every face block.

- [ ] Replace the independent NumPy identity with a small-system run executing both monolithic and hierarchical routes.
- [ ] Run it; expect failure until production checkpoint/operator outputs are available.
- [ ] Perform one global metric-unitary alignment and compare occupied projector, overlap, `H_kin+DG+V_ion`, WW/WP and every face block at `wannier_sawf_equivalence_tolerance`.
- [ ] Run all Task-0 focused tests and `cmake --build build-mpi-eigenexa -j 2`; expect zero failures.
- [ ] Request findings-first review of the Task-0 commit range and fix all P0/P1 findings.
- [ ] Commit with `git commit -m "test: validate hierarchical SAWF operator equivalence"` and stop before Task 1 for the required report.
