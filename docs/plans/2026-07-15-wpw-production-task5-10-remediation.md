# WPW Production Task 5-10 Remediation Implementation Plan

> **Historical/removed:** This document describes an obsolete experimental DG route
> removed on 2026-07-31. It is retained only as an implementation record and is
> not executable guidance.

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Connect the distributed `windowed_kg` basis and canonical-face assembly to production DG-DC and midpoint Exp, ending with a provenance-locked Si64 2x2x2-fragment smoke test.

**Architecture:** Build one rank-local sparse operator pipeline consumed by both GS and RT matrix-free callbacks.  Serialize the converged DG-DC state through a shared checkpoint and require static DC/RT identity plus field-off stationarity before enabling the accepted length-gauge field-on route.

**Tech Stack:** Fortran 2008, MPI, SALMON LDA/Hartree infrastructure, LAPACK/EigenExa bounded Rayleigh-Ritz, Python source contracts and linked fixtures, CMake.

---

### Task 5A: Make missing production integration fail visibly

**Files:**
- Create: `tests/dg/check_dg_wpw_production_consumer.py`
- Create: `tests/dg/test_dg_wpw_production_operator_mpi.f90`
- Create: `tests/dg/run_dg_wpw_production_operator_mpi.py`
- Inspect: `src/gs/main_dft.f90`
- Inspect: `src/rt/dg/rt_dg_fragment_ops.f90`

1. Require a public GS-neutral production context initializer, rank-local quadrature builder, canonical-face scanner call, sparse builder call, H/S callback binder, and bounded occupied-subspace consumer reachable from `main_dft` under an explicit namelist route.
2. Reject `run_dg_wpw_fixed_scf` in `main_dft`, global dense H/S/WP/PP allocations, `do ifrag=1,n_frag` in scanners, hidden MPI in the scanner, G-only production columns, and PP face candidates.
3. Run `python3 tests/dg/check_dg_wpw_production_consumer.py`; expect RED because no production consumer exists.
4. Add a two-rank linked fixture that compares the assembled rank-local operator action against the existing dense mathematical oracle, including one canonical face.
5. Run the fixture; expect compile/link RED.

### Task 5B: Connect rank-local quadrature to sparse blocks

**Files:**
- Create: `src/common/dg_wpw_production_context.f90`
- Create: `src/rt/dg/rt_dg_wpw_trace_halo_provider.f90`
- Create: `src/rt/dg/rt_dg_wpw_production_builder.f90`
- Modify: `src/common/CMakeLists.txt`
- Modify: `src/rt/CMakeLists.txt`
- Modify: `tests/dg/test_dg_wpw_production_operator_mpi.f90`

1. Define plain distributed layout, owned W/P ids, support ids, sparse blocks, trace-halo state, and validity/fingerprint fields without embedding global owner tables.
2. Bind a persistent provider context whose callback reads prepared local/halo traces only; reject stale or missing halo epochs before scanning.
3. Traverse only locally owned volume boxes and the rank-local canonical face list.  Pack candidates immediately and reject non-owned output candidates.
4. Build rank-local H/S blocks with `build_windowed_sparse_wpw_operators`; leave PP face structurally absent.
5. Run the source contract and MPI fixture until GREEN, then run all `check_dg_wpw_*.py` contracts and `git diff --check`.
6. Review the scoped diff.  Do not commit because the user requested no additional commits.

### Task 5C: Implement bounded matrix-free DG-DC consumer

**Files:**
- Create: `src/gs/dc/dg_wpw_matrix_free_scf.f90`
- Modify: `src/gs/dc/CMakeLists.txt`
- Modify: `src/gs/main_dft.f90`
- Modify: `src/io/salmon_global.f90`
- Modify: `src/io/inputoutput.f90`
- Create: `tests/dg/test_dg_wpw_matrix_free_scf.py`

1. Add namelist controls for the production route, extra states, gap threshold, metric cutoff, mixing, residual tolerances, and maximum iterations.
2. Write a linked two-rank RED fixture requiring retained dimension `n_occ+extra`, S-orthonormality, charge, gap, residual, density/potential/energy/projector convergence, and one-map fixed-point reproduction.
3. Implement block subspace iteration using only batched H/S callbacks and bounded Rayleigh-Ritz.  Diagnose rank-local invalid blocks before collective failure.
4. Reconstruct rank-local mixed density and call the existing LDA/global-Hartree adapter; never allocate a full density matrix.
5. Route `main_dft` to this solver only for the explicit production namelist value; keep the dense solver unreachable.
6. Run the fixture, source contracts, full build, and a fragment-count scaling fixture; require GREEN and bounded per-rank storage.
7. Review the scoped diff; do not commit.

### Task 6: Add the shared DG-DC checkpoint

**Files:**
- Create: `src/common/dg_wpw_checkpoint.f90`
- Modify: `src/common/CMakeLists.txt`
- Modify: `src/gs/dc/dg_wpw_matrix_free_scf.f90`
- Create: `tests/dg/test_dg_wpw_checkpoint_roundtrip.py`

1. Write RED cases for valid roundtrip, corrupt checksum, truncation, nonfinite values, excessive dimensions, fingerprint mismatch, operator-convention mismatch, and interrupted publication.
2. Store distributed layout/fingerprints, occupations, metric metadata, retained eigenpairs, potential, fixed sparse blocks, face convention, and checksums in a versioned bounded schema.
3. Write a temporary file, flush/close, and atomically publish only after all ranks validate; delete no prior valid checkpoint on failure.
4. Reload with bounded allocation and mandatory checksum/finiteness/convention validation.
5. Run roundtrip tests, Task 5 fixtures, full build, and `git diff --check`; review without committing.

### Task 7: Prove DG-DC to RT identity

**Files:**
- Modify: `src/rt/dg/rt_dg_fragment.f90`
- Modify: `src/rt/dg/rt_dg_fragment_ops.f90`
- Create: `tests/dg/check_dg_dc_rt_handoff.py`
- Create: `tests/dg/test_dg_dc_rt_handoff_mpi.f90`

1. Write RED tests requiring checkpoint-backed RT initialization with no conventional-orbital projection and no metric reselection.
2. Load the shared checkpoint into the same distributed layout and rebuild rank-local callbacks.
3. Before reductions, compare fixed H/S blocks and report the detecting rank, block, ids, and residual; then check occupied generalized residual and S orthonormality.
4. Run the MPI handoff fixture, checkpoint suite, full build, and `git diff --check`; review without committing.

### Task 8: Implement production midpoint Exp

**Files:**
- Create: `src/rt/dg/rt_dg_wpw_exp_production.f90`
- Modify: `src/rt/CMakeLists.txt`
- Modify: `src/rt/dg/rt_dg_integrator_expdiag.f90`
- Modify: `src/io/salmon_global.f90`
- Modify: `src/io/inputoutput.f90`
- Create: `tests/dg/check_wpw_exp_production_route.py`
- Create: `tests/dg/test_wpw_exp_midpoint.py`

1. Audit WPW/BPW/MIXED-Z scientific environment variables; move result-changing controls into namelist or reject the production route when they are set.
2. Write RED tests for exact S-unitarity, saved-`C_n` midpoint restart, cumulative-corrector rejection, nonconvergence, and S-norm failure.
3. Implement one namelist production route using Task 5 callbacks and the shared generalized algebra.  Iterate midpoint density/potential from saved `C_n` and publish `C_(n+1)` only on convergence.
4. Run a checkpoint-started field-off fixture and require stationary density, occupied projector, energy, S norm, and `Delta_Pz`.
5. Run Tasks 2/3/5/7/8 tests, full build, and `git diff --check`; review without committing.

### Task 9: Connect the accepted length-gauge observable

**Files:**
- Modify: `src/rt/dg/rt_dg_fragment_ops.f90`
- Modify: `src/rt/dg/rt_dg_observables.f90`
- Modify: `src/io/write.f90`
- Create: `tests/dg/check_wpw_length_gauge_observable.py`
- Create: `tests/dg/test_wpw_length_gauge_observable.py`

1. Write RED tests for the accepted WW/WP/PP volume blocks, absence of PP face terms, raw metric-Hermiticity residual, field-sign oddness, discrete commutator consistency, continuous branch tracking, and `Jz=d Delta_Pz/dt`.
2. Build the length-gauge action from the same distributed layout and trace data as H/S; add no unproven interface correction.
3. Publish `Pz`, `Delta_Pz`, and `Jz` with explicit volume/branch metadata.
4. Rerun the field-off gate, a short ±field fixture, full build, and `git diff --check`; review without committing.

### Task 10: Create and run the Si64 2x2x2 checkpoint

**Files:**
- Create: `samples/exercise_dg_fragment_rt/diamond64_dc_flux_mac/inputfile_stage2d_dg_dc_wpw_smoke`
- Create: `samples/exercise_dg_fragment_rt/diamond64_dc_flux_mac/inputfile_stage2d_wpw_fieldoff_smoke`
- Create: `samples/exercise_dg_fragment_rt/diamond64_dc_flux_mac/inputfile_stage2d_wpw_laser_smoke`
- Create: `samples/exercise_dg_fragment_rt/diamond64_dc_flux_mac/stage2d_wpw_exp_manifest.tsv`
- Create: `tests/dg/check_stage2d_wpw_exp_inputs.py`

1. Derive all three inputs from one verified Si64 2x2x2 geometry and converged conventional-GS provenance.  Record pseudopotential hashes, basis/checkpoint ids, fragment/grid layout, occupation, metric, `dt=2`, propagator, observable, branch, and volume normalization.
2. Write the manifest consistency test first and confirm RED while files are absent.
3. Run DG-DC to convergence and validate its checkpoint; stop on any convergence or provenance failure.
4. Run checkpoint-backed field-off RT and require the Task 8 stationarity gate.
5. Run a 20-step Sin2 laser smoke using the explicitly recorded laser parameters; require finite observables, bounded S-norm/charge drift, and successful midpoint convergence at every step.  Do not interpret an HHG spectrum yet.
6. Run all Task 5-10 focused tests, all WPW source contracts, full build, grid/face/MPI fixtures, and `git diff --check`.

### Task 10 review checkpoint

1. Review all changes from the clean `125bdb3` baseline findings-first: production reachability, storage scaling, face uniqueness, halo epochs, callback lifetime, checkpoint publication, DC/RT identity, midpoint restart, length-gauge definition, provenance, and failure behavior.
2. For every valid finding, add a minimal RED regression before changing production code.
3. Repeat review/fix cycles until there are no Critical or Important findings.
4. Rerun the complete Task 5-10 verification bundle and confirm the worktree contains only the intended uncommitted Task 5-10 changes.
5. Do not commit, push, or create a pull request unless the user changes that instruction.
