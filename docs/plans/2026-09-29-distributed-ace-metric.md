# Distributed ACE Metric Implementation Plan

> 2026-10-06：本記録中の`developer_tests`は当時のローカル開発検証です。GitHub配布から除外しました。通常の回帰試験は`testsuites`を使用します。

> Execute inline in the existing hybrid feature checkout using executing-plans.

**Goal:** Remove replicated ACE metric/factor matrices from the native Gamma spatial/orbital path with ScaLAPACK.

**Architecture:** Stream metric columns into two-dimensional block-cyclic tiles. Diagonalize tiles with ScaLAPACK, retain distributed factors for sparse W, and use vector reductions for rotations/application. Preserve serial/no-ScaLAPACK fallback and factorized application, thresholds and physics.

**Tech Stack:** Fortran, communication wrappers, BLACS/ScaLAPACK, MPI synthetic tests.

## Tasks
1. Add `developer_tests/651_hybrid_exchange/metric/driver.f90` and runner. Compare old/new dense and packed actions over spatial/orbital/mixed layouts, including empty ranks, arbitrary targets and error cases. First compile must fail on the missing optional combined communicator argument.
2. Add `src/xc/exx_distributed_metric.f90`, primitive tile metadata in `s_exx_ace`, optional `comm_matrix` in `orbital_ace_build`, native `icomm_ro` wiring and CMake source. Keep BLACS ownership temporary so copied ACE states are safe. No global N-square scratch.
3. Update manual test object dependencies. Build GNU ScaLAPACK; run synthetic MPI and native short regressions. Compile/run no-ScaLAPACK fallback. Check coding rules and scalar IEEE inquiries.
4. Record tested tile memory and remaining replicated MLWF/hermitian-correction matrices. Review complete patch; no production timing or Fugaku validation claims without measurements.

## Constraints and progress
- Existing feature checkout is reused; production executable is a separate saved copy.
- No new physical truncation, user memory cap, or CTest registration redesign.
- Scope is the ACE construction/application step; MLWF and LCFO temporaries remain subsequent work.

## Completion record
- Missing `comm_matrix` API compile failure observed before implementation.
- Distributed build/apply, CMake/native integration, fallback and script dependencies implemented.
- MPI 1/2/3/4/8, N=2/7/65 and MPI 1/2/4/8 N=128 passed; source/fallback logs are in workspace `work/distributed-ace-*`.
- Native HSE06/PBE0/PBEh40 RT/guard tests: 7 passed. MPI and non-MPI HSE builds passed.
- Independent read-only review found no important defects; retain distinction between factor allocation and whole-process RSS.
- Additional projected-seed regression exposed a stale expectation of the old identity saddle. Reproduced identical failure using HEAD ACE modules, then updated the assertion to require Gamma Jacobi localization; no localization algorithm change.
- No CTest registration redesign, production benchmark rerun, or push in this task.
