# WPW Physical-Cell Wannier Normalization Implementation Plan

> **Historical/removed:** This document describes an obsolete experimental DG route
> removed on 2026-07-31. It is retained only as an implementation record and is
> not executable guidance.

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Normalize each projected occupied Wannier on the unique canonical points available in prepared P and apply one identical column scale to its core values, P values, and analytic gradients.

**Architecture:** Add a partial-P canonical norm assembler and a checked column-normalization helper that derives `1/sqrt(norm)` and scales the projected transform and already materialized fields. Production globally collects owner norms, validates them before mutation, applies one replicated scale, and for full coverage rebuilds canonical links and retains the existing unit-norm diagnostic gate.

**Tech Stack:** Fortran 2008, MPI collective synchronization, existing DG source-contract and focused test drivers.

---

### Task 1: Checked column normalization helper

**Files:**
- Modify: `src/gs/dc/lcfo_wannier_sawf_seed.f90`
- Modify: `tests/dg/test_dg_wpw_core_wannier_seed_mpi.f90`

1. Add failing tests for `normalize_sawf_projected_wannier_columns`: exact complex scaling of `polar_transform`, core values, and P values; unit input idempotence; shape mismatch; non-positive and non-finite norm rejection.
2. Run `python3 tests/dg/run_dg_wpw_core_wannier_seed_mpi.py`; expect RED because the API is absent.
3. Implement `scale=1/sqrt(norm)` with transactional validation before mutation. Apply each scale to the matching W column of all three arrays.
4. Re-run the focused test; expect PASS.
5. Request findings-first review and fix all Critical/Important findings.

### Task 2: Unique canonical-point norm assembler

**Files:**
- Modify: `src/gs/dc/dg_wpw_occupied_w_basis.f90`
- Modify: `tests/dg/test_dg_wpw_occupied_w_basis.f90`

1. Add failing B=5 partial-P and B=10 aliased-P tests for
   `assemble_dg_wpw_canonical_buffer_norm`. Require each canonical point to be
   counted once, with missing points contributing zero, and use a non-unit
   quadrature weight to prove the result is an integral rather than a bare sum.
2. Run `python3 tests/dg/run_dg_wpw_occupied_w_basis.py`; expect RED because
   the API is absent.
3. Implement the same zero-based-origin canonical mapping and representative
   rule as periodization, selecting/counting one present representative per
   canonical point and multiplying by the positive finite quadrature weight.
4. Re-run the focused test; expect PASS, then request review.

### Task 3: Two-pass production normalization

**Files:**
- Modify: `src/gs/dc/lcfo_flux.f90`
- Modify: `tests/dg/check_dg_wpw_fixed_h_relaxation.py`

1. Add failing source-contract assertions requiring: initial unique-canonical
   norm assembly; owner placement into the global W vector; collection over
   `dc%icomm_tot`; identical local-slice extraction on every rank; collective
   validation before scaling; scaled `polar_transform`, `source_core`, and
   `buffer_values`; explicit update of the `source_values` local slice;
   clearing and rebuilding of `occupied_w_p`, canonical values, and every
   norm/link local/partial/global array; a second unique-canonical norm
   assembly/global collection with `abs(N_w-1)<=1e-8` for every B; analytic gradients after scaling; and
   unchanged `require_unit_norm=.true.` for full coverage.
2. Run the source-contract test; expect RED.
3. For every B>4, assemble first-pass unique-canonical norms. Owner roots place
   them in the existing global W vector and collect it over `dc%icomm_tot`.
   Every rank extracts the identical local fragment slice, validates it via a
   named collective stage before mutation, and applies the helper locally.
   Update `source_values(:,local_slice)=source_core`, rebuild `occupied_w_p`,
   and clear all canonical/norm/link arrays. Reassemble and globally collect
   the unique-canonical norm for every B and enforce the unchanged unit-norm
   tolerance. For globally uniform full coverage, additionally re-extract and
   reassemble links and run the existing spread diagnostic. Partial P
   continues to gradients and descriptor publication without spread.
   Do not change raw `spsi`, IDs, or gates.
4. Run serial/MPI focused tests and `cmake --build build-mpi-eigenexa -j4`; expect PASS.
5. Request findings-first review and fix all Critical/Important findings.

### Task 4: Fresh B=5/B=6/B=10 convergence validation

**Files:**
- Modify: `docs/plans/2026-07-22-wpw-physical-periodic-p.md`
- Modify: `docs/plans/2026-07-19-wpw-density-carrying-seed-session-handoff.md`

1. Run fresh B=5, B=6, and B=10 in new directories with the same executable and eight MPI ranks.
2. Require post-normalization norm error at or below `1e-8` for all 128 W at
   every B. Require no canonical mismatch for full-cell B=6/B=10; require
   successful unique-canonical norm assembly for partial B=5.
3. Require B=5 normalized descriptor construction without a larger-buffer request. Compare B=6/B=10 stable-ID centers modulo the cell, Ω, widths, and outer-shell ratios; record differences as B-convergence data.
4. Run focused tests, full build, `git diff --check`, and final findings-first review.
