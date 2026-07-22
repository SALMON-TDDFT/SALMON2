# WPW Physical-Periodic P Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Build the WPW bootstrap P from one deterministic physical-cell representative per canonical point before Wannier projection and differentiation.

**Architecture:** A pure helper selects nearest-core representatives in unwrapped P and broadcasts each representative to all physical-period aliases. Production calls it on the assembled occupied-orbital P before projected-W values and analytic gradients are formed. Existing image extraction remains an independent postcondition.

**Tech Stack:** Fortran 2008, MPI stage synchronization, existing DG test drivers and MPI/EigenExa build.

---

### Task 1: Physical-periodic P helper

**Files:**
- Modify: `src/gs/dc/dg_wpw_occupied_w_basis.f90`
- Modify: `tests/dg/test_dg_wpw_occupied_w_basis.f90`

1. Write failing B=5/B=6/B=10, absolute origin/seam, core-preservation,
   idempotence, tie-break, partial-coverage, unsupported-geometry, and
   non-finite tests for
   `periodize_dg_wpw_fragment_buffer`.
2. Run `python3 tests/dg/run_dg_wpw_occupied_w_basis.py`; expect RED because
   the API is absent.
3. Implement nearest-core representative selection and exact alias fill for
   aliases present in P. Accept partial physical-cell coverage; validate
   `nxyz_domain <= total_shape` componentwise. For one-based P array index
   `i`, use logical grid coordinate `g=i-buffer` and the explicit mapping
   `modulo(fragment_origin + g - 1,total_shape)+1`, equivalently
   `modulo(fragment_origin + i - buffer - 1,total_shape)+1`.
4. Run the driver; expect PASS.
5. Request findings-first review and fix Critical/Important findings.

### Task 2: Production placement before projection and derivatives

**Files:**
- Modify: `src/gs/dc/lcfo_flux.f90`
- Modify: `tests/dg/check_dg_wpw_fixed_h_relaxation.py`

1. Add failing source-contract assertions requiring periodization after
   `unwrapped_buffer_order` and before `occupied_buffer` construction,
   projected-W transformation, and analytic gradients.
2. Run the contract test and verify RED.
3. Periodize `occupied_stencil` collectively into the WPW copy; add the
   `physical_periodic_p` synchronized stage. Do not change raw `spsi`.
4. Run focused serial/MPI tests and the full build. Check equality of aliases
   present in P for every B>4. Run full canonical-cell extraction and spread
   only when P has full physical-cell coverage; B=5 must construct P and D
   successfully without requesting a larger buffer.
5. Request review and fix Critical/Important findings.

### Task 3: Fresh B=6/B=10 validation

**Files:**
- Modify: `docs/plans/2026-07-19-wpw-density-carrying-seed-session-handoff.md`
- Modify: `docs/plans/2026-07-22-wpw-physical-periodic-p.md`

1. Run fresh B=6 and B=10 in new directories.
2. Require canonical extraction success and identical spread distributions
   within `1e-10 Angstrom` centers modulo the cell and
   `1e-10 Angstrom^2` spread.
3. Record outer-shell results without changing its tolerance.
4. Run final focused tests, full build, diff check, and findings-first review.
