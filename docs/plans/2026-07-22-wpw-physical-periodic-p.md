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

### Task 3: Fresh B=5/B=6/B=10 validation

**Files:**
- Modify: `docs/plans/2026-07-19-wpw-density-carrying-seed-session-handoff.md`
- Modify: `docs/plans/2026-07-22-wpw-physical-periodic-p.md`

1. Run fresh B=5, B=6, and B=10 in new directories.
2. Require post-normalization norm error below `1e-8`; require canonical
   extraction for full-cell B=6/B=10. Compare spreads as B-convergence data;
   do not assume fragment SCFs at different B are numerically identical.
3. Record outer-shell results without changing its tolerance.
4. Run final focused tests, full build, diff check, and findings-first review.

#### 2026-07-22 normalized validation result

- B=5 (`20260722_task11_normalized_b5`): the input is formally compatible
  with the radius-4 derivative when the independent production window remains
  2. The Si64 fragment SCF became unstable after iteration 75 and remained at
  `diff=1.2466E-01` at iteration 167, so it was stopped before WPW bootstrap.
  B=5 is therefore not a converged Si64 parameter point, and its end-to-end
  WPW path remains untested because of this upstream SCF failure. Partial-P
  construction without requesting a larger buffer is verified by the focused
  helper and source-contract tests, not by this runtime.
- B=6 (`20260722_task11_normalized_b6/run.log`): SCF converged in 87
  iterations (`diff=9.8498E-10`). All 128 W had post-normalization maximum
  norm error `1.11022E-15`, no canonical mismatch, and width range
  `1.26877--1.29094 A` with mean `1.28540 A`.
- B=10 (`20260722_task11_normalized_b10/run.log`): SCF converged in 85
  iterations (`diff=9.4464E-10`). All 128 W had post-normalization maximum
  norm error `1.15463E-14`, no canonical mismatch, and width range
  `1.23411--1.26600 A` with mean `1.25543 A`.

At printed precision, the largest stable-ID center-component difference
between B=6 and B=10 is `7.531E-03 A`, the largest Omega difference is
`8.676E-02 A^2`, and the largest width difference is `3.466E-02 A`. This is a
small but resolved buffer dependence, not `1E-10` identity. Both runs retain
the unchanged outer-shell failure: B=6 ratio `7.2686E-03`, B=10 ratio
`2.2432E-02`, versus tolerance `1E-08`.
