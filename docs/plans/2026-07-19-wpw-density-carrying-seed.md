# WPW Density-Carrying Seed Implementation Plan

> **Historical/removed:** This document describes an obsolete experimental DG route
> removed on 2026-07-31. It is retained only as an implementation record and is
> not executable guidance.

> **For Codex:** REQUIRED SUB-SKILL: Use `executing-plans` to implement this plan task-by-task.

**Goal:** Build a qualified W+P fixed-H checkpoint from deterministic core-owned projected Wannier functions and their communicated tails without running the Flux-inclusive dense LCFO eigensolver.

**Architecture:** Construct deterministic bond-centered projected Wannier functions, select them by unique periodic core ownership, communicate their bounded tails, and qualify their collective reconstruction of the converged DC core density. Construct both W and P overlaps, solve `S C=B`, transform the non-diagonal occupation matrix through S-normalization, then solve the frozen H0 problem and continue the complete interface operator to one.

**Tech Stack:** Fortran 2008, MPI, BLAS/LAPACK, SALMON DC-LCFO/WPW bounded operators, Python source-contract tests.

**Repository constraints:** Work in the existing dirty worktree. Do not reset, checkout, recreate Tasks 0–10, commit, push, or open a PR. New runs must use dedicated directories and must not overwrite sample results.

---

### Task 0: Re-establish the protected workspace

**Files:**
- Read: `docs/plans/2026-07-19-wpw-density-carrying-seed-design.md`
- Read: `docs/plans/2026-07-19-wpw-density-carrying-seed-session-handoff.md`

1. Run `pwd`, `git branch --show-current`, `git rev-parse HEAD`, and
   `git status --short`.
2. Require branch `codex/singlescale-vortex-observables` and HEAD
   `87435fb88f8242c828753c931b2578e9dba6a47f`, or a documented descendant that
   contains commits `192f3af` and `87435fb`. Record the actual starting SHA.
3. Confirm that the worktree remains intentionally dirty. Do not clean it.
4. Run `git diff --check`.
5. Run `cmake --build build-mpi-eigenexa -j4`.
6. Run the focused tests listed in the handoff document and save their exact
   results before editing.

### Task 1: Repair and qualify the frozen-H snapshot

**Files:**
- Modify: `src/gs/dc/lcfo_flux.f90`
- Test: `tests/dg/check_dg_wpw_fixed_h_relaxation.py`
- Test: `tests/dg/run_dg_wpw_production_operator_mpi.py`

1. Extend the source contract so frozen validation requires WW values and IDs,
   WP/PP volume/nonlocal/face values, bounded H0/interface caches, layout IDs,
   and provenance.
2. Run the contract and confirm RED.
3. Audit the partially added `wpw_frozen_*` arrays from the previous session.
   Do not assume that grouped deallocation or deep assignment is transactional.
4. Add `stat=` allocation paths and collective failure synchronization. Ensure
   every root enters the same production-communicator collectives.
5. Make validation shape-safe: test allocation and shape before comparing
   array values.
6. Add MPI mutation oracles for WW projector-local and projector-cross values,
   one WP/PP nonlocal value, one H0 cache, one interface cache, sparse/cache
   index arrays, ownership/support IDs, and one routing/callback binding
   record. Each value, shape, or ID mutation must be detected collectively
   without deadlock.
7. Verify rollback leaves the last valid operator and callbacks intact.

### Task 2: Replace the false projection oracle with `S C=B`

**Files:**
- Modify: `src/gs/dc/dg_wpw_matrix_free_scf.f90`
- Modify: `tests/dg/test_dg_wpw_matrix_free_scf_mpi.f90`
- Test: `tests/dg/test_dg_wpw_matrix_free_scf.py`

1. Replace the current oracle that normalizes coefficient guesses.
2. Build a two-rank fixture with nonzero W and P rows and a known Hermitian
   positive-definite W+P metric containing WP coupling.
3. Choose known coefficients `C_ref`, form `B=S C_ref`, and require the new
   projection routine to recover the same S-metric occupied projector.
4. Rotate `C_ref` by an occupied unitary matrix and require projector
   invariance.
5. Add a rank-deficient RHS case and require fail-closed behavior.
6. Add a nonfinite RHS on one rank and require collective failure.
7. Use separate input/output buffers in callback tests; never alias an
   `intent(in)` and `intent(out)` actual argument.
8. Run the test and confirm RED.
9. Implement a bounded distributed metric solve for `S C=B`. Report the true
   relative equation residual and every RHS residual. A fixed iteration cap or
   collective stagnation condition must return failure, not a partially
   converged result.
10. Apply rank-revealing S-orthonormalization only after the solve.
11. Add positive diagonal Jacobi preconditioning from the locally owned
    diagonal of `S`. Reject nonfinite or nonpositive diagonal entries
    collectively; report global diagonal spread and per-RHS residual history.
    Forbid global dense fallback and all-state gathers. Do not introduce a
    block-local preconditioner unless a focused oracle first proves diagonal
    Jacobi insufficient and a separate design review approves the change.
12. Add a fixture in which one RHS converges before the others and require the
    remaining active RHS columns to continue without zero-denominator failure.
13. Run the MPI oracle and expect PASS.

### Task 2B: Add diagnostic-first fragment-local coupled W+P block Jacobi

**Files:**
- Modify: `src/common/dg_wpw_bounded_operator.f90`
- Modify: `src/gs/dc/dg_wpw_matrix_free_scf.f90`
- Modify: `src/gs/dc/lcfo_flux.f90`
- Modify: `tests/dg/test_dg_wpw_matrix_free_scf_mpi.f90`
- Modify: `tests/dg/test_dg_wpw_production_operator_mpi.f90`
- Modify: `tests/dg/check_dg_wpw_fixed_h_relaxation.py`

1. Preserve the Si64 evidence that diagonal Jacobi is insufficient: true and
   recursive residual agreement, Hermitian defect, 256-step residual history,
   and PCG/Lanczos Ritz estimate. Remove diagnostic-only extra work from the
   normal solver path after replacement qualification.
2. Write a failing two-rank MPI oracle with an SPD W+P metric made ill
   conditioned by fragment-local WP coupling, for which coupled block Jacobi
   converges within the fixed cap. Include nonzero W and P rows on both ranks.
3. Show in a comparison oracle that W/P-separated blocks fail or improve
   materially less than the coupled block.
4. Add fail-closed oracles for a non-Hermitian block, nonfinite values, a
   nonpositive eigenvalue, and an eigenvalue at or below the local relative
   cutoff. Verify collective failure and exact rollback.
5. Run the focused tests and confirm RED.
6. Build each principal `S_f` only from W/P IDs canonically owned by one
   fragment. Locate owned IDs in the required-ID arrays and extract WW from
   `ww_s_dense(required_w,required_w)`, WP from
   `wp_s_dense(required_w,owned_p)`, and PP from
   `pp_s_dense(owned_p,required_p)`. Check the assembled owned principal block
   for Hermitian consistency. Reuse bounded IDs and schedules to validate the
   layout, but forbid global dense assembly and all-state gathers.
7. Qualify `S_f` by Hermitian eigendecomposition. Record dimension, Hermitian
   defect, eigenvalue extrema, and condition number. Never silently truncate.
8. Store spectral factors and the inverse action transactionally in a distinct
   preconditioner object tagged with the bounded operator epoch and layout
   fingerprint. Build a candidate, synchronize validation, and replace the
   valid factors only on collective success. Failed rebuild leaves the last
   valid factorization unchanged; failed apply returns zero output and does not
   mutate factors, epochs, IDs, the bounded operator, or frozen values.
9. Apply the coupled inverse to the canonical-owned residual W/P rows. This is
   rank-local in the current one-fragment-owner-per-production-rank layout and
   must add no owner exchange. Synchronize validation failures before any rank
   enters the next `apply_s` collective. If a future layout owns several
   fragments on one rank, require one block descriptor per fragment.
10. Require all focused RHS to converge, including early-inactive columns, and
    verify the result by explicitly recomputing `B-SC`.
11. Run a dedicated Si64 preflight comparing diagonal and block Jacobi on the
    same seed. Require all RHS below `dg_wpw_metric_cutoff` within the cap,
    agreement with explicitly recomputed `B-SC`, and at least a one-decade
    reduction from the diagonal-Jacobi condition estimate (`1.06e5` baseline,
    block-Jacobi target at most `1.06e4`).
12. If step 11 fails, stop for a separate overlap-1 additive Schwarz design.
    Do not tune cutoff, rescale occupations, add a global dense fallback, or
    proceed to production.
13. Only after step 11 passes, resume Task 3 density/captured-norm and fixed-H
    gates. This task alone does not qualify the DG+TDDFT initial state.

**Recorded Si64 result:** The corrected full-cell shell check passed with zero
outer-shell norm. Coupled block Jacobi reduced the estimated condition number
from `1.06e5` to `2.31e3`, but after 256 iterations the aggregate and worst-RHS
residuals were respectively `1.41e-7` and `3.40e-7`, above the `1e-10` cutoff.
Step 11 therefore fails its strict convergence clause. Continue only through
Task 2C's explicitly nonpublishable physical diagnostic route before deciding
whether overlap-1 Schwarz is required. Do not change the iteration cap or
cutoff during that diagnostic.

### Task 2C: Continue a qualified best metric iterate into physical diagnostics

**Files:**
- Modify: `src/gs/dc/dg_wpw_matrix_free_scf.f90`
- Modify: `src/gs/dc/lcfo_flux.f90`
- Modify: `tests/dg/test_dg_wpw_matrix_free_scf_mpi.f90`
- Modify: `tests/dg/check_dg_wpw_checkpoint_publication.py`

1. Add a RED MPI oracle in which block-preconditioned PCG remains above its
   strict tolerance at the fixed cap while every RHS improves by at least one
   decade. Require the solver to return each RHS's best iterate, its explicit
   residual, and a `diagnostic_continuation` flag only when explicitly enabled.
2. Add RED cases proving that nonfinite values, nonpositive curvature,
   collective callback failure, or less than one-decade improvement still fail
   and return zero coefficients.
3. Keep the default solver behavior fail-closed so existing callers do not
   acquire diagnostic continuation implicitly.
4. Track the best coefficient column and best explicit residual independently
   for every RHS. At the cap, collectively qualify finiteness, positive
   curvature history, recursive/explicit agreement, and the one-decade
   reduction rule before returning the best block.
5. In the fixed-H route, record the diagnostic-only flag in run state and log
   the strict target plus returned per-RHS residual range.
6. Continue through captured norm, source/projected/normalized density and
   charge, S-rank/orthogonality, occupied projector, zero-interface fixed-H,
   and interface continuation without weakening any of those gates.
7. Make checkpoint publication fail before creating rank files or a manifest
   whenever the diagnostic-only flag is set. Do not start field-on or
   production RT from this route.
8. Run focused MPI tests, checkpoint source-contract tests, full MPI build, and
   the same Si64 diagnostic directory with a new unique run name.
9. Reconsider the `1e-10` metric target versus overlap-1 Schwarz only after the
   physical diagnostic results are recorded.

**Recorded Task 2C result:** Diagnostic continuation returned the per-RHS best
iterate with aggregate/worst residuals `1.39e-7`/`2.73e-7`. The subsequent
physical gate failed before fixed-H: captured norm `8.49e5`, projected-density
relative error `1.36e6`, projected charge `2.717e7` versus source charge 32,
and source-to-DC density residual `4.38e-2`. Normalization-density residual
`1.42e-14` localizes the failure before Löwdin normalization. Do not relax the
metric gate. Before coding overlap-1 Schwarz, add focused diagnostics that
separate unresolved original-metric extreme modes from an `S`/routed-`B`/
real-space reconstruction contract mismatch.

### Task 2D: Reconstruct physical W values from the full support coefficients

**Files:**
- Modify: `src/gs/dc/dg_wpw_core_w_provider.f90`
- Modify: `src/gs/dc/lcfo_flux.f90`
- Modify: `tests/dg/test_dg_wpw_core_w_provider.f90`
- Modify: `tests/dg/check_dg_wpw_fixed_h_relaxation.py`

1. Add a RED fixture with a nonzero coefficient on a halo-only W row. Require
   the reconstructed core value to contain that neighboring Wannier tail; an
   owned-W-only contraction must give a different result.
2. Add a source-contract RED check requiring the density-carrying seed to use
   the support-W reconstruction for both `C_raw` and the normalized occupied
   coefficients.
3. Reuse `evaluate_dg_wpw_core_w_support` at every core point, then contract
   all `wpw_support_w_ids` with `density_rw_support`. Do not extract only
   `wpw_owned_w_ids` into a smaller coefficient block.
4. Apply the same reconstruction to the raw physical diagnostics and the
   normalized density. Leave the already identity-consistent P reconstruction
   unchanged.
5. Rerun the focused provider fixture, source-contract test, full MPI build,
   and a new unique Si64 physical-diagnostic run. Require routed-W and direct-W
   capture to agree within numerical tolerance before interpreting the density
   and charge gates.
6. Continue on the current route through the physical gates. Reconsider the
   all-RHS `1e-10` criterion or overlap-1 Schwarz only after those physical
   results are recorded; this diagnostic-only route remains nonpublishable.

**Recorded diagnosis motivating Task 2D:** The routed and direct P captures
agree (`-1.2911e3`), while routed W is `8.5035e5` and direct W is `1.2937e3`.
Overlap assembly evaluates all support W rows, but physical reconstruction
currently discards nonowned support coefficients. This violates the `B` versus
real-space `A C` identity whenever a neighboring Wannier tail enters the core.

**Recorded Task 2D result:** The support-W reconstruction restores the intended
identity: routed/direct total capture are both `8.4905e5`, routed/direct W are
both `8.5035e5`, routed/direct P are both `-1.2911e3`, and the relative defect
is `1.10e-15`. The physical gate still fails: projected/normalized charge is
`2.6624e14` for source charge 32, while S-orthogonality is `5.25e-12` and the
normalization-density residual is `1.30e-15`. Therefore neither stricter PCG
convergence nor overlap-1 Schwarz addresses the current failure. Diagnose the
contract between the assembled metric `S` and the real-space Gram
`A^dagger A`, split into W-W, W-P, and P-P contributions, before changing the
preconditioner or entering fixed-H.

### Task 2E: Compare assembled S with the real-space Gram

**Files:**
- Modify: `src/gs/dc/lcfo_flux.f90`
- Modify: `tests/dg/check_dg_wpw_fixed_h_relaxation.py`

1. Add a RED source contract requiring an explicit application of assembled
   `S` to `C_raw` and a `[DG-WPW-METRIC-REALSPACE-GRAM]` diagnostic before
   Löwdin normalization.
2. Compute the occupation-weighted quadratic form `C_raw^dagger S C_raw` on
   canonical fragment owners using the production `wpw_apply_s` callback.
3. From the already reconstructed raw fields, compute the same real-space norm
   and split it into W-W, W-P (including both conjugate terms), and P-P.
4. Normalize all reported values by the same occupation-weighted source norm.
   Synchronize callback and finiteness failures before later collectives.
5. Run focused contracts, full MPI build, and a new unique Si64 diagnostic.
   If assembled S disagrees with the real-space total, inspect the dominant
   split component and its metric assembly; do not modify PCG.

**Recorded Task 2E result:** For the same `C_raw`, assembled S gives
`8.4905e5`, while the real-space Gram gives `8.3200e12`. The real-space split
is WW `8.3277e12`, WP `-7.7479e9`, and PP `4.6702e7`; WW overwhelmingly
dominates the discrepancy. The production adapter currently publishes the WW
metric under the `orthonormal_ww` convention. Before changing that convention,
apply assembled S to W-only and P-only copies of the same `C_raw` to obtain an
assembled WW/WP/PP split, and probe real-space W norms/overlaps including
cross-fragment tails. This must distinguish an invalid orthonormal-W metric
assumption from a scaling or duplication error in the W value/halo path.

### Task 2F: Split assembled S and probe global W norms

**Files:**
- Modify: `src/gs/dc/lcfo_flux.f90`
- Modify: `tests/dg/check_dg_wpw_fixed_h_relaxation.py`

1. Add RED contracts for `[DG-WPW-ASSEMBLED-METRIC-SPLIT]` and
   `[DG-WPW-W-REALSPACE-NORM]` diagnostics.
2. Apply production S separately to `(C_W,0)` and `(0,C_P)`. Report the
   occupation-weighted assembled WW, combined Hermitian WP/PW, and PP terms
   with the same source-norm denominator used by Task 2E.
3. During full-core support-W evaluation, integrate `|W_i|^2` for every
   support row, reduce spatial/orbital partials, and route each scalar to its
   canonical W owner using the existing owner schedule.
4. Report global min/max/mean owned-W norm and require finite positive values.
   Do not alter S or the halo path in this diagnostic task.
5. If W diagonals are near one but assembled/real-space WW disagree, add a
   focused cross-fragment off-diagonal overlap diagnostic before replacing
   `orthonormal_ww`. If W diagonals are not near one, localize scaling or
   duplicate integration first.

**Recorded Task 2F result:** Global W norms range from `1.0024` to
`2.5096e6` with mean `5.9761e4` over 1320 canonical W IDs, proving a severe W
value/tail scaling error rather than a small neglected cross-fragment overlap.
The assembled split is WW `4.7554e7`, WP `-9.3406e7`, PP `4.6702e7`; PP agrees
exactly with real-space PP `4.6702e7`, while real-space WW is `8.3277e12` and
WP is `-7.7479e9`. Continue only by separating each W norm into owner-core and
halo-tail contributions.

### Task 2G: Localize abnormal W norms to owner values or halo tails

1. Extend the routed W norm diagnostic with total, owner-core, and halo-tail
   columns for every canonical W ID.
2. On each canonical owner report the maximum-norm W ID and its two component
   norms. Require `total = owner_core + halo_tail` within roundoff.
3. If halo-tail dominates, inspect periodic image mapping and buffer packing
   for that ID. If owner-core dominates, inspect the underlying
   `local_basis_value` normalization/index mapping. Do not modify the metric.

**Recorded Task 2G result:** Every fragment's worst row is its last local W
ID (`165, 330, ..., 1320`). Each has owner-core norm `1.0000` and halo-tail
norm `2.5096e6`; the total identity holds. Thus core normalization and owner
routing are correct, while the abnormal norm is created entirely by the
buffer/tail values before or during halo packing. Inspect whether the active
path is `sawf_explicit_buffer` or the `basis_transform*spsi%rwf` continuation,
then record the maximum tail value with source local index, buffer coordinate,
target fragment, and periodic image before proposing a fix. In particular,
check whether a core-only near-null closure direction is normalized and then
applied without a buffer-norm bound.

### Task 2H: Trace the maximum W tail through halo packing

1. Add a RED contract for `[DG-WPW-W-HALO-PACK-MAX]`.
2. Without changing payloads, scan every packed send on each source fragment
   and report the largest absolute W value, W ID, canonical grid coordinate,
   destination rank/fragment, route image, periodic image, and active source
   path (`explicit_buffer` or `transformed_spsi`).
3. Verify the packed maximum equals the corresponding pre-pack buffer value.
4. Use the result to inspect the exact buffer-basis column and transform; do
   not introduce a clamp or discard the last W row.

**Recorded Task 2H result:** All ranks use `transformed_spsi`, not
`explicit_buffer`. The largest packed amplitudes are `2.91e2`--`2.97e2` at
buffer coordinates outside the owner core. Every packed value equals its
pre-pack buffer value exactly (`pack_defect=0`), including nontrivial route and
periodic images. Halo packing and image canonicalization therefore transmit,
but do not create, the blow-up. The active construction applies
`basis_transform(:,ibasis)` obtained by core-only overlap diagonalization and
Gram-Schmidt to buffered `spsi%rwf`; cancellation/normalization guaranteed on
the core is not bounded on the buffer.

### Task 2I: Restore the fragment occupied-W basis contract

The present WPW row layout incorrectly reuses all 165 retained Flux-LCFO core
basis directions per fragment as W rows. This conflicts with the approved
fragment contract: W rows are the core-owned occupied projected Wannier
functions (16 per fragment for Si64), while P supplies the augmentation and
extra-state capacity. The pathological rows are precisely high-index
core-orthonormalized LCFO directions whose buffer continuation is unbounded.

1. Before implementation, trace every consumer of `wpw_owned_w_ids`, WW
   Hamiltonian components, face traces, and W halo values and separate the
   occupied-W representation from the legacy 165-row LCFO basis.
2. Define stable W IDs from the already constructed bond/image Wannier IDs and
   require 16 owner rows per Si64 fragment and 128 globally.
3. Build owner-core and communicated buffer values directly from those
   projected Wannier functions, not from `basis_transform*spsi%rwf` columns.
4. Assemble the WW metric from the actual tail-carrying W fields; do not label
   it `orthonormal_ww` unless the global core-tiling Gram proves that property.
5. Preserve P rows and the W+P coupled solve. Re-run W norm, routed/direct,
   metric/real-space Gram, density, and charge gates before fixed-H.
6. Do not repair the legacy 165-row continuation with a tail clamp, an
   arbitrary eigenvalue cutoff, or deletion of only the last row.

**Task 2I dependency review:** The current order is `calc_basis` -> legacy W
row layout/halo -> W/P volume accumulator -> WW/WP/PP operator -> fixed-H seed;
the occupied projected-Wannier ensemble is constructed only inside that final
seed. Therefore changing the row count in `build_wpw_w_row_layout` alone would
misindex all downstream data. Refactor the bootstrap in this order:

1. Build support geometry, canonical faces, G modes, and a W-independent core
   grid/P quadrature descriptor.
2. Construct the core-owned occupied projected Wannier functions and their
   stable IDs, including buffer tails, before W row-layout/context creation.
3. Initialize W owner/support layouts from those IDs and publish their core,
   halo, and trace values.
4. Assemble WW/WP/PP metric and H components from the new W representation,
   then bind callbacks and enter the existing density-carrying solve.

Affected consumer groups are projector overlaps/nonlocal WW, bounded-operator
WW/WP rows, fixed-H snapshots and checkpoint IDs, density/potential
reconstruction, canonical face traces, and dynamic potential updates. Each
must receive the same stable occupied-W IDs. The legacy `calc_basis` may remain
for the non-WPW LCFO path but must not define WPW W rows.

### Task 3: Build both W and P overlap blocks from the DC ensemble

**Files:**
- Modify: `src/gs/dc/lcfo_flux.f90`
- Modify or extract pure helpers from: `src/gs/dc/lcfo_wannier_sawf_seed.f90`
- Reuse bond-center construction from: `src/gs/dc/lcfo_flux.f90`
- Test: `tests/dg/check_dg_wpw_fixed_h_relaxation.py`
- Create: `tests/dg/test_dg_wpw_core_wannier_seed_mpi.f90`
- Create: `tests/dg/run_dg_wpw_core_wannier_seed_mpi.py`

1. Write RED pure oracles for canonical fractional wrapping, exact-upper-bound
   wrapping, boundary values at `boundary +/- center_tolerance`, face/edge/
   corner ownership, translations by multiple cells, canonical image triples,
   stable bond IDs, and ID collision rejection.
2. Write a RED projected-Wannier oracle using deterministic bond-center/SAWF
   trials. Require full trial rank, bounded Gram condition, invariant density
   under rotation of the input occupied eigenvectors followed by reconstruction,
   and deterministic phase/permutation under repeated execution.
3. Write a RED cross-fragment fixture whose selected source Gram is nonidentity.
   Require the test to show that diagonal occupations after Löwdin
   normalization change density, while transformed
   `F_Q=T^{-1} F_src T^{-dagger}` preserves it.
4. Write RED halo fixtures for a required cross-core tail, missing and duplicate
   records, wrong periodic image, asymmetric send/receive schedules, multiple
   periodic images, and an insufficient-buffer case. The independent reference
   must evaluate origin-buffer values before halo packing.
5. Write a RED routed-overlap fixture with nonzero W and P rows in which a row
   owned on rank 0 requires a Wannier tail supplied by rank 1. Omitting either
   W, P, or the tail must fail before `S C=B`. Include a full-rank but
   deliberately incomplete W+P basis: `S C=B` must converge, while the
   source-to-projected physical density and captured-deficit gates fail.
6. Run `python3 tests/dg/run_dg_wpw_core_wannier_seed_mpi.py`; require RED for
   missing production helpers, not compilation errors in the fixture.
7. Extract or reuse pure bond-center/SAWF trial construction. Project the
   trials into the converged fragment occupied `spsi` subspace and apply the
   polar/Löwdin factor. Do not invoke global `diag_eigenexa`, import a global
   `coef_wf`, or launch an external iterative Wannier90 optimization.
8. Implement canonical center ownership and stable source IDs from ordered bond
   endpoint atom IDs plus canonical image triple. Apply the deterministic phase
   and permutation rules from the design.
9. Select core-owned centers. For Si64 require the stoichiometric oracle of 16
   centers per fragment and 128 globally; independently require reconstructed
   global charge 256. Permit nonuniform local counts in general fixtures.
10. Implement bounded tail discovery over every periodic image represented in
    the fragment buffer. Apply `dg_wpw_wannier_tail_norm_tolerance`, accumulate
    intentionally omitted available-tail norm and charge bounds, and measure
    the outer-buffer-shell norm. Fail if the shell norm or independent DC
    density reconstruction cannot qualify finite-buffer sufficiency.
11. Build the symmetric halo schedule keyed by source ID, destination core,
    image triple, and point/support ID. Validate counts and uniqueness before
    exchanging values; synchronize every allocation/validation failure before
    any peer enters communication.
12. Reconstruct core density from local functions and received tails. Compare
    it with both the pre-pack origin-buffer reference and the converged DC
    density. Record relative L2 error, charge error, omitted-tail bounds, image
    range, record/value counts, and peak storage.
13. Accumulate and route both `W^dagger phi` and `P^dagger phi` by stable source
    and basis-row IDs, then solve `S C_raw=B`. Record the raw metric residual and
    captured norm before normalization. Reconstruct `density[A C_raw,F_src]`
    on the global core tiling and require its relative L2 error, charge error,
    and `abs(1-captured_norm)` to satisfy
    `dg_wpw_wannier_density_tolerance` against both source and DC densities.
14. Form `G_C=C_raw^dagger S C_raw`, compute the rank-revealing transform `T`,
    set `Q=C_raw T`, and carry the Hermitian occupation matrix
    `F_Q=T^{-1}F_src T^{-dagger}`. Verify the density identity before fixed-H
    relaxation by reconstructing `density[A Q,F_Q]` and comparing it with the
    accepted `density[A C_raw,F_src]`. Reject rank loss, nonfinite values,
    non-Hermiticity, excessive condition number, or density/charge mismatch.
15. Restrict per-Wannier identity tests to phase, permutation, and periodic
    image equivalence. For an arbitrary occupied eigenvector rotation, rerun
    projected-Wannier construction and compare the stable center-ID set,
    density matrix, captured norm, and S-projector rather than individual
    columns.
16. Rename output and checkpoint provenance to
    `core_owned_projected_wannier_density_seed` and include localization method,
    center/phase/ID conventions, source-Gram condition, occupation-matrix
    identity residual, source-to-W+P density and charge residuals,
    captured-deficit gate, tail tolerance, omitted-tail bounds, and halo-layout
    fingerprint.
17. Add the named tolerance profile from the design and RED sensitivity oracles
    at one decade above and below the center/tail defaults. Do not accept a
    tolerance choice that changes ownership, rank, or qualified density.
18. Run `python3 tests/dg/run_dg_wpw_core_wannier_seed_mpi.py` and the existing
    fixed-H source contract; require GREEN.

### Task 4: Complete deterministic extra states without polluting occupancy

**Files:**
- Modify: `src/gs/dc/dg_wpw_matrix_free_scf.f90`
- Modify: `src/gs/dc/lcfo_flux.f90`
- Test: `tests/dg/test_dg_wpw_matrix_free_scf_mpi.f90`

1. Add a failing test with an extra seed containing a deliberate occupied
   component.
2. Require `Q_occ^dagger S Q_extra=0` after completion.
3. Require occupied projector invariance before and after extra completion.
4. Rank-qualify only the extra block after occupied projection.
5. Do not run a final full-block normalization that can rotate occupied and
   extra columns together unless the oracle proves occupied-projector
   invariance.
6. Reject insufficient extra rank; do not replace occupied columns randomly.

### Task 5: Implement zero-interface fixed-H relaxation

**Files:**
- Modify: `src/gs/dc/lcfo_flux.f90`
- Modify: `src/gs/dc/dg_wpw_matrix_free_scf.f90`
- Test: `tests/dg/test_dg_wpw_matrix_free_scf_mpi.f90`

1. Set the bounded operator to `lambda_interface=0` transactionally.
2. Validate the frozen state before and after every solver stage.
3. Start block CG/LOBPCG from the qualified density-carrying seed.
4. Keep density, Hartree, XC, local potential, and all H0 blocks fixed.
5. Forbid every call to `wpw_potential_step` in this route.
6. Use occupied trace, maximum generalized residual, S-orthogonality, and
   occupied-projector change as acceptance metrics.
7. Reconstruct density only for diagnostics and establish the zero-interface
   representation baseline.
8. Fail closed at the iteration cap.

### Task 6: Continue the complete interface operator to one

**Files:**
- Modify: `src/common/dg_wpw_bounded_operator.f90`
- Modify: `src/gs/dc/lcfo_flux.f90`
- Test: `tests/dg/test_dg_wpw_production_operator_mpi.f90`

1. Preserve the existing λ=0, 0.5, 1 dense oracle.
2. Add rejected-step and rollback cases.
3. Snapshot accepted coefficients, lambda, solver merit, and frozen invariant
   state before every trial.
4. On rejection, restore the last accepted state and reduce the lambda step.
5. Scale only `ww_face_self`, `ww_cross_value`, and `wp_h_face`.
6. Prove H0 caches and transport IDs are bitwise unchanged across all trials.
7. Require convergence at `lambda_interface=1`.

### Task 7: Qualify and publish the checkpoint

**Files:**
- Modify: `src/common/dg_wpw_checkpoint.f90`
- Modify: `src/gs/dc/lcfo_flux.f90`
- Modify: `tests/dg/check_dg_wpw_checkpoint_publication.py`

1. Extend checkpoint metadata with fixed-H mode, projected-Wannier method and
   stable-ID convention, center/tail tolerances, source-Gram condition,
   transformed-occupation density-identity residual, source-to-W+P projected
   density/charge residuals, captured-deficit acceptance, omitted-tail norm/charge,
   halo-layout fingerprint, metric residual, captured norm, projection rank,
   projection charge, final interface lambda, tolerance profile, and frozen
   fingerprints.
2. Confirm the publication contract is RED.
3. Publish only after lambda one, converged fixed-H residuals, bounded density
   diagnostics, and exact frozen-state validation.
4. Ensure failure leaves no manifest that could be mistaken for a valid
   checkpoint.
5. Run checkpoint read-back and identity validation.

### Task 8: Review and preflight before Si64

1. Run all focused tests listed in the handoff document.
2. Run `cmake --build build-mpi-eigenexa -j4` and `git diff --check`.
3. Perform a findings-first review focused on collective ordering, metric
   projection mathematics, support-row overlap routing, raw-versus-normalized
   diagnostics, frozen-state coverage, rollback, and provenance.
4. Fix every P0/P1 finding and rerun the focused suite.
5. Run only a short low-cost preflight. Confirm MPI ranks, OpenMP threads,
   estimated memory, output path, and free capacity first.
6. Require the preflight to report per-fragment core electron count,
   core-owned Wannier count, global source count, unique-center ownership, tail
   halo record counts, image range, omitted-tail bounds, source-Gram condition,
   transformed-occupation density identity, peak halo storage, and
   communication-enabled source and W+P-projected density reconstruction errors.
7. Do not start the full Si64 production checkpoint until the preflight proves
   that every metric RHS converges without stagnation, routed overlap assembly
   passes its identity checks, projection and zero-interface solve converge,
   and at least one continuation trial finishes with finite diagnostics.

### Task 9: Resume the production sequence

After the fixed-H checkpoint passes review:

1. Run the Si64 2x2x2 DG-DC GS/checkpoint calculation in a dedicated run
   directory.
2. Stop if GS convergence, charge, energy, provenance, or checkpoint identity
   fails.
3. Run checkpoint-backed field-off RT and review finite observables, midpoint
   convergence, S-norm/charge drift, and field-off Pz/Jz drift.
4. Only after field-off acceptance, run the 20-step Sin2 laser smoke from the
   same checkpoint.
5. Do not run or interpret long-time HHG.
