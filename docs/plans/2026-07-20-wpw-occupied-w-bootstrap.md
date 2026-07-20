# Occupied-W WPW Bootstrap Implementation Plan

> **For Codex:** REQUIRED SUB-SKILL: Use `executing-plans` and
> `test-driven-development` to implement this plan task-by-task.

**Goal:** Replace the 165-row core-orthonormalized Flux-LCFO W representation
in the WPW production route with the core-owned occupied projected-Wannier
representation (16 W rows per fragment and 128 globally for Si64).

**Architecture:** Split the current mixed W/P volume accumulator into a
W-independent core-grid/P bootstrap followed by an occupied-W descriptor.
Construct deterministic projected Wannier IDs, core values, buffer values, and
gradients before creating the W owner schedule. Assemble every WW/WP metric and
Hamiltonian consumer from that descriptor; keep the legacy LCFO basis only for
the non-WPW route.

**Tech Stack:** Fortran 2008, MPI owner/halo schedules, existing SAWF projected
Wannier helpers, CMake, Python source-contract tests, focused MPI fixtures.

**Repository rule:** Preserve the dirty worktree and generated run directories.
Do not commit, push, or open a PR unless the user explicitly requests it.

---

### Task 1: Extract a W-independent core-grid/P bootstrap

**Files:**
- Modify: `src/gs/dc/dg_wpw_rank_local_quadrature.f90`
- Modify: `src/gs/dc/lcfo_flux.f90`
- Test: `tests/dg/test_dg_wpw_rank_local_quadrature.f90`
- Test: `tests/dg/check_dg_wpw_fixed_h_relaxation.py`

1. Add a RED fixture constructing a core-grid/P descriptor with zero W rows.
   Require positive point count, stable grid IDs, weights, density/potential,
   P values/gradients, and PP matrices.
2. Add `s_dg_wpw_core_p_accumulator` (or an equivalently named type) whose
   allocation and point insertion do not require `nw>0`.
3. Move grid IDs, weights, density/potential and P storage into this bootstrap.
   Keep the existing mixed accumulator API temporarily as an adapter so
   unrelated tests remain green.
4. Change projected-Wannier source construction to consume the bootstrap
   `npoint`, `grid_ids`, and `weights`, not `wpw_volume_accumulator%w_points`.
5. Run the focused quadrature fixture and source contract. Expected: PASS.

**Completed:** Added the transactional W-independent core/P accumulator and a
focused zero-W behavioral fixture. PP integration, grid IDs, weights,
density/potential storage, P values, and P gradients pass independently of W.
`lcfo_flux` now populates both accumulators during the transition, checks local
PP matrices and grid IDs for equality, and the projected-Wannier builder uses
only the W-independent bootstrap's `npoint` and `grid_ids`. The builder source
contract forbids references to the legacy mixed accumulator.

### Task 2: Publish a transactional occupied-W descriptor

**Files:**
- Create: `src/gs/dc/dg_wpw_occupied_w_basis.f90`
- Modify: `src/gs/dc/CMakeLists.txt`
- Modify: `src/gs/dc/lcfo_flux.f90`
- Create: `tests/dg/test_dg_wpw_occupied_w_basis.f90`
- Create: `tests/dg/run_dg_wpw_occupied_w_basis.py`

1. Add a RED fixture with two fragments and nontrivial bond/image keys.
   Require lexicographically deterministic global IDs, unique ownership,
   strictly increasing owner/support ID lists, and rollback on duplicate keys,
   nonfinite values, inconsistent shapes, or wrong local/global counts.
2. Define a descriptor carrying stable keys/IDs, owner fragment, core grid IDs,
   owner core values, rectangular buffered values, gradients, buffer bounds,
   source condition, epoch, fingerprint, and validity.
3. Extend `build_core_owned_projected_wannier_density_seed` into a builder that
   returns this descriptor before operator/context initialization. Reuse the
   existing polar transform and bond/image canonicalization; do not recompute a
   second Wannier ensemble later.
4. Require Si64 diagnostics `local_count=16` and `global_count=128` before W
   layout creation.
5. Run the new fixture, source contracts, and `git diff --check`.

**Progress:** Added `s_dg_wpw_occupied_w_basis` with lexicographically stable
global bond/image IDs, owner-count and local-order validation, core/buffer
values and gradients, epoch/source-condition metadata, and a content
fingerprint. The focused fixture proves deterministic IDs and transactional
rollback for duplicate global keys, inconsistent buffer shape, and invalid
source condition. Remaining Task 2 work is to make the projected-Wannier
builder populate this descriptor collectively and enforce the Si64 16/128
count oracle before W layout creation.

The descriptor now also has a transactional MPI collective initializer. It
gathers variable local key counts, assigns the same lexicographic global IDs
on every fragment owner, enforces an expected global count, and commits only
after all ranks validate. A two-rank fixture proves reversed key/owner order
and collective rollback on a global-count mismatch. Production builder wiring
remains pending because the current source values are spatially distributed;
they must be gathered or represented without duplicating stable keys before
calling the fragment-root collective.

The projected-Wannier reconstruction now exposes a validated reusable
`apply_sawf_projected_wannier_transform` operation. The core projection alone
computes the polar transform; buffer samples reuse that exact transform rather
than performing a second polar decomposition. The same operation is the
required path for the three occupied-state buffer-gradient components, so
core values, buffer values, and gradients cannot acquire different numerical
Wannier gauges. A focused MPI oracle covers successful transformation and
fail-closed shape validation. The full MPI build passes after this change.

### Task 3: Build occupied-W owner layout and halo directly from the descriptor

**Files:**
- Modify: `src/gs/dc/lcfo_flux.f90`
- Modify: `src/rt/dg/rt_dg_wpw_volume_halo_provider.f90` only if stable-key
  metadata cannot be retained by the existing integer-ID payload
- Test: `tests/dg/test_dg_wpw_core_wannier_seed_mpi.f90`
- Test: `tests/dg/test_dg_wpw_core_w_provider.f90`

1. Add RED MPI assertions that each Si64-like fragment owns only its projected
   occupied W rows and that a neighboring core receives each required tail
   exactly once with the correct periodic image.
2. Replace `build_complete_wpw_w_row_layout(...n_basis,index_basis...)` in the
   WPW route with descriptor IDs. Do not change the legacy LCFO layout.
3. Pack owner buffer values/gradients from the descriptor, not
   `basis_transform*spsi%rwf`.
4. Reuse the existing canonical halo exchange and support evaluator.
5. Require global W norms to be finite and bounded near the descriptor's
   independently integrated norms; remove Task 2H tracing only after this
   passes.

### Task 4: Assemble the metric from the actual occupied-W fields

**Files:**
- Modify: `src/gs/dc/dg_wpw_rank_local_quadrature.f90`
- Modify: `src/gs/dc/lcfo_flux.f90`
- Modify: `src/common/dg_wpw_production_context.f90`
- Modify: `src/common/dg_wpw_bounded_operator.f90`
- Test: `tests/dg/test_dg_wpw_production_operator_mpi.f90`

1. Add a RED MPI fixture with nonzero cross-fragment W overlap. Require WW,
   WP, and PP quadratic forms to equal independent real-space integration.
2. Accumulate support-W WW and WP overlap blocks on every core, route rows to
   canonical owners, and preserve Hermitian symmetry collectively.
3. Publish a metric convention describing the actual tail-carrying Gram.
   Remove `orthonormal_ww` from the WPW occupied-W route.
4. Keep the coupled W+P block preconditioner and fail-closed rank checks.
5. Run production-operator MPI tests and the metric/real-space identity tests.

### Task 5: Rebuild occupied-W Hamiltonian consumers

**Files:**
- Modify: `src/gs/dc/lcfo_flux.f90`
- Modify: projector/nonlocal helpers selected by the existing WPW route
- Test: `tests/dg/test_dg_wpw_production_operator_mpi.f90`
- Test: `tests/dg/check_dg_wpw_fixed_h_relaxation.py`

1. Add RED oracles for occupied-W kinetic/local-potential, nonlocal projector,
   canonical-face WW self/cross, and WP face terms.
2. Assemble these terms from descriptor core/buffer values and gradients.
   Do not import legacy `mat_H_*` arrays by positional LCFO basis index.
3. Preserve fixed-H H0/interface separation, frozen fingerprints, collective
   validation, and rollback.
4. Update dynamic potential refresh and checkpoint ID mapping to use stable
   occupied-W IDs.
5. Run focused operator, frozen-state, checkpoint, and rollback tests.

### Task 6: Reconnect the density-carrying solve and physical gates

**Files:**
- Modify: `src/gs/dc/lcfo_flux.f90`
- Modify: `docs/plans/2026-07-19-wpw-density-carrying-seed*.md`

1. Remove the second projected-Wannier construction from the seed; use the
   single bootstrap descriptor as both source ensemble and W basis provenance.
2. Solve the same coupled `S C=B` contract and retain diagnostic continuation
   only under its existing explicit nonpublishable flag.
3. Re-run focused tests, full MPI build, and a new Si64 run.
4. Require: 16 W/fragment, 128 global W, bounded W norms, routed/direct capture
   identity, assembled/real-space WW/WP/PP identity, physical charge/density
   gates, occupied rank, and frozen invariants.
5. Stop before fixed-H or publication on any failed gate. Reconsider the
   all-RHS `1e-10` target only after the physical route reaches fixed-H.
