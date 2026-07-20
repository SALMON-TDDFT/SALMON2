# New Session Instructions and Review Points

## Objective

Continue the fixed-H WPW checkpoint work from the existing dirty worktree.
Repair Task 3 so that the initial state is a true nonorthogonal W+P projection
of deterministic core-owned projected Wannier functions and communicated
tails, then implement zero-interface
relaxation and interface continuation. Do not begin long Si64 production until
the focused implementation and review gates pass.

## Workspace identity

```text
worktree: /Users/otobetoshihito/SALMON-dev/SALMON2_RTDG
branch:   codex/singlescale-vortex-observables
HEAD:     87435fb88f8242c828753c931b2578e9dba6a47f
```

This is a local worktree. The unrelated OneDrive workspace shown in some app
environment metadata is not the worktree for this task.

## Non-negotiable constraints

- Preserve all existing uncommitted Task 5–10 and review-fix changes.
- Do not use checkout, reset, clean, or destructive restoration.
- Do not redo Tasks 0–10 from the beginning.
- Do not commit, push, or open a PR.
- Do not overwrite historical sample results.
- Use dedicated directories for every run.
- Do not proceed from a failed or provenance-ambiguous checkpoint.
- Do not run or interpret long-time HHG.

## Documents to read first

1. `docs/plans/2026-07-19-wpw-density-carrying-seed-design.md`
2. `docs/plans/2026-07-19-wpw-density-carrying-seed.md`
3. `docs/plans/2026-07-19-wpw-fixed-hamiltonian-basis-flux-relaxation-design.md`
4. `docs/plans/2026-07-19-wpw-fixed-hamiltonian-basis-flux-relaxation.md`
5. The earlier Task 5–10 remediation and ownership/adapter plans listed in the
   original session handoff.

## State inherited from the previous session

Implemented and previously passing:

- explicit opt-in `yn_dg_wpw_fixed_h_relaxation` control;
- fixed-H snapshot/validation route, initially fail-closed;
- H0/interface decomposition in the bounded operator;
- `lambda_interface` cache rebuild at λ=0, 0.5, and 1;
- dense MPI oracle proving the interface mapping and unchanged transport IDs;
- an S-orthonormalization helper and deterministic extra-state completion.

Important: Task 3 was reopened again after the Si64 preflight. The current
`build_wpw_density_carrying_fragment_seed` counts every occupied eigenvector in
each buffered fragment, producing 1024 source columns instead of the required
128, and fails collectively before projection. Replace it with the
core-owned projected-Wannier construction, tail halo, W/P overlap assembly,
transformed occupation matrix, and `S C=B` route in the revised Task 3.

The previous session also began strengthening frozen-value snapshots. Treat
those newest edits as unverified until compilation, focused tests, and review
are rerun.

## First commands

```bash
pwd
git branch --show-current
git rev-parse HEAD
git status --short
git diff --check
cmake --build build-mpi-eigenexa -j4
python3 tests/dg/check_dg_wpw_fixed_h_relaxation.py
python3 tests/dg/test_dg_wpw_matrix_free_scf.py
python3 tests/dg/run_dg_wpw_production_operator_mpi.py
python3 tests/dg/run_dg_wpw_gs_bounded_apply_mpi.py
python3 tests/dg/check_dg_wpw_bounded_index_cache.py
```

If a command fails, inspect the rank-local diagnostic that first failed. Do not
infer a method failure from only the reduced aggregate result.

## Required corrections before proceeding

### P1: false nonorthogonal projection

Current behavior:

```text
q_W = W^dagger phi
q_P = 0
normalize q in S
```

Required behavior:

```text
b_W = W^dagger phi
b_P = P^dagger phi
solve S c = (b_W,b_P)
rank-qualify and S-normalize c
```

The MPI oracle must use nonzero W/P coupling so the old implementation fails.

### P1: incomplete distributed overlap right-hand side

`B` is owned by W/P basis-row ID, not by source fragment. Each fragment core
must evaluate every support W/P row that is nonzero there, then route partial
overlaps to the canonical row owner. A local fragment root may not populate
only its locally owned rows and columns. Add a two-fragment oracle where a W
row owned on one rank requires a nonzero overlap contribution from the other
fragment core.

Compute captured norm from the raw metric solution before occupied
S-orthonormalization:

```text
C_raw = solve(S,B)
captured_norm = Tr(F B^dagger C_raw)
source_norm = Tr(F Phi_src^dagger Phi_src)
```

Here `F` contains the converged source occupations. Do not compute this
diagnostic with the normalized occupied block. A gauge oracle rotates the input
occupied fragment eigenvectors and repeats deterministic projected-Wannier
construction; it compares ensemble density, stable center-ID set, and projector,
not individual Wannier columns.

### P1: nonorthogonal source density across normalization

Projected Wannier functions from independent buffered fragments are not
assumed mutually orthogonal. If `Q=C_raw T` and
`T^dagger(C_raw^dagger S C_raw)T=I`, carry
`F_Q=T^{-1}F_src T^{-dagger}` and verify
`C_raw F_src C_raw^dagger = Q F_Q Q^dagger`. Replacing `F_Q` by a diagonal
occupation vector before this identity check is forbidden. Report source-Gram
condition and density-identity residual collectively.

Separately qualify projection loss: reconstruct `density[A C_raw,F_src]` on the
global core tiling and compare it with both the selected source density and the
converged DC density. Require relative L2 error, charge error, and
`abs(1-captured_norm)` below `dg_wpw_wannier_density_tolerance`. Repeat with
`density[A Q,F_Q]`. A converged metric equation and full rank do not waive this
publication gate.

### P1: source selection, naming, and provenance

The source is the direct sum of core-owned projected Wannier functions. For the
Si64 route, construct them deterministically by projecting existing
bond-center/SAWF trials into the converged fragment occupied subspace and
applying the polar/Löwdin factor. It is
not the direct sum of every occupied eigenvector in each buffered fragment
calculation and is not the global Flux-LCFO eigenvector set `coef_wf`. Rotate
each converged fragment occupied subspace to Wannier functions, verify density
invariance, wrap their centers periodically, and select centers with the
half-open core ownership rule `[lower,upper)` in every axis. Do not select the
first `N` fragment eigenvectors by energy order.

For non-spin-polarized Si64 2x2x2, require 32 core electrons and 16 doubly
occupied Wannier functions per fragment, hence 128 source columns and 256
electrons globally. The current `count(system%rocc>0)` result of 128 columns per
fragment and 1024 columns globally is a failing oracle, not a valid source.
Every selected stable bond-center ID must have exactly one fragment owner. Do
not run `diag_eigenexa`, import global `coef_wf`, or launch an external
iterative Wannier90 optimization merely to generate this seed. Record the
projected-Wannier method, bond/image ID, phase convention, ownership rule, and
source provenance in diagnostics and checkpoint metadata.

Do not infer that center selection preserves density merely because the input
occupied projector is gauge-invariant. Communicate each selected Wannier tail to every
neighboring core intersected by its buffered support, then reconstruct the
selected ensemble on the global core tiling and compare it with the converged
DC density. Halo records must identify stable source-Wannier ID, source and
destination fragments, periodic image, and destination support. Missing,
duplicate, wrong-image, or asymmetrically scheduled records fail collectively.
Route overlaps separately by stable source ID and basis-row ID. Representation
failure is a stop condition and occupations must not be rescaled to conceal it.
Enumerate every periodic image represented by the fragment buffer, bound
deliberately discarded available tails using the named tail tolerance, measure
the outer-buffer-shell norm, and stop unless the independent DC density
comparison qualifies the finite buffer. Run the specified one-decade tolerance
sensitivity oracle; do not tune tolerances solely to make Si64 pass.

### P1: incomplete frozen invariant

Compare actual values and shapes for every frozen component and bounded cache.
A stored fingerprint alone is insufficient because an accidental mutation may
not recompute it. Include WP/PP nonlocal and face contributions and transport
IDs.

### P2: allocation and oracle safety

- Add `stat=` and collective failure handling to new large allocations.
- Avoid grouped deallocation unless all members are known allocated.
- Avoid derived-type deep assignment as the only transactional mechanism.
- Remove callback argument aliasing in the MPI projection test.
- Test both nonzero W and P rows.
- Report per-RHS metric residuals, detect collective stagnation, and exercise a
  case where one RHS converges before the remaining active columns.
- Require positive diagonal Jacobi preconditioning from owned `S` diagonal
  entries before Si64. Report global diagonal spread and per-RHS residual
  history; do not use a global dense fallback. The focused failure and separate
  design review are now complete: Si64 gives a preconditioned Ritz condition
  estimate near `1.06e5`, no cutoff-level near-null modes, and no convergence
  after 256 PCG iterations. Implement Task 2B's canonical-fragment coupled W+P
  block Jacobi as a diagnostic first stage. Extract the owned principal block
  using the actual required-W/owned-P/required-P dense-store shapes; the
  preconditioner action itself is local on canonical-owned Krylov rows and must
  not add an owner exchange. Qualify every local block spectrally and fail
  rather than silently truncate. Require every RHS to converge and the
  estimated condition number to improve by at least one decade. If Si64
  improvement is insufficient, stop for a separate overlap-1 additive Schwarz
  design review.

The dedicated `20260720_block_jacobi_preflight` reached the strict solver
target stop condition.
The corrected full-cell buffer shell had zero norm. Coupled block Jacobi
reduced the estimated condition number from `1.06e5` to `2.31e3`, but the
256-step aggregate/worst-RHS residuals were `1.41e-7`/`3.40e-7`, above the
`1e-10` cutoff. Recursive and explicit residuals agreed and no near-null mode
was found. Implement Task 2C: return per-RHS best iterates only under an
explicit diagnostic-continuation option, run the physical representation and
fixed-H diagnostics, and prohibit all checkpoint/manifest publication. Do not
silently increase the cap or relax the cutoff. Decide whether overlap-1
Schwarz is required only after those physical diagnostics are reviewed.

Task 2C physical diagnostics are now available from
`20260720_metric_physical_diagnostic2`. The best metric iterate had aggregate/
worst residuals `1.39e-7`/`2.73e-7`, but captured norm was `8.49e5`, projected
density error `1.36e6`, and projected charge `2.717e7` versus source charge 32.
The normalization-density residual was `1.42e-14`, so the failure precedes
Löwdin normalization. Stop before fixed-H. Diagnose original-metric extreme
mode amplification versus an `S`/routed-`B`/real-space reconstruction mismatch
before assuming overlap-1 Schwarz alone is sufficient.

Task 2D corrected that reconstruction mismatch in
`20260720_support_w_reconstruction`: raw and normalized W fields now contract
the full support coefficient block, including neighboring Wannier tails.
Routed/direct total, W, and P captures agree with relative defect `1.10e-15`.
The common value is nevertheless unphysical (`8.4905e5`), and the
projected/normalized charge is `2.6624e14` versus source charge 32 despite
S-orthogonality `5.25e-12`. The next diagnostic must compare assembled `S`
against the real-space `A^dagger A` Gram and split W-W, W-P, and P-P. Do not
return to PCG tuning or overlap-1 Schwarz until this metric/operator contract
is resolved.

Task 2E (`20260720_metric_realspace_gram`) applied the production S explicitly
to `C_raw`. The normalized quadratic forms were assembled S `8.4905e5` versus
real-space total `8.3200e12`, split as WW `8.3277e12`, WP `-7.7479e9`, PP
`4.6702e7`. The discrepancy is therefore WW-dominated. The WW adapter labels
its metric `orthonormal_ww`; verify that assumption against the actual
tail-carrying W functions with an assembled block split and real-space W
norm/overlap probes before implementing a replacement metric.

Task 2G localized the W norm blow-up completely to halo tails. On all eight
fragments the worst W is the last local row (`165, 330, ..., 1320`), with
owner-core norm exactly `1.0000` and halo-tail norm `2.5096e6`. The next step
is to identify the active buffer-value path and log the worst source index,
buffer coordinate, destination/image, and pre-pack value. A likely mechanism
to test is amplification of a core-near-null closure direction when its
core-orthonormalizing transform is applied to the buffer, but this is not yet
proven and must not be fixed by an arbitrary tail cutoff.

Task 2H proves packing is innocent. The active path is `transformed_spsi`;
pre-pack and packed maxima are identical (`2.91e2`--`2.97e2`, defect zero) at
outside-core buffer coordinates across valid periodic images. The blow-up is
already present when a core-only `basis_transform` is applied to buffered
fragment states. Return to the approved occupied-W contract: Si64 must expose
16 core-owned projected Wannier W rows per fragment, not all 165 retained
Flux-LCFO core directions. Plan the representation change across WW, traces,
halos, and IDs before editing production code.

## Findings-first review checklist

Review in this order:

1. **Projection mathematics:** Does the code construct both overlap blocks and
   route every support-row contribution to the canonical owner before solving
   `S C=B`? Is the reported residual the equation residual for `C_raw`?
2. **Source provenance:** Are the orbitals exactly the core-owned occupied
   projected-Wannier ensemble? Is projected construction deterministic and
   invariant to input occupied gauge? Does every wrapped
   center have exactly one half-open core owner? Are local/global counts and
   charge consistent collectively, including Si64's 16x8=128 contract? Does
   the selected ensemble, including explicitly communicated cross-core tails,
   reproduce the DC core density independently of the input-gauge
   check? Is every tail received exactly once with the correct periodic image?
   Are omitted-tail norm/charge and finite-buffer sufficiency proved?
3. **Occupied projector:** Is occupied rank preserved? Is the S-projector
   invariant after deterministic reconstruction from a rotated input gauge and
   after extra-state completion? Is captured
   norm computed from `C_raw`, not the normalized block? Is the non-diagonal
   occupation matrix transformed through S-normalization with a verified
   density identity? Does the physical W+P-projected density match source and
   DC densities before and after normalization?
4. **Frozen H0:** Can any density, potential, H0 value, cache, or ID change
   without detection?
5. **Interface map:** Are only WW self/cross face and WP face terms scaled?
6. **MPI ordering:** Can one rank return before peers enter a collective? Are
   allocation and validation failures synchronized before communication?
7. **Transactional behavior:** Do failed solves, lambda trials, publications,
   and callback bindings restore the last valid state?
8. **Finite and shape checks:** Are allocation and shape tested before array
   comparison or indexing? Are all derived diagnostics finite?
9. **Performance:** Are operations rank-local plus bounded owner/halo exchange?
   Is any global dense matrix or all-state gather reintroduced?
10. **Provenance/publication:** Can an incomplete or diagnostic-only run leave a
    checkpoint that RT might accept?

## Stop conditions

Stop and report before Si64 production if:

- the metric projection needs an unbounded global dense solve;
- any support-row overlap contribution is missing or ownership-ambiguous;
- captured norm is not computed from the raw metric solution;
- any metric RHS stagnates or lacks a bounded preconditioner for Si64;
- occupied rank is lost;
- source occupation count/charge is inconsistent;
- any Wannier center is multiply owned, unowned, or boundary-ambiguous;
- projected-Wannier construction is not invariant to the input occupied gauge
  within tolerance;
- the selected ensemble does not reproduce the DC density on the global core
  tiling within the representation tolerance;
- any required neighboring Wannier tail is missing, duplicated, assigned the
  wrong periodic image, or uses an asymmetric halo schedule;
- omitted tail norm/charge exceeds tolerance or the finite buffer cannot bound
  the missing support;
- source-Gram conditioning or transformed-occupation density identity fails;
- source-to-W+P projected density/charge or captured-norm deficit exceeds the
  named density tolerance;
- frozen-state mutation is not detected collectively;
- zero-interface convergence is not demonstrated in a small preflight;
- interface rollback is not exact;
- memory/rank/thread/storage estimates are missing;
- the checkpoint manifest does not distinguish the density-carrying seed from
  global LCFO eigenvectors.

## Expected handoff after implementation

Report:

- files changed in the new session;
- exact focused tests and build result;
- projection residual, captured norm, rank, S-orthogonality, and charge from
  the preflight;
- per-fragment core electron count, core-owned Wannier count, and collective
  unique-center ownership result;
- Wannier-tail halo send/receive counts, periodic-image validation, duplicate
  detection, image range, omitted-tail bounds, peak storage, and density
  reconstruction error with communication enabled;
- source-Gram condition and transformed-occupation density-identity residual;
- source-to-W+P projected density and charge residuals before and after
  normalization, plus captured-norm deficit;
- zero-interface solver iterations and residuals;
- accepted/rejected continuation steps and final lambda;
- frozen-state validation result;
- findings-first review, including “no P0/P1 findings” only if supported by a
  fresh independent review.
