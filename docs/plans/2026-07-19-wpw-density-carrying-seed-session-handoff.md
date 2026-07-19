# New Session Instructions and Review Points

## Objective

Continue the fixed-H WPW checkpoint work from the existing dirty worktree.
Repair Task 3 so that the initial state is a true nonorthogonal W+P projection
of the DC density-carrying fragment ensemble, then implement zero-interface
relaxation and interface continuation. Do not begin long Si64 production until
the focused implementation and review gates pass.

## Workspace identity

```text
worktree: /Users/otobetoshihito/SALMON-dev/SALMON2_RTDG
branch:   codex/singlescale-vortex-observables
HEAD:     125bdb30d1cd1e49ab09835888797c015f75aa6d
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

Important: Task 3 was reopened after review. The current implementation in
`build_wpw_projected_occupied_seed` is not acceptable as final code because it
uses W overlaps as coefficients and explicitly zeros P coefficients. Replace
that path with W/P overlap construction plus `S C=B`.

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

Here `F` contains the converged DC occupations. Do not compute this diagnostic
with the normalized occupied block. A unitary-rotation density-invariance test
must be restricted to equal-occupation blocks or transform the full occupation
matrix as `F -> U^dagger F U`.

### P1: source naming and provenance

The source is the direct sum of occupied DC fragment orbitals that carries the
converged DC density. It is not the global Flux-LCFO eigenvector set `coef_wf`.
Do not run `diag_eigenexa` merely to generate a seed. Record this distinction in
diagnostics and checkpoint metadata.

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
  history; do not use a global dense fallback. Block-local preconditioning
  requires a focused failure oracle and a separate design review.

## Findings-first review checklist

Review in this order:

1. **Projection mathematics:** Does the code construct both overlap blocks and
   route every support-row contribution to the canonical owner before solving
   `S C=B`? Is the reported residual the equation residual for `C_raw`?
2. **Source provenance:** Are the orbitals exactly the DC density-carrying
   fragment ensemble? Are occupations and charge consistent collectively?
3. **Occupied projector:** Is occupied rank preserved? Is the S-projector
   invariant under unitary rotations and extra-state completion? Is captured
   norm computed from `C_raw`, not the normalized block?
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
- zero-interface solver iterations and residuals;
- accepted/rejected continuation steps and final lambda;
- frozen-state validation result;
- findings-first review, including “no P0/P1 findings” only if supported by a
  fresh independent review.
