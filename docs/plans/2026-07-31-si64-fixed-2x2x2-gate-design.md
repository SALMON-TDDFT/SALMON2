# Si64 Fixed 2x2x2 Gate Design

## Decision

Task 9 uses one physical DC decomposition, `2x2x2`, on eight MPI ranks
and fixes buffer 5 as the initial production profile.
The immediate objective is to prove that the complete ground-state route
runs normally and is insensitive to the candidate window within buffer 5.
Buffer 6 provides an individual normal-operation and sensitivity
diagnostic, not a cross-buffer acceptance threshold. MPI-placement
invariance is deferred.

## Rationale

On the fixed 32-cubed Si64 grid, `4x2x1` makes the x core only eight grid
points wide. A five- or six-point buffer on each side then dominates the
fragment box. Changing that physical decomposition changes the DC
approximation rather than merely changing its MPI placement. It is
therefore not an appropriate equivalent-decomposition gate for this small
system.

The failed `4x2x1` evidence remains preserved as a diagnostic. It is not
accepted as evidence against the fixed `2x2x2` production route.

## Acceptance Matrix

The two production acceptance rows are:

- `2x2x2`, buffer `(5,5,5)`, candidate 192, target 48;
- `2x2x2`, buffer `(5,5,5)`, candidate 256, target 48;

The following two diagnostic rows remain required:

- `2x2x2`, buffer `(6,6,6)`, candidate 192, target 48;
- `2x2x2`, buffer `(6,6,6)`, candidate 256, target 48.

Every row uses eight MPI ranks and one OpenMP thread. Each must pass the
existing route, numerical-quality, memory-evidence, forbidden-promotion,
checkpoint-publication, and checkpoint-reuse gates.

For each recorded convergence quantity, the two buffer-5 production rows
must agree within the existing relative tolerance of `1e-4`. Both
buffer-6 diagnostic rows must individually pass all route, numerical, and
checkpoint gates, but their difference from buffer 5 is recorded rather
than used as an acceptance threshold. The normal reference remains
required. No comparison with `4x2x1` is part of Task 9 acceptance.

## Safety and Scope

This change modifies only the Si64 Task 9 runner, checker, and results
documentation. It does not change conventional DC, LCFO, EigenExa,
direct real-space DG, WPW, normal checkpoint publication, or RT route
selection. Task 10 remains blocked until the revised Task 9 gate passes.

MPI-count and large-system physical-decomposition tests remain future
validation work.
