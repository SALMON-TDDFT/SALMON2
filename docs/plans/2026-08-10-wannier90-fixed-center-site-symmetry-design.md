# Wannier90 Fixed-Center Site-Symmetry Design

## Problem

The full 1536-operation affine proof passes for the 384-state Si64 LCFO
subspace, and the Wannier90 transform is unitary to about `7e-14`.  Standard
Wannier90 nevertheless reaches its 200-iteration limit with a large
gauge-dependent spread and produces centers that are not closed under the
full-system symmetry.  The adapter currently mistakes reduction from the
initial spread for convergence.

## Architecture

Keep the two proof layers separate.

1. The complete affine group remains a streamed LCFO-subspace proof.  It is
   not materialized as representation matrices and is not recomputed after a
   unitary MLWF gauge change.
2. The maximal full-system fixed-center point subgroup, of order at most 48,
   supplies Wannier90 symmetry constraints and later operator projection.
   Its center is derived from the full atomic catalog, never from fragments.

Before Wannier90 setup, select the fixed-center subgroup.  For each selected
operation, assemble one distributed 384 by 384 representation matrix in the
orthonormal LCFO basis.  Use the same matrix for `D_band` and `D_wann`; with
the orthonormal LCFO basis as the projection gauge, `A=I` and the DMN
covariance is exact.  Send one operation at a time to the existing DMN writer,
which validates unitarity and group relations while spilling matrices to
scratch.  No `48*Nwann**2` tensor is retained on every rank.

Wannier90 runs with `site_symmetry = .true.` and a strict `symmetrize_eps`.
The resulting transform must remain Gamma-real and unitary, the optimization
must report convergence before its iteration limit, and the recomputed center
set must close under the fixed-center subgroup.  Translation and cocycle
evidence remains in the full-affine V3 proof.

## Memory and Parallelism

Spatial-core ranks assemble row-owned overlaps in orbital tiles.  At most one
dense point-group matrix is gathered for the DMN writer at a time.  Rank 0
owns the writer scratch files; all ranks retain only their core wavefunctions,
one row block, and communication tiles.  V3 records the measured peak rather
than an estimate alone.

## Failure Policy

Reject missing or stale DMN files, nonclosed fixed-center operations,
nonunitary representations, covariance failure, Wannier90 iteration-limit
exhaustion, nonfinite spreads, noncanonical Gamma gauges, and center-orbit
failure.  Do not fall back to unconstrained Wannier90 or the custom localizer.

## Verification

Use genuine RED tests for DMN publication, site-symmetry setup, bounded
one-operation workspace, convergence rejection, and center closure.  Compare
small streamed matrices and DMN receipts with dense references on MPI
1/2/4/8.  Then run the clean full-feature overlay and ideal undisplaced Si64
GS before LR/HHG acceptance.
