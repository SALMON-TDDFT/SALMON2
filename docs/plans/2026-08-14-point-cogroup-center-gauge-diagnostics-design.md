# Point-Cogroup Center-Gauge Diagnostics Design

## Goal

Identify whether the post-Wannier90 point-cogroup representation preserves the
localized-center blocks, without weakening the existing affine center-orbit
gate.

## Evidence and scope

For the Si64 production case, affine operation 2 is the nonsymmorphic point
operation

```text
R = [ 0 -1  0 ]    tau = (1/8, 1/8, 1/8)
    [ 1  0  0 ]
    [ 0  0  1 ]
```

It is not a pure translation.  The final Wannier90 centers fail operation 2
with a bottleneck residual of about `2.03e-1`.  The post-character gauge lowers
the first reported residual to `4.77e-2`, but the required tolerance is
`3.16e-3`.  Reversing the affine action does not close the Wannier90 centers.
The periodic moment magnitudes are moderate, so this is not explained solely
by undefined centers from vanishing first moments.

This change is diagnostic only.  It does not relax the tolerance, skip an
operation, or accept a non-closed center orbit.

## Architecture

Add a distributed diagnostic primitive that consumes the final row-owned
Wannier basis, its spatial weights, one affine point map, and the measured
periodic centers.  Reuse the existing streamed symmetry-overlap assembly to
obtain the retained-space representation for that operation.  Match affine
mapped centers to their nearest measured center blocks and report two defects:

1. **Monomial defect:** how far the absolute-squared representation is from a
   one-to-one orbital permutation.
2. **Center-block leakage:** the maximum representation weight connecting an
   orbital to targets outside the tolerance-compatible mapped-center block.

The routine must remain row distributed.  It may retain one `N x N/P`
representation and `O(N)` matching metadata, but must not create a replicated
`N x N` matrix.  All shape, allocation, MPI-status, and workspace arithmetic
checks follow the existing factored-proof contracts.

## Data flow

The production path already has `ow_core_values`, `ow_core_weights`, affine
maps, and measured centers immediately before the failing center gate.  When
the center gate rejects, or when diagnostics are explicitly requested by the
call site, it evaluates the first failing operation and emits one line with:

```text
operation, center_bottleneck_residual, monomial_defect,
center_block_leakage, representation_unitarity_defect,
workspace_peak_bytes
```

The existing actionable center mismatch remains the authoritative failure.
The new line explains whether failure comes from internal point-gauge mixing
or from center data inconsistent with an otherwise monomial representation.

## Error handling

The diagnostic primitive is collective.  Rank-disagreeing metadata,
duplicate/missing row ownership, allocation failure, MPI failure, nonfinite
input, unsafe extent arithmetic, or a nonunitary retained representation cause
a collective diagnostic failure.  Production still stops at the original
center-orbit gate and prints both messages when the diagnostic succeeds.

## Testing

Add MPI fixtures for:

- an exactly permuted two-center representation with zero leakage;
- a unitary internal rotation between different center blocks, which preserves
  the retained subspace but produces nonzero center-block leakage;
- an allowed unitary rotation inside a repeated-center block, which has zero
  center-block leakage even though it is not elementwise monomial;
- duplicate row ownership and rank-disagreeing operation metadata rejection;
- rank-independent receipts on MPI 1/2/4/8.

Run the focused construction MPI fixture, the production route checker, the
release build, and `git diff --check`.  Only after the diagnostic identifies
the failed contract will a separate design choose between a point-gauge repair
and a revised mathematical center contract.
