# Point-Closed Periodic-Center Gauge Design

## Evidence

The Si64 final retained representation is unitary to `1.17e-14`, but affine
operation 2 has monomial defect `0.742`, center-block leakage `1`, and center
orbit residual `0.175`.  The failure is therefore an internal gauge failure,
not loss of the retained subspace.  The current joint periodic-center routine
uses one identity-start Jacobi solve and never consumes the point-cogroup
representation.  It also accepts a stationary point without bounding the
reported off-diagonal objective.

## Construction

Build the compact reference-sector point representation

```
D_p = A^H U_p A
```

by the existing row-owned streamed overlap assembler.  Gather only the small
`m x m x |P|` result.  Extend the position objective from six matrices `H_q`
to their full point orbit

```
H_(p,q) = D_p^H H_q D_p.
```

Jointly diagonalize this closed set.  Since left multiplication by any `D_g`
only permutes the point index, the objective is point-cogroup invariant.  Gate
both the Jacobi update and the normalized residual objective.  Validate every
`D_p` as unitary and require the final representation to be block-monomial
between tolerance-equal center blocks.

The added persistent data are `O(|P| m^2)`.  No global dense spatial operator,
all-character tensor, or replicated retained-rank projector is introduced.

## Production flow

Before reference-sector canonicalization, assemble the compact point matrices
from `translation_reference_spatial`, `ow_core_weights`, and the 48 local point
maps.  Pass them into the canonicalizer, bind them to its fingerprint, then
continue the existing translation-character/Gamma schedule unchanged.

## Acceptance

- Synthetic point-swap and repeated-center blocks pass on MPI 1/2/4/8.
- Corrupt, nonunitary, and rank-disagreeing point matrices reject collectively.
- The route contract requires production to provide point representations.
- Si64 operation 2 must reduce center-block leakage and center residual below
  the unchanged tolerance; otherwise retain the diagnostic stop.

