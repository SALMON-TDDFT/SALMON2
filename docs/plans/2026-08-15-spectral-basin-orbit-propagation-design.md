# Spectral basin orbit propagation design

## Goal

Diagonalize only one canonical basin per symmetry orbit and generate every
symmetry-related trial block without retaining or diagonalizing the remaining
basin operators.

## Design

The smallest canonical basin index is the representative of each orbit.  Its
selected eigenvectors are supplied in retained-state coordinates.  A
deterministic breadth-first traversal of the validated basin generator maps
records one parent generator for each target basin.

For one target at a time, reconstruct its generator word and apply the
row-owned retained-state generator matrices.  Each application forms owned
output rows locally and performs one `Nstate x block_rank` collective sum.
Only two full block buffers are live; completed channels are stored directly in
the final row-owned `local_state_rows x Nstate` result.

The routine verifies exact state-row ownership, equal selected rank along every
basin-generator edge, complete basin reachability, finite propagated values,
and the final distributed Gram matrix.  Failure to obtain an orthonormal full
retained frame is collective.  The output fingerprint binds catalog and
representation provenance plus quantized row-owned coefficients.

This changes symmetry-related basin diagonalization from one dense eigensolve
per basin to one per orbit.  Communication is proportional to the propagated
selected blocks, not to dense basin operators.

