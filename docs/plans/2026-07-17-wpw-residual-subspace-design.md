# WPW Residual-Augmented Matrix-Free Eigensolver Design

## Problem

The production WPW smoke reaches the matrix-free algebra callback but the
fixed-step Richardson update exhausts 160 inner iterations with `info=40`.
The measured metric orthogonality is about `1e-15`, while the generalized
eigen residual remains about `1e1`.  The failure is therefore the update rule,
not loss of metric orthogonality or a communicator mismatch.

## Accepted update

Replace the fixed scalar Richardson step with a locally optimal block
preconditioned conjugate-gradient (LOBPCG-style) Rayleigh--Ritz update.  For
the current retained block `Q`, residual `R`, and previous search block `P`,
compute

```text
R = H Q - S Q Lambda
Z = [Q, R, P]
```

and solve the reduced generalized problem in `span(Z)`.  Keep its lowest
`nretain` Ritz vectors as the next `Q`, and form the next `P` from the `R` and
old-`P` coefficient blocks.  The initial `P` is zero and is removed by the
rank-revealing metric solve.  All `H`, `S`, and Gram operations
continue through the existing distributed callbacks.  No global WPW operator,
density matrix, or fragment scan is formed.

## Rank-revealing reduced solve

`Z^H S Z` can be rank deficient: the initial search block is zero, converged
residuals vanish, search columns can be dependent, and `3*nretain` can exceed the global basis dimension in a
small fixture.  The reduced solver therefore diagonalizes the reduced metric,
retains eigenvectors with eigenvalue greater than
`metric_cutoff * max(metric_eigenvalue)`, and fails only when fewer than
`nretain` directions remain.  It then solves the Hermitian Hamiltonian in that
effective metric-orthonormal space and returns the lowest requested Ritz
pairs.

The rank-revealing operation is confined to reduced matrices of dimension at
most `3*nretain`; it is not a dense production fallback.

## Failure behavior

- Non-Hermitian reduced H or S fails closed.
- A non-positive reduced metric fails closed.
- Effective metric rank below `nretain` fails closed.
- Callback and MPI failures retain their existing propagation.
- Exhausting the bounded inner iteration count still returns `info=40`.

## Verification

Extend the two-rank matrix-free fixture with a diagonal generalized problem
having a wide spectral range that the Richardson update cannot converge.
Require the expected lowest eigenvalues, generalized residual, and metric
orthogonality.  Then run the existing matrix-free SCF fixture, focused WPW
contracts, the full MPI/EigenExa build, and the two-rank production smoke.

No commit, push, or pull request is part of this session.
