# WPW Rank-Local Dense Apply Design

## Problem

The Si64 2x2x2 production preflight reaches the first matrix-free LOBPCG step, but bounded H/S application treats rank-local WW/WP/PP blocks as scalar sparse entries.  For 160 retained states (480 vectors in the expanded subspace), this performs hundreds of millions of small vector updates per application even though WP and PP are complete bounded rectangles.

## Design

Keep the existing sparse arrays as the authoritative checkpoint and provenance representation.  During transactional bounded-operator initialization, also pack H and S into rank-local dense caches indexed by the already validated required/owned endpoint indices:

- WW: `required_w x required_w`
- WP: `required_w x owned_p`
- PP: `owned_p x required_p`

The caches contain only the current fragment root's bounded support.  No global matrix is allocated and owner-targeted coefficient exchange is unchanged.  Application fetches support coefficients, evaluates the three blocks with BLAS-backed `matmul`, and reduces the W partial back to its owners exactly as before.

Duplicate sparse coordinates are accumulated while packing, preserving the current additive sparse semantics.  Allocation, index, and finite-value failures invalidate the candidate collectively; the previously published operator remains untouched.

## Verification

- Source contract forbids scalar sparse-entry loops in the production apply path.
- MPI dense oracle compares cached and expected H/S results.
- A multi-vector benchmark verifies the 480-vector path completes within a bounded preflight time.
- Si64 `dg_wpw_scf_max_iter=1` must enter and exit the first algebra step fail-closed rather than remaining in H/S application.

Checkpoint schema, sparse provenance, physical conventions, and RT interfaces do not change.
