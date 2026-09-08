# WF+PW Fragment-Compatible Symmetry Design

## Context

The divided WF+PW production basis receives compact affine generators whose
spatial maps act on distributed physical grid points.  A `2 x 2 x 2` DC
partition need not be a block system for any of those nonidentity generators.
Filtering only the supplied generators can therefore leave no operation, even
though the identity always preserves every DC fragment.

## Design

The production PW adapter shall prepend an explicit distributed identity map
and the identity reciprocal rotation.  It shall then test every supplied
nonidentity generator with the existing collective window-distribution
contract and retain only operations that map each complete DC fragment to one
complete DC fragment.  Spatial maps and reciprocal rotations shall use the
same selected column indices.

The existing permutation, fragment covariance, reciprocal cutoff closure, and
window covariance checks remain authoritative.  An incompatible operation is
not treated as valid; it is omitted from the symmetry packet generators.  The
identity guarantees a nonempty compatible set without weakening a check.

With identity only, the catalog still contains every reciprocal vector inside
the cutoff and every fragment window.  Thus the spanned WF+PW space is not
reduced; only symmetry-based packet grouping is reduced.  Compatible
nonidentity generators remain available when the DC decomposition admits
them.

## Scope and invariants

- Change only the divided WF+PW production-basis path.
- Do not change conventional DC+LCFO or Wannier90 behavior.
- Do not add a post-LCFO density SCF or convergence gate.
- Preserve exactly one final distributed LCFO diagonalization.
- Select operations collectively and deterministically for all MPI layouts.

## Verification

Add an MPI fixture whose supplied nonidentity operation splits a fragment.
The fixture must first fail without identity insertion, then pass on 1, 2, 4,
and 8 ranks while retaining a valid production catalog.  Rebuild SALMON, run
the divided and legacy route contracts, and finally rerun Si64 with eight MPI
ranks, `OMP_NUM_THREADS=1`, and no time cutoff.
