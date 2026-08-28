# DG Production Selection Boundary Design

## Purpose

The production DG continuation route must close the requested WF blocks and
PW packets under the authoritative physical group before it freezes or
materializes a basis.  The current production PW helper cannot satisfy this
contract because it silently drops operations that do not map whole fragments,
constructs every packet immediately, and discards the action maps.

## Boundary

Split production basis preparation into three explicit phases.

1. Analyze the supplied physical operations against the spatial partition.
   Return an immutable analysis receipt and the accepted spatial, fragment,
   reciprocal, and packet actions.  A supplied nonidentity operation that is
   not fragment covariant is an error.  Identity-only operation is valid only
   when represented by one explicit identity action and a successful receipt.
2. Construct the finite PW packet universe inside the requested cutoff.  Give
   every packet a stable ID and return the complete packet action.  This phase
   does not select or materialize columns.
3. Accept externally closed effective packet IDs, validate that they are a
   duplicate-free invariant subset of the universe, and freeze only those
   packets into the production catalog.  Ownership and fingerprints are
   recomputed from the effective catalog.

The existing convenience entry point remains available for protected routes.
It uses the complete packet universe and therefore preserves their observable
behavior.  Only the new continuation route uses the phased API.

WF symmetry is kept at block/subspace level.  A general complex Wannier
representation is not converted into a fictitious permutation of individual
orbitals.  Task 7b supplies the authoritative WF-block action obtained from
the accepted overlapping-Wannier analysis and closes block IDs separately
from PW packet IDs.

## Failure and provenance

All ranks must agree on analysis inputs and results.  The receipt records
successful analysis, identity-only status, operation count, and a nonzero
fingerprint over the complete action payload.  Empty operation lists,
inconsistent distributed inputs, non-covariant fragmentations, invalid
effective IDs, and non-closed effective selections fail before catalog
freezing.

## Verification

MPI tests run on 1, 2, 4, and 8 ranks and cover nontrivial symmetry,
identity-only symmetry, rejection of a non-covariant operation, deterministic
packet actions, closure of a split requested orbit, and rejection of a
non-closed effective selection.  The existing full-universe convenience API
is also exercised as a protected-route regression.
