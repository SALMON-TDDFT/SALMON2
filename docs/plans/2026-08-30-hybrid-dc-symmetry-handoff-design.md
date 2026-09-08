# Hybrid DC Symmetry Handoff Design

## Context

The Si64 Hybrid continuation acceptance uses the same `2x2x2` decomposition
for DC and DG on eight MPI ranks.  Wannier construction already certifies the
retained Wannier space against the complete physical affine group.  The
subsequent Hybrid production-basis preparation currently repeats a stronger,
fragment-level condition: every physical operation must map each complete DC
fragment to one complete DG fragment.  The Si64 affine translations do not
preserve that rectangular fragment partition, so this redundant condition
rejects an otherwise symmetry-capable basis.

The final whole-system LCFO solve is responsible for recovering the physical
symmetry from the certified Wannier basis.  Fragment covariance is therefore
not a required property of the intermediate DC representation.

## Design

Keep DC and DG on exactly the same decomposition.  Validate the fragment
count, fragment origins, extents, and physical point ownership before Hybrid
initialization, and reject any mismatch.

Remove only the Hybrid DC production-basis requirement that the complete
physical symmetry action map whole fragments.  Do not silently replace the
physical group with an identity-only group and do not manufacture a
fragment-compatible subgroup receipt.

Reuse the authoritative symmetry provenance already produced by Wannier
construction.  The Hybrid path must retain the complete certified Wannier
basis and its fingerprint.  Adding PW functions is allowed because it does
not remove the certified Wannier subspace.  A missing or mutated Wannier
certificate remains fatal.

Run the final whole-system LCFO solve on the retained WF+PW space.  Preserve
the existing post-LCFO checks of electron count, residuals, and covariance of
the complete occupied projector under the actual physical group.  These are
the final symmetry acceptance checks; no individual-state symmetry is
required.

## Error handling

- Reject DC/DG decomposition mismatches before basis construction.
- Reject missing, failed, or fingerprint-inconsistent Wannier symmetry
  provenance.
- Reject failure of the final LCFO occupied-projector covariance check.
- Do not reject solely because a physical operation cuts across DC fragment
  boundaries.

## Verification

Use TDD to add a production fixture in which a known physical operation cuts
across fragments while the authoritative Wannier certificate is complete.
Require Hybrid preparation to continue without relabeling the group as
identity-only.  Retain negative fixtures for decomposition mismatch,
certificate mutation, and failed final occupied-projector covariance.

Run the focused Hybrid production-basis and continuation suites, all protected
symmetry and RT routes, a fresh build, and finally the eight-rank Si64
GS-to-zero-field-RT acceptance without a timeout.
