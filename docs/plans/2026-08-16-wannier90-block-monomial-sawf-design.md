# Wannier90 Block-Monomial SAWF Design

## Goal

Make Wannier90 preserve the full affine symmetry of individual Wannier-center
orbits, while retaining Wannier90 as the localization engine and avoiding a
second dense post-localization gauge optimization.

## Constraint exposed by the Si64 failure

At Gamma, Wannier90 requires

`D_band A = A D_wann`.

The current full-rank trial matrix makes `D_wann=A^H D_band A`, which is a
general dense representation.  It constrains the retained subspace but does
not require individual Wannier functions to transform between center blocks.
The unchanged Si64 centers are therefore a valid result of the current SAWF
input, not an affine-action-direction bug.

## Design

Before writing the DMN file, compute periodic centers and center-confidence
receipts for the pre-Wannier trial frame.  For each affine generator, analyze
the existing dense band action in that frame and require a block-monomial
structure:

1. every source center block maps to exactly one target center block;
2. leakage outside that target block is below tolerance;
3. the internal source-to-target block is unitary;
4. the generator blocks satisfy the supplied affine group relations.

When these gates pass, use the measured block-monomial action as `D_wann` and
the matching pre-Wannier frame as `A`.  Do not threshold or invent a target
representation when the gates fail.  A failed gate identifies the seed frame,
not Wannier90, as the component requiring redesign.

After Wannier90 converges, retain the full center-orbit gate and additionally
measure SAWF generator covariance from the returned transform.

## Memory and scaling

Generator actions remain streamed one at a time.  Persistent storage is the
distributed trial frame, center metadata, and one `N x N` generator matrix on
the coordinator, matching the existing DMN writer contract.  No all-generator
`G*N^2` tensor and no second full-space localization pass are introduced.

## Testing

- Synthetic repeated-center fixture with a center permutation and a nontrivial
  internal `p`-like unitary block.
- Leakage and ambiguous-center REDs.
- Covariance RED proving a block-monomial `D_wann` cannot be paired with a
  mismatched `A`.
- Si64 pre-Wannier feasibility receipt, followed by the full center-orbit gate.

