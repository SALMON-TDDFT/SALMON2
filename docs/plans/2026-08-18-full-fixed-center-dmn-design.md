# Full Fixed-Center DMN Design

## Problem

The Si64 production path constructs a closed fixed-center site-symmetry group of order 12, but writes only the affine identity and five generators to the Wannier90 DMN file. Wannier90 applies a one-pass average over the supplied operations. A one-pass Reynolds average is a projector onto the invariant subspace only when the supplied set is the complete closed group. Averaging identity plus generators is generally neither idempotent nor invariant under the generators.

The writer also finishes the transaction with the first six fixed-center operation descriptors even though the six matrices were assembled from the global affine identity and global affine generators. Closure validation is disabled. Consequently, the DMN payload, its operation metadata, and the post-Wannier covariance contract do not describe one common group action.

## Design

Publish every operation in `fixed_center_operations` in its canonical fixed-center order. Assemble each operation directly from the corresponding column of `fixed_center_symmetry_map`; do not translate through `global_affine_generators` or assume that two independently constructed catalogs share indices.

Set the DMN operation count to `fixed_center_group_order`, pass the identical complete `fixed_center_operations` array to `finish_sawf_dmn`, and retain the default closed-group validation. The identity flag is determined from `fixed_center_identity_operation`, not from loop position.

This keeps memory bounded: one fixed-center representation is assembled, gathered, appended, and released at a time. No `N^2 * |G|` tensor is introduced.

## Error handling

All existing collective assembly, allocation, writer abort, and MPI agreement paths remain in use. Failure of full-group closure becomes a fatal pre-Wannier diagnostic rather than being deferred to the post-Wannier covariance gate.

## Verification

1. Add a route regression that rejects generator-count DMN publication and requires use of the full fixed-center map and closed-group finish.
2. Add a numerical regression showing that a one-pass identity-plus-generator average for `Z3` is not invariant/idempotent, while the complete-group average is.
3. Run the DMN format test and route checker in the RED state.
4. Make the minimal production change.
5. Run focused DMN, route, W90 MPI 1/2/4/8, SAWF/DMN, build, and diff checks.
6. Re-run Si64 with MPI 8 / OMP 1 and confirm that Wannier90 completes and the post-Wannier covariance gate is reached and passes, or record the next measured failure without weakening a gate.

