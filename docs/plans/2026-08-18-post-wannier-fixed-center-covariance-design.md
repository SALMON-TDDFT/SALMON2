# Post-Wannier Fixed-Center Covariance Design

## Evidence

The complete 12-operation fixed-center DMN completes Wannier90 normally for Si64. Reading the produced checkpoint and DMN independently gives a maximum fixed-center covariance defect of `2.78510548e-12`, below the configured `1e-10` tolerance. Production nevertheless rejects the transform because its immediate post-Wannier loop tests the five generators of the 1536-element global affine group, which were deliberately not supplied to Wannier90.

## Decision

The immediate post-Wannier gate validates exactly the same complete fixed-center group and target representations that were supplied through DMN. It streams all 12 operations from `fixed_center_symmetry_map` in canonical order and retains the maximum defect and workspace receipt.

This does not weaken the full physical symmetry contract. Translation covariance is constructed and validated by the later character-sector gauge. The 48 point-cogroup representatives and translation cocycle then prove the final full-affine closure. The immediate gate answers only whether Wannier90 honored its own DMN contract.

## Rejected alternatives

- Supplying all 1536 affine operations to Wannier90 confuses site symmetry with space-group orbit transport and violates the bounded-memory design.
- Removing the immediate covariance gate loses a useful component-boundary check.
- Retaining the global-affine generator gate before the translation correction demands a property that the architecture explicitly constructs only later.

## Verification

- Route RED requires the post-Wannier loop to use `fixed_center_group_order` and `fixed_center_symmetry_map`, and rejects `global_affine_generators` inside that loop.
- Focused covariance, W90 MPI 1/2/4/8, DMN, route, fragment symmetry, build, and diff checks pass.
- Si64 MPI 8 / OMP 1 passes the immediate fixed-center covariance gate and advances into the translation-character correction. Subsequent failures are measured without relaxing their gates.

