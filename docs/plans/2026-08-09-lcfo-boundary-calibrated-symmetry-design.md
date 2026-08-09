# LCFO Boundary-Calibrated Occupied Symmetry Design

## Purpose

Preserve a rank-fixed full-system occupied space without mistaking the known DC+LCFO fragment-boundary stitching defect for physical bulk symmetry breaking.

## Error model

The converged DC fragments and the LCFO reconstruction are not pointwise smooth across fragment boundaries.  Applying an exact atomic symmetry to a reconstructed LCFO state therefore creates small components outside the numerical rank-128 occupied projector even when the underlying material is inversion symmetric.  Counting every such component as a new state makes the symmetry closure rank depend on the stitching error and can inflate it toward the complete 128-by-group-order orbit.

The grid is partitioned into an interior and a boundary layer.  The boundary layer contains unique-core points within one SALMON derivative-stencil radius of a fragment face.  It is not removed from normalization or physical observables.  It is used only to classify the source of the covariance residual.

## Measurements and gates

LCFO/EigenExa returns at least 384 eigenvector coefficient columns even when the conventional DFT state count is smaller.  The complete retained candidate is reconstructed from the same expression on every fragment,

`fragment LCFO basis(core+buffer) * global LCFO coefficients`,

and contains the first 384 LCFO states.  No independently constructed fragment-local complement or direct sum is admitted.  Only these final 384 distributed real-space rows are materialized; a full-system real-space wavefunction array is never gathered.

For the metric-orthonormalized LCFO rows, measure each full-system atomic operation before any correction:

1. the best-fit rank-128 representation from distributed overlaps;
2. total projector leakage;
3. leakage norm on the interior and boundary layer separately;
4. the LCFO value and gradient jump scale across matching fragment faces;
5. density covariance and the inversion-odd density norm on the interior and boundary separately.

The boundary leakage allowance is calibrated from the measured LCFO face value/gradient mismatch, grid spacing, and stencil radius.  No free material-specific tolerance is introduced.  Interior leakage must satisfy the strict symmetry tolerance plus the propagated roundoff part of that baseline.  Boundary leakage may be accepted only when it is bounded by the measured stitching baseline.  A residual of the same magnitude placed in the interior must be rejected.

## Rank-preserving correction

When the calibrated gates pass, retain exactly 384 LCFO states, with the first 128 defining the occupied projector.  Build the measured group action inside the 384-dimensional space, metric-polar project each operation to a unitary matrix, and synchronize the matrices to the product table.  Symmetry-average within rank 384; do not append numerical boundary-residual directions.  Reject the transaction if synchronization changes the occupied projector, density, or occupied energy beyond the calibrated boundary error, or if a complete symmetry block would be cut.

Atomic projectors may be used only as a localizer inside the LCFO rank-384 space.  They do not supply additional basis rows.

## Evidence

Synthetic MPI tests place identical perturbations either in the boundary layer or interior.  Boundary-calibrated perturbations pass; interior perturbations fail.  Genuine ideal-Si64 evidence reports raw and corrected total/interior/boundary leakage, face value and gradient mismatch, density inversion-odd norms, occupied rank, correction norm, and energy change.  Subsequent V3, field-off, linear-response, and long-pulse HHG calculations continue to use polarization as the primary observable.
