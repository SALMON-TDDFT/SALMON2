# Spatial Representation Convention Design

## Problem

The Si64 post-Wannier Hamiltonian has order-one covariance defects in its kinetic, local, and nonlocal components. The current code then group-averages that matrix, producing a formally symmetric matrix and an artificial occupied-boundary degeneracy. A prior change to the row-action formula did not remove the underlying defect.

The representation is projected from spatially permuted basis functions, while the basis diagnostic and operator symmetrizer separately assume a matrix-action convention. Existing tests inject a representation matrix directly and therefore do not prove that representation construction, spatial action, Cartesian gradient action, and operator covariance use one convention.

## Design

Add one minimal MPI fixture that starts from a point permutation and complex nonsymmetric basis functions. The fixture must use the production representation builder, spatial exchange, gradient transformation, weak-operator assembly, and distributed operator covariance measurement. It will compare the four possible transpose/conjugation conventions and require exactly the convention derived from the spatial action.

The fixture will first be run against the current implementation and must fail for the expected convention mismatch. Production code will then receive the smallest convention correction needed to make the end-to-end fixture pass. Tests that construct a representation matrix directly are insufficient for this decision and will only remain as lower-level algebra tests.

After the convention is fixed, a failed raw gradient or Hamiltonian covariance check becomes fatal before group averaging. Group averaging may remove roundoff-scale residuals only; it must not convert an order-one noncovariant operator into the published Hamiltonian.

## Verification

- Focused MPI fixture on 1, 2, 4, and 8 ranks.
- Existing fragment-symmetry and construction MPI fixtures.
- Route checker and `git diff --check`.
- Si64 production rerun only after focused tests establish the convention.

## Non-goals

- No change to DC density construction, Wannier90 iteration policy, or eigensolver thresholds.
- No additional symmetry repair layer.
- No relaxation of the occupied-boundary cluster gate.
