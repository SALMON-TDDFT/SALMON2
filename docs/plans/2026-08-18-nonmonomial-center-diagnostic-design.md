# Nonmonomial Wannier-Center Diagnostic Design

## Evidence

The corrected Si64 route passes the fixed-center Wannier covariance gate (`5.09e-13`), the factored 48-element point-cogroup/cocycle proof (`2.84e-12`), and the periodic-position Gram/canonical gates. It then rejects operation 2 solely because individual first-moment centers do not form a permutation orbit. The measured retained representation remains unitary (`5.03e-12`) but has monomial defect `0.904` and center-block leakage `1.0`.

This is a dense unitary action within the retained space, not a broken symmetry action. Individual Wannier centers are subsequently used only to choose storage owners. Hamiltonian, overlap, position, velocity, transition matrices, density, and checkpoint receipts are rebuilt from the actual basis and have independent numerical and symmetry gates.

## Decision

Keep the center-orbit test and point-center gauge analysis as diagnostics. If a center orbit does not close but the diagnostic successfully proves that the retained action is unitary, report the failed operation, center mismatch, monomial defect, leakage, and workspace, then continue to deterministic per-orbital center ownership.

Fail closed if the diagnostic itself cannot establish a valid unitary retained representation. Do not relax center matching tolerance, alter the basis, or add another localization pass.

## Physical acceptance

Acceptance remains controlled by the already mandatory full-affine subspace/cocycle proof and by the later overlap, Hamiltonian, generalized-eigenproblem, density, Hermiticity, observable, exact-group, and checkpoint gates. The center diagnostic is not used as evidence for those physical contracts.

## Verification

- Route RED forbids an unconditional `localized Wannier center orbit failed` stop and requires continuation to center-owner assignment after the diagnostic.
- A focused synthetic dense-unitary center fixture remains diagnostic evidence rather than an acceptance failure.
- Existing invalid/nonunitary diagnostic REDs remain fail closed.
- Run construction MPI 1/2/4/8, route, W90 MPI 1/2/4/8, build, and diff checks.
- Rerun Si64 through center redistribution and the subsequent stitched operator/physical gates.

