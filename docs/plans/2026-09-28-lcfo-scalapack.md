# Distributed complex LCFO full diagonalization

Goal: eliminate replicated global Hamiltonians/eigenvectors in a selectable full eigensolver for complex LCFO. The existing CheFSI and LAPACK paths remain unchanged. No fixed memory ceiling is introduced.

Architecture: stream fragment diagonal and halo blocks into a two-dimensional block-cyclic matrix, diagonalize with PZHEEV, verify Hermiticity/orthogonality/residuals through distributed PBLAS, and reduce only fragment coefficient rows to their output representatives. Eigenvalues and small fragment blocks remain replicated; full matrices never gather on a root. The solver still has quadratic total storage and cubic full-spectrum work.

1. Add a distributed dense matrix utility and a deterministic complex-Hermitian MPI comparison against LAPACK, including small dimensions and empty local tiles.
2. Wire lcfo_eigensolver='scalapack' for complex LCFO and reject unavailable builds/real LCFO explicitly. Preserve halo conjugation and repeated periodic-image accumulation.
3. Verify HSE/PBEh DC energies/eigenvalues and reconstructed RT, plus multi-k complex Si. Build with/without ScaLAPACK. Document process layouts and remaining replicated fragment storage.

## Verification

The bounds-checked matrix oracle passed 24 cases (six matrix sizes and four process counts), with reversed communicator ordering, empty local tiles, partial-row recovery and non-Hermitian rejection. A scalar analytic path handles the observed installed PZHEEV scalar shortcut returning a nonunit vector.

HSE/PBEh DC+reconstructed impulse tests passed on 2/8 ranks. Maximum compared eigenvalue difference was 1.95e-11 Ha and maximum RT energy history difference 7.55e-12 in output units. The ScaLAPACK build also passed the 17 existing Ehrenfest tests. The unmodified complex Si verification passed all four k points and 192 eigenvalue records. The conventional LCFO diagnostic prefix is preserved with an appended LCFO_SCALAPACK marker.

Independent review found no blocking issue. Remaining replication is at fragment-block/coefficient and eigenvalue level; this full-spectrum solver has quadratic total dense storage and cubic work. The default LAPACK route is retained, including its original replicated storage.

Final compatibility checks: 13 existing HSE/PBEh SCF/DC/RT tests passed with ScaLAPACK disabled; the new option is explicitly rejected by that build. HSE-disabled and ScaLAPACK-enabled builds also passed.
