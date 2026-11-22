# Small-eta screening stabilization

## Accepted design
Preserve the original Hermitian Hamiltonian selected by the existing lower-triangle LAPACK convention. Use the tridiagonal factorization as a preconditioner, then refine against that original Hamiltonian. Do not replace it with approximate projected eigenvalues or change projection thresholds. Retain fractional populations, both eta channels, and the existing metallic head construction. RT propagates the same occupied orbitals.

Compensated residual products use fixed-point slices and BLAS; high+low solution storage and compensated Hermitian Gram accumulation avoid losing the small real response beside a large imaginary pole. Adapt precision and refinement iterations; explicitly reject inadequate precision or convergence. Python dense-LU refinement is an independent validation implementation, not a production dependency.

## Validation
- 90-digit original-H oracle with imperfect tridiagonal representation, complex off-diagonal responses, eta 5e-4 through 1e-13.
- Dense 2048-dimensional analytical case: require agreement through eta 1e-9; explicit nonconvergence at more extreme eta is recorded as a limitation.
- Finite values and NaN/Infinity rejection; lower-triangle convention.
- MPI full-W pipeline and partially occupied metallic regression with unchanged tolerances.
- Si eta=0.001 update, independent Python W, external-W current, restart current, and measurable updated-W action using a delayed verification probe.
- Mandatory source and application-object audits before new executable use; recorded SHA256 on reuse.
- Fugaku production compiler flags and 4 MPI x 12 OMP execution. Preserve all failed evidence and existing production binaries/jobs.

## Reporting
Record accuracy, extra memory and runtime costs separately. Do not claim production convergence or acceleration relative to the former unrefined kernel. Store only logs and analysis in canonical calc2026; no GS/checkpoint payload collection.

## Prior hypothesis
Spectral anchoring was rejected in independent review because it changes H for imperfect eigenstates. An initially proposed exact-zero q=0 test was also rejected as an oracle when supplied levels differ from finite-H eigenvalues. Their diagnostic records remain preserved; acceptance uses the original-H reference above.
