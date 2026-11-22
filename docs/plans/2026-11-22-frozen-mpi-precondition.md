# Frozen W k/q MPI and preconditioner

Use streamed source/target frames, never replicated all-k density matrices. Keep legacy service as reference.
1. Add failing native MPI layout/action-failure and preconditioner residual tests.
2. Distribute q groups and source k batches; each worker solves locally after two frame transfers.
3. Import single-k HSEX action and rebuild nonlocal PP phases on scratch slots. Retain A and snapshot potential.
4. Add optional kinetic inverse basis expansion. Certify the original equation with unchanged residual thresholds.
5. Audit source/objects, test layouts and failures, benchmark small pumped Si on Mac at 8 MPI, 1 OMP/BLAS.
6. Record W differences, solver counts, measured time and memory; preserve existing jobs and GS data.

q response assembly and the existing Gamma closure remain unchanged to isolate numerical effects.

## Completed validation

- Added streamed k/q workers, source/target PP/HSEX frames, optional kinetic basis expansion, and one-pass H-basis orthogonalization. Original residual checks remain unchanged.
- Legacy expansion order remains the default reference. Native dense-oracle, coherent, finite/NaN/Inf,1/4MPI layout/failure, Gamma and OMP checks passed with OpenBLAS; Accelerate model checks also passed. Nr12³ nonorthogonal/nonzero-k+A inverse-kinetic plane waves match their analytic oracle.
- Full HSE-on and HSE-off builds passed with source/application-object prohibition checks retained. All final real-Si runs used one immutable audited SHA at8MPI×1OMP/BLAS and completed800steps with one1fs update and4occupied orbitals. No GS/restart payload was collected.
- Independent review findings (undefined failed-preconditioner output, unvalidated+U frame bounds, general Nr support, line length) were fixed and rechecked.

- final3-legacy: W 368.501265s, relative difference 0.000e+00, max 0.000e+00Ha, residual bound 8.020e-09, peak rank RSS 625.7MiB.
- final3-q1: W 36.957427s, relative difference 0.000e+00, max 0.000e+00Ha, residual bound 8.020e-09, peak rank RSS 643.5MiB.
- final3-q2: W 37.316282s, relative difference 0.000e+00, max 0.000e+00Ha, residual bound 8.020e-09, peak rank RSS 647.4MiB.
- final3-kinetic: W 40.035808s, relative difference 5.717e-14, max 7.938e-15Ha, residual bound 8.020e-09, peak rank RSS 634.5MiB.
- final3-q2-repeat: W 39.374695s, relative difference 0.000e+00, max 0.000e+00Ha, residual bound 8.020e-09, peak rank RSS 636.6MiB.

Kinetic preconditioning does not reduce the full512-dimensional spaces in this fixture and is left opt-in. Distributed solves plus one-pass expansion give the benefit. Root q-matrix storage/assembly and larger-grid memory/cost remain limitations. The small real-Si fixture uses identical finite-shell Gamma rules solely for numerical/performance comparison; physical optical-limit convergence is not established.

Results: canonical calc2026/si-frozen-mpi-precondition-20261122. No commit requested.
