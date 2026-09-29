# Distributed Gamma MLWF Implementation Plan

**Goal:** Extend the approved dense-matrix memory reduction to native Gamma MLWF, retaining its transport, seed and Jacobi physics.

**Architecture:** Two-dimensional block-cyclic gauge/link/SVD arrays over the same spatial-orbital communicator used by ACE. Temporary BLACS context per refresh; copied states own arrays and mapping only. Stream grid columns to build overlaps and to rotate mesh functions in either direction. Distributed Jacobi uses the same pair order and 3x3 rotation rule, communicating pairs of rows/columns without any N-square gather. ScaLAPACK-free and single-rank paths remain available.

**Alternatives:** Root-only localization reduces replication but retains one N-square bottleneck. Removing only temporary copies reduces constants but leaves replicated scaling. Use distributed tiles as already authorized; measure numerical invariants, not bitwise eigenvector identity.

## Tasks
1. Write a Gamma localization probe comparing spread, gradient, source/projector, phase-aligned transported orbitals and exchange invariants across spatial/orbital/mixed layouts. Demonstrate missing distributed API before implementing.
2. Add `exx_distributed_gauge.f90`: layout, overlap, SVD transport, projected seed, Jacobi minimization, distributed functional and forward/adjoint rotation. Scalars only for IEEE inquiries.
3. Integrate optional combined communicator into `exx_spatial` and native calls, including canonical cleanup and adaptive inverse gauge action. Preserve fallback and copied state lifetime.
4. Build and test MPI 1/2/4 and no-ScaLAPACK; existing source-ACE native RT and Gamma regressions. Record tile allocation versus RSS distinction. Review numerical and collective edge cases.

No new memory limits, localization approximation, user namelist, or production spectral runs.

## Verification record
- Missing `comm_matrix` compile failure observed before implementation.
- Implemented tiled overlap/SVD/seed/Jacobi/forward and inverse rotations; native and manual-probe dependencies updated.
- 8 MPI layouts and the same 8 no-ScaLAPACK layouts pass, including transport, reset, retained gauge, deep copy and non-finite previous state.
- Empty-rank test exposed PZGESVD INFO=N+1 for near-degenerate singular values. Diagnostic logs established 1e-15 process-grid differences; added distributed Newton-Schulz retry with unchanged rank-loss threshold and explicit orthogonality acceptance. Failed test then passed.
- Review identified communication as a performance concern. Small 2x2 blocks now screen out rotations before row/column traffic; rank 0 broadcasts one rotation decision for collective consistency. Independent review found no important defect in both changes.
- MPI/ScaLAPACK and non-MPI HSE builds pass. Native source-ACE/pair-screening RT tests: 10 passed.
- No claim of whole-process RSS reduction or production speedup; no new production spectral run or push.
