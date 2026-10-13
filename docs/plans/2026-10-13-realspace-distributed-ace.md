# Real-space distributed ACE implementation plan

**Goal:** Keep real-space orbitals and ACE factors distributed, including the orbital metric and learned-factor alignment. Never project RT dynamics onto an LCFO basis.

**Architecture:** Reuse the existing combined spatial/orbital communicator and ScaLAPACK block-cyclic metric. Extend distributed matrix operations to factor alignment, then remove the orbital-distribution restriction only after complete end-to-end validation.

**Tech stack:** Fortran, SALMON communication wrappers, MPI, BLAS, ScaLAPACK.

## Current implementation verified by static inspection

- `src/xc/exx_distributed_metric.f90`: constructs local metric tiles from orbital columns; PZHEEV supplies the inverse-square-root factor. The metric and eigenvectors are distributed. Replicated eigenvalue and reduction buffers are only O(Nocc).
- `src/xc/exx_orbitals.f90`: `orbital_ace_build_core` selects this route when `comm_matrix` is supplied and ScaLAPACK is enabled with multiple ranks. `distributed_metric_rotate` constructs local real-space factor columns.
- Dense `orbital_ace_apply` already streams one factor column and local target overlaps. It does not allocate the full overlap matrix. Communication is currently per column.
- `src/xc/exx_native.f90`: learned-factor configuration explicitly requires `info%isize_o == 1`. Removing that check alone would be incorrect.
- `src/xc/exx_factor_history.f90`: `factor_align` replicates dense overlap, SVD vectors, rotation, and an 8*Nocc**2 complex workspace. This is the principal learned-ACE matrix-memory obstacle.
- Existing `orbital_hermitian_action` and gauge processing must also be traced in each supported route. Distributed metric construction alone does not prove bounded memory for the entire reconstruction.

## Chosen design

Keep Phi, W, and Xi as local grid rows and local orbital columns. Build C = X^dagger Xref directly into block-cyclic tiles. Use distributed complex SVD C = U Sigma V^dagger and form R = U V^dagger with distributed GEMM. Stream factor panels and rotation tiles to form X R without gathering C, U, V, R, or Xi. Preserve the existing singular-value rejection threshold.

A polar decomposition is an alternative, but squaring the overlap to form C^dagger C worsens conditioning. A root-only SVD removes replication but leaves a root memory bottleneck. Neither is the primary design.

Keep three local history factors. Reduce the six regression scalars across the combined spatial/orbital communicator once, avoiding duplicated spatial contributions. Synchronize readiness and failure across all relevant groups. Prediction remains an elementwise operation on local Xi. Semilocal potentials and RT propagation remain real-space operations.

Use bounded panels for ACE application. A panel width b requires O(b*Nr_local + b*Ntarget_local) buffers. It must not scale as the full Nocc**2 overlap. Uneven local orbital counts and zero tile extents must be supported.

## Implementation sequence and acceptance criteria

1. Add a standalone development MPI action check outside testsuites. Compare distributed metric construction and application to a dense reference on deterministic complex orbitals. Cover 1, 2, and 4 ranks, spatial-only, orbital-only, combined layouts, and uneven orbital counts. Check action relative error < 1e-11 and Hermiticity.
2. Introduce a reusable distributed matrix-layout helper with explicit ownership and BLACS lifetime. Reuse it for metric and alignment instead of duplicating descriptor logic. Test ownership uniqueness and complete coverage.
3. Implement distributed Procrustes alignment. Test arbitrary unitary factor rotations, singular overlap rejection, and non-divisible tile sizes. Compare aligned operators and residuals with the serial SVD reference. Query ScaLAPACK workspaces and record their actual sizes before claiming memory scaling.
4. Add optional alignment and reduction callbacks to history acceptance. Keep existing serial interfaces valid. Change teacher reductions to the combined communicator when orbital distribution is active.
5. Audit materialization, expectation values, readiness, and checkpoint metadata for local orbital ranges. Only then allow orbital MPI in factor_configure. Validate uninterrupted and resumed results outside testsuites, preserving the user requirement not to add restart tests there.
6. Batch ACE application into bounded panels. Verify numerical action and collective-count reduction against the original column algorithm. Keep reduction/accumulation order where practical.
7. Record actual route, global/local orbital counts, local matrix sizes, queried solver workspace and peak RSS in isolated benchmarks. Run HSE and non-HSE builds, existing hybrid regression tests, and a short real-space laser comparison. Do not overlap heavy verification with the active production calculation.

## Memory targets

For fixed orbital count and process count P, matrix storage should be O(Nocc**2/P), plus solver workspace measured explicitly. Orbital/factor frames remain O(Nr*Nocc/P). History retains three additional local frames. Replicated O(Nocc) spectra and scalar regression state are acceptable. No claim of end-to-end scaling is made until all temporary allocations and communication buffers are measured.

## Preservation

Keep the pending fixed-LCFO RT removal intact. Do not modify frozen production binaries or existing results. No new long calculation, commit, or push is implied by this plan. This document records the original plan. Gamma orbital-distributed learned ACE has since been implemented and validated with small MPI and real-space RT tests; see docs/changes/2026-10-13-distributed-learned-ace-ja.md. Large-system RSS and timing remain unmeasured.
