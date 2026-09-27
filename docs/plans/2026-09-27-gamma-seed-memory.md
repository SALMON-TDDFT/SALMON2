# Gamma seed peak memory implementation plan

Goal: remove LCFO root's rank-changing coefficient copy and avoid overlap between QR and SVD matrices, without a new physical approximation.

Architecture: add a Gamma-only 2D coefficient seed API. Keep the original general k-point API as the numerical reference. Use exactly the same conjugate-transposed coefficients, ZGEQP3 pivoting, selected rows, ZGESVD options/workspace size and singular-value threshold. Allocate QR arrays first; release QR columns/tau after obtaining pivots; only then allocate SVD matrices. LCFO no longer constructs a zero position array or reshape of the full coefficients. Root gather, global pivoted QR, and dense U remain scalability limitations.

1. Add test comparing new API with zero-k old API for complex full-rank, square, rectangular, singular and undersized coefficient sets; assert input preservation, U equivalence and unitarity. Verify missing API failure.
2. Implement new API in src/xc/hse_wannier_gauge.f90 and use it in src/xc/lcfo_rt_wannier.f90. Preserve fallback and binary snapshot format.
3. Run unit checks and rebuild. Measure isolated old/new seed RSS in separate sequential processes with the same generated coefficients and BLAS1. Retain original general API and original caller reshape in baseline.
4. Run native MPI direct-WF regression and a sequential C128 comparison with frozen before/after binaries and same GS. Compare initial coefficients/U/centers and RT current/density/energy.
5. Review, update development notes and 3D limitations, commit and push. Do not infer 8³ readiness from a small isolated RSS measurement or run large jobs.
