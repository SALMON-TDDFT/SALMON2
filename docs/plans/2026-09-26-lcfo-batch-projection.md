# Batched exchange FFT and compact projection

Goal: reduce FFT call overhead and zero-valued work in LCFO exchange projection without new physical approximations.
Architecture: each worker compacts the nonzero pair densities in a small target tile and executes a cached FFTW plan for the actual count, including tails. Each result column retains source accumulation order. Cache nonzero rows/columns and the conjugate-scaled core basis for projection, then scatter the rectangular product and retain the original Hermitian average.
Tech stack: Fortran, FFTW plan_many, OpenMP, BLAS, existing MPI.

The user approved the two proposed optimizations and sequential diamond benchmarking. Reuse the branch. Multi-node work remains deferred. Conservative alternatives are scalar FFT and dense projection; keep scalar selection available for validation and measure each component before claiming speedup.

1. Extend the direct-convolution test for batch sizes1/2/4/8, empty tiles, compact tails, repeated plan rebuilding, tiny complex products, OMP1/2/4. Add failing complex nonorthogonal compact-projection test against the old expression, including zero rows/columns and empty basis.
2. Implement cached per-worker per-count FFT plans, exact-zero compaction and counters. Default general Wannier path to scalar initially; choose LCFO batch size only after standalone timing. Add compact projection module and integrate fixed core plan.
3. Run direct tests, existing K-point/fractional Wannier and native MPI2/MPI4 regressions. Review races, normalization and Hermiticity. Measure batch1/4/8 on representative 32x16x16/60source/192target data sequentially.
4. Paired C64/8MPI and C128/16MPI R6 RT against beccb272, same GS, dt=.02, nt16, ACE1, OMP1/BLAS1. Use detailed timers for FFT action and projection to distinguish effects. One numerical job at a time.
5. Record timings/errors/limits, update both notes, commit and publish after verification.
