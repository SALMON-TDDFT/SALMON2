# Exact zero-pair FFT elimination

Goal: avoid screened-exchange FFTs whose pair density is exactly zero, with no new threshold or physical approximation.
Architecture: after forming conjugate(source)*target in each thread-local FFT buffer, skip only all-zero products. Replace redundant per-pair target magnitude scans. Count total/executed FFT pairs with 64-bit counters and an OpenMP reduction, and expose the counts in LCFO refresh diagnostics. Preserve source accumulation order and existing kernel/normalization.

User authorization: continue other performance optimization; multi-node experiments are explicitly deferred. Reuse the existing development branch. Batch FFT and orbital-distributed exchange construction are alternatives, but defer them until this smaller exact optimization is measured.

1. Add a failing compiled test with disjoint, overlapping, zero, and tiny nonzero complex source/target support. Compare action against a direct periodic convolution independently reconstructed from the reciprocal kernel; assert retained pair count and Hermiticity. Test OMP1/2/4.
2. Add exact product-zero screening, counters and LCFO diagnostics.
3. Run native Wannier multi-k/fractional/ACE tests, native LCFO MPI2/4 regression and build. Review for OpenMP races and hidden thresholds.
4. Sequential paired C64/8MPI and C128/16MPI R6 benchmarks against 35d1d38e, same GS, dt=.02, nt16, ACE1, OMP1/BLAS1. Compare action diagnostics, current, density, energy, timings. No concurrent numerical jobs.
5. Record results in both notes, commit/push after validation.
