# Measured exchange FFT planning implementation plan

Goal: test cached FFTW_MEASURE plans for repeated exchange transforms without changing grid, pairs, support or integrator.
Architecture: add opt-in worker planning mode; include it in cache identity. Plan only serially on disposable worker scratch before pair densities are filled. Keep ESTIMATE default and use the existing MPI root validation/broadcast pattern for SALMON_LCFO_RT_FFT_MEASURE=0/1.
Alternatives: batching previously failed to improve synthetic action timing; sparse propagation needs a separate radius/gauge design. Measured planning is bounded and preserves the operator, but adds initial latency and machine-dependent choices.

1. Extend exact-pair convolution regression to both planning modes and toggles at fixed worker/batch sizes; observe missing-field compile failure.
2. Implement worker mode/cache key and collective RT setting. Verify OMP1/2/4, tiny/empty pairs, tails and translated multi-k reference.
3. Build, compare identical GS C128 with mode0/1 sequentially, measure preparation separately from repeated exchange and RT, check current/density/energy.
4. Review, document both costs and benefits, retain opt-in unless wider evidence warrants a default change. Push within existing user authorization.

Ruling: continue inline in the existing user-authorized development checkout; keep default behavior and numerical tolerances. No new physical approximation or multiple simultaneous SALMON jobs.

Validation ledger: missing fft_measure field produced the expected RED compile failure. Both planning modes passed independent convolution with OMP1/2/4. Native MPI2/4 measured-mode parity and invalid flags passed. Read-only review found no blocker and requested explicit cache-key and translated multi-k coverage; both were added and passed. Default remains0. Performance runs use the identical binary for mode0 and mode1 and alternate ordering by system size.
