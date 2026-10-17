# Dense-k production source promotion (2026-10-17)

Promote the source used by the Fugaku dense-k runs onto `hse-symmetry-fft-reuse`, based on f69069552. Production inputs and binaries are unchanged. Source files were retrieved through authenticated Firefox Open OnDemand and checked against remote SHA256 values.

## Behavior

- Permit learned ACE with a symmetry-reduced representative mesh for k-only MPI, retaining all q points and the existing physical restrictions. Validate representative ranges and history/factor dimensions.
- Check finite positive normalized weights on `wtk(:system%nk)`. Symmetry leaves unused zero-padded tail storage; it must not be treated as an active k point. Numerical thresholds are unchanged.
- Scope skinny density-product BLAS to one thread and restore the parent thread configuration afterwards.
- Emit SVD reconstruction residual, left/right unitarity, and singular-value extrema; reject nonfinite or failed quality checks. Residual limit 1e-10, unitarity limit 1e-10 sqrt(n), singular-value floor smax*1e-12.
- Emit `STEP_WALL rank/step/seconds` on every rank around the time-evolution step. No added barrier; checkpoint and timing-output costs lie outside this interval. The sum of per-step rank maxima is a reporting metric, not an exact globally synchronized wall clock.

GS does not enable learned history. FFT reuse OFF still retains crystallographic symmetry. Learned runs use the established K8/nteach8/history3 automatic-warmup policy in their input/environment.

## Provenance and verification

`source-provenance.json` records the exact three source hashes, production binary hashes, and fresh verification. The diagnostic-only production artifact and the whole-step-timed artifact share the same exchange/history implementation; only the latter contains the timer additions. This commit incorporates both.

All three changed Fortran sources compiled with GNU MPI/OpenMP against existing module dependencies. The existing 2-rank history test passed spatial, orbital, and combined layouts; shifted/unshifted FFT reuse oracle tests passed. These are small code checks, not additional physical production calculations. No full application rebuild, GPU build, or OpenMP-disabled build was run during promotion. Existing unrelated LCFO working-tree edits were excluded from this commit.
