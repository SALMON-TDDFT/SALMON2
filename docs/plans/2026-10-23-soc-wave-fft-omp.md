# SOC wave-worker FFT OpenMP implementation plan

**Goal:** Distribute independent spatial-pair k-mesh FFTs over OpenMP, as approved by the user on 2026-10-23.

**Architecture:** Retain the existing four-spin density tile and its 64 MiB cap. Partition its interleaved pair lanes into disjoint batches, with serial FFTW plans created and destroyed outside OpenMP. Execute forward FFT, kernel multiplication and backward FFT for independent spin/batch tasks; MPI and contractions remain outside this region. No per-thread density replica and no threaded FFTW library. Preserve a one-thread route through the same code.

**Tech stack:** Fortran, FFTW plan-many, OpenMP, existing MPI workers.

Alternatives considered: threading each small FFT adds overhead; copying whole tiles per thread multiplies memory. Disjoint plan-many batches retain bounded storage and expose the independent work directly.

1. Add a focused production-module probe: compare against explicit multidimensional DFT, test partial and zero rows/columns, 1/2/4 threads, and repeated executions. First confirm missing feature fails.
2. Implement plan lifecycle and disjoint lane execution in `src/xc/exx_spinor_wave_fft.f90`; connect `exx_spinor_wave_workers.f90` without changing MPI placement, numerical thresholds or tile sizes.
3. Register guarded CMake regression, source/object audit, build with bounds checks and without OpenMP, run focused tests and existing owner/helper regression when feasible.
4. Record timing and measured memory separately from scientific/full-scale performance acceptance. No existing jobs or binaries changed; no remote submission implied.
