# Collinear independent FFT OpenMP implementation plan

**Goal:** Apply the approved bounded independent-pair FFT distribution to collinear exchange.
**Architecture:** A serial-plan/thread-safe-execute helper partitions existing storage without replicas. Cache full and Hermitian tail plans in the kernel. Preserve contiguous and strided layouts and reduced-q residue classes. Route reference, crystal-orbit and Hermitian-pair FFTs through it; MPI and BLAS remain outside application OMP regions.
**Tech stack:** Fortran, serial FFTW, OpenMP, existing MPI.

1. Add a failing direct-DFT regression for both layouts, reduced-q groups, partial batch, padding preservation, and OMP1/2/4.
2. Add reusable bounded FFT plan lifecycle; connect existing kernel/full-q orbit/pair routes with unchanged kernel thresholds and counters.
3. Build after source audit, audit objects, test frozen action OMP equivalence and reference/orbit/Hermitian agreement including empty ranks and finite rejection. Check SOC regressions.
4. Measure isolated tile timings and process RSS; save non-GS evidence to canonical Calc2026. No remote production changes or submissions.
