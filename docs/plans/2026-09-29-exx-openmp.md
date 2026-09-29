# EXX OpenMP implementation plan

User approved parallelizing FFT batches and independent exchange loops.
Keep MPI and FFTW planning/destruction on the caller thread. Split independent
1D transforms into balanced contiguous chunks; each worker has its own serial
FFTW plan and disjoint work slice. Cache includes requested worker count.
No new FFTW thread-library dependency; OpenMP-disabled build retains one worker.
Parallelize FFT packing/copy and exchange pointwise loops with private indexes.
Do not parallelize source loop containing MPI or use nested FFTW threads.

1. Add FFT regression for thread-count changes, partial batches, forward reference,
   inverse normalization, both layouts and MPI1/2/4. Demonstrate serial baseline gate.
2. Implement chunked plans and independent loops. Build MPI and no-OMP variants.
3. Run FFT/rVV10 and hybrid source-ACE regression, compare OMP1/2/4 outputs.
4. Time same Si short RT with OMP1/2/4, reuse GS, then resume laser comparison.

Long laser run stopped before modifying/building to avoid competing timing jobs.
Preserve completed PBE and partial old HSE output; do not overwrite them.

## Final decision and verification

Batch16/64 reduced some synchronization cost but regressed whole Si RT. Retained
batch4 and the original low-overhead serial path for small FFTs. Larger transforms
use a persistent OMP team, >=65536 complex points per worker, and primary-thread
MPI calls. Pointwise work uses vector expressions with bounded teams for>=262144
components. No physics approximation or new input knob.

Large FFT MPI4 OMP1→4 roundtrip14.096→8.841ms, local FFT5.129→1.671ms.
Small Si RT remains approximately1.1s/step; no demonstrated OMP speedup there.
FFT/large exchange/empty targets/rVV10 and13 hybrid RT regressions pass.
See docs/exx-openmp-ja.md and docs/reports/exx-openmp for full measurements,
failed tuning experiments and the pre-existing HSE GS test convergence caveat.
