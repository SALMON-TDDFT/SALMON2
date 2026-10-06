# FFT packing and streamed ACE OpenMP implementation plan

**Goal:** Reduce measured serial overhead in Cartesian FFT and streamed ACE action.
**Architecture:** Retain Cartesian ownership, FUNNELED MPI and factor summation order. Use existing work buffers; bypass send/recv for an unpartitioned axis. Parallelize independent lines/target orbitals, with no per-thread mesh allocation.
**Tech Stack:** Fortran, OpenMP, FFTW, MPI.

Approved design: user approved FFT data movement reduction and ACE inner-loop threading on 2026-09-30. Baseline Si4 MPI2 OMP2 sampling found ~65% main-thread stack occupancy in mesh_transform, ~17% in orbital_ace_apply and ~90% worker idle time.

1. Modify src/xc/fftw_blocks.f90: unpartitioned axes pack directly to work and unpack from work; parallelize line pack/unpack and distributed reorder loops. Preserve padding correctness and MPI only on the master thread, with barriers publishing received data. Do not allocate a full-mesh replica.
2. Modify src/xc/exx_orbitals.f90: parallelize streamed dense ACE overlaps across target orbitals, scalar grid summation per orbital; parallelize action updates across target orbitals. Retain outer factor order and collective order.
3. Build work/hybrid-rules-build. Run 657_unit_fftw_blocks roundtrip/reference suite for MPI1/2/4/8 and OMP1/2/4, and distributed metric tests for MPI1/2/4 and OMP1/2/4. Include no-OMP build. Existing numerical tests should pass before and after; the baseline performance deficiency is measured, not a correctness failure.
4. Snapshot new binary and run identical Si4 MPI2 short RT for OMP1/2/4, separately from production load. Compare current/energy/norm with immutable baseline, log peak RSS and timing. Sampling overhead must not be presented as clean baseline runtime.
5. Record measured outcomes and limitations. Do not push or replace the running production binary.

## Verification and measured outcome

Final FFT team retains only explicit barriers after MPI; omp do end barriers preserve array dependencies. Timer updates execute only on master. Two intermediate variants were inspected; their timings are not the final comparison.

GNU MPI build and no-MPI build succeeded. Final FFT regression:20 passes (MPI1/2/4/8, OMP1/2/4, padded line counts and cache changes); ACE metric/action:324 passes (MPI1/2/4, OMP1/2/4). Separate compilation without enabling OpenMP passed FFT MPI1/2/4 and ACE MPI1/2 probes. Diff whitespace/132-column checks passed.

Si4 diagnostic:32Si, mesh64x16x16, MPI2 spatial2x1x1, dt0.08au,16steps. GS reused; same binary baseline a4dd8744. Runs isolated by suspending production ranks, no profiling/build during final timed cases. Single measurements, not statistical speedup estimates.

|Version|OMP|RT s/step|Peak MiB/rank|
|---|---:|---:|---:|
|baseline-clean|2|3.301687|273.58|
|final1|1|2.877875|268.20|
|final2|2|2.901625|269.97|
|final4|4|3.555625|270.83|

Final OMP2 RT time is12.1% lower than baseline OMP2, but OMP1 andOMP2 are effectively tied; OMP4 regresses in this small case. Do not claim good thread scaling. The largest gain is data movement/work reduction. Maximum current difference across final runs is7.96e-17au; energy column difference1.01e-12 (output units), printed norm differences zero. Not a long-time laser accuracy check. No additional mesh array allocation was introduced.

Raw evidence and reproduction: work/si-hse-hotspot-mpi2/verify-final.py, no-omp-check.py, final-comparison.json and final1/2/4 folders. Production8-cell binary unchanged, ranks resumed, excluded pause/build windows in external-load-events.json.

Final diagnostic binary SHA256: 5a31e64a99110d32edcda6741637e6d934286dab411980fa1713a06bed68a73e

## Follow-up: remove per-line peer division

Split packing into complete owner/slot rectangles and a short tail handled by
master. Both packing and restoration derive ell=1+owner+peers*slot; the source
no longer divides or takes a remainder by peers for each line. Compute
full_slots=lines/peers once per axis transformation. Padding is separate and
empty when lines is divisible by peers. collapse(2) preserves work distribution
when the number of peers is smaller than the OMP team. No additional barrier,
parallel region or mesh buffer is introduced. Compiler-generated collapse
indexing is not asserted to be division-free machine code.

Build,21 MPI/OMP FFT cases and3 no-OMP cases passed; added a small cross-section
case whose single-field cache test has fewer lines than peers. Raw logs:
work/si-hse-hotspot-mpi2/owner-{build,fft,no-omp}.log. Performance has not been
remeasured for this follow-up. Production ranks resumed with their old binary.

Quotient reuse follow-up: all seven ell-to-band/line index conversions now
compute the quotient once and derive the remainder by integer multiply/subtract.
The six coordinate conversions use the same pattern. Temporary q is private in
both OMP regions. MPI/OMP21 cases and no-OMP3 cases passed again; see
quotient-{build,fft,no-omp}.log. This is source-level common-subexpression
elimination; no additional performance gain is claimed without measurement.

## Approved follow-up: fused local lines and packed ACE overlaps

Implemented independent gather/FFT/scatter per unpartitioned line and parallel packed ACE overlaps. No change to overlapping sparse action accumulation. Regression and benchmark results are in docs/reports/si4-fused-fft-2026-09-30/report-ja.md. Compared FFT batches and ESTIMATE/MEASURE in a standalone probe; retain production ESTIMATE until end-to-end batch validation.

## End-to-end local FFT batching

Integrated plan_many for unpartitioned strided axes using existing workspace, scalar tails and cached forward/backward plans. Compared widths 4 and 16 with MPI2/OMP2 and OMP4. Keep width 4 and ESTIMATE: 37.082 s per 16 RT steps versus prior fused 38.466 s at OMP2; single-run evidence only. Width 16 gave 38.527 s. FFT regressions 21 and no-OMP/bounds-check regressions 3 passed. Four Si cases completed with current differences below 8e-17 au and no energy difference at output precision. Production resumed unchanged. See docs/reports/si4-batched-fft-2026-09-30/report-ja.md.
