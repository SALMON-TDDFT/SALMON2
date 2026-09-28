# H2 2x2x2 fragment-block scaling implementation plan

**Goal:** Remeasure PBEh40 GS preparation and 16-step impulse native RT using eight H2 per fragment core.
**Architecture:** Buffered periodic DC-SCF at 300 K with canonical full-support exchange; LCFO export reconstructed onto the global real-space grid. Compare fraction=1 and .999 on identical seeds, using the same spatial EXX backend. Production solver unchanged.
**Tech stack:** Existing Python MPI benchmark helpers, ScaLAPACK-enabled SALMON, rank wait4 RSS measurement.

## Approved design
User approved this workflow on 2026-09-28. Core size 16 cubed bohr, spacing .5 bohr, eight H2 (bond 1.4 bohr). Split directions use 4 bohr buffers; unsplit directions have zero buffer as required by the input validator. Fragment occupied capacity follows padded volume with four additional states for finite-temperature filling. PBEh cutoff remains 4 bohr to isolate geometry changes, not a physical convergence claim. DC MLWF disabled. PBE pre-SCF threshold 1e-4; target residual 1e-10. RT has fixed occupations, dt=.02, impulse=1e-4, 16 steps. No assigned RT temperature.

Weak arrays: 1x1x1, 2x1x1, 4x1x1, 8x1x1, 16x1x1, 2x2x1, 4x4x1, 2x2x2; one rank per fragment for GS, same total ranks for native spatial RT. Strong native RT: 4x4x1 blocks (128 H2), MPI1/2/4/8/16. Three sequential paired repeats, reverse mode order on alternate repetitions, prioritize minimum times and retain median/range. GS measured separately once per shape; do not imply repeated GS timings.

## Tasks
1. Create testsuites/benchmark_h2_blocks/run.py reusing existing geometry, RSS wrapper and RT parser. Store generated inputs, binary/source/seed hashes and build configuration. Resume only verified completed runs. Validate atom/state/grid counts and finite-buffer geometry before MPI.
2. Run 1- and 2-fragment pilot through DC-SCF, LCFO, and paired RT. Require convergence, correct electron count, complete outputs and normal rank exits. Record support radius, norm loss, local FFT counts/points, localization fallback, current/energy differences. Do not label approximate RT numerically validated solely because it completed.
3. Run the full weak/strong matrix, summarize times/memory and compare the old narrow cells. Produce standalone plots and Japanese note with raw reproducibility archive. Review measured facts vs claims, commit without pushing.

## Progress / rulings
- Reuse existing feature checkout pbeh40-rvv10-water-md, clean at c5728525.
- Benchmark input/script changes are reversible instrumentation; verify by actual pilot and geometry assertions rather than implementation-mirroring unit tests.
- Fixed R remains static-only; RT comparison uses automatic .999 support as approved.

- Pilot ruling: gauss10 initialization failed the actual core-weighted state-capacity check before SCF in the buffered two-fragment cell. Gaussian centers concentrate in part of the periodic fragment; use random initial orbitals with broadly uniform core weights for all GS seeds. Preserve the failed pilot in work/h2-block-pilot. No occupation guard or solver code is relaxed.
