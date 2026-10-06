# H2 2x2x2 fragment-block scaling implementation plan

> 2026-10-06：以下の`benchmarks`は当時のローカル性能測定です。マージ対象から除外し、測定コードはブランチ外へ保存しています。

**Goal:** Remeasure PBEh40 GS preparation and 16-step impulse native RT using eight H2 per fragment core.
**Architecture:** Buffered periodic DC-SCF at 300 K with canonical full-support exchange; LCFO export reconstructed onto the global real-space grid. Compare fraction=1 and .999 on identical seeds, using the same spatial EXX backend. Production solver unchanged.
**Tech stack:** Existing Python MPI benchmark helpers, ScaLAPACK-enabled SALMON, rank wait4 RSS measurement.

## Approved design
User approved this workflow on 2026-09-28. Core size 16 cubed bohr, spacing .5 bohr, eight H2 (bond 1.4 bohr). Split directions use 4 bohr buffers; unsplit directions have zero buffer as required by the input validator. Fragment occupied capacity follows padded volume with four additional states for finite-temperature filling. PBEh cutoff remains 4 bohr to isolate geometry changes, not a physical convergence claim. DC MLWF disabled. PBE pre-SCF threshold 1e-4; target residual 1e-10. RT has fixed occupations, dt=.02, impulse=1e-4, 16 steps. No assigned RT temperature.

Weak arrays: 1x1x1, 2x1x1, 4x1x1, 8x1x1, 16x1x1, 2x2x1, 4x4x1, 2x2x2; one rank per fragment for GS, same total ranks for native spatial RT. Strong native RT: 4x4x1 blocks (128 H2), MPI1/2/4/8/16. Three sequential paired repeats, reverse mode order on alternate repetitions, prioritize minimum times and retain median/range. GS measured separately once per shape; do not imply repeated GS timings.

## Tasks
1. Create benchmarks/benchmark_h2_blocks/run.py reusing existing geometry, RSS wrapper and RT parser. Store generated inputs, binary/source/seed hashes and build configuration. Resume only verified completed runs. Validate atom/state/grid counts and finite-buffer geometry before MPI.
2. Run 1- and 2-fragment pilot through DC-SCF, LCFO, and paired RT. Require convergence, correct electron count, complete outputs and normal rank exits. Record support radius, norm loss, local FFT counts/points, localization fallback, current/energy differences. Do not label approximate RT numerically validated solely because it completed.
3. Run the full weak/strong matrix, summarize times/memory and compare the old narrow cells. Produce standalone plots and Japanese note with raw reproducibility archive. Review measured facts vs claims, commit without pushing.

## Progress / rulings
- Reuse existing feature checkout pbeh40-rvv10-water-md, clean at c5728525.
- Benchmark input/script changes are reversible instrumentation; verify by actual pilot and geometry assertions rather than implementation-mirroring unit tests.
- Fixed R remains static-only; RT comparison uses automatic .999 support as approved.

- Pilot ruling: gauss10 initialization failed the actual core-weighted state-capacity check before SCF in the buffered two-fragment cell. Gaussian centers concentrate in part of the periodic fragment; use random initial orbitals with broadly uniform core weights for all GS seeds. Preserve the failed pilot in work/h2-block-pilot. No occupation guard or solver code is relaxed.

- Initial 32-H2 adaptive RT failed to obtain a converged gauge at tolerance 1e-7, including a diagnostic rerun with maxiter=1000. A controlled rerun at the existing default tolerance 1e-6 and maxiter=1000 completed. Ruling: use 5/1000/1e-6 for every final RT case/mode; retain failed 1e-7 experiment separately, reuse only input/hash-identical GS preparations. No solver convergence guard is bypassed.
- Independent script review: fixed resume validation to include launcher, MPI version, pseudopotential and thread settings; current comparisons use columns 14–16 only.
- Task 1 complete: geometry generator, provenance/resume validation, shared-seed MPI runner, time/RSS/support/error analyzer implemented and reviewed.
- Task 2 complete: one- and two-fragment pilots converged; 32-H2 localization controls validated independently before final measurements.
- Task 3 running in work/h2-block-scaling-final; original 1e-7 results and localization probes preserved separately. Do not combine their RT times with final 1e-6 runs.

- User runtime steering: 64/128-H2 conditions get one repetition initially; <=32-H2 conditions retain three. All cases retain 16 impulse steps. Label large cases single observations, not three-run minima. Import only same-binary, same-input, same-seed completed RT rows; retain source metadata/hashes and rank data.

## Current execution handoff (after Jacobi correction)

- Authoritative active measurement: parent-workspace `work/h2-block-scaling-measured/results.json` and `.log` (log is `work/h2-block-scaling-measured.log`). Runner uses `--repeat 3 --large-repeat 1`, 40 RT jobs total, weak shapes first then strong MPI8/4/2/1. Each RT remains16 steps. Twenty completed same-binary samples were imported with input/seed/parser validation from `work/h2-block-scaling-jacobi`. Never merge the generic-localizer partial runs.
- Active executable: `work/pbeh40-scalapack-build/salmon`, includes solver commit f7c473ba. Metadata stores actual binary SHA256. Do not rebuild/replace it or overlap numerical jobs while the benchmark is running.
- Seven canonical DC GS payloads imported from `work/h2-block-gs-all`; original executable hash is retained per preparation. The eighth, 2x2x2, will be prepared in a fresh folder by the runner. The old incomplete 2x2x2 folder was intentionally interrupted and is not a valid seed. Canonical DC bypasses the modified localizer; 20 regressions and 27 exchange configurations passed.
- Latest complete same-binary weak pairs: 8H2/MPI1 full2.3447/adaptive2.7256 s;16H2/MPI2 7.9622/6.4952 s;32H2/MPI4 37.510/26.982 s (all minima of3);64H2/MPI8 220.51/109.03 s (one run). Max-rank RSS medians/single:171.0/170.9,298.0/291.5,555.7/565.2,1130.4/1135.5 MiB respectively. These are intermediate, not the full result.
- `benchmarks/benchmark_h2_blocks/analyze.py RESULTS OUTPUT --partial` explicitly makes an incomplete interim report and lists missing cases; without --partial it requires complete=true. Latest interim: `work/h2-block-interim`.
- Remaining deliverables after complete=true: run analyzer without --partial; verify all40 rows, seed hashes and normal completion. Review energy-width/current/MPI differences and localization/support counters. Finish `docs/h2-block-scaling-ja.md` (currently an uncommitted draft), including GS preparation table, weak/strong timings and RSS, speedups, limits of single measurements and finite-buffer initial state. Use actual measured values only.
- Explain the measured/local-code bottlenecks: source localization reduces per-pair FFT volume, but all Nocc^2 pairs remain; compact action uses batches of4 targets, limiting active spatial FFT workers. Native RT stores all occupied columns per grid rank in this nproc_ob=1 experiment. Do not claim linear weak scaling or memory savings without evidence.
- Archive results, source/build hashes, generated inputs, stdout, energy/current and per-rank RSS under `docs/benchmarks/2026-09-28-h2-blocks`; avoid copying large seed wavefunctions into Git. Preserve local seeds and reference their hashes/original paths. Include a manifest for every archived file. Render/inspect final scientific plot.
- Review the final note against data, commit completed artifacts (no push requested this turn), and report completion with main comparisons. If a job fails, preserve partial data, identify the failure and report it; never label the matrix complete or fabricate missing rows.


## User-directed stop before pair optimization

Stopped at 32 completed RT runs / 8 GS preparations. Weak scaling complete; strong only MPI16 complete. MPI8 full run interrupted and excluded. Reference archive: `docs/benchmarks/2026-09-28-h2-blocks-before-pairs/`. Frozen old executable: `../work/salmon-before-pair-generation`. Do not resume this old matrix; implement `2026-09-28-exx-pair-candidates.md`, validate, then remeasure.
