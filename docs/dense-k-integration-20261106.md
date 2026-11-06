# Dense-k HSE integration and FFT consensus, 2026-11-06

The upstream merge retains the established dense-k EXX worker implementation (automatic blocked layout, symmetry-reduced rows, bounded density tiles and balanced strips). This follow-up fixes automatic FFT startup and accepts uneven strip decompositions. No additional algorithm-selection switches are introduced. Remote jobs and executables are unchanged; no push is performed.

## Selection based on Fugaku results

The saved original `exx-worker-performance-20261101/final-53272164/comparison.json` compares baseline with workspace/plan reuse + packed FFT + direct target transforms at identical Nr16 eight-step conditions. Numerical error is at most 1.021405182655144e-14, but maximum-rank worker-call median increases from 32.04582759994082 to 32.70011614006944 seconds (`median_speedup=0.985310390599675`, `faster=false`). It does not justify replacing the baseline.

The original monitoring record for job 53550165 reports EXIT0/PJM0 and four passing numerical gates. Relative to the reuse version, the later batched contraction has ABBA worker speedups 1.01094/1.01871, below its 5% requirement, and contraction speedups 0.972976/0.975822 (slower). Full metadata-inclusive independent analysis and final recovery remain pending; no metadata-inclusive speedup is claimed. These are separate experiments; their absolute times must not be compared across allocations.

Local synthetic 512/1000-k scalar/SOC comparisons favored the batched candidate on average (eight-case initial geometric speedup 1.145, final ABBA 1.171), but some pairs regressed. Grid 4^3, states 4/16, MPI4/OMP2; whole `apply_wave_workers` call, excluding outer service metadata and full RT. Direct-operator errors were below 1e-10. Those results do not override the Fugaku acceptance failure. Both candidate implementations and benchmark evidence are retained outside production source. The baseline is retained as the supported path; no universal fastest implementation is claimed.

## FFT fix

Rank-local trial timings previously selected inconsistent layouts. Before creating layout-specific plans, the native k communicator now collectively rejects local initialization/trial errors, checks uniform automatic/manual mode and explicit layout flags, and takes the maximum rank time for each layout. The existing 1e-5 second noise floor and 5% margin remain. A callback keeps the standalone serial kernel independent of MPI. Deferred-plan and early-error paths participate exactly once; later failures retain the outer collective gate.

## Uneven block work

Balanced strips assign quotient/remainder rank fanouts to rows, then quotient/remainder columns within each row. Fanouts differ by at most one, and within a row strip sizes differ by at most one. Across rows, differing fanouts can necessarily produce larger size differences: 164x4096 pairs on1280 ranks gives512..586 pairs/rank; 3x10 on7 ranks gives3..5. These are accepted. Every pair remains covered exactly once, shared column-group members retain identical start/count, and peak allocation retains the64MiB density cap. Invalid dimensions, insufficient ranks for the strip scheme, empty strips and over-cap allocations still reject; the existing blocked scheme supports other resource shapes.

## Validation and limits

Deterministic MPI2/4/8 auto/manual/error fixtures compare actual FFT operators to the direct oracle. Six real MPI4 application starts with OMP1/2/4 exercise automatic choices. Serial and MPI GNU15 HSE builds and source/application-object prohibition checks are required. Uneven MPI17/25 action probes and independent integer coverage tests verify the strip change; numerical, NaN/Inf and mismatch rejection thresholds remain unchanged. Final test results are recorded in the canonical evidence package. Fujitsu/GPU execution of this follow-up is not certified.

Final local results: 40/40 focused CTest cases passed, including MPI2/4/8 FFT consensus, integer layout/coverage oracles and MPI17/25 uneven action checks. The additional MPI9 uneven action matched the direct operator to1.2412670766236366e-16. Six real automatic-layout application starts passed. GNU15 serial/MPI builds and all application-object subtree audits passed. The final64-step HSE RT result is recorded in `validation/rt64.log` in the evidence package. No Fugaku speedup or full metadata analysis is inferred from these local tests.
