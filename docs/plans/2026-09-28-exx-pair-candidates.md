# Bounded spatial EXX pair candidate generation

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Reduce expensive pair generation, pair FFTs and temporary memory for localized HSE/PBEh native mesh exchange, then repeat the approved H2 benchmark.

**Architecture:** Use localized targets without truncating them further. A sparse catalogue of spatial-block target amplitude maxima supplies candidate targets for each finite-support source. Certified absent pairs consume the existing exchange-action error budget; candidate pairs retain the existing tighter bounds, Hermitian correction and unscreened ACE fallback. Process only selected columns in compact FFT work.

**Tech Stack:** Fortran/MPI/FFTW, existing ScaLAPACK build, Python regression drivers.

## Approved scope / reference

User approved stopping the current measurement and completing pair-generation optimization before remeasurement. The preceding design called for localized targets, locality-based pair candidates and compact work arrays; dense RT/ACE representation changes are a later phase. Current 32 completed RT runs and 8 GS preparations are archived under `docs/benchmarks/2026-09-28-h2-blocks-before-pairs/`. Interrupted MPI8 is excluded. Frozen executable: `../work/salmon-before-pair-generation`.

Existing isolated worktree/branch `SALMON2-pbeh40-rvv10` / `pbeh40-rvv10-water-md` is reused. No push requested.

## Mathematical contract

For source q, convolution norm lambda=max|Fourier multiplier|, and target t bounded by a throughout the nonzero support of q:

`|| q * C(conjg(q)*t) ||_2 <= max|q| * lambda * ||q||_2 * a`.

Use physical mesh-weighted L2 norms and a conservatively rounded threshold. Each omitted ordered pair receives at most `budget=tolerance/(Nsource*sqrt(Ntarget))`; summing per target and then in Frobenius norm stays within tolerance. No 99.9%-norm argument substitutes for this bound. Zero budget excludes only exactly zero pairs. Nonfinite or underflow-ambiguous bounds retain work.

Index blocks are a computational partition, not a new physical cutoff. Build a sparse block catalogue retaining target maxima above the smallest source threshold. Source bounding boxes cover every nonzero source point (including periodic wrap; a loose box is safe). The target lookup unions only entries in intersecting blocks. Avoid a dense Nsource×Ntarget candidate table and avoid evaluating grid products for noncandidates. Exact native target orbitals are reconstructed by the unitary gauge; fractional-occupation cases stay on existing supported paths.

## Task 1: Sparse candidate catalogue and independent tests

Files: create `src/xc/exx_pair_candidates.f90`, `testsuites/unit_hse_wannier/pair_candidates_probe.f90`, and Python driver; add module to `src/xc/CMakeLists.txt`.

Interface: build catalogue collectively for local target columns over comm_r; query a source bounding box and amplitude threshold to return unique target indices. Store CSR block entries and per-target marks, never an Nsource×Ntarget matrix. Keep counters for catalogue entries and candidate examinations. Metadata may grow dense for physically delocalized input; do not impose a memory cap or silently discard amplitudes.

1. Write and run RED tests against independent brute-force support overlap: separated supports, boundary wrap, complex tails, zeros, unequal/empty orbital partitions, MPI1/2/4.
2. Add catalogue construction from block maxima with sparse collective metadata exchange and conservative source bounding boxes. No full-grid gathers.
3. Verify candidate completeness, uniqueness and linear candidate count on a replicated localized chain; verify exact-zero and finite-tail bounds. Expected: all tests pass across ranks.
4. Commit module/tests and record evidence.

## Task 2: Integrate bounded candidates and PBEh

Files: `src/xc/hse_spatial.f90`, `src/xc/hse_native.f90`, `src/io/inputoutput.f90`, `src/io/salmon_global.f90`, existing pair-screen tests.

Consumes Task1 catalogue. Produces screened action with unchanged error-budget/ACE acceptance semantics.

1. Extend RED pair action tests to Coulomb omega=0; assert bounded action error, compact/global parity, zero tolerance, tiny finite values, diagnostics parity and MPI/orbital layout parity. Add an independent chain case proving distant pairs avoid full product generation.
2. Build/source-select catalogue only in screening modes. Use broad-phase absent bounds plus existing fine pair bounds. Keep `diagnose` exact and `off` unchanged. Generalize HSE-only input guard and localized action omega to the actual functional. Canonical DC fragments remain unlocalized.
3. Allocate compact target/action work only for surviving columns; distribute compact pair batches across available spatial workers without a fixed four-worker bottleneck.
4. Run bounded-pair and exchange/ACE suites. Expected: numerical bounds and ACE interpolation hold; no regression in off/diagnose.
5. Commit integration/tests.

## Task 3: Native RT validation and benchmark configuration

Files: `testsuites/unit_pbeh_rvv10/test_pair_screen.py`, `testsuites/benchmark_h2_blocks/run.py`, analysis/README, `docs/inputs/exx-mlwf.md`.

1. Add RED same-seed PBEh native RT on/off/diagnose, full/.999 supports, MPI1/2 and orbital split validation; test fixed-radius guards remain unchanged.
2. Verify actual small H2 seed pair reductions and accepted bound/fallback behavior before launching the long matrix. Compare energy/current, not merely successful exit.
3. Extend benchmark to explicit screening controls and preserve them in reuse checks/metadata. Compare full unscreened reference vs .999 screened with tolerance=1e-6 initially, reporting support error separately from screened-vs-unscreened .999 error. If no useful reduction occurs, inspect bounds and locality; never loosen the budget silently.
4. Run ScaLAPACK build and pertinent regression suites; fresh read-only whole-change review, fix material findings with RED/GREEN tests.
5. Remeasure same approved geometry, 16 steps, small cases3 repeats, >=64H2 one; reuse validated canonical DC seeds, never old RT timings across binaries. Save full raw evidence and Japanese report, clearly distinguishing residual dense RT memory from reduced exchange work.

## Execution ledger

- Planning: old benchmark stopped by user; all child processes exited; old binary and completed data frozen.
- Ruling: preserve dense native RT wavefunctions and ACE factors for this phase. Pair/workspace optimization is measurable independently; total-memory linear scaling is not claimed.
- Ruling: block maxima bound candidates instead of treating nonoverlapping 99.9% spheres as exact zero. This retains explicit error control without changing target wavefunctions.

## Approved amendment: local-support ACE before remeasurement

The bounded untruncated-target pilot at64 H2 retained4096/4096 pairs, with32,768 catalogue entries: off10.70s vs on11.389s for2 steps, peak1139.67 vs1176.78MiB, no fallback, current difference1.02e-20 and energy difference9.95e-14Ha. Synthetic sparse tests pass, but this does not provide a useful production speedup. Do not launch the long benchmark yet.

User explicitly approved: use the99.9%-masked MLWFs as ACE construction vectors and exactly exclude nonoverlapping supports; retain real-space RT wavefunctions and compare against the existing method.

### Task 4: Source-support ACE (before Task3 remeasurement)

Add opt-in `exx_ace_support='source'`, default `'occupied'`. Native source-support ACE forms S from the already masked source, W=K_S S, then ACE from(S,W). Nonorthogonal S is allowed by the existing Gram-based ACE construction; retain its strict Hermitian/positive/conditioning checks. Use zero-budget candidate screening to remove exact disjoint supports only. Do not silently reinterpret the finite pair-action error budget: reject combining source-support ACE with explicit pair screening initially.

Apply the resulting ACE to original occupied mesh orbitals for cached action and exchange energy; time propagation remains on the mesh. Alias S into the three-dimensional ACE interface where possible instead of adding another full mesh copy. If the source ACE is inadmissible, recompute the existing occupied-vector ACE and report fallback. Initially restrict source-support mode to fully occupied native fixed-ion RT (existing DC preparations remain unchanged). Full support must reproduce occupied-vector ACE.

Tests: RED input/control and full-support equivalence against old executable; independent nonorthogonal ACE interpolation; finite-source overlap pairs with complex orbitals; MPI/orbital layouts; same-seed16-step native RT. Compare .999 source-support against .999 occupied ACE and full support, report approximation differences separately from exact-zero pair pruning. Rerun short64 H2 pilot and inspect pair count/time/RSS before full matrix.

Benchmark: add `--ace-support source` for adaptive runs only; full reference stays occupied. Store controls in metadata/reuse checks and report local-support ACE acceptance/fallback. Native dense RT and ACE factors remain a later storage-optimization phase; this amendment must not claim linear total memory.

## User-directed reference reuse

User subsequently instructed reusing existing Full results instead of rerunning
them. This supersedes the earlier same-binary reference requirement. Reuse the
16 completed Full runs from the stopped baseline, retain original executable and
measurement provenance, and check exact inputs/seeds/output. Import completed
new-source-ACE adaptive runs without rerunning them. Only missing Full strong
cases (MPI1/2/4/8) require new calculations. Report cross-executable/time speedups
as historical comparisons, without presenting them as same-binary controlled runs.
