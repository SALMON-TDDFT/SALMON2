# H₂ PBEh(40): impulse RT for 16 steps

Uses the same shapes, geometry, grid, cutoff and y/z spatial layouts as
[the GS scaling harness](../benchmark_h2_pbeh/README.md), but measures real-time
propagation. Production SALMON code is unchanged.

```sh
python3 testsuites/benchmark_h2_pbeh_rt/run_rt.py \
  --binary /absolute/path/to/salmon \
  --output /absolute/path/to/new-results --repeat 3
```

## Initial states and propagation

The previous GS measurements did not save wavefunctions. The supported native
PBEh RT reader uses DC-LCFO data. Therefore each shape is prepared once with a
**single fragment covering the entire cell**, zero buffer, and all occupied
states retained. This avoids a truncated-fragment approximation. Preparation
uses at most eight spatial ranks and the same gauss10 / alpha_mb=0.1 settings
as the final GS benchmark, with zero-temperature occupations and convergence
threshold 1e-10. Preparation time and its process memory are not RT measurements.

The resulting state is reconstructed onto the real-space mesh by the native
streaming reader; RT does not retain a projected LCFO basis. Every repetition
and every strong-scaling rank count reads the identical shape-specific payload.
Hashes before and after RT verify the payload is unchanged. The payloads are
archived to permit identical initial states in future measurements.

- PBEh(40), no rVV10; fixed ions, fixed occupations, Gamma point.
- H₂ bond 1.4 bohr, basic cell 8³ bohr, mesh spacing 0.5 bohr.
- Coulomb interaction cutoff fixed at 4 bohr for every shape.
- Full EXX source support; MLWF interval/maxiter/tolerance 5/100/1e-7.
- Default native Taylor4 with ACE and predictor/corrector.
- x-directed impulse `e_impulse=1e-4` a.u.; `dt=0.02` a.u.; **16 steps**
  (total time 0.32 a.u.). Energy and current output every step.
- Fresh processes; three sequential repeats, no intentional overlap of jobs.
- OpenMP and BLAS each one thread. Default OpenMPI launcher `--bind-to none`.

Weak shapes: 1×1×1, 2×1×1, 4×1×1, 8×1×1, 16×1×1, 2×2×1, 4×4×1, 2×2×2,
with one MPI rank per H₂. Strong: 4×4×1 at MPI1/2/4/8/16. The common MPI16
case is shared: 12 distinct cases, 36 RT jobs, plus eight GS preparations.
MPI1/2/4/8/16 uses `nproc_rgrid=1,1,1 / 1,2,1 / 1,2,2 / 1,4,2 / 1,4,4`
and one orbital group. All ranks retain all occupied columns, so weak scaling
holds owned grid size, not orbital-pair work or wavefunction memory, constant.
MPI1 uses the legacy single-rank Wannier backend; other ranks use spatial EXX.

## Measurements and checks

Primary time: **maximum rank `rt iterations` timer for all 16 steps**, including
per-step potential/density work and current/energy output. The independent
`time propagation` timer isolates the instrumented propagation portion. Root
total calculation time includes initialization; launcher wall includes process
startup/teardown. Avoid adding nested timers together.

The user reports competition with other applications, particularly at MPI16.
Use the **minimum** observed time among the three runs as the main comparison,
and retain medians and ranges. Since every run has exactly 16 steps, loop-time
and per-step speedups are identical. Minima do not prove absence of contention.

Memory is the **median of the maximum-rank lifetime peak RSS of each RT job**,
measured from the SALMON child with `wait4`. It includes initial-state loading,
mesh reconstruction and library state, but excludes the separate GS process,
Python wrapper and MPI launcher. It is not a phase-isolated RT-loop peak and
not simultaneous whole-node RSS. No memory ceiling is imposed.

Require 16 output steps, correct time grid and impulse, finite data, successful
rank exits and normal SALMON completion. Check post-impulse energy variation
(excluding pre-kick t=0) below 1e-7 Ha. For the common 4×4×1 state, compare all
16 current rows and 17 energy rows across ranks: current tolerance 1e-10 a.u.,
energy tolerance 1e-8 Ha. This short-run check does not establish long-time
accuracy or scientific convergence of the grid/cutoff/impulse size.

Raw data includes source snapshots, build cache, binary and pseudo hashes,
seed payload hashes, every rank RSS, timings, currents, energies and MLWF
status lines. Preserve the output directory until it has been archived.

Recorded results: [Japanese RT measurement note](../../docs/h2-pbeh-rt-scaling-ja.md).

## Paired adaptive-support measurements

Reuse the archived GS payloads and original `results.json` in a seed directory:

```sh
python3 testsuites/benchmark_h2_pbeh_rt/run_adaptive.py \
  --binary /absolute/path/to/salmon \
  --seeds /absolute/path/to/original-measurements \
  --output /absolute/path/to/new-adaptive-measurements --repeat 3
python3 testsuites/benchmark_h2_pbeh_rt/analyze_adaptive.py \
  /absolute/path/to/new-adaptive-measurements/results.json \
  /absolute/path/to/analysis
```

Runs 72 jobs: the same 12 cases, three repeats, norm fractions 1.0 and 0.999.
Both use the spatial EXX backend even at MPI1; 1.0 is a same-backend full-support
reference. Modes alternate within each repetition, with reversed order on even
repetitions. Initial seeds are identical and checked against archived hashes.
`exx_local_fft='auto'`; all other physical conditions remain the same.

Unlike the original full-support harness, the paired driver records approximate
support energy widths even when they exceed 1e-7 Ha, and sets an explicit flag.
It retains the original parser's completion, finiteness, rank, and step checks.
The analyzer records current and energy differences across modes and layouts;
`complete=true` means all measurements finished, **not numerical validation**.
Support counters are emitted only while adaptive masking is active. Their pair
counts exclude full-support fallback updates. The FFT volume ratio describes
local pair transforms only; it excludes setup and communication costs.

Results: [99.9% support measurement](../../docs/h2-pbeh-adaptive-rt-scaling-ja.md).
