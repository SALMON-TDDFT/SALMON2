# Native Gamma matrix memory/time benchmark

This compares **current-source routing variants**, not historical commits. Only
`exx_native.f90` is recompiled; all other objects and compiler/link settings are
identical. The replicated variant calls the existing replicated gauge and ACE
paths. The ACE variant is the adopted production route; the all variant additionally
enables the experimental distributed MLWF route. Production source/build objects are not edited.

Build provenance includes routing sources, compile/link commands, executable and
shared-object SHA256 hashes, working-tree diff, and new distributed module sources.

```sh
python3 benchmarks/benchmark_matrix_memory/build_variants.py \
  --build /path/to/mpi-scalapack-build --output /path/to/results/variants
python3 benchmarks/benchmark_matrix_memory/run.py \
  --work /path/to/results/mpi8-omp1 --variants /path/to/results/variants \
  --reference-root /path/to/h2-exx-memory-final --states 64 128 --repeats 1
```

The reference tree must contain the completed `adaptive-8x1x1-mpi8-r1` and/or
`adaptive-16x1x1-mpi16-r1` input files and `data_dcdft` seed directories. Seed data
are reused read-only. Runs use impulse 16 steps, source ACE, 0.999 norm support,
and the template's other physical/numerical controls. The runner defaults to
MPI 8/OMP 1 (BLAS 1), and also supports MPI 4/OMP 2. It runs one job at a time,
reversing variant order on alternate repetitions. Choose a new output directory.

Timing is SALMON's maximum-rank RT-iteration and propagation timers plus launcher
wall time. Rank memory is each SALMON child process's lifetime high-water RSS
from `wait4`; it includes LCFO reconstruction and initialization. A 0.5-second
RSS sampler adds context but cannot resolve short allocation peaks. RSS sums
include shared pages more than once and are not unique physical node memory.
Current and total-energy arrays are compared with the first replicated run.
No zero-field calculation or current subtraction is performed.

Aggregate one or more completed batches (identical binaries, inputs and thread
configuration are required):

```sh
python3 benchmarks/benchmark_matrix_memory/summarize.py \
  /path/to/results/mpi8-omp1/results.json \
  /path/to/results/mpi8-omp1-repeat/results.json \
  --output /path/to/report
```

This recomputes observable differences across all batches, rejects incomplete
batches and numerical mismatches, and writes `summary.json` and
`time-memory.png`. Bars show medians; error bars show the observed range, not
confidence intervals. A single sample cannot establish timing variability.

For a single adopted-route update, use `--modes ace --states 128 --repeats 1`
with a manifest pointing to the updated executable. Compare with a saved batch:

```sh
python3 benchmarks/benchmark_matrix_memory/compare_update.py \
  --before /path/to/baseline/results.json --after /path/to/update/results.json \
  --output /path/to/comparison.json
```

The comparison checks input hashes, seed paths, execution configuration and
current/energy agreement. It does not claim statistical significance from one run.
