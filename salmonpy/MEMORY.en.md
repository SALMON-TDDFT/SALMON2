# Memory estimates and advisory budgets

Budgets generate warnings only. Excess does not stop computation/checkpointing,
set OS limits or change scheduler settings. A below-budget estimate is not a
feasibility guarantee. One GiB is2**30 bytes.

## Python entry

```python
from salmonpy.memory import MemoryEstimate, process_memory, memory_advice

live, basis = process_memory()
# Replace this illustrative512MiB with your backend's EXTRA candidate storage.
estimate = MemoryEstimate(live_bytes=live, candidate_bytes=512 * 2**20,
                          reserve_bytes=128 * 2**20, basis=basis)
report = memory_advice(estimate, budget_bytes=4 * 2**30)
print(report.to_dict())
```

user-written calculation implementations supply measured/array-derived byte counts and a basis description.
Use None for unavailable values and unknown for omitted components, never zero.
known_peak_bytes sums known estimates; estimated_peak_bytes is None if any part
is unknown. exceeds_budget=False means only that the supplied estimate fits,
not that execution is safe. True also reports excess_bytes; unknown verdicts
are None.

With MPI, pass an existing comm and call on every rank. The budget covers the
SUM across that communicator, not each rank. Supply node-local communicators
for per-node budgets. Different rank budgets are reported as unknown. Only
rank0 warns. If calculation implementation estimation fails, handle it locally and pass estimate=None
or an unknown estimate; never skip the collective on that rank.

process_memory lazily reads process peak RSS on Linux/macOS, not current use or
child memory. Shared pages may be counted per process. Probe failure/unsupported
platform returns None with a reason. Windows users can supply estimates, but no
RSS probe is inferred. Import performs no measurement/MPI initialization/build.
Warnings promoted to exceptions are caught in warning_error, not scientific
failure. Invalid API arguments still raise InputError.

## Ordinary mixed sample

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 mpiexec -n 4 \
  python -m salmonpy.examples.mpi_mixed_native BUILD_JSON RUNTIME GS_BUNDLE NEW_OUTPUT \
  --memory-budget-gib 4 --rank-budget 4
```

Four ranks with one thread each use four cores. When configured,
NEW_OUTPUT/memory-advice.json records:

- Before initialization: grid/LD-pole, SBE matrix and wavefunction-input heuristics,
  with unknown initialization workspaces identified.
- After initialization and at segment end: rank process peak RSS plus128MiB reserve per rank.
- Before protected restore: peak RSS plus4 times solver payload sizes,32 times
  JSON sizes, extra snapshot buffers and reserve.

This includes simultaneous original/candidate/read/copy storage. The multipliers
are conservative heuristics, not proved heap bounds; already-freed allocations
and shared pages can be counted again. Child worker memory remains unknown.
Inspect rank breakdowns, basis and unknowns. Budget is runtime advice and may
change at resume without changing physical inputs. Measurement, warning delivery
and advice-file write failures are diagnostic only and do not set scientific
failed flags; distinguish them from real invalid input/numerical/MPI failures.
The sample disables advice when no budget is specified. Automatic estimates here
cover the ordinary CPU mixed sample, not every DFT/ASE/external calculation implementation.

## Expanded grid and parallel checks on a larger Mac

Keep the4GiB/4-core beginner recipe as its own verified configuration. Ask for actual available resources and reserve headroom for the OS and other work. For the supplied64GB/18-core Mac, this campaign uses a48GiB advisory budget and at most15 MPI ranks with one thread each. These are campaign settings, not universal defaults or hard limits.

`mpi_mixed_native --field-shape 96 96 96` changes Maxwell grid counts: each axis must be even and >=40 and match its MPI decomposition. Default40³ and2-bohr spacing are unchanged. Enlarging the grid adds outer vacuum, so compare decompositions only at the same grid/input. Default material-point physical coordinates are preserved. Electronic GS grid/k sampling is unchanged: this is not a96³ TDDFT test. The3D fixture tests interfaces;1D remains recommended for production TDDFT coupling.

Field2+SBE1+TDDFT2 uses `mpiexec -n 5 ... --field-grid 1 2 1 --tddft-ranks 2`. Field12+SBE1+TDDFT2 uses `mpiexec -n 15 ... --field-grid 3 2 2 --field-shape 96 96 96 --tddft-ranks 2`. Replace `...` with the existing Python module invocation and BUILD,RUNTIME,GS,new-output arguments. Set `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1`.

Source-distribution tool `tools/test_mixed_scaling.py BUILD RUNTIME GS NEW_OUTPUT --rank-budget 18 --memory-budget-gib 48` uses the installed wheel and runs4/5/7/11/15 ranks sequentially twice on identical96³ inputs. A/E/j histories are compared with the first case. Different grid shapes are rejected. Supply MPI-vendor-specific flags with `--launcher-arg`; avoid other heavy work during timing.

`coupled_loop_seconds_max` excludes initialization, output and serial reference replay. Communication metrics include waits, not pure link bandwidth; FDTD called directly from Python halo exchanges are not separately decomposed. Memory is the sum of per-rank peak RSS values, not simultaneous aggregate RSS or whole-system usage. Resource/timer measurement failures are distinct from calculation failures. These short repeats do not establish universal scalability or physical accuracy.

## Explicit field observation cost measurements (development)

Source-checkout tool `tools/test_diagonal_field_energy_mpi.py` varies field size/ranks and observation interval without GS data. This command starts an actual local calculation; it is not an input-only management-node check. The user must select compute nodes and the site's launcher. Use a new results directory outside the code and a field library matching the installed wheel.

```sh
mpiexec -n 5 python -I tools/test_diagonal_field_energy_mpi.py \
  BUILD.json RESULTS/new-profile --profile tensor-background \
  --grid-count 72 --field-ranks 4 --observe --observation-stride 8
```

This uses four field ranks plus one bookkeeping rank. Grid spacing/dt stay fixed; larger grid counts extend vacuum, not a grid-spacing convergence test. The stride applies to native point output, Python history and individual-step energy queries; the final step is always retained in history. Sparse energy observations do not constitute a cumulative energy balance.

`report.json` separates `calculation_status` from `resource_measurement_status`. It records per-rank setup, advance, backup capture/release, observation and final-checkpoint wall times, plus peak RSS from `process_memory()`. Missing RSS stays None with its reason and does not fail the calculation. `--fail-resource-probe` explicitly tests that policy. Summed process high-water RSS is not a simultaneous job peak. File inventory includes checkpoints/static inputs, with point-stream bytes separate; the report itself and external launcher log are excluded.

Timings include MPI waits; pure communication time/traffic is not isolated. Observation communication stays within the field group. Measurement records are gathered application-wide once after computation. This fixture captures old state only for observed updates; it is not an optimization that removes per-step failure recovery from protected coupled calculations. Stored observation maxima depend on sampling interval and must be distinguished from physical state when comparing checkpoint files.


## Measuring many-point protected mixed calculations

The manual measurement driver `tools/measure_mixed_performance.py` launches four
separate jobs, sequentially: checked/performance at 32 points × 512 steps and
64 points × 2048 steps. It requires exactly nine local ranks: four for a 64³ E/H
grid containing two LD materials, two two-rank TDDFT groups, and one SBE rank.
One quarter of the macro points use independent TDDFT states from shared models;
the others use independent SBE states. All numerical updates use Fortran kernels.
Every update explicitly retains the state needed for protected restoration in
both modes. This does not measure unprotected performance-mode updates.

```sh
python tools/measure_mixed_performance.py BUILD.json EXISTING_SYNTHETIC_GS NEW_RESULTS
```

Prepare the package, LD/SBE kernels and embedded MPI library beforehand. The
command neither builds nor downloads them, runs no SCF, and uses test-only
read-only GS links. The existing GS must contain `rt.inp` and `data_for_restart`
with the synthetic periodic PZ fixture used by this test. This is a developer
reproduction tool, not an arbitrary-material input generator or scheduler script.
Run it on allocated compute resources; it starts MPI jobs immediately. Do not run
it on a site's login/management node. Site launchers and multi-node layouts need
separate tests; this driver assumes `mpiexec -n 9`.

Each job writes `report.json` and complete macro-point A/E/current arrays in
`traces.npz`. The driver compares every stored step and both LD-region E/H point
outputs. Observation stride is 1 for the smaller case and 16 for the larger case;
macro-point histories are always retained at every step. Changing stride does
not change the physical update. This pair is not an isolated stride benchmark.

Report interpretation:

- `loop_seconds_max` is the slowest rank's elapsed update-loop time, including
  waiting, protected backups and Python history collection. Setup and total
  launcher time are separate. Communication instrumentation itself has overhead.
- MPI counts/times cover Python orchestration calls only; Fortran-internal MPI is
  excluded. Blocking times include waiting and cannot be summed across ranks as
  network time. Logical A/J bytes omit metadata, protection, repetition and wire
  protocol overhead. The final measurement gather is outside the timed loop.
- RSS values are per-process high-water marks. Their sum is not simultaneous
  job memory, and can count shared pages repeatedly. It is not a hard memory cap.
- I/O reports regular-file counts and logical sizes, with observation bytes and
  launcher-log bytes separate. GS symlink targets are excluded. These numbers
  are neither physical disk traffic nor filesystem throughput. OS block counters
  are reported in platform block units without guessing a byte conversion.
- Resource probe failures are recorded separately from calculation failures;
  unavailable memory remains null. The test-only environment variable
  `SALMONPY_TEST_RESOURCE_FAILURE=1` exercises missing resource probes without
  changing numerical work. An actual calculation error still fails the run.

One measured run per mode/size is evidence of reproducibility in that scope,
not a statistically established speedup or physical-convergence result. The
small synthetic TDDFT/SBE material models are not production material inputs.
See the active state-isolation checklist for measured numbers and remaining
multi-node/platform work.
