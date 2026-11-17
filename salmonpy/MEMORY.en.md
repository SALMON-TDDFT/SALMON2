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

Custom backends supply measured/array-derived byte counts and a basis description.
Use None for unavailable values and unknown for omitted components, never zero.
known_peak_bytes sums known estimates; estimated_peak_bytes is None if any part
is unknown. exceeds_budget=False means only that the supplied estimate fits,
not that execution is safe. True also reports excess_bytes; unknown verdicts
are None.

With MPI, pass an existing comm and call on every rank. The budget covers the
SUM across that communicator, not each rank. Supply node-local communicators
for per-node budgets. Different rank budgets are reported as unknown. Only
rank0 warns. If backend estimation fails, handle it locally and pass estimate=None
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
- Before protected restore: peak RSS plus4 times native payload sizes,32 times
  JSON sizes, extra snapshot buffers and reserve.

This includes simultaneous original/candidate/read/copy storage. The multipliers
are conservative heuristics, not proved heap bounds; already-freed allocations
and shared pages can be counted again. Child worker memory remains unknown.
Inspect rank breakdowns, basis and unknowns. Budget is runtime advice and may
change at resume without changing physical inputs. Measurement, warning delivery
and advice-file write failures are diagnostic only and do not set scientific
failed flags; distinguish them from real invalid input/numerical/MPI failures.
The sample disables advice when no budget is specified. Automatic estimates here
cover the ordinary CPU mixed sample, not every DFT/ASE/external backend.
