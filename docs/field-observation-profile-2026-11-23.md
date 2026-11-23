# Local field observation profile — 2026-11-23

This measures the existing field/LD observer, not an optimization of propagation. The maintained calculation libraries are unchanged from the validated Ohmic stage (runtime SHA256 ef886367d07c2df94a8984d2d5122755ec645dd489dd9b8f6a7af91609a7a33a). No source/object audit is repeated at job startup; the recorded library SHA is verified when loading. Newly built lightweight wheel artifacts retain their build-time audits.

## Measured scope

Single Mac; user reports64GB memory and18 cores. At most five MPI processes, one thread per rank: one/four field ranks plus one bookkeeping rank. Two repeats of40³/72³ grids, no/every-step/every8-step observation give24 runs, each64 updates. Same reciprocal tensor background material, source and point geometry as the previous field fixture; h=2 au and dt=.003 au stay fixed. Larger grids extend vacuum, not the deferred material-grid convergence campaign. No GS/SCF, remote jobs or manual CI dispatch.

The development tool now accepts --grid-count and --observation-stride and records setup, advance, backup capture/release, observation, loop, checkpoint and close times per rank. Imports no MPI until after argument validation. It reuses process_memory() and records missing resource measurements separately from completed calculation. Rank reports are gathered once after computation; no profiling collective is added inside the loop.

## Results at72³

Times are medians of two repeats of per-rank maxima; loop ranges show both observed values. Different phase maxima can come from different ranks, so do not sum them as a precise critical-path decomposition. Memory is the sum of individual process high-water RSS, including the bookkeeping rank. It is not a simultaneous aggregate peak, a job memory bound or whole-system memory usage.

|Field ranks|Observation stride|Loop seconds median [min,max]|Energy query seconds|Backup capture seconds|Summed peak RSS MiB|Native point bytes|
|---|---|---|---|---|---|---|
|1|off|20.65 [20.49, 20.81]|0.000|0.000|750.0|10520|
|1|1|31.62 [23.00, 40.24]|2.733|0.624|1007.0|10520|
|1|8|23.48 [16.20, 30.77]|0.361|0.080|1007.6|1448|
|4|off|13.15 [6.68, 19.61]|0.000|0.000|1540.7|10520|
|4|1|15.80 [7.72, 23.87]|1.738|0.434|2295.7|10520|
|4|8|8.34 [5.58, 11.09]|0.172|0.051|2175.1|1448|

The observation-query time decreases with fewer queries. In the single-field-rank72³ case, measured query time is2.733s versus0.361s. Point output falls from10,520 to1,448 bytes (64 versus8 records plus headers). File count remains constant for a given rank layout:8 before report.json at one field rank,11 at four ranks; final rank checkpoints account for the difference. This fixture has one output point and does not prove metadata scaling for thousands of material points.

Peak memory does not fall in proportion to observation frequency: a sampled update still needs a complete old-state backup and cross-operator scratch. For the single-field-rank case, summed peaks stay about1007MiB. Extra MPI ranks add duplicated process state and halos. These data do not certify operation under a4GiB hard limit.

Loop timings have large repeat spread. Background activity/thermal/scheduling effects were not isolated. Do not infer optimal rank count, reliable speedup or multi-node scalability. Stride changes both native point sampling and Python/energy sampling; total loop differences cannot be attributed solely to the energy query. Advance includes native point diagnostics and MPI waits. Pure network time/traffic is not measured.

The current tensor-background query invokes three field-group reduction helpers and twelve halo helpers per field rank; those counts come from source inspection, not a network trace. Backup/advance operations have their own communication and are excluded from those counts. No application-wide per-step profiler agreement is added.

## Integrity and failure checks

- Default instrumentation preserves the previous final checkpoint bytes. Within the same grid/rank layout, physical payload and all other checkpoint bytes match across repeats/observation intervals, except the two explicitly parsed sampled E/H maxima. Current GNU int/logical32, real64, little-endian tensor format is the only format interpreted by the validation script. Original files are not rewritten.
- Sparse A/E/history and energy queries equal matching dense frames exactly. Native sparse point E/H output equals dense rows8,16,...64 exactly. Sparse queries cover the individual immediately preceding update only, not all intervening updates; they cannot be summed into a full energy balance.
- Across one/four field ranks, mapped A/E and energy budgets agree at rtol2e-10, atol2e-12. This is not a full-field repartition comparison.
- Injected RSS-probe failure produces calculation_status=completed and resource_measurement_status=partial, RSS=None with reason. History and full final checkpoint bytes equal the successful-probe control. File-inventory errors are also reported diagnostically. Actual calculation failures remain failures.
- The fixture captures an old state only for queried updates. This is a development cost measurement, not removal of per-step failure recovery from protected coupled production workflows. Existing protected workflows can reuse their existing guard.
- Source and installed-wheel suites pass696 tests each, one optional skip each. The latest wheel reproduces sparse observations and final state bytes exactly. Wheel docs/AI JSON match the checkout; the sdist contains the exact fixture source. Dependencies are supplied by the validation environment, not an offline installation test.

Raw runs, per-rank metrics, checkpoint hashes, medians.csv, exact fixture and runner/analysis scripts, audited runtime and final wheel: canonical calc2026/salmonpy-field-observation-profile-20261123. Local OneDrive storage is confirmed; cloud synchronization is not checked.

## Follow-up

Keep observer interval configurable by the calling program. Always retain per-step work/flux accounting when a cumulative energy balance is required; do not substitute sparse samples. Before optimizing, separately measure real protected mixed sessions with many points and long runs, including file metadata on the target shared filesystem. Facility MPI/multi-node and other-platform runtime evidence remain open.
