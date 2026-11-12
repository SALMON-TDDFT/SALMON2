# Validated optimization candidates

These alternative snapshots each define wannier_neighbor_owner. Compile one
snapshot in place of the parent module, using the other parent modules.

The parent directory remains the adopted descriptor/task-cache snapshot.
These candidates preserve subsequent optimization work for review. Adoption and
native SALMON integration require a separate decision.

| Snapshot | Change | Fugaku comparison job | Baseline / candidate seconds |
| --- | --- | --- | --- |
| preparation.f90 | Outer-loop preparation OpenMP | 53989239 | 3.77476541 / 3.29398716 |
| source_plan.f90 | Input communication plan cache | 53989489 | 3.77767210 / 2.39104698 |
| input_pipeline.f90 | Input communication/computation overlap | 53990460 | 2.82208035 / 2.34051710 |

Each job compared eight Si512 tiles with 1024 MPI ranks, one rank per node and
48 OpenMP threads. Both paths returned status zero and bitwise identical actions,
with 4198400 evaluated pairs, 448862256000 source bytes and 8581996544 return bytes.
Timings are a single fixed-order comparison. Baselines differ between stages,
so speedups must not be combined. These timings cover the sampled tiles.
Local bounds-checked comparisons include empty ownership, larger tiles and
OpenMP-disabled builds. Candidate-specific tests also cover multiple windows,
heterogeneous teams and disabled prefetch where recorded.

SHA256 values pin the full-tile source snapshots copied without edits. The
Fugaku test harness additionally limits execution to the first eight tiles.
Diagnostic-only locality instrumentation is maintained separately.
