# Cost versus accuracy: native HSE propagation

Si8, 12³ grid, 4³ k points, Apple M5 Pro, MPI8, one OpenMP/BLAS thread
per rank. Re-measured each case for five steps twice, reversing case order
on the second pass. Runs execute sequentially; the separate Python MPI8
production reference remains active. These are shared-machine short pilots,
not isolated-node or long-trajectory benchmarks. Timing is the maximum-rank
RT iteration timer divided by five: includes RT output and final checkpoint,
excludes initialization. Prior paused/oversubscribed batch times are not used.

| Propagator | dt (au) | s/step | min per simulated fs (estimate) | 20 fs (estimated hours) | Short current error | Short induced-density error |
|---|---:|---:|---:|---:|---:|---:|
| Taylor4 + ACE | 0.16 | 1.53 | 6.59 | 2.20 | 0.01057% | 0.04232% |
| Taylor4 + ACE | 0.08 | 1.51 | 13.00 | 4.33 | 0.00215% | 0.00911% |
| PT-CN + ACE | 0.32 | 3.19 | 6.86 | 2.29 | 0.25421% | 1.72768% |
| PT-CN + ACE | 0.16 | 2.74 | 11.79 | 3.93 | 0.06530% | 0.62835% |
| PT-CN + ACE | 0.08 | 1.75 | 15.05 | 5.02 | 0.01656% | 0.17050% |

Errors are from the independently completed 0.1935-fs convergence runs,
relative to Taylor dt=.04; they are not errors measured at 20 fs. Time
projections assume constant per-step cost, and exclude initial GS convergence.
At 0.774 fs the current errors versus Taylor dt=.08 are 0.473% (PT dt=.32)
and 0.149% (PT dt=.16); density errors are 4.31% and 2.93%. No resolved
exciton-spectrum accuracy or long-time optimal step has been established.

## Decision from these measurements

- Taylor dt=.16 is about 2.28 times cheaper per physical time than PT dt=.08,
  while giving smaller current and induced-density errors in the short test.
- PT dt=.32 halves the step count relative to Taylor dt=.16 but costs about
  2.08 times as much per step. Net cost is approximately equal (PT 4% higher,
  within practical shared-machine variability), with worse measured accuracy.
- Taylor dt=.08 is about 9% more expensive per physical time than PT dt=.16,
  but much more accurate in the short test. The small timing difference is
  not a robust performance win for either method under this background load.
- PT-CN has excellent norm/energy behavior and permits dt=.32 where the
  tested Taylor+ACE path diverges. Stability alone does not establish accuracy
  or cost-effectiveness. Retain PT-CN as an independent reference/option;
  Taylor+ACE is the current favorable baseline for this Si case.

## ACE, native implementation, and remaining cost

Full-Fock Taylor dt=.16 takes 7.41 s/step versus 1.53 s/step with ACE: **4.84x** speedup at the same dt and MPI8. The separate 0.1935-fs full/ACE audit gives 0.007109% current RMS difference and 0.003575% induced-density difference. Full-Fock timing includes the audit path's ACE refresh overhead; this is a measured implementation comparison, not an optimized full-Fock lower bound.

The earlier matched three-step PT-CN benchmark found native Fortran MPI8
4.68x faster than Python MPI8, and native MPI8 4.95x faster than MPI1.
These factors concern different comparisons and must not be multiplied.
See ../si-hse-native/README.md for definitions and raw evidence.

In that native PT-CN pilot, full exchange dominated (~1.75–1.88 s/step);
ACE construction itself was ~0.012 s and ACE applications ~0.14 s. Thus
reducing full-exchange refresh cost is more promising than tuning ACE builds.
The existing PT-CN memory snapshot is ~2.31 GiB summed rank RSS (not unique
or peak memory); no matched Taylor memory measurement was made. Sources
are still replicated and no real-space support cutoff is used, so these
results establish neither a Wannier-locality speedup nor linear k scaling.

Raw repeated timing records: costs.json. Accuracy details: README.md.
