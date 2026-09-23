# Si HSE06: ACE implementation and short RT validation

The existing independent HSE reference now supports adaptively compressed exchange (ACE), constructed from the full MLWF exchange action. Native SALMON hybrid support remains unchanged. Si8, 16 occupied orbitals, 4³ shifted k mesh, 12³ primitive grid, all orbital pairs, HSE06 fraction 0.25 and omega 0.11 bohr^-1. No k convergence study.

## Algorithm and correctness

Each k block uses the grid-weighted negative exchange metric to construct a negative Hermitian low-rank operator. No eigenvalue/rank truncation is applied; singular or non-Hermitian metrics are rejected. On the actual Si construction orbitals, the full action is reproduced to relative **1.95e-15**. Tests also cover complex gauge covariance, arbitrary-target Hermiticity, sign, invalid input, and interpolation for nonorthogonal construction vectors.

The RT integration is **self-consistent implicit midpoint in the original Bloch gauge**, not PT-CN. Inner iterations rebuild the Hartree/semilocal terms and reuse ACE. Outer iterations recompute untruncated exchange and accept only if the full relative midpoint-equation residual is below 1e-10 (inner tolerance 1e-12). No fixed exchange freezing interval is introduced. Subsequent steps seed the inner solve with the previous accepted midpoint's ACE, without bypassing the full-residual acceptance gate.

This isolates ACE from a change of orbital gauge. The implementation does not yet realize the larger-time-step benefit of parallel-transport gauge. Midpoint is second order; it cannot be declared equivalent to fourth-order RK4 at the same dt merely because exchange is converged. Full exchange is also used for initial/final energy checks.

## Measured application cost

Apple M5 Pro, single-thread FFTW and OPENBLAS_NUM_THREADS=1. In the first positive-impulse pilot:

| Operation | Wall time |
|---|---:|
| Full MLWF exchange action | 11.08 s |
| ACE compression after full action | 0.00785 s |
| ACE application, median of 7 warmed calls | 0.00371 s |
| Complete first midpoint step, including 3 full builds | 35.14 s |
| Earlier full RK4 step, 4 full evaluations | 46.02 s |

ACE application is about 3,000 times cheaper than one full exchange action for this case, but **this is not a 3,000-fold propagation speedup**. Building ACE dominates total step cost. The first-step cost ratio is only 1.31, and the two time integrators have different accuracy. These are short observed step timings, not repeated long-job throughput benchmarks. Tests/metadata collection lasting less than a second were performed during parts of the exploratory pilot sequence; no second simulation job ran concurrently.

## Time-step comparison

At the same final time 0.08 a.u. = 0.001935 fs, positive impulse A_z=1e-4 a.u.:

| Midpoint integration | Relative current difference from full RK4 (2×0.04) | Relative density difference |
|---|---:|---:|
| 1×0.08 | 2.93e-5 | 4.06e-8 |
| 2×0.04 | 7.55e-6 | 1.03e-8 |

The decrease is consistent with second-order convergence in this short pilot. Energy changes were 2.44e-12 Ha and 2.49e-14 Ha respectively. The first dt=0.08 full residual was 8.80e-11; electron count was 31.9999999999997 and maximum orbital Gram error 3.2e-14. Good norm preservation does not establish nonlinear energy conservation or long-time stability.

Raw result, trajectory, construction, and comparison data are collected in `validation.json`. Large states are preserved locally in `calculations/si_hse_reference/ace_*` and ignored by Git. `collect_results.py` regenerates the compact report data from those checkpoints. `samples/hse_mlwf_reference/README_ACE.md` documents execution and limitations.

## Scope

No HSE optical spectrum, exciton shift, pump–probe result, or long-time stability claim follows from these pilots. Next substantial algorithmic work is PT-gauge propagation and a matched-accuracy longer-duration benchmark. The old full RK4 implementation is retained as a reference.

References: [Lin, ACE (2016)](https://arxiv.org/abs/1601.07159); [Jia and Lin, PT+ACE real-time hybrid calculations (2019)](https://arxiv.org/abs/1809.09609). The latter motivates subsequent PT work; it is not the integrator implemented here.

## Reusing the previous midpoint operator

The seeded 3-step run (dt=0.08, A_z=1e-4) required **3, 2, 2 full exchange builds**, with complete step times **35.82, 24.28, 24.25 seconds**. Steps 2–3 each used 16 cheap inner applications and passed full residuals 7.34e-11. Endpoint electron number stayed within 3e-13 of 32 and maximum Gram error was below 9e-14. The last two step times are about 1.90 times shorter than the earlier 46.02-second RK4 step, **a same-dt cost comparison, not a matched-accuracy speedup**. Total simulated time is only 0.24 a.u. (0.00581 fs).

## Final verification

All 40 Python tests passed, including the existing HSE functional derivative, native projector/current fixture, new ACE tests and midpoint failure gates. Independent mathematical/code review found no blocker for this documented scope. The zero-impulse dt=0.08 pilot gave energy change 2.42e-12 Ha, electron number 32.0, maximum Gram error 1.05e-14, and full residual 8.00e-11. Its small residual current (~1e-12 a.u.) is recorded, not treated as a physical signal. No native Fortran source changed in the ACE extension.
