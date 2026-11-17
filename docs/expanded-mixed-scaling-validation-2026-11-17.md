# Expanded mixed MPI/grid validation — 2026-11-17

Host capacity supplied by the user:64GB RAM,18 cores. Campaign:48GiB advisory budget, OMP/BLAS1, sequential jobs, at most15 MPI ranks. No OS hard cap installed. The4GiB/4-core beginner profile is unchanged. Existing audited binaries were reused; no native numerical source was changed.

## Scope and measurements

Maxwell grid96³ =884,736 cells,13.824 times the original40³. Spacing2 bohr and default material-point physical positions retained; larger outer vacuum changes the physical domain. Comparisons are only between identical96³ inputs with different MPI decompositions. TDDFT uses the same small16³,2³-k Libxc-PBE Si fixture, two orbital ranks. This is not a96³ electronic calculation or a physical/grid-convergence study.

Each case contains two local LD materials, one SBE point and one TDDFT point.128 steps; two repetitions per layout. All ten runs pass serial-response checks. Every A/E/j array is exactly equal across all layouts/repetitions; maximum recorded difference0.

| FDTD ranks | Total ranks | Median loop seconds | Range seconds (N=2) | Speedup | Maximum summed rank peaks GiB |
|---:|---:|---:|---:|---:|---:|
| 1 | 4 | 30.36 | 28.37–32.36 | 1.00× | 1.017 |
| 2 | 5 | 17.39 | 17.34–17.44 | 1.75× | 1.089 |
| 4 | 7 | 12.40 | 12.14–12.66 | 2.45× | 1.284 |
| 8 | 11 | 7.96 | 7.51–8.42 | 3.81× | 1.532 |
| 12 | 15 | 7.42 | 7.30–7.53 | 4.09× | 1.822 |

Loop time excludes initialization, array output and serial reference replay. End-to-end times and all individual rank/communication measurements are retained in scaling.json. No affinity, frequency or background-load controls were imposed. These are regression/smoke measurements, not publication-quality strong-scaling evidence. Two samples do not establish an optimal rank allocation. Gains beyond8 field ranks are small here; fixed electronic work and native field communication were not separately profiled. Response-interface communication includes waits, and does not separately measure native FDTD halo exchanges.

RSS is the sum of independently observed per-rank peaks (before the later serial reference child), not simultaneous aggregate RSS or system-wide memory. It must not be confused with heuristic initial-plan memory allowances.

## Protected restart and regression

- Newly tested40³ field2+SBE1+TDDFT2, total5 ranks: late staging failure preserves every original native byte and Python state; original remains usable; retry and following3 frames exactly match uninterrupted execution; zero replayed steps.
- Same5-rank topology with native fast-PR HSE PT-CN k2 (40³ field), SBE and two LD materials: identical acceptance passes. HSE error recovery beyond previously supported restore candidates is not expanded.
-96³ field12(3×2×2)+SBE1+TDDFT2, total15 ranks: same protected-failure, independent continuation and exact retry pass, including ten pre-load identity rejection/agreement cases. Step64 saved, failure injected at66, future3 frames checked.
-15-rank failure-injection segment ends with summed rank peak RSS about5.071GiB; normal runs reach at most1.822GiB. Restore advisory known allowance reaches about9.505GiB, with uncertainty retained; the initial heuristic is larger. These are not universal upper bounds.
- Source and installed wheel **524 tests passed, one skipped**, sdist→wheel build/install pass, synchronized JP/EN/AI guides and sample/controller included. Input tests reject invalid budgets and incompatible shape before launching MPI; diagnostic clock failure returns unknown. Read-only review found no blocker. No Actions run.

## Reproduce

Source tool: `tools/test_mixed_scaling.py BUILD RUNTIME GS NEW_OUTPUT --rank-budget 18 --memory-budget-gib 48`. It executes the installed wheel, never compiles/import-builds, runs layouts sequentially, and compares every case with the first. MPI-specific launcher flags are explicit. Example OpenMPI timeout: `--launcher-arg=--timeout --launcher-arg=900`.

Persistent artifacts: `/Users/otobetoshihito/.codex/outputs/salmonpy-expanded-mixed-20261117/`. Includes all normal runs/logs, full scaling/summary records,5-rank and15-rank protected restore evidence, source/installed test logs and distribution hashes. Local paths are provenance, not user installation requirements.

## Remaining

General native post-update recovery, unrestricted GS/TDDFT ownership/handle ABI, many material points/larger electronic grids, archive migration/topology redistribution, and Linux/Intel Mac/runtime delivery remain open. Molecular work remains paused. Physical convergence is user-owned.
