# WPW Metric-aware Correction Implementation Plan

> **Historical/removed:** This document describes an obsolete experimental DG route
> removed on 2026-07-31. It is retained only as an implementation record and is
> not executable guidance.

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Compare the current unpreconditioned LOBPCG residual correction with a fragment-local `S_block^{-1}r` correction.

**Architecture:** Add a default-off, diagonal-mutually-exclusive input control and route the existing validated metric block-Jacobi inverse through a collective production adapter for fixed-H and continuation only. Preserve the normal production route's legacy no-callback behavior and every solver/publication decision outside correction-vector construction.

**Tech Stack:** Fortran, MPI, Python source contracts, CMake.

---

### Task 1: Add the metric correction comparison route

**Files:**
- Modify: `tests/dg/test_dg_wpw_matrix_free_scf.py`
- Modify: `tests/dg/check_dg_wpw_production_consumer.py`
- Modify: `tests/dg/test_dg_wpw_production_operator_mpi.f90`
- Modify: `src/common/dg_wpw_bounded_operator.f90`
- Modify: `src/io/salmon_global.f90`
- Modify: `src/io/inputoutput.f90`
- Modify: `src/gs/dc/lcfo_flux.f90`

1. Add RED source assertions for declaration, namelist membership, default `n`,
   broadcast, variables log, y/n validation, and rejection when both diagonal
   and metric preconditioners are `y`.
2. Add RED assertions that fixed-H and continuation select exactly one of
   diagonal, metric block, or no callback, while the normal production algebra
   route remains callback-free under default controls. Preserve retained-search
   propagation in every branch.
3. Add a RED production-adapter contract requiring collective metadata validation
   on the bounded-operator communicator before delegation on every rank.
4. Extend the production-operator two-rank MPI fixture to invoke the real adapter
   on its non-identity SPD block, verify the exact local block solve, verify
   continued validity after `set_dg_wpw_interface_lambda`, and verify rejection
   after a genuine operator epoch/layout fingerprint change.
5. Run the matrix-free and production-operator tests; expect failure because
   the control and route do not exist.
6. Add `yn_dg_wpw_metric_preconditioner='n'` to global/input plumbing and enforce
   mutual exclusion with `yn_dg_wpw_preconditioner` during argument validation.
7. Add a public common adapter with the eigensolver callback signature.  Reduce
   eigenvalue/RHS validity collectively on `op%w_schedule%comm`, require finite
   eigenvalues and matching RHS count, then make every rank delegate to
   `apply_dg_wpw_fragment_block_preconditioner`; eigenvalue values do not enter
   the block inverse.
8. Add the LCFO wrapper that calls this adapter.  Branch fixed-H and continuation
   through metric / diagonal / none.  Leave `wpw_algebra_step` callback-free and
   retain its search-history keyword.
9. Log the effective correction mode and test all configurations: default
   diagonal=y/metric=n, none n/n, metric n/y, and invalid y/y.
10. Run the matrix-free and production-operator tests; expect PASS.
11. Run `git diff --check` and commit only intended implementation/test files.

### Task 2: Verify, review, and compare B=6

**Files:**
- Modify: `docs/plans/2026-07-23-wpw-metric-correction.md`
- Modify: `docs/plans/2026-07-19-wpw-density-carrying-seed-session-handoff.md`

1. Run matrix-free, fixed-H, generalized-algebra, occupied-W MPI tests, full
   MPI/EigenExa build, and `git diff --check`.
2. Request code review for mutual exclusion, callback lifecycle/fingerprint,
   MPI collective ordering, route-specific legacy behavior, and default
   compatibility.
3. Resolve every Critical/Important finding and rerun verification.
4. Commit only intended source/test changes.
5. Create a fresh Task-16-equivalent B=6 run with
   `yn_dg_wpw_preconditioner='n'`,
   `yn_dg_wpw_metric_preconditioner='y'`, and retained search history.
6. Run 8 MPI ranks through the unchanged fixed-H boundary.
7. Record the fragment block dimension/condition and compare inner 32/96/160
   occupied/extra residuals, correction ratios, effective rank, Ritz boundary
   defects, final `info`, and publication state with Task 16. Interpret failure
   only as failure of the local block-Jacobi metric approximation.
8. Record the evidence and next hypothesis in this plan and the session handoff.
9. Run `git diff --check` and commit only the two result documents.

## 2026-07-23 実施結果

Task 1は`803f6df`で実装した。`yn_dg_wpw_metric_preconditioner`はdefault `n`で、
既存diagonal controlとの同時`y`を入力検証で拒否する。fixed-H/continuationだけを
metric block / diagonal / noneの3経路に分け、通常のproduction algebraはcallback-freeの
まま維持した。共通adapterは固有値数・有限性とRHSをbounded-operator communicator上で
collectiveに検証してから、既存のfragment-local block inverseへ全rankで委譲する。

RED確認後、matrix-free 2-rank、real production-operator 2-rank、fixed-H source contract、
generalized algebra、occupied-W basis/row-layout MPI、MPI/EigenExa全体buildがPASSした。
real fixtureは非identity SPD blockのexact solve、interface lambda変更後の因子有効性、
operator epoch/layout fingerprint変更後のstale因子拒否を確認した。コードレビューは
Critical/Importantなし、`git diff --check`もPASSした。

fresh B=6比較は
`stage2d_wpw_runs/20260723_task19_metric_block_b6/run.log`で8 rank実行した。
Task 16と同じnormalized seed、diagonal=`n`、metric block=`y`、search history=`y`である。
各fragmentのblock dimensionは138、固有値範囲は`4.1210E-06`--`1.8573E+00`、
condition numberは`4.5069E+05`だった。

| inner | Task 16 occupied / extra | metric block occupied / extra | effective rank Task 16 / metric |
|---:|---:|---:|---:|
| 32 | `1.9667E-03 / 4.0432E-03` | `7.6809E-05 / 6.2909E-03` | `437 / 337` |
| 96 | `3.6079E-04 / 1.1093E-03` | `2.4120E-06 / 1.2688E-03` | `292 / 210` |
| 160 | `2.5607E-04 / 2.0289E-04` | `1.9816E-06 / 9.5958E-04` | `213 / 175` |

occupied residualはinner 32/96/160で約25.6/149.6/129.2倍改善した。しかしextra residualは
約1.56/1.14/4.73倍悪化し、最終max residualもTask 16の`2.5607E-04`から
`9.5958E-04`へ約3.75倍悪化した。correction norm ratioはinner 32で
occupied/extra=`7.7504E+02/5.5529E+02`、inner 96で
`3.6288E+02/7.0644E+02`、inner 160で`1.5280E+02/8.4597E+02`であり、
conditionの大きい局所逆が特にextra correctionを強く増幅した。

31→32、95→96、159→160のRitz post/direct相対defectは最大でもoccupied
`5.96E-10`、extra`1.81E-09`、post metric orthogonalityは`1.97E-12`から
`1.08E-13`で、更新整合性は維持された。runはfixed-H `info=40`、publicationなしで
正常終了した。

結論は「metric-aware correction一般が無効」ではなく、今回のfragment-local
`S_block^{-1}r` block-Jacobi近似がoccupiedを大幅改善する一方、extra statesを悪化させて
全体収束を改善しなかった、である。基底表現不足よりcorrection生成手続きが主対象という
従来の判断は強まった。次は局所逆をそのまま全stateへ適用せず、extra側の過増幅を抑える
regularization/state partition、または`H-epsilon S`を近似するmetric-consistent correctionを
同じ基底・publication gateで比較する。
