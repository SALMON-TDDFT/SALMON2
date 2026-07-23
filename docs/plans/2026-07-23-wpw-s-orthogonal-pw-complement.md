# WPW S-Orthogonal PW Complement Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Make selected PW modes an exact S-orthogonal complement of the occupied-Wannier span without changing the retained WPW space or suppressing Hamiltonian W--PW coupling.

**Architecture:** Build a distributed transform `T=[[I,-A],[0,I]]`, where `S_WW A=S_WP`, from the completed bounded metric. Wrap existing H/S actions with `T` and `T^dagger`, map seed/density coefficients consistently, and bind transform validity to a metric-only fingerprint rather than Hamiltonian epoch. First expose this only as a default-off fixed-H comparison; promote it to normal SCF/checkpoint use only after the B=6 gate.

**Tech Stack:** Fortran 2008, MPI, LAPACK test oracle, Python source contracts, CMake.

---

### Task 1: Add the distributed complement transform and dense oracle fixture

**Files:**
- Create: `src/common/dg_wpw_s_orthogonal_complement.f90`
- Modify: `src/common/CMakeLists.txt`
- Create: `tests/dg/test_dg_wpw_s_orthogonal_complement_mpi.f90`
- Create: `tests/dg/run_dg_wpw_s_orthogonal_complement_mpi.py`
- Modify: `tests/dg/check_dg_wpw_production_consumer.py`

1. Add a RED two-rank fixture with a nonidentity SPD `S_WW`, nonzero `S_WP`, and nonidentity `S_PP`. Require initialization of a public `s_dg_wpw_s_orthogonal_complement` and verify failure because the module/API does not exist.
2. In the fixture, form the dense oracle
   `A=solve(S_WW,S_WP)`, `T=[[I,-A],[0,I]]`, `S'=T^H S T`, and require relative `maxabs(S'_WP)<=1d-11`.
3. Require coefficient maps
   `to_original: (cW,cP)=(cW'-A*cP',cP')` and
   `from_original: (cW',cP')=(cW+A*cP,cP)` to round-trip within `1d-11`.
4. Require wrapped S and H actions to equal dense `T^H*S*T*x` and `T^H*H*T*x` within `1d-11`; use nonzero `H_WP` to prove physical Hamiltonian mixing is retained.
5. Require collective rejection for NaN metric input, rank-disagreeing dimensions, singular/cutoff-failing `S_WW`, and stale metric fingerprint.
6. Add a source contract that forbids global `S_WW`/`A` allocation on a single production rank and requires owner-schedule fetch/reduce plus collective validation before every distributed solve/action.
7. Run the new runner and contract; expect FAIL for the missing module.
8. Implement the module with this public API:
   - `initialize_dg_wpw_s_orthogonal_complement(op,relative_cutoff,transform,info)`
   - `apply_h_dg_wpw_s_orthogonal_complement(op,transform,xw,xp,yw,yp,info)`
   - `apply_s_dg_wpw_s_orthogonal_complement(...)`
   - `map_dg_wpw_complement_to_original(...)`
   - `map_dg_wpw_original_to_complement(...)`
   - `validate_dg_wpw_s_orthogonal_complement(op,transform,info)`
   - `release_dg_wpw_s_orthogonal_complement(transform)`
9. Store distributed `A` rows following owned-W ownership and canonical required-P IDs. Compute `S_WP` RHS from the installed real-space metric and solve `S_WW A=S_WP` collectively with block PCG/global Gram; do not use the Task-19 fragment-local inverse.
10. Build a metric-only 64-bit fingerprint from W/P IDs, ownership, metric convention, and all distributed WW/WP/PP S entries. Reduce it collectively; do not include H entries, interface lambda, or operator epoch.
11. Report solve residual, `||W^H S P_perp||`, numerical P-complement rank, cutoff-low metric weight, minimum Rayleigh estimate, condition estimate, and metric fingerprint. Do not delete P columns in this comparison. Reject initialization unless solve residual and cross-metric defect are `<=1d-11` relative.
12. Run the new two-rank fixture and source contract; expect PASS.
13. Run `git diff --check` and commit only the new module, CMake entry, fixture, runner, and contract hunk.

### Task 2: Integrate a default-off fixed-H comparison route

**Files:**
- Modify: `src/io/salmon_global.f90`
- Modify: `src/io/inputoutput.f90`
- Modify: `src/gs/dc/lcfo_flux.f90`
- Modify: `tests/dg/test_dg_wpw_matrix_free_scf.py`
- Modify: `tests/dg/check_dg_wpw_fixed_h_relaxation.py`
- Modify: `tests/dg/test_dg_wpw_production_operator_mpi.f90`

1. Add RED assertions for `yn_dg_wpw_s_orthogonal_pw`: declaration, namelist, default `n`, broadcast, variables log, and y/n validation.
2. Require comparison mode `y` to be accepted only with `yn_dg_wpw_fixed_h_relaxation='y'` in the first implementation. This prevents silently changing normal SCF/checkpoint semantics before the physical gate.
3. Add RED route assertions that default `n` uses the legacy H/S callbacks and coefficient flow byte-for-byte, while `y` installs complement-wrapped H/S callbacks for fixed-H and continuation.
4. Extend the real production-operator fixture to verify that `set_dg_wpw_interface_lambda` preserves the complement fingerprint and wrapped action, an H-only potential rebuild preserves transform validity, and a genuine S/ID change rejects it collectively.
5. Add RED coefficient-flow checks: convert the density-carrying seed from original to complement coordinates before the first solve; convert occupied Ritz coefficients back to original coordinates before density diagnosis/publication-facing reconstruction.
6. Run matrix-free, fixed-H, and production-operator tests; expect FAIL for missing control/integration.
7. Add the input plumbing and initialize the transform after the completed real-space metric/bounded operator is installed and before the density-carrying fixed-H seed is solved.
8. Add LCFO wrappers `wpw_apply_h_complement` and `wpw_apply_s_complement`. Keep correction/preconditioner selection independent: the first physical run uses neither diagonal nor metric-block correction.
9. Keep solver state in complement coordinates. At every boundary to existing physical density/projector/checkpoint-facing code, call `map_dg_wpw_complement_to_original`; never mix coordinate systems in one Gram or residual calculation.
10. Log `[DG-WPW-PW-COMPLEMENT] mode=... fingerprint=... solve_residual=... cross_metric_defect=... numerical_p_rank=...` once per fragment root.
11. Run the focused tests and real two-rank fixture; expect PASS.
12. Run the full MPI/EigenExa build, `git diff --check`, and request code review focused on collective ordering, coordinate-map boundaries, metric fingerprint lifecycle, and default-off compatibility.
13. Resolve all Critical/Important findings, rerun verification, and commit only intended source/test changes.

### Task 3: Run the B=6 physical gate

**Files:**
- Modify: `docs/plans/2026-07-23-wpw-s-orthogonal-pw-complement.md`
- Modify: `docs/plans/2026-07-19-wpw-density-carrying-seed-session-handoff.md`

1. Clone the Task-16 B=6 restart into a fresh run directory.
2. Set `yn_dg_wpw_s_orthogonal_pw='y'`, both preconditioners=`n`, and search history=`y`; keep every basis, cutoff, tolerance, continuation, and publication parameter unchanged.
3. Run 8 MPI ranks through the existing fixed-H boundary.
4. Record complement solve residual, W--P cross-metric defect, unchanged PW column count, numerical complement rank, transformed metric condition/effective rank, and coefficient round-trip defect.
5. Compare inner 32/96/160 occupied/extra residuals, Ritz boundary defects, final `info`, and publication state with Task 16 and Task 19.
6. Compare selected generalized eigenvalues with the untransformed run; require differences consistent with solver residual rather than treating nonconverged values as exact invariants.
7. Interpret improvement as evidence for coordinate/metric overlap; interpret no improvement as evidence that a span-preserving coordinate transform is insufficient.
8. Record the result and next hypothesis in this plan and the session handoff.
9. Run `git diff --check` and commit only the two result documents.

#### 2026-07-23 B=6 physical gate result

Task 16のrestartを`/tmp/20260723_task3_s_orthogonal_b6`へ複製し、8 MPI rankで
diagonal/metric preconditioner=`n`、search history=`y`、S-orthogonal complement=`y`
として実行した。basis、cutoff、tolerance、continuation、publication設定は変更していない。
DC-SCFは87 iteration、最終差`9.8498230E-10`で既存seedを再現した。

補空間初期化のsolve residualとW--P cross-metric defectはいずれも
`6.4725E-13`、metric fingerprintは`-8374564770650107427`だった。raw PW列数は
8 fragment x 128=`1024`のまま、cutoff `1.0E-10`でnumerical complement rankは`976`
（低metric-weight 48方向）である。production runに不足した診断値だけを、solver stateを
変更しない一時logging buildで再実行したところ、original→complement→originalの最大係数
round-trip defectは`4.2352E-22`だった。

inner 32/96/160のoccupied/extra residualは
`1.8708E-03/3.4680E-03`、`2.9849E-04/8.7240E-04`、
`2.0653E-04/1.5936E-04`で、search metric effective rankは`426/284/202`だった。
retained最小固有値比は`1.0433E-10/1.0063E-10/1.1730E-10`
（condition proxy約`9.59E9/9.94E9/8.53E9`）。Ritz post/direct residual-vectorの
occupied/extra相対defectは順に`2.96E-12/1.26E-12`、
`1.35E-11/2.69E-11`、`4.94E-11/1.29E-10`、metric orthogonality defectは
`6.80E-12/1.29E-12/3.50E-13`だった。

Task 16のresidual `1.9667E-03/4.0432E-03`、
`3.6079E-04/1.1093E-03`、`2.5607E-04/2.0289E-04`に対する改善は約
1.05--1.27倍に留まった。Task 19 metric-block correctionに対してoccupiedは約
24--124倍悪い一方、extraは約1.45--6.02倍良い。全runともfixed-Hは`info=40`で、
checkpoint/publicationは行われなかった。

同じ一時logging buildによるuntransformed比較では、selected generalized eigenvalue
`(e1,e_occ,e_extra1,e_last)`の絶対差はそれぞれ
`1.13E-08/8.33E-08/2.66E-08/1.61E-06`だった。これは非収束runの
`O(1E-4)` residualより十分小さく、座標変換による固有値の有意な変化は認めないが、
非収束値を厳密な不変量とは扱わない。

従ってS-orthogonal化はmetric overlapを正しく除去し、同一retained spanを保つものの、
plateauを質的には改善しない。座標/metric overlapだけが主因ではなく、span-preservingな
座標変換だけでは不十分である。default-off診断routeは保持し、normal outer SCF、
checkpoint、RTへは昇格させない。次の仮説はextra stateの過増幅を抑えた
state-partitioned/regularized correction、または`H-epsilon S`に基づく
metric-consistent correctionである。

### Task 4: Promotion decision checkpoint

Do not integrate complement coordinates into normal outer SCF, checkpoint schema, or RT handoff automatically.

1. Request user review of the B=6 gate.
2. If accepted, write a separate design/plan for normal SCF and checkpoint/RT promotion using the metric fingerprint and explicit coordinate provenance.
3. If rejected, leave the default-off diagnostic route available and proceed to `H-epsilon S` or state-partitioned correction design.
