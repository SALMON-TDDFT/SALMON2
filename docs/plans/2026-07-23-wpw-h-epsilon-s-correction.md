# WPW H-epsilon S Correction Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add a default-off fixed-H comparison that generates WPW residual corrections from a regularized fragment-local approximation to `(H-epsilon*S)^{-1}`.

**Architecture:** Factor each fragment's bounded generalized H/S block once, then apply a signed-floor spectral inverse separately for every Ritz value. Route the correction through the existing eigensolver callback, project it out of the current Ritz span in the global S metric, and preserve all normal-production, search-history, convergence, and publication behavior.

**Tech Stack:** Fortran 2008, MPI, LAPACK, Python source contracts, CMake, ScaLAPACK/EigenExa integration build.

---

### Task 1: Add the fragment-local generalized spectral correction

**Files:**
- Modify: `src/common/dg_wpw_bounded_operator.f90`
- Modify: `tests/dg/test_dg_wpw_production_operator_mpi.f90`
- Modify: `tests/dg/check_dg_wpw_production_consumer.py`

1. Add a RED two-rank fixture using a nonidentity SPD S block and a noncommuting Hermitian H block. Call a new public factor initializer and correction application API; require compilation failure because they do not exist.
2. In the fixture, solve the dense oracle with LAPACK, normalize `U^dagger*S*U=I`, and form
   `z=-U*g(theta-epsilon)*U^dagger*r`. Include Ritz values below, inside, and above the local spectrum, plus one denominator below the relative floor.
3. Require relative max-absolute agreement with the dense oracle within `1d-11`, finite bounded corrections, deterministic positive flooring at exact zero, and the same rule for columns labelled occupied and extra.
4. Require diagnostics for raw/floored minimum denominator, floored mode count, floored residual weight, and maximum inverse magnitude.
5. Require collective rejection of nonfinite H/S/Ritz/RHS input, non-SPD S, invalid regularization, rank-disagreeing shapes, generalized eigensolver failure, and stale operator epoch/layout fingerprint.
6. Extend the source contract to require collective validation before each factorization/application and prohibit silent fallback or single-rank global production allocation.
7. Run:
   `python3 tests/dg/run_dg_wpw_production_operator_mpi.py`
   and
   `python3 tests/dg/check_dg_wpw_production_consumer.py`.
   Expect failure for the missing API.
8. Add `s_dg_wpw_h_epsilon_s_factor` storing epoch, layout fingerprint, W/P dimensions, relative floor, S conditioning, generalized eigenvalues/eigenvectors, and validity.
9. Add:
   `initialize_dg_wpw_h_epsilon_s_factor(op,relative_floor,factor,info)`,
   `apply_dg_wpw_h_epsilon_s_correction(op,factor,eigenvalues,rw,rp,zw,zp,diagnostics,info)`,
   and `release_dg_wpw_h_epsilon_s_factor(factor)`.
10. Assemble local owned H/S blocks in the same canonical W+P ordering used by the existing metric factor. Use LAPACK `zhegv` or an equivalent Cholesky-reduced Hermitian solve; reject cutoff-failing S instead of deleting modes.
11. For each state use `scale=max(1d0,maxval(abs(theta)),abs(epsilon_i))`,
    `floor=relative_floor*scale`, and
    `denominator=sign(max(abs(theta-epsilon_i),floor),theta-epsilon_i)`,
    with exact zero mapped to positive `floor`. Apply the leading minus sign required by the correction equation.
12. Compute diagnostics before overwriting spectral residuals and validate the final correction collectively.
13. Re-run the runner and source contract; expect PASS.
14. Run `git diff --check`.
15. Commit only the common module and its fixture/contract hunks:
    `git commit -m "feat(dg-wpw): add H-epsilon-S correction factor"`.

### Task 2: Integrate the default-off fixed-H route and global S projection

**Files:**
- Modify: `src/io/salmon_global.f90`
- Modify: `src/io/inputoutput.f90`
- Modify: `src/common/dg_wpw_production_context.f90`
- Modify: `src/gs/dc/dg_wpw_matrix_free_scf.f90`
- Modify: `src/gs/dc/lcfo_flux.f90`
- Modify: `tests/dg/test_dg_wpw_matrix_free_scf.py`
- Modify: `tests/dg/test_dg_wpw_matrix_free_scf_mpi.f90`
- Modify: `tests/dg/check_dg_wpw_fixed_h_relaxation.py`
- Modify: `tests/dg/test_dg_wpw_production_operator_mpi.f90`

1. Add RED source assertions for a default-`n` control named
   `yn_dg_wpw_h_epsilon_s_correction`: declaration, namelist, broadcast,
   variables log, y/n validation, and mutual exclusion with both existing
   preconditioner controls.
2. Add RED route assertions that the new callback is used only by fixed-H and
   continuation, that normal production algebra stays callback-free, that
   search-history propagation is unchanged, and that S-orthogonal comparison
   remains an independent control.
3. Add a RED matrix-free MPI fixture for
   `z <- z-Q*(Q^dagger*S*z)`. Use a nonidentity S action and require
   `maxabs(Q^dagger*S*z)<=1d-11`, while preserving a correction already
   S-orthogonal to Q.
4. Add RED lifecycle checks: initialize the factor after final fixed-H H/S
   installation; preserve it during read-only solves; reject it after H/S
   epoch or layout changes; release it with the production context.
5. Run:
   `python3 tests/dg/test_dg_wpw_matrix_free_scf.py`,
   `python3 tests/dg/check_dg_wpw_fixed_h_relaxation.py`,
   `python3 tests/dg/run_dg_wpw_matrix_free_scf_mpi.py`, and
   `python3 tests/dg/run_dg_wpw_production_operator_mpi.py`.
   Expect failure for missing plumbing/projection.
6. Add the input control with default `n`. Reject any configuration selecting
   more than one of diagonal, metric-block, or H-epsilon-S correction.
7. Store the new factor in `s_dg_wpw_production_context`; initialize it only
   when the comparison mode is selected and the fixed-H bounded operator is
   complete.
8. Add a matrix-free helper that globally computes `Q^dagger*S*z` through the
   existing Gram callback and subtracts `Q` times that block. Validate all
   shapes and finite values collectively.
9. Add `wpw_h_epsilon_s_precondition` in LCFO. Delegate to the common spectral
   correction, then invoke the global S projection before returning the
   correction to the three-block search.
10. Preserve the solver's current Ritz values, residual construction,
    retained-search recurrence, Ritz diagnostics, convergence, and publication
    gates byte-for-byte outside callback selection.
11. Log factor metadata once and selected-iteration correction diagnostics at
    inner 32/96/160.
12. Re-run all four focused commands; expect PASS.
13. Run the generalized-algebra, occupied-W basis, and row-layout MPI runners
    referenced by `tests/dg/CMakeLists.txt`; expect PASS.
14. Run the full MPI/EigenExa build with the repository's existing configured
    build procedure; require the `salmon` target to link successfully.
15. Run `git diff --check`.
16. Request code review focused on generalized-factor mathematics, sign and
    floor policy, collective ordering, S-projection, lifecycle, mutual
    exclusion, and default-off compatibility.
17. Resolve every Critical/Important finding and repeat Steps 12--15.
18. Commit only intended source/test changes:
    `git commit -m "feat(dg-wpw): integrate H-epsilon-S fixed-H route"`.

### Task 3: Run the B=6 physical gate

**Files:**
- Modify: `docs/plans/2026-07-23-wpw-h-epsilon-s-correction.md`
- Modify: `docs/plans/2026-07-19-wpw-density-carrying-seed-session-handoff.md`

1. Clone the Task 16 B=6 restart into a fresh run directory.
2. Keep the Task 16 basis, seed, cutoff, tolerance, continuation, history, and
   publication settings unchanged.
3. Set diagonal=`n`, metric-block=`n`, S-orthogonal complement=`n`, and
   H-epsilon-S correction=`y`.
4. Run 8 MPI ranks through the existing fixed-H boundary.
5. Record factor dimension, S condition, generalized spectral range,
   regularization, and epoch/fingerprint.
6. At inner 32/96/160 record occupied/extra raw and corrected residuals,
   amplification ratios, denominator floors, floored mode counts/weights,
   projected-away correction fraction, search-metric effective rank, and Ritz
   post/direct defects.
7. Compare occupied/extra residual histories, final `info`, and publication
   state with Task 16, Task 19, and the S-orthogonal gate.
8. Interpret improvement only if both occupied and extra residuals improve
   without Ritz inconsistency, severe effective-rank loss, or publication
   safety regression.
9. If occupied improves while extra worsens, reject only the fragment-local
   spectral approximation and use the recorded spectral weights to choose
   global iterative or explicitly state-partitioned follow-up.
10. Record the evidence and next hypothesis in this plan and the session
    handoff.
11. Run `git diff --check`.
12. Commit only the two result documents:
    `git commit -m "docs(dg-wpw): record H-epsilon-S B=6 gate"`.

#### 2026-07-23 B=6 physical gate result

Task 16 restartを読み取り専用の比較元として
`/tmp/20260723_task3_h_epsilon_s_b6`へcloneし、Task 2検証で作成した一時overlay
source/buildの`salmon`を使って8 MPI rankで実行した。このbinaryは親作業ツリーの
未コミット前提実装に依存するが、その前提実装はこのbranchへ取り込んでいない。
Task 16から変更した入力はdiagonal=`n`、metric-block=`n`、
S-orthogonal complement=`n`、H-epsilon-S=`y`、history=`y`の明示だけであり、
basis、seed、cutoff、tolerance、continuation、publication設定は維持した。
DC-SCFは同じ87反復でfixed-H境界へ到達した。

rank 0 fragment-local factor recordはdimension `138`、S spectrum
`[4.1210E-06,1.8573E+00]`、condition `4.5069E+05`、generalized
`H-epsilon S` spectrum `[-1.4491E-01,1.5978E+01]`、relative floor
`1.0000E-10`、epoch `1`、fingerprint `-6581501611198156709`だった。
selected iterationの証拠は次の通り。

| inner | occupied raw | extra raw | occupied correction norm/ratio | extra correction norm/ratio | denominator min occupied/extra | floored modes/weight | projected fraction | effective rank |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 32 | `2.2451E-03` | `1.7313E-02` | `3.9594E-01` / `1.7635E+02` | `6.4026E+00` / `3.6981E+02` | `1.3677E-05` / `8.4375E-06` | `0/0`, weight `0/0` | `2.1826E-01` | 205 |
| 96 | `2.2449E-03` | `1.7318E-02` | `3.9541E-01` / `1.7614E+02` | `6.4071E+00` / `3.6996E+02` | `1.3675E-05` / `8.5677E-06` | `0/0`, weight `0/0` | `2.1824E-01` | 205 |
| 160 | `2.2445E-03` | `1.7322E-02` | `3.9544E-01` / `1.7618E+02` | `6.4088E+00` / `3.6997E+02` | `1.3672E-05` / `8.7591E-06` | `0/0`, weight `0/0` | `2.1821E-01` | 206 |

最大inverse magnitudeはinner 32/96/160で`1.1852E+05`/
`1.1672E+05`/`1.1417E+05`だった。floorは一度も発火せず、約21.8%を
global S projectionが除去してもcorrectionはoccupiedを約176倍、extraを約370倍へ
増幅した。search metricのdiscarded residual fractionはinner 160でoccupied
`6.4268E-01`、extra `2.9668E-01`まで増えた。

Task 16比の残差比（H-epsilon-S / Task 16）はinner 32/96/160でoccupied
`1.142/6.222/8.765`、extra `4.282/15.612/85.376`であり、両state群が悪化した。
Task 19比ではinner 160がoccupied約`1.13E+03`倍、extra約`18.1`倍悪く、
S-orthogonal gate比でもoccupied約`10.9`倍、extra約`109`倍悪い。
Ritz post/direct relative defectの最大は`2.21E-12`、post metric orthogonalityの
最大は`1.85E-14`で、stale imageやRitz更新不整合は見られない。一方effective rankは
Task 16の`437/292/213`に対し`205/205/206`となり、重いrank lossが生じた。

fixed-Hは`info=40`で終了し、checkpoint/manifest/RT publicationは行われなかった。
従ってこのfragment-local dense spectral approximationはrejectする。floor未発火にも
かかわらず強い過増幅とsearch-space collapseが起きたため、次の仮説は単なるfloor調整では
なく、global iterative correction solve、またはoccupied/extraを明示分割してextra側の
inverse amplificationを制限するstate-partitioned correctionである。default-off routeは
比較用に維持し、normal SCF/checkpoint/RTへは昇格させない。

### Task 4: Decision checkpoint

Do not promote this correction into normal outer SCF, checkpoint schema,
publication, or RT handoff automatically.

1. Request user review of the B=6 gate.
2. If accepted, write a separate promotion design with explicit correction
   provenance and restart semantics.
3. If rejected, keep the default-off route available and design either the
   global iterative correction solve or a state-partitioned regularization
   based on Task 3 diagnostics.
