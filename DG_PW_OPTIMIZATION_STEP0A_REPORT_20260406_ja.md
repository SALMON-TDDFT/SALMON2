# DG+PW 最適化 Step0/A 実施結果

作成日: 2026-04-06
ブランチ: codex/dg-phi-box-cache-phase-a
基準コミット: ac9606ae

## 1. 実施条件

- MPI: 4 ranks
- OMP (性能比較): 10
- OMP (安定性確認): 2
- 入力: /tmp/h2_periodic20_rt2frag_2x2_smoke_clean_20260406/inputfile_pw_on
- 実行: /usr/bin/time -p mpirun --bind-to none -np 4 <repo>/build/salmon

## 2. Step 0 (baseline)

- run_baseline_1: exit=0, end SALMON=yes, fatal/signal=no, real=0.42
- run_baseline_2: exit=0, end SALMON=yes, fatal/signal=no, real=0.41
- run_baseline_3: exit=0, end SALMON=yes, fatal/signal=no, real=0.40
- median_baseline = 0.41 s

## 3. Step 1 (施策Aのみ)

変更対象:
- src/rt/dg/rt_dg_density_reconstruct.f90

変更要点:
- 共有 rho_send への OpenMP atomic 更新を削減
- 並列区間では remote 寄与をブロック局所配列へ格納
- 並列外で一括集約して rho_send へ反映

評価:
- run_A_1: exit=0, end SALMON=yes, fatal/signal=no, real=1.14
- run_A_2: exit=0, end SALMON=yes, fatal/signal=no, real=0.36
- run_A_3: exit=0, end SALMON=yes, fatal/signal=no, real=0.31
- median_A = 0.36 s
- improvement_A = (0.41 - 0.36) / 0.41 * 100 = 12.20%

判定:
- 妥当性チェック全通過
- 改善率 >= 3%
- 採用候補

## 4. Step 3 (最終確認: 採用構成=Aのみ)

- run_final_1: exit=0, end SALMON=yes, fatal/signal=no, real=0.41
- run_final_2: exit=0, end SALMON=yes, fatal/signal=no, real=0.31
- run_final_3: exit=0, end SALMON=yes, fatal/signal=no, real=0.38
- median_final = 0.38 s
- improvement_final = (0.41 - 0.38) / 0.41 * 100 = 7.32%

OMP=2 安定性:
- run_final_omp2_1: exit=0, end SALMON=yes, fatal/signal=no, real=0.27

## 5. 結論

- 施策Aは採用可能。
- baseline 比で中央値改善を確認 (7.32%)。
- 次段として施策B (array temporary 警告経路の個別解消) へ進む。
