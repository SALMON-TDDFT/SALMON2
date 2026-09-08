# NOTE_PLUSU

Status: Current (index note)

このノートは DFT+U（+U）関連ドキュメントの統合入口です。

## 1. 対象スコープ
- RT-TDDFT 中の +U 密度行列更新
- +U ポテンシャル更新と電流評価への寄与

## 2. 一次情報（コード実体）
- `src/rt/time_evolution_step.f90`
- `src/common/density_matrix_and_energy_plusU_sub.f90`
- `src/common/calc_jxyz_plusu_sub.f90`

## 3. 現在の扱い
- +Uの主要ノートはバグ報告起点のため、実運用時はコード実体で確認する

## 4. 既知の重要論点
- RT中に `calc_density_matrix_and_energy_plusU` を更新呼び出しすること
- 電流計算に +U 寄与を含めること

## 5. 旧資料リンク（Reference）
- 旧+Uノートは統合に伴い削除済み（2026-03-05）

## 6. 更新ルール
- +U修正を入れたら、本ノートの「既知の重要論点」と対象コードを更新する
