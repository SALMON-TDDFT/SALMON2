# NOTE_HSE

Status: Current (index note)

このノートは HSE（Plan A / Plan C / CD-RI）関連ドキュメントの統合入口です。

## 1. 対象スコープ
- HSE交換項のDG-Fragment RT統合
- RI/DF（Plan C）および CD-RI による高速化・省メモリ化

## 2. 一次情報（コード実体）
- `src/xc/xc_hse.f90`
- `src/xc/xc_hse_ri.f90`
- `src/xc/auxiliary_basis.f90`
- `src/rt/rt_dg_fragment.f90`

## 3. 現在の扱い
- HSEノートは「実装報告」「提案」「移行メモ」が混在している
- 現行判断はコード実体 + 統合完了レポートを優先する

## 4. 推奨参照順
1. `RT_INTEGRATION_COMPLETE.md`
2. `HSE_IMPLEMENTATION_CHECK_REPORT.md`
3. コード実体（上記）

## 5. 旧資料リンク（Reference/Archive/Proposal）
- `HSE_IMPLEMENTATION_REPORT.md`
- `HSE_RI_IMPLEMENTATION_README.md`
- `RT_DG_HSE_COMPLETE_IMPLEMENTATION_v2.2.2.md`
- `CD_RI_IMPLEMENTATION.md`
- `DISTANCE_SCREENING_IMPLEMENTATION.md`
- `HSE_RI_DF_METHOD_PLAN_C.md`
- `HSE_ERI_OPTIMIZATION_PROPOSAL.md`
- `GPU_IMPLEMENTATION_GUIDE.md`
- `GPU_URGENCY_PLAN_Si64.md`
- `OPENACC_GPU_IMPLEMENTATION.md`

## 6. 更新ルール
- HSE仕様変更時は、本ノートの「推奨参照順」と「一次情報」を先に更新する

## 7. 最新実装メモ（ACE / RT系）
- ACEは RT 主経路で `src/rt/time_evolution_step.f90` から `src/xc/xc_ace_update_manager.f90` を呼び出して動作する
- 有効化は入力ファイル項目ではなく環境変数 `SALMON_ACE_RT`（`1/y/true/on`）
- 更新判定は2段階:
	- 第1段階: 波動関数変化量 `d_max`（`ace_first_stage_update_trigger`）
	- 第2段階: サンプル残差 `R_max`（`ace_stage2_residual_check`）
- `R_max > eps_R` または強制更新条件で ACE 再構築を実行し、成功時のみ適用を継続
- 交換作用素適用は `apply_exchange_ACE` を通して `phi_curr` に直接加える（RTステップ内）

### 7.1 主要ACE環境変数
- ON/OFF: `SALMON_ACE_RT`
- しきい値: `SALMON_ACE_EPSD`, `SALMON_ACE_EPSR`
- 更新間隔: `SALMON_ACE_NMIN`, `SALMON_ACE_NMAX`
- サンプリング: `SALMON_ACE_NV`, `SALMON_ACE_TOPM`, `SALMON_ACE_NRAND`, `SALMON_ACE_SEED`
- ログ: `SALMON_ACE_LOG_EX`, `SALMON_ACE_EX_STRIDE`
- 数値安定化: `SALMON_ACE_DIAG_SHIFT`, `SALMON_ACE_KSR_MIN`

### 7.2 DG-Fragmentとの関係
- 現在のDG-Fragment側HSEは `src/rt/rt_dg_fragment.f90` の RI/DF（`yn_hse_ri`）経路が主で、ACE管理器はDG専用ループには接続していない
- DGで `yn_hse='y'` かつ `yn_hse_ri='y'` でも、`phi_frag`（実空間基底）が無い場合は RI/DF を自動無効化する
