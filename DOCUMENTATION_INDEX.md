# SALMON v2.2.2 Documentation Index (整理版)

このファイルは、散在している Markdown を「今読むべきもの」と「履歴・提案」に分離するためのハブです。

---

## 1) まずここから（現行導線）

1. `README.md`  
   - プロジェクト本体（上流 SALMON）の概要

2. `README_JP.md`  
   - 日本語向けの入口

3. `doc/NOTE_DG.md`  
   - DG系統合ノート（入口）

4. `doc/NOTE_PLUSU.md`  
   - +U系統合ノート（入口）

5. `doc/NOTE_HSE.md`  
   - HSE系統合ノート（入口）

6. `doc/DOCUMENT_STATUS_MATRIX.md`  
   - 文書の現行/履歴/提案の分類表

7. `DG_DC_TO_DG_RT_USAGE_JA.md`
   - DC から DG-RT へ buffered seed を渡す現行手順

注: DG関連の説明文書は更新遅延があるため、現時点では参照専用として扱います。

---

## 2) テスト実行ドキュメント（参照）

- `samples/exercise_dg_rt_hse_test/README.md`
- `samples/exercise_dg_rt_hse_test/START_HERE.md`
- `samples/exercise_dg_rt_hse_test/EXECUTION_GUIDE.md`
- `samples/exercise_dg_rt_hse_test/DG_FRAGMENT_WORKFLOW.md`

補助資料:
- `samples/exercise_dg_rt_hse_test/LDA_BASELINE_TEST.md`
- `samples/exercise_dg_rt_hse_test/TEST_PREPARATION.md`
- `samples/exercise_dg_rt_hse_test/WORKFLOW_CORRECTED.md`

---

## 3) 履歴・セッション記録（参照専用）

旧DG/+Uノート群は統合に伴い削除済み（2026-03-05）。
残る履歴資料は以下です。

- `SALMON_v2.2.2_COMPREHENSIVE_EVALUATION.md`

---

## 4) 提案・設計案（未確定/条件付き）

以下は「実装済み仕様」ではなく、提案・比較・計画文書です。

- `HSE_RI_DF_METHOD_PLAN_C.md`
- `HSE_ERI_OPTIMIZATION_PROPOSAL.md`
- `GPU_IMPLEMENTATION_GUIDE.md`
- `GPU_URGENCY_PLAN_Si64.md`
- `OPENACC_GPU_IMPLEMENTATION.md`
- `ADAPTIVE_BASIS_IMPLEMENTATION.md`
- `ADAPTIVE_BASIS_IMPLEMENTATION_ja.md`
- `DUAL_RHO_VH_ONLY_IMPLEMENTATION_SPEC_JA.md`
- `dg_fragment_rk4_adaptive_notes.md`

---

## 5) HSE関連の扱い（重要）

HSE RI/DF 関連には、作成時点の異なる文書が混在しています。

最新実装の要点（2026-03時点）:
- ACE は標準RT経路（`src/rt/time_evolution_step.f90` → `src/xc/xc_ace_update_manager.f90`）で利用され、`SALMON_ACE_RT` 系環境変数で制御
- RT-DG の PW/adaptive basis は DG統合ノート（`doc/NOTE_DG.md`）を一次参照
- DG側HSEは現状 RI/DF 経路が主で、実空間基底が無い場合は自動無効化

- 現行側の基準:
   - `doc/NOTE_HSE.md`
  - `RT_INTEGRATION_COMPLETE.md`
  - コード実体: `src/xc/xc_hse_ri.f90`, `src/rt/rt_dg_fragment.f90`

- 履歴/移行メモとして扱う:
  - `HSE_RI_IMPLEMENTATION_README.md`（「Partial Implementation」記載を含む）

---

## 6) 分類表

全 Markdown の分類一覧は以下を参照してください。

- `doc/DOCUMENT_STATUS_MATRIX.md`

---

## 7) メンテナンスルール（今後）

1. 新規文書は、必ず「現行 / 履歴 / 提案」のいずれかを冒頭に明記する。  
2. 実装完了時は、提案文書に「Superseded by: ...」を追記する。  
3. `DOCUMENTATION_INDEX.md` と `doc/DOCUMENT_STATUS_MATRIX.md` を同時更新する。

