# DOCUMENT STATUS MATRIX

この一覧は、SALMON v2.2.2 ワークスペース内の Markdown を実運用目線で分類したものです。

- **Current**: 現行の参照先として推奨
- **Reference**: 参照価値はあるが、最新仕様の一次情報ではない
- **Archive**: 履歴・記録（再現や経緯確認向け）
- **Proposal**: 提案・計画（未確定）

---

## Root (`./`)

| File | Status | Notes |
|---|---|---|
| README.md | Current | 上流 SALMON 全体の入口 |
| README_JP.md | Current | 日本語入口 |
| DOCUMENTATION_INDEX.md | Current | 本リポジトリ向け索引 |
| doc/NOTE_DG.md | Current | DG系統合ノート |
| doc/NOTE_PLUSU.md | Current | +U系統合ノート |
| doc/NOTE_HSE.md | Current | HSE系統合ノート |
| RT_INTEGRATION_COMPLETE.md | Current | Plan C 統合完了の要約 |
| HSE_IMPLEMENTATION_CHECK_REPORT.md | Reference | HSE 実装点検レポート |
| HSE_IMPLEMENTATION_REPORT.md | Reference | Plan A 実装レポート |
| RT_DG_HSE_COMPLETE_IMPLEMENTATION_v2.2.2.md | Reference | 統合ガイド（網羅） |
| CD_RI_IMPLEMENTATION.md | Reference | CD-RI 実装報告 |
| DISTANCE_SCREENING_IMPLEMENTATION.md | Reference | 距離スクリーニング報告 |
| GITHUB_BUGFIX_FOR_DEVELOPERS.md | Reference | CMake/開発者向け修正ガイド |
| HSE_RI_IMPLEMENTATION_README.md | Archive | "Partial" 記載を含む移行期メモ |
| SALMON_v2.2.2_COMPREHENSIVE_EVALUATION.md | Archive | 総合評価（時点依存） |
| ADAPTIVE_BASIS_IMPLEMENTATION.md | Proposal | 適応基底の提案/実装メモ |
| ADAPTIVE_BASIS_IMPLEMENTATION_ja.md | Proposal | 同上（日本語） |
| dg_fragment_rk4_adaptive_notes.md | Proposal | RK4/基底更新ノート |
| HSE_RI_DF_METHOD_PLAN_C.md | Proposal | Plan C 提案ドキュメント |
| HSE_ERI_OPTIMIZATION_PROPOSAL.md | Proposal | HSE 最適化提案 |
| GPU_IMPLEMENTATION_GUIDE.md | Proposal | GPU 実装方針 |
| GPU_URGENCY_PLAN_Si64.md | Proposal | GPU 導入優先度整理 |
| OPENACC_GPU_IMPLEMENTATION.md | Proposal | OpenACC 方針 |

---

## `doc/`

| File | Status | Notes |
|---|---|---|
| doc/NOTE_DG.md | Current | DG系統合ノート |
| doc/NOTE_PLUSU.md | Current | +U系統合ノート |
| doc/NOTE_HSE.md | Current | HSE系統合ノート |
| doc/DOCUMENT_STATUS_MATRIX.md | Current | 本分類表 |

---

## `samples/exercise_dg_rt_hse_test/`

| File | Status | Notes |
|---|---|---|
| samples/exercise_dg_rt_hse_test/README.md | Reference | テストケース総覧（実行条件は要再確認） |
| samples/exercise_dg_rt_hse_test/START_HERE.md | Reference | 実行導線（更新遅延の可能性） |
| samples/exercise_dg_rt_hse_test/EXECUTION_GUIDE.md | Reference | 実行手順（更新遅延の可能性） |
| samples/exercise_dg_rt_hse_test/DG_FRAGMENT_WORKFLOW.md | Reference | GS→RT ワークフロー（参照用） |
| samples/exercise_dg_rt_hse_test/LDA_BASELINE_TEST.md | Reference | ベースライン検証 |
| samples/exercise_dg_rt_hse_test/TEST_PREPARATION.md | Archive | 準備完了記録 |
| samples/exercise_dg_rt_hse_test/WORKFLOW_CORRECTED.md | Archive | 修正履歴 |

---

## Other Markdown (upstream/tooling)

| File | Status | Notes |
|---|---|---|
| src/ttm/README.md | Current | サブモジュール説明 |
| utility/ttm/README.md | Current | ユーティリティ説明 |
| utility/linear_response_from_transition_moments/README.md | Current | ツール説明 |
| utility/linear_response_from_transition_moments/example/README.md | Current | 例題説明 |
| samples/exercise_dg_fragment_C2H2/README.md | Current | サンプル説明 |
| samples/exercise_dg_fragment_rt/README.md | Current | サンプル説明 |

---

## 運用ルール

1. 新規作成時に `Status: Current/Reference/Archive/Proposal` を冒頭に追記する。  
2. Proposal が実装完了したら、Proposal 側に `Superseded by: ...` を追記する。  
3. 索引更新時は `DOCUMENTATION_INDEX.md` と本ファイルを同時更新する。
