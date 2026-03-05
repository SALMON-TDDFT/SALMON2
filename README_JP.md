# 日本語ドキュメントガイド

SALMON v2.2.2 ワークスペース向けの、現在有効な日本語導線です。

## まず読む
1. `DOCUMENTATION_INDEX.md`
2. `doc/DOCUMENT_STATUS_MATRIX.md`

## カテゴリ別統合ノート
- DG: `doc/NOTE_DG.md`
- +U: `doc/NOTE_PLUSU.md`
- HSE: `doc/NOTE_HSE.md`

## 実装の一次参照（コード）
- DG: `src/rt/rt_dg_fragment.f90`, `src/rt/main_tddft.f90`
- +U: `src/rt/time_evolution_step.f90`, `src/common/density_matrix_and_energy_plusU_sub.f90`
- HSE: `src/xc/xc_hse.f90`, `src/xc/xc_hse_ri.f90`, `src/rt/rt_dg_fragment.f90`

## 最新実装メモ（要点）
- HSEのACEはRT主経路（`time_evolution_step` → `xc_ace_update_manager`）で動作し、`SALMON_ACE_RT` 系環境変数で有効化・更新条件を制御
- RT-DGのPW混合は `yn_plane_wave_basis`（`n_plane_waves_dg`, `k_cutoff_plane_wave`）で有効化され、`coef` と `coef_pw` を同時に時間発展
- adaptive basisは `yn_adaptive_basis` と `basis_update_threshold` で制御し、トリガ時は基底更新処理へ遷移
- DG側HSEは現状RI/DF経路が主で、実空間基底が無い場合は自動で無効化

詳細は以下を参照:
- DG: `doc/NOTE_DG.md`
- HSE: `doc/NOTE_HSE.md`

## 補足
- 旧DG/+Uノート群は統合に伴い削除済みです（2026-03-05）。
- 提案文書は `Proposal` として `doc/DOCUMENT_STATUS_MATRIX.md` に整理されています。
