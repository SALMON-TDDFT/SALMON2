# NOTE_DG

Status: Current (index note)

このノートは DG-Fragment RT 関連ドキュメントの統合入口です。

## 1. 対象スコープ
- DG-Fragment RT の実装・運用・検証
- 断片基底の読み込み、時間発展、ハミルトニアン再構成

## 2. 一次情報（コード実体）
- `src/rt/rt_dg_fragment.f90`
- `src/rt/main_tddft.f90`

## 3. 現在の扱い
- DG関連Markdownは作成時点が古いものが多いため、**参照専用**として扱う
- 仕様判断はコード実体を優先する

## 4. 主要パラメータ（参照）
- `yn_dg_fragment_rt`
- `yn_dg_fragment_from_dcdft`

## 5. 旧資料リンク（Archive/Reference）
- 旧DGノート群は統合に伴い削除済み（2026-03-05）
- 以後の一次参照は `src/rt/rt_dg_fragment.f90` と `src/rt/main_tddft.f90`

## 6. 更新ルール
- DG仕様変更時は、まず本ノートの「一次情報」と「主要パラメータ」を更新する
- 詳細ノートを追加する場合は、冒頭に `Status:` を明記する

## 7. 最新実装メモ（PW / adaptive basis）

### 7.1 RT-DGでのPW（plane wave）扱い
- 入力は `src/io/salmon_global.f90` の `yn_plane_wave_basis`, `n_plane_waves_dg`, `k_cutoff_plane_wave`
- 初期化は `src/rt/rt_dg_plane_wave.f90:init_plane_wave_basis`
	- `yn_plane_wave_basis='y'` のときのみ有効
	- カットオフ球内のk点（Γ除外）を作成し、`n_plane_waves_dg` で上限をかける
	- `coef_pw`, `k_pw`, 混合ハミルトニアン/重なり行列（`H_mat_mixed`, `S_mat_mixed_prop`）を確保
- 伝播は fragment係数 `coef` と PW係数 `coef_pw` を同時更新（RK系実装）
- 基底更新時（mixedモード）は `diagonalize_mixed_basis_pw` + `project_coefficients_mixed_state` で状態を新しい混合基底へ射影する

### 7.2 RT-DGでのadaptive basis扱い
- 入力は `yn_adaptive_basis`, `basis_update_threshold`
- 判定対象は「場なし成分のハミルトニアン変化」で、`A(t)` 依存成分は更新トリガ指標から外す方針（PW拡張側で吸収）
- トリガ時は当該ステップの `H_mat` 反映を一旦戻してから `trigger_basis_update` を実行
- 非PWモードでは DC-CG更新→全体系対角化→重なり行列健全性検査→係数射影
- PW混合モードでは「伝播状態（fragment+PW係数）」更新を行い、実空間 `phi_frag` 自体は置換しない

### 7.3 補足（HSEとの接続）
- DGループのHSEは基本的に RI/DF 経路（`init_hse_ri_data`）で、実空間基底が無い場合は自動でOFF
