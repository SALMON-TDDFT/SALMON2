# 引き継ぎ仕様書: DG transition-probe 検証 (2026-04-07)

## 1. 目的
- DG/DG+PW の電流を `diag/offdiag` に分解して、励起時コヒーレンス寄与を時系列で確認できるようにする。
- dense 行列化に依存せず、`momentum_blocks` を直接たどって集計する。

## 2. 変更内容 (コード)
- 変更ファイル:
  - `src/rt/dg/rt_dg_observables.f90`
- 実装ポイント:
  - `SALMON_DG_TRANSITION_PROBE` / `SALMON_DG_TRANSITION_PROBE_STRIDE` で出力制御。
  - `momentum_blocks` の各ブロック `(ifrag_row, ifrag_col)` を直接縮約して
    `Jdiag`, `Joffdiag`, `Jtotal` を計算。
  - `Jdiag` は同一ブロック内の対角要素 `(ifrag_row==ifrag_col かつ ib==jb)` のみを加算。
  - `Joffdiag = Jtotal - Jdiag`。
- 注意:
  - 以前の dense 化経路は使わない設計に変更済み。

## 3. ビルド状態
- ワークツリー:
  - `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/.worktrees/codex/dg-phi-box-cache-phase-a`
- ビルド:
  - `cmake --build build -j 8` 成功済み
- 重要:
  - MPI 実行を使うビルドでは `--enable-mpi` が必須。
  - セッションを変えた直後に `init_communicator_dft` で停止した場合は、コード変更より先に「そのビルドが `--enable-mpi` 付きか」を確認すること。
  - この停止は入力不整合に見えても、実際には MPI 無効ビルドの再使用が原因のことがある。

## 4. 実行と主要結果

### 4.1 成功ケース (パルス終端まで)
- 実行ディレクトリ:
  - `/tmp/h2_2mol_rt_dg_from_fresh_20260405`
- 入力:
  - `inputfile_nt1` をコピーし、`nt=240` のみ変更した `inputfile_nt240_from_nt1`
- 実行コマンド:
  - `OMP_NUM_THREADS=10 SALMON_DG_TRANSITION_PROBE=1 SALMON_DG_TRANSITION_PROBE_STRIDE=10 mpirun --bind-to none -np 4 <worktree>/build/salmon < inputfile_nt240_from_nt1 > run_nt240_transition_probe.log 2> run_nt240_transition_probe.err`
- 判定:
  - exit code: 0
  - `end SALMON`: あり
- 時間条件:
  - `nt=240`, `dt=0.05`, `tw1=6.0`
  - `T_final = nt*dt = 12.0 = 2*tw1` (パルス終端相当まで到達)

### 4.2 transition-probe 出力
- ログ:
  - `/tmp/h2_2mol_rt_dg_from_fresh_20260405/run_nt240_transition_probe.log`
- 出力行数:
  - `transition-probe:` 行は 24 点 (itt=1,11,...,231)
- 主要所見:
  - `Jdiag` は全サンプルで 0 (数値上ゼロ)
  - `Joffdiag` は有意に時間振動
  - 比率 `|Jdiag|/|Joffdiag|` の中央値 (x/y/z) はすべて 0

## 5. 失敗した試行と原因
- `/tmp/h2_2mol_rt_dg_from_fresh_20260405/rerun_dg` で `nt=240` 試行時:
  - `data_dcdft lg mismatch` により停止
  - 典型メッセージ: `STOP data_dcdft: input mismatch (lg)`
- DC-LCFO データを読む設定 (`yn_dg_fragment_from_dcdft='y'`) では、格子サイズ不一致時に失敗する。

## 6. 未完タスク (次担当)
- ユーザー要求の未完了項目:
  1. Full 電流と DG/DG+PW 電流の振幅差を定量比較
  2. 外場 `A(t)` と `J(t)` の位相差を Full vs DG/DG+PW で比較
- データ在庫監査結果 (`/tmp/h2_2mol_rt_dg_from_fresh_20260405`):
  - 比較に使える既存結果は一部あるが、未完タスクをそのまま完結できるだけの対応データは不足。
  - 2mol 実データ `H2_periodic_2mol_dc_2frag_short_rt.data` と `H2_periodic_2mol_dc_2frag_short_pulse.data` は存在。
  - `variables.log` 上で `yn_dg_fragment_rt=y` のため、上記は DG 側データ。
  - 2mol 同条件の Full 側波形が見当たらない。
  - 20mol 側は `H2_periodic_20_dc_2frag_smoke_rt.data` があるが、現行ファイルはヘッダのみで比較波形として不十分。
  - `rerun_full` / `rerun_dg` には入力・ログはあるが、対応 `rt.data` が残っていない。
- 推奨手順:
  1. 2mol で Full を再実行し、DG と同条件・同時間軸の `rt.data` を再生成する。
  2. 20mol は Full/DG の両側を再出力して、ヘッダのみファイルを置換する。
  3. 同一窓で振幅指標 (peak-to-peak / RMS / 基本波振幅) を計算する。
  4. 位相差はクロススペクトルまたは基本周波数フィットで算出する。
  5. 結果を比較表 (Full基準差分) でまとめる。

## 7. 現在の git 状態 (worktree)
- `M src/rt/dg/rt_dg_observables.f90`
- `?? build_check_bounds/`
- `?? build_dbg/`

## 8. 参照ログ/ファイル
- コード:
  - `src/rt/dg/rt_dg_observables.f90`
- 成功ログ:
  - `/tmp/h2_2mol_rt_dg_from_fresh_20260405/run_nt240_transition_probe.log`
  - `/tmp/h2_2mol_rt_dg_from_fresh_20260405/run_nt240_transition_probe.err`
- 入力:
  - `/tmp/h2_2mol_rt_dg_from_fresh_20260405/inputfile_nt240_from_nt1`
