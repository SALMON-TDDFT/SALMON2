# DG Buffered DC-RT Audit Checklist (JA)

## 1. 目的

本チェックリストは、DC で生成した buffered 基底・状態ファイルと、RT での読み込み・利用が最初から整合しているかを段階的に監査するための実務用チェックリストである。

対象は次の 5 段階である。

1. DC 側の buffered 基底計算
2. DC 側の buffered 出力
3. RT 側の読み込み
4. RT 側の演算子構築と初期状態整合
5. 実行時検証

本書は次の既存仕様を補完する。

- [DG_BUFFERED_RT_SEED_CONTRACT_SPEC_JA.md](/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/DG_BUFFERED_RT_SEED_CONTRACT_SPEC_JA.md)

## 2. 監査の基本原則

- buffered 基底と buffered 状態系列は、同じ state labeling を共有しなければならない
- DC 出力と RT 読み込みは、metadata だけでなく基底表現の意味も一致しなければならない
- startup Lowdin / startup projection は補正であり、契約不整合の代替ではない
- no-kick / fixed-H で drift が出る場合は、まず DC→RT 契約不整合を疑う

## 3. Stage A: DC 側 buffered 基底計算

### A-1. buffered basis の定義

- [ ] `basis_functions_buffered.bin` は `basis_functions.bin` と同じ `lambda` 基底定義を持つ
- [ ] 違うのは実空間 support のみである
- [ ] DC 側の `basis_transform` が buffered / non-buffered の両方で同じ基底 ordering を保っている
- [ ] buffered box の再配置規約が RT 側の `phi_frag` 配置と一致する

確認対象コード:

- [lcfo.f90](/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/.worktrees/codex/dg-phi-box-cache-phase-a/src/gs/dc/lcfo.f90)

確認観点:

- `wrk_buffer -> buffered_basis` の並べ替え
- `basis_transform(jo,io,ispin)` の適用方向
- core / pos_buf / neg_buf の並べ替え規約

### A-2. local Hamiltonian と基底の対応

- [ ] `mat_H_local` が buffered basis の state ordering と矛盾していない
- [ ] fragment-local diagonalization に使う行列と出力基底の ordering が一致する

## 4. Stage B: DC 側 buffered 出力

### B-1. 必須ファイルの同時生成

- [ ] `basis_functions_buffered.bin`
- [ ] `wavefunctions_buffered.bin`
- [ ] `occupations_buffered.bin`
- [ ] `eigenvalues_buffered.bin`

が同一 run で同時生成されている

### B-2. metadata 一致

各ファイルについて次が一致すること:

- [ ] `n_frag`
- [ ] `nspin`
- [ ] `nstate_frag`
- [ ] `nstate_tot`
- [ ] `n_mat`
- [ ] `n_basis`
- [ ] `index_basis`

### B-3. state lineage 一致

- [ ] `coef(:,i)` と `esp(i)` が同じ generalized EVP 由来
- [ ] `occ(i)` も同じ state labeling を参照
- [ ] `esp=0` の未設定状態が残っていない

### B-4. buffered / non-buffered 混在禁止

- [ ] buffered basis を使う state set に対して `wavefunctions.bin` を混ぜていない
- [ ] non-buffered basis を使う state set に対して `wavefunctions_buffered.bin` を混ぜていない

### B-5. 実空間成分ダンプ診断

- [x] `src/gs/dc/lcfo.f90` の `hpsi_basis` に、代表 1 状態の `psi/Tpsi/Vlocpsi/Vnlpsi/Hpsi` を出す診断経路を追加した
- [x] `Tpsi` は `hpsi(..., ttpsi)` optional 出力から取得する実装になっている
- [x] 選択状態のみを 3 次元配列で集約する形にして `comm_summation` の generic 解決エラーを回避した
- [x] コンパイル確認は通過した
- [x] `&dc` を含む checker 整合入力で `dc_rs_component_if1_sp1_io1_fragroot.dat` の run-time 生成を確認した
- [x] サンプル直下 run では `hubbard_param.dat` 自動検出により無関係の +U 初期化例外へ入るため、検証 run は `/tmp` から実施した
- [x] MPI 複数 rank での `SIGTRAP` は、rank 選択条件内で collective を呼んでいたことが原因と特定した
- [x] `src/gs/dc/lcfo.f90` のダンプ部を修正し、全 rank が `comm_summation` に参加する実装へ変更した
- [x] 修正後、C64 条件の MPI 8 実行で `[DC-RS-COMP]` ログと dat 生成、`end SALMON` を確認した

## 5. Stage C: RT 側読み込み

### C-1. strict 読み込み

- [ ] buffered mode で `eigenvalues_buffered.bin` を strict 読み込みする
- [ ] `SALMON_DG_ALLOW_MISSING_BUFFERED_ESP` に依存せず正常動作する
- [ ] `SALMON_DG_ENFORCE_STATE_TRUNCATION=1` で正常動作する

### C-2. cutoff 一致

- [ ] `nstate_keep = min(runtime,file)` が `coef/occ/esp` に共通適用される
- [ ] cutoff を超える `index_basis` 参照が無い
- [ ] truncation 後の state 集合に対して metadata が still valid

### C-3. helper consistency

- [ ] `SALMON_DG_HELPER_CONSISTENCY_CHECK=1` で FATAL が出ない
- [ ] fragment 間で `n_mat` / `n_basis` / `index_basis` が一致する

確認ログ:

- `DG helper consistency check: ON`
- `SEED-CONSISTENCY-OK`

## 6. Stage D: RT 側利用整合

### D-1. overlap 整合

- [ ] `C^† S C` の対角が 1 に近い
- [ ] `offdiag_rms(C^† S C)` が小さい
- [ ] DC 想定 overlap と RT の `apply_overlap_operator_batch` が一致する

確認ログ候補:

- `HS-EIGENSTATE`
- `HS-OVERLAP-COMPARE`

### D-2. frozen H 整合

- [ ] `Hc - eSc` 残差が小さい
- [ ] `Rayleigh-esp` 差が小さい
- [ ] startup-pre と startup-post の差が理解可能

確認ログ候補:

- `HS-EIGENSTATE`
- `DG-FIXED-H`
- `DG-FIXED-H-RAY`

### D-3. core / nonlocal 契約

- [ ] `H_mat` は core-only であることを理解している
- [ ] buffered canonical full reference は `H_mat(core) + H_nl_cache`
- [ ] buffered canonical Vnl は self-only

確認ログ候補:

- `HS-HC-COMPARE-NOTE`
- `HS-HC-CONTRACT`

### D-4. nonlocal projector 整合

- [ ] `uVphi_self` と direct overlap が一致
- [ ] `proj_local`, `proj_global`, `proj_weight`, `contrib`, `y` が canonical 定義と一致
- [ ] deprecated な dense local/global 指標を正式判定に使っていない

確認ログ候補:

- `VNL-STAGE-TRACE`
- `HS-VNL-COMPARE`
- `HS-VNL-ROUTE-SPLIT`

### D-5. m 項の位相生成

- [ ] fragment basis 側 `m` 項に `-i` が正しく入っている
- [ ] `dcoef_dt_m(im2) > 0` が kick 時に立つ
- [ ] no-kick では `dcoef_dt_m = 0`

確認ログ候補:

- `prop-im-trace`
- `rk-im-trace`
- `coef-state-trace`
- `dg-para-pair-imtrace`

## 7. Stage E: 実行時物理チェック

### E-1. no-kick 基線

- [ ] no-kick で `J_dia = 0`
- [ ] no-kick で `J_total = J_para`
- [ ] no-kick で `J_para` が 0 近傍に留まる

もし FAIL なら:

- [ ] initial state と frozen `H,S` の不整合を再確認
- [ ] startup projection の有無で drift が減るか確認

### E-2. kick 応答

- [ ] kick/no-kick で `J_para` 差分が出る
- [ ] `J_dia` が定常 offset として立つ
- [ ] 時間発展で `J_para` が `J_dia` を打ち消す方向へ成長する

### E-3. fixed-H sanity

- [ ] fixed-H / no-kick で `J_para` drift が極小
- [ ] fixed-H / kick で expected な初期応答が出る
- [ ] `itt=0 -> 1` で residual が急悪化しない

## 8. 推奨監査順

次の順に潰す。

1. Stage B metadata/state-lineage
2. Stage C strict 読み込み
3. Stage D overlap/frozen-H 整合
4. Stage D nonlocal / m-term 整合
5. Stage E no-kick 基線
6. Stage E kick 応答

## 9. PASS/FAIL 判定メモ

### PASS と見てよい例

- `SEED-CONSISTENCY-OK`
- `HS-OVERLAP-COMPARE rel_norm_diff ~ 0`
- `HS-TVLOC-COMPARE route-dense ~ 0`
- `VNL-STAGE-TRACE LOCAL-CHECK status=PASS`
- kick で `coef im_norm > 0`, no-kick で `coef im_norm = 0`

### FAIL と見てよい例

- strict mode で metadata/cutoff FATAL
- `C^†SC` の対角が 1 から大きく外れる
- `Rayleigh-esp` が大きくずれる
- no-kick で `J_para` が大きく drift
- kick/no-kick で `J_para` 差分が立たない

## 10. 実行メモ欄

監査 run ごとに次を書き残す。

- [ ] input file
- [ ] data directory
- [ ] binary path
- [ ] commit or patch state
- [ ] enabled `SALMON_DG_*` env vars
- [ ] relevant log files
- [ ] PASS/FAIL summary

### 2026-04-23 監査実行ログ（codex/dg-phi-box-cache-phase-a）

対象:

- input file: `.namelist_dg_rt_c64_dt002_empties32_nopw_nokick.tmp`（派生 `nt=1`, `nt=200`）
- data directory: `data_dcdft` -> `data_dcdft_step1_new`（symlink）
- binary path: `.worktrees/codex/dg-phi-box-cache-phase-a/build/salmon`
- env vars（主）:
	- `SALMON_DG_ENABLE_PROBES=1`
	- `SALMON_DG_HELPER_CONSISTENCY_CHECK=1`
	- `SALMON_DG_ENFORCE_STATE_TRUNCATION=1`
	- `SALMON_DG_ALLOW_MISSING_BUFFERED_ESP=0`
	- `SALMON_DG_CURRENT_COMPONENT_TRACE=1`
	- `SALMON_DG_STARTUP_LOWDIN=1`
	- `SALMON_DG_STARTUP_STATIONARY_PROJECTION=0/1`（A/B）
	- `SALMON_DG_ALLOW_BUFFERED_STARTUP_STATIONARY_PROJECTION=1`
- relevant logs:
	- `run_diag_nokick_nt1_strict_read.log`
	- `run_diag_nokick_nt1_proj0_trace.log`
	- `run_diag_nokick_nt1_proj1_trace.log`
	- `run_diag_nokick_nt200_l1p0.log`
	- `run_diag_nokick_nt200_l1p1.log`

今回の判定（実測反映）:

- Stage B-1 必須ファイル同時生成
	- [x] `basis_functions_buffered.bin`
	- [x] `wavefunctions_buffered.bin`
	- [x] `occupations_buffered.bin`
	- [x] `eigenvalues_buffered.bin`
	- 判定根拠: `data_dcdft_step1_new/fragments/000001..000008` の全 fragment で存在を確認

- Stage C-1 strict 読み込み
	- [x] buffered mode で strict 読み込みが完走（`EXIT_STRICT:0`）
	- [x] `SALMON_DG_ALLOW_MISSING_BUFFERED_ESP=0` 依存で FATAL なし
	- [x] `SALMON_DG_ENFORCE_STATE_TRUNCATION=1` で完走

- Stage C-3 helper consistency
	- [x] `DG helper consistency check: ON` を確認
	- [x] `SEED-CONSISTENCY-OK` の明示ログ確認（`run_diag_nokick_nt1_seedcheck.log`）
	- [ ] 注意: 同 run は終了コード 133（起動後の既知不安定）。seed metadata ログ自体は出力済み

- Stage E-1 no-kick 基線（Lowdin=1, projection ON/OFF 比較）
	- [x] no-kick で `J_dia = 0`
	- [x] no-kick で `J_total = J_para`
	- [ ] no-kick で `J_para` が 0 近傍維持（未達）
	- 定量（nt=200）:
		- projection OFF: `|J_para,z|_max = 6.619228e-3`, `RMS = 2.994857e-3`
		- projection ON : `|J_para,z|_max = 5.244387e-3`, `RMS = 2.078411e-3`
		- 改善率（ON/OFF）: max 0.792, RMS 0.694
	- 補足: startup projection は drift を低減するが根治ではない

- Stage D-2 frozen H 整合（nt=1, fixed-H residual trace）
	- [ ] `Hc - eSc` 残差が小さい（未達）
	- [ ] `itt=0 -> 1` で residual 急悪化なし（未達）
	- 定量（`DG-FIXED-H resid_abs_max`）:
		- Lowdin=1, projection OFF: `1.1071e+0 -> 5.6203e+0`
		- Lowdin=1, projection ON : `3.2812e-1 -> 1.1218e+1`
	- 判定: startup projection で `itt=0` は改善するが、`itt=1` で急増するケースがあり frozen-H 一貫性は未解決

- Stage D-1 overlap 整合（HS診断）
	- [ ] `C^† S C` 対角/offdiag の定量確認（未取得）
	- [ ] `HS-OVERLAP-COMPARE` の rel/max 差分取得（未取得）
	- 実行メモ: `run_diag_nokick_nt1_hscheck.log` は Hamiltonian 準備完了まで進んだが、その後停滞。`HS-EIGENSTATE` 系ログ未出力のため再実行条件の見直しが必要

- Stage D-5 m 項位相生成（既存ログ再監査）
	- [x] kick 時に `prop-im-trace m(re2/im2)` の虚部チャネルが立つ
	- [x] no-kick では `prop-im-trace m(re2/im2)` がゼロ
	- [x] kick 時に `coef-state-trace` / `dg-para-pair-imtrace` の虚部が非ゼロ
	- 判定根拠ログ:
		- kick: `run_kick_nt2_force_i_mterm.log`
		- no-kick: `run_nokick_nt2_propim_trace.log`, `run_nokick_nt2_force_i_mterm.log`

PASS/FAIL summary:

- PASS: Stage B-1, Stage C-1（strict read smoke）, Stage E-1 の構造条件（`J_dia=0`, `J_total=J_para`）
- FAIL/継続課題: Stage E-1 の「`J_para≈0`維持」、seedcheck run の終了安定性（exit 133）

## 11. 現時点の既知論点

現時点の切り分けで、少なくとも次は既知論点である。

- startup projection は no-kick 基線 drift を低減するが、根治ではない
- buffered/self-only canonical Vnl 契約は妥当
- nonlocal projector route は canonical 定義と整合
- fragment basis 側 `m` 項の `i` 欠落は実際に kick 応答を殺していた
- それでも最終的な `J_total` が 0 中心へ戻るかは、まだ no-kick 基線の汚染に影響される

以上を踏まえ、今後の監査では「まず no-kick 基線を潰す」ことを優先する。
