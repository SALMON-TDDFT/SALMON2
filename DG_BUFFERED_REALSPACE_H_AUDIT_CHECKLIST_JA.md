# DG Buffered 実空間 Hamiltonian 監査チェックリスト

## 0. 目的
- DC 側フラグメント計算と DG/RT 側フラグメント計算について、まず基底行列ではなく実空間作用 `Hψ(r)` の定義が一致しているかを確認する。
- 比較対象は、同じフラグメント、同じ buffered box、同じ実空間波動関数 `ψ(r)` に対する作用結果 `Hψ(r)` とする。
- 本チェックは、`<φ_i|H|φ_j>` 比較や startup projection の前段に置く。

## 0.1 core 弱形式の設計原則（必須）

この監査の最終目標は、以下の原則を満たす canonical self Hamiltonian を構成すること。

- [ ] 基底自由度（行列添字 `i,j`）は core 領域上の DOF に対応している
- [ ] 積分領域は core のみである
	- `S_ij = ∫_core φ_i φ_j`
	- `H_ij = ∫_core φ_i (Hφ_j)`
- [ ] `Hφ_j` の評価では buffer を参照してよいが、最終積分点は core のみである
- [ ] 近傍値が不足する場合は「フラグメント周期」を使って補う（全系周期ではない）
- [ ] 非局所 PP は self-only 契約を守り、`⟨φ_i|Vnl|φ_j⟩` を core 弱形式で閉じる
- [ ] 固有値問題は generalized EVP `Hc = εSc` として解く（`S` を省略しない）

補足（実装判断）:

- buffer 上で固有関数でも、core 積分で定義した `H,S` の固有関数とは限らない
- よって最終的に使う定義が core 弱形式なら、固有状態もその `H_core,S_core` で求める
- startup で再対角化を行う場合は、global state の破壊的再ラベリングを禁止し、既存 state set の対応を維持する

推奨フロー（安全側）:

1. DC で buffered データ生成
2. RT/DG で `S_core,H_core` を構築
3. `H_core c = ε S_core c` を解く
4. 得られた状態を初期状態に使う

## 1. 前提固定

### 1.1 比較対象フラグメント
- [ ] 比較対象 `ifrag` を 1 つ固定した
- [ ] `ispin` を固定した
- [ ] まずは self-only 条件で比較する方針にした

### 1.2 比較対象波動関数
- [ ] LCFO 後の係数ではなく、まずはフラグメント基底関数 `|lambda>` そのもの、または 1 本の代表 `ψ(r)` を使う
- [ ] 同じ `ψ(r)` を DC 側と DG 側で同一 box 上に置く
- [ ] 必要なら `cg_orbitals_buffered.bin` / `basis_functions_buffered.bin` のどちらを使うか明記した

### 1.3 比較対象作用素
- [ ] まず canonical fragment Hamiltonian を `H_frag = T + Vloc + Vnl` と定義した
- [ ] DG 面項、jump/average、時間依存項 `A(t)`、`A^2/2` は最初の比較から除外した
- [ ] `Vnl` は self-only projector / self-only cache のどちらを正とするか固定した

## 2. 実空間 box と座標の一致

### 2.1 box metadata
- [ ] `nxyz_domain`
- [ ] `nxyz_buffer`
- [ ] `nxyz_box`
- [ ] `lgnum_total`
- [ ] `jxyz_tot`
を DC 側と DG 側で並べて確認した

### 2.2 座標原点と並び順
- [ ] buffered box の並び順が DC と DG で一致していることを確認した
- [ ] `core, pos_buf, neg_buf` と `neg_buf, core, pos_buf` の再配置規約を明記した
- [ ] 実際の `ix,iy,iz -> x,y,z` 対応を 2,3 点サンプルして一致確認した

### 2.3 周期条件
- [ ] フラグメント自身を周期境界として扱うことを前提にした
- [ ] DC 側も DG 側も同じ periodic wrap 規約で比較している
- [ ] 全系周期とフラグメント周期が混ざっていない

## 3. 入力波動関数 `ψ(r)` の一致

### 3.1 生波動関数の一致
- [ ] DC 側 `ψ(r)` を buffered box 上でダンプした
- [ ] DG 側 `ψ(r)` を同じ box 上へ再配置した
- [ ] `max|ψ_dc(r)-ψ_dg(r)|` を確認した
- [ ] `sum |ψ|^2 hvol` が双方で一致することを確認した

### 3.2 基底関数の場合
- [ ] `basis_functions_buffered.bin` の代表 basis について z profile / xy-integrated profile を比較した
- [ ] `|lambda>` の符号差だけでなく振幅差がないことを確認した

## 4. `Tψ(r)` の一致

### 4.1 微分作用
- [ ] DC 側の kinetic 作用 `Tψ(r)` を取得した
- [ ] DG 側の kinetic 作用 `Tψ(r)` を取得した
- [ ] stencil / wrap / halo 依存を切った self-only 条件で比較した

### 4.2 判定
- [ ] `max|Tψ_dc - Tψ_dg|` を評価した
- [ ] `L2` ノルム差を評価した
- [ ] 境界近傍と内部点を分けて誤差を見た

## 5. `Vloc ψ(r)` の一致

### 5.1 local potential の一致
- [ ] DC 側で使った `Vloc(r)` を固定した
- [ ] DG 側で使った `Vloc(r)` を固定した
- [ ] same box 上で `Vloc_dc(r)` と `Vloc_dg(r)` を比較した

### 5.2 作用結果
- [ ] `Vloc(r) ψ(r)` を DC 側と DG 側で比較した
- [ ] `max|Vlocψ_dc - Vlocψ_dg|` を評価した

## 6. `Vnl ψ(r)` の一致

### 6.1 projector 定義
- [ ] projector `β_a(r)` の定義域が同じであることを確認した
- [ ] buffered/self-only 契約を明示した
- [ ] halo 混合 projector を canonical 比較から外した

### 6.2 実空間作用
- [ ] DC 側で `Vnl ψ(r)` を取得した
- [ ] DG 側で `Vnl ψ(r)` を取得した
- [ ] `uVphi_self`, `proj_global`, `proj_weight`, `contrib`, `y` の経路が canonical 定義と一致することを別途確認した

### 6.3 判定
- [ ] `max|Vnlψ_dc - Vnlψ_dg|` を評価した
- [ ] `L2` ノルム差を評価した

## 7. 合成 `Hψ(r)` の一致

### 7.1 段階比較
- [ ] `Tψ`
- [ ] `Vlocψ`
- [ ] `Vnlψ`
を個別に比較した

### 7.2 合成比較
- [ ] `Hψ = Tψ + Vlocψ + Vnlψ` を DC/DG 双方で合成した
- [ ] `max|Hψ_dc - Hψ_dg|` を評価した
- [ ] `L2` ノルム差を評価した

### 7.3 点ごとの可視化
- [ ] 代表断面（x, y, z のいずれか）で `Hψ(r)` を可視化した
- [ ] core 内部と buffer 境界近傍で差の局在を確認した

## 8. 比較条件の固定

### 8.1 時間依存項を切る
- [ ] `A=0`
- [ ] `m` 項なし
- [ ] `A^2/2` なし
- [ ] 固定 `Vh/Vxc`
でまず比較した

### 8.2 DG 固有項を切る
- [ ] 面項なし比較を先に行った
- [ ] canonical self Hamiltonian と DG coupling term を分けて扱った

## 9. 合格基準

### 9.1 最低限
- [ ] `ψ(r)` 自体が一致
- [ ] `Tψ`, `Vlocψ`, `Vnlψ` の各段で差分源が特定できている
- [ ] `Hψ` 合成差が単一項の差分へ還元できる

### 9.2 基底行列比較へ進んでよい条件
- [ ] `Hψ(r)` の差が許容誤差内
- [ ] 周期条件・box 並び順・projector 契約が明文化済み
- [ ] DC 側と DG 側で「同じ実空間作用素を見ている」と言える

## 10. 実装メモ
- DC 側では `hamiltonian_local.bin` を書いているが、RT はそれを直接読まず、基底から `H_mat` / `H_nl_cache` を再構成する。
- したがって、まず比較すべきは行列ではなく、同じ `ψ(r)` に対する実空間作用 `Hψ(r)` である。
- `basis_functions_buffered.bin` は「同じ lambda 基底を、support だけ拡張したもの」という契約である。

## 11. 関連資料
- [DG_BUFFERED_RT_SEED_CONTRACT_SPEC_JA.md](/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/DG_BUFFERED_RT_SEED_CONTRACT_SPEC_JA.md)
- [DG_BUFFERED_DC_RT_AUDIT_CHECKLIST_JA.md](/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/DG_BUFFERED_DC_RT_AUDIT_CHECKLIST_JA.md)

## 12. 2026-04-23 監査進捗（実測）

対象:

- 実行ディレクトリ: `.worktrees/codex/dg-phi-box-cache-phase-a`
- データ: `data_dcdft -> data_dcdft_step1_new`
- 主ログ:
	- `run_diag_nokick_nt1_strict_read.log`
	- `run_diag_nokick_nt1_seedcheck.log`
	- `run_diag_nokick_nt1_proj0_trace.log`
	- `run_diag_nokick_nt1_proj1_trace.log`
	- `run_diag_nokick_nt1_energy_probe.log`
	- `run_nokick_nt2_propim_trace.log`
	- `run_kick_nt2_force_i_mterm.log`

今回埋めた項目（実測あり）:

- 2.1 box metadata（DG側）
	- [x] `nxyz_domain = 12 x 12 x 12` を確認
	- [x] `nxyz_buffer = 4 4 4`（`Buffered basis half-width`）を確認
	- [x] `nxyz_box = 20 20 20`（`Buffered basis box size`）を確認
	- [x] `nproc_rgrid = 2 2 2`、`Number of fragments = 8` を確認（`run_diag_nokick_nt1_energy_probe.log`）
	- [ ] `lgnum_total = 24 x 24 x 24` は上記から推定可能だが、DC 側との明示並列表は未作成
	- [ ] `lgnum_total` / `jxyz_tot` の DC-DG 並列表は未作成

- 8.2 DG 固有項を切る（条件確認）
	- [x] buffered 基準で halo exchange 無効（self寄り条件）を確認
	- [ ] 面項なし比較として DC 側 `Hψ` と 1:1 比較は未実施

- 8.1 時間依存項を切る（部分）
	- [x] no-kick で `prop-im-trace: m(re2/im2)=0` を確認
	- [x] kick で `prop-im-trace: m(re2/im2)` の虚部チャネルが立つことを確認
	- [ ] `A=0`, `A^2/2` 完全固定での DC-DG `Hψ` 比較は未実施
	- [ ] 固定 `Vh/Vxc` 条件での `Hψ` 直接比較は未実施

- 関連 sanity（前段整合）
	- [x] `SEED-CONSISTENCY-OK` を確認（metadata 整合）
	- [x] fixed-H 残差の初期値を取得（`DG-FIXED-H`）

- 4/5 節に向けた追加実測（energy probe, no-kick nt=1）
	- [x] `energy-global compare` を取得
	- [x] `kin_mat = kin_rs = 9.677861E+02`（itt=1）を確認
	- [ ] `static_mat` と `static_rs`、`vloc_mat` と `vloc_rs` に大きな差（itt=1）が残存
	- [ ] ただし現状は空間積分量の比較のみで、`max|Δ|` / `L2(Δ)` は未取得

- 4/5/6/7 節に向けた追加実装（DG 側の実空間成分ダンプ経路）
	- [x] `compute_realspace_energy_probe` に代表状態の成分ダンプを追加
	- [x] 環境変数で `ifrag/ispin/io/itt` を固定できるようにした
	- [x] 出力成分は `psi, Tpsi, Vlocpsi, Vnlpsi, Hpsi(=H0psi+Vnlpsi)`
	- [x] rank 単位で同一フォーマットの実空間ダンプを出力
	- [x] `maxabs` / `L2` をローカル指標としてログ出力

実装で使う環境変数（DG 側）:

- `SALMON_DG_RS_COMPONENT_DUMP=1`
- `SALMON_DG_RS_COMPONENT_DUMP_IFRAG=<ifrag>`
- `SALMON_DG_RS_COMPONENT_DUMP_ISPIN=<ispin>`
- `SALMON_DG_RS_COMPONENT_DUMP_IO=<state index>`
- `SALMON_DG_RS_COMPONENT_DUMP_ITT=<itt>`（`-1` で毎ステップ）
- `SALMON_DG_RS_COMPONENT_DUMP_PREFIX=<prefix>`（任意）

検証 run（nt=1, no-kick）:

- [x] run log: `run_diag_nokick_nt1_rs_component.log`
- [x] 成分ファイル: `dg_rs_component_itt1_if1_sp1_io1_rank0.dat`
- [x] ログ例（`ifrag=1, ispin=1, io=1, itt=1`）
	- `maxabs-local`: `T=2.900669E-01`, `Vloc=9.111734E-01`, `Vnl=5.936859E-01`, `H=4.598584E-01`
	- `l2-local`: `T=8.081856E-01`, `Vloc=2.926811E+00`, `Vnl=1.134314E+00`, `H=2.757164E+00`

- 4/5/6/7 節に向けた追加実装（DC 側の実空間成分ダンプ経路）
	- [x] `src/gs/dc/lcfo.f90` の `hpsi_basis` に環境変数ゲート付きダンプを追加
	- [x] 出力成分は `psi, Tpsi, Vlocpsi, Vnlpsi, Hpsi`
	- [x] `Tpsi` は `hpsi(..., ttpsi)` の optional kinetic 出力から取得
	- [x] `Vnlpsi = Hpsi - Tpsi - Vlocpsi` で再構成
	- [x] 選択 1 状態のみを 3 次元配列で `comm_summation` する形へ修正し、`comm_summation` generic 解決エラーを回避
	- [x] `cmake --build build -j4 --target salmon` でコンパイル通過
	- [x] checker 整合入力で `dc_rs_component_if1_sp1_io1_fragroot.dat` の run-time 生成を確認

実装で使う環境変数（DC 側）:

- `SALMON_DC_RS_COMPONENT_DUMP=1`
- `SALMON_DC_RS_COMPONENT_DUMP_IFRAG=<ifrag>`
- `SALMON_DC_RS_COMPONENT_DUMP_ISPIN=<ispin>`
- `SALMON_DC_RS_COMPONENT_DUMP_IO=<state index>`
- `SALMON_DC_RS_COMPONENT_DUMP_PREFIX=<prefix>`

実行確認メモ:

- 旧小型 DC-LCFO 入力 (`samples/exercise_dg_rt_hse_test/H2/inputfile_dc_lcfo_gs`) は `yn_dc_lcfo` / `nstate_frag` を `&parallel` に置いており、現行 input checker では `&dc` 不足として `1 error(s) in input.` で停止した
- `&dc` に `yn_dc_lcfo`, `yn_dc_lcfo_diag`, `num_fragment`, `nproc_rgrid_tot`, `nstate_frag` を置き、`dl` を外した checker 整合入力では DC ダンプ経路が起動した
- `samples/exercise_dg_rt_hse_test/H2` 直下からの実行は `hubbard_param.dat` 自動検出で無関係の +U 初期化例外へ入るため、検証 run は `/tmp` から実行した
- `/tmp/run_dc_lcfo_rs_component.log` で以下を確認した
	- `[DC-RS-COMP] config prefix=dc_rs_component ifrag=1 ispin=1 io=1`
	- `[DC-RS-COMP] select-local: ifrag=1 ispin=1 io=1 norm=1.000000E+00`
	- `end SALMON`

追加進捗（C64/DC-DG 直接比較の前段）:

- [x] C64 条件 (`.namelist_dc_gs.tmp` 派生、`nscf=3`) で DC 側ダンプを MPI 8 並列実行し、`dc_rs_component_if1_sp1_io1_fragroot.dat` を生成
- [x] MPI 8 実行で発生していた `SIGTRAP` は、`src/gs/dc/lcfo.f90` の DC ダンプ部で rank 選択条件の内側に collective (`comm_summation`) が置かれていたことが原因
- [x] 修正: 全 rank が collective に参加し、対象 fragment のみ非ゼロ寄与する形へ変更して再ビルド
- [x] 修正後ログで `[DC-RS-COMP] config/select-local/maxabs-local/l2-local` と `end SALMON` を確認

DC-DG 差分（ifrag=1, ispin=1, io=1, 実測）:

- 比較ファイル
	- DG: `.worktrees/codex/dg-phi-box-cache-phase-a/dg_rs_component_itt1_if1_sp1_io1_rank0.dat`
	- DC: `.worktrees/codex/dg-phi-box-cache-phase-a/dc_rs_component_if1_sp1_io1_fragroot.dat`
- 格子範囲
	- DG: `1..12` 各軸（1728 点）
	- DC: `1..14` 各軸（2744 点）
	- 重なり領域 `1..12` 各軸（1728 点）で差分評価
- 差分指標（重なり領域）
	- `psi`: `maxabs=5.963512100E-01`, `L2=3.150199980E+00`
	- `Tpsi`: `maxabs=2.900668513E-01`, `L2=1.928541916E+00`
	- `Vlocpsi`: `maxabs=2.029547375E+00`, `L2=7.573884007E+00`
	- `Vnlpsi`: `maxabs=5.052786068E+00`, `L2=1.660556887E+01`
	- `Hpsi`: `maxabs=5.479693128E+00`, `L2=1.597034920E+01`

未完の主要ブロッカー:

- 本チェックリストの中核差分（`max/L2`）は重なり領域で取得済み。ただし DC 側 box が `14^3`、DG 側 box が `12^3` で完全同一 box ではないため、最終判定には box 契約を揃えた再比較が必要。

次アクション（このチェックリストの続き）:

1. DC 側の dump box を `12^3` に合わせるか、DG 側を `14^3` へ再配置して完全同一 box 比較へ移行
2. 同一 box で `max|Δ|` と `L2(Δ)` を再計測し、今回の重なり領域結果と差を評価
3. 必要に応じて断面プロファイル（x/y/z）を追加して差の局在を確定

## 13. 2026-04-24 Poisson/BC 経路差分の切り分け（今回の到達点）

確認済み（コード経路）:

- [x] DG-RT 側の Hartree 更新は `hartree_dg_distributed` を直接呼ぶ
	- `src/rt/dg/rt_dg_density_hamiltonian_update.f90`
	- `src/rt/dg/rt_dg_integrator_stage_update.f90`
- [x] `hartree_dg_distributed` は `poisson_ft/poisson_ffte`（または HSE-SR 版）へ直行し、`hartree_sub::hartree` の `iperiodic/method_poisson` 分岐を通らない
	- `src/poisson/poisson_dg_distributed.f90`
- [x] DC/GS 側は `hartree_sub::hartree` を経由し、`iperiodic` と `method_poisson` による分岐を持つ
	- `src/poisson/hartree.f90`

今回の入力条件（確認）:

- [x] DG RT 入力 / DC GS 入力ともに `yn_periodic='y'`（= `iperiodic=3`）
- [x] どちらの入力にも `yn_ffte`, `yn_fftw`, `method_poisson` の明示指定なし
- [x] 既定値は `yn_ffte='n'`, `yn_fftw='n'`, `method_poisson='cg'`（ただし `iperiodic=3` では `method_poisson` は実質未使用）
	- `src/io/inputoutput.f90`
- [x] 実行時ログ (`variables.log`) でも `yn_ffte=n`, `yn_fftw=n`, `yn_put_wall_z_boundary=n`, `method_poisson=cg` を確認

解釈（現時点）:

- [x] 「rho を FFT 後に高周期だけ別途除去する処理」が DC と DG で二重実装されている証拠は未検出
- [x] 乖離主因候補は、FFTフィルタ差分よりも Poisson 実行経路差（分岐入口、境界項、拘束条件、box 契約）

次に埋めるチェック（Poisson/BC 4軸）:

- [x] 境界条件: wall 付加 (`yn_put_wall_z_boundary`) は現行比較条件で `n`（variables.log）
- [x] 平均値拘束の観測量として出力 Vh の平均値差を定量化
	- DC vs DG(Vh isolation from DC rho): `N=8000`, `relL2=1.370724956`, `centered_relL2=1.232155264`, `corr=0.035101658`, `mean_offset(dg-dc)=+0.721727776`
- [x] box 定義: DC 側 Vh ダンプと DG isolation 出力の index 範囲は一致（ともに `ix,iy,iz=1..20`, 8000点）
- [x] 単純シフト否定: z-shift 全探索でも `best corr=0.26451 (shift=16)`、限定 xyz 探索でも `best corr=0.29254 (3,3,16)`
- [ ] FFT/solver 条件: `yn_ffte/yn_fftw` と HSE-SR Hartree (`SALMON_HSE_SR_HARTREE`) を両側で固定して A/B 実験

次の実験セット（最短）:

1. `SALMON_HSE_SR_HARTREE=0/1` を DC/DG 両側で同期させた A/B を実施し、上記 5 指標を再計測
2. `yn_ffte=n/y` を両側で固定して A/B を実施（FFTE ビルド有効時）
3. `run_dc_rs_component_*.log` と `run_vh_iso*.log` の同一条件ペアごとに、指標を1表へ集約

13.1 HSE-SR Hartree A/B 実測（2026-04-24）:

- 条件:
	- 同一入力: `.namelist_dg_rt_c64_nt1_rsdump.tmp`
	- 同一 DC rho 注入: `dc_rs_component_rho_if1_sp1_fragroot.dat`
	- 比較先: `dc_rs_component_vloc_if1_sp1_fragroot.dat` の `Vh` 列
	- 切替: `SALMON_HSE_SR_HARTREE=0/1`
- 実行結果:
	- [x] `flag=0`: run 完走 (`ec=0`), isolation 出力生成
	- [x] `flag=1`: run 完走 (`ec=0`), isolation 出力生成
- 指標比較（DC Vh との差）:
	- `flag=0`: `relL2=1.370724956`, `centered_relL2=1.232155264`, `corr=0.035101658`, `mean_offset=+0.721727776`
	- `flag=1`: `relL2=1.364188535`, `centered_relL2=1.226919882`, `corr=0.037509153`, `mean_offset=+0.716705408`
- 解釈:
	- [x] HSE-SR Hartree を ON にしても改善は軽微（相関・誤差とも依然として大きく不一致）
	- [x] 主因は HSE-SR の有無単独では説明しにくく、Poisson/BC の他要素（拘束・境界・定義差）を優先

13.2 Poisson 契約トレース追加（2026-04-24）:

- 目的:
	- [x] DC 経路 / DG 経路で同一フォーマットの `rho in -> solver(kernel) -> Vh out` を段別に取得し、gauge/拘束/定義域差を直接比較する
- 有効化:
	- [x] 環境変数 `SALMON_POISSON_CONTRACT_TRACE=1`
- 追加した計測点:
	- [x] `src/poisson/hartree.f90`
		- `DC-HARTREE-rho-in`
		- `DC-HARTREE-vh-out`
	- [x] `src/poisson/poisson_dg_distributed.f90`
		- `DG-HARTREE-rho-in`
		- `DG-HARTREE-vh-out-prewall`
		- `DG-HARTREE-vh-out-postwall`
	- [x] `src/poisson/poisson_periodic.f90`
		- `POISSON-FT-IN`, `POISSON-FT-KERNEL`, `POISSON-FT-OUT`
		- `POISSON-FFTE-IN`, `POISSON-FFTE-KERNEL`, `POISSON-FFTE-OUT`
- 出力される主指標:
	- [x] 実空間スカラー統計: `sum`, `mean`, `l2`, `min`, `max`, `lg_num`, `mg_is`
	- [x] カーネル統計: `g0_coef`, `g0_rhoG(re,im)`, `g0_kernel_out(re,im)`, `||rhoG||_2`, `||kernel_out||_2`
- 実行時の確認キーワード:
	- [x] `\[POISSON-CONTRACT\]`
	- [x] `POISSON-FT-KERNEL|POISSON-FFTE-KERNEL|DC-HARTREE|DG-HARTREE`
- 次アクション:
	1. 同一 `dc_rs_component_rho_if1_sp1_fragroot.dat` 注入条件で DC / DG を別々に実行し、`[POISSON-CONTRACT]` 行だけ抽出して 1 表へ並べる
	2. `mean(Vh_out)-mean(rho_in)` と `g0_kernel_out` の DC-DG 差を first-pass 指標にする
	3. その後 `yn_ffte=n/y` A/B でも同じ表を再計測し、solver 実装依存か契約依存かを切り分ける

13.3 Poisson 契約トレース実測（HSE-SR A/B, 2026-04-24）:

- 実行条件:
	- 実行バイナリ: `.worktrees/codex/dg-phi-box-cache-phase-a/build_mpi_enabled_relaxed/salmon`
	- 入力: `.namelist_dg_rt_c64_nt1_rsdump.tmp`
	- 共通 env:
		- `SALMON_POISSON_CONTRACT_TRACE=1`
		- `SALMON_DG_VH_ISOLATION_FROM_DC_RHO=1`
		- `SALMON_DG_VH_ISOLATION_IFRAG=1`
		- `SALMON_DG_VH_ISOLATION_ITT=1`
		- `SALMON_DG_ENFORCE_STATE_TRUNCATION=0`（互換運用）
	- 切替: `SALMON_HSE_SR_HARTREE=0/1`
- ログ:
	- `run_vh_iso_contract_hsesr_0_newbin.log`
	- `run_vh_iso_contract_hsesr_1_newbin.log`
	- いずれも `end SALMON` まで完走

- first-pass 比較（タグ初回出現値）:
	- DC-HARTREE-rho-in (mean):
		- `flag=0`: `1.054489256019869E-01`
		- `flag=1`: `1.054489256019869E-01`
	- DC-HARTREE-vh-out (mean):
		- `flag=0`: `-6.705669969914791E-16`
		- `flag=1`: `-6.397917175472980E-16`
	- DG-HARTREE-rho-in (mean):
		- `flag=0`: `1.054489256019868E-01`
		- `flag=1`: `1.054489256019868E-01`
	- DG-HARTREE-vh-out-prewall (mean):
		- `flag=0`: `-9.621932880084689E-16`
		- `flag=1`: `-9.868649107779170E-16`
	- POISSON-FT-KERNEL (`g0_kernel_out(re,im)`):
		- `flag=0`: `(0, 0)`
		- `flag=1`: `(0, 0)`

- 追加所見:
	- `POISSON-FT-KERNEL` の `||rhoG||_2` で `flag=1` に異常スパイクを観測
		- `max ||rhoG||_2 (flag=0) = 8.249453689432553E-01`
		- `max ||rhoG||_2 (flag=1) = 1.722965298945772E+79`
	- `g0_kernel_out` は両ケースともゼロのままのため、異常は G=0 ではなく高G側またはノルム評価経路の破綻が疑わしい

- 判断:
	- [x] Poisson 契約トレース機構は実測で有効化できた
	- [x] DC/DG の `rho-in` と `Vh-out` の mean は first-pass では同程度（gauge 平均値差は小さい）
	- [ ] `SALMON_HSE_SR_HARTREE=1` での `||rhoG||_2` 異常値は追加切り分けが必要（数値/実装バグ候補）

13.4 監査方針の確定（2026-04-24）:

- 結論整理:
	- [x] 本件の主要不一致（DC vs DG の `Vh` 不一致）は `HSE-SR=0/1` 切替でほぼ残存
	- [x] A/B 差は軽微であり、現時点で主因を HSE-SR 分岐へ帰属しない
	- [x] `HSE-SR=1` で観測した `POISSON-FT-KERNEL` の巨大スパイク（過去実測 `~1.7E+79`）は、独立バグ疑いとして分離して扱う

- 本件（本筋）の扱い:
	- [x] 以降の Poisson/境界条件/定義差の比較・監査は `SALMON_HSE_SR_HARTREE=0` を基準固定で実施する
	- [x] 比較表・監査指標（`mean/l2/g0`）は `HSE-SR=0` 系列を正規トラックとする

- 別件（独立バグ）の扱い:
	- [x] `HSE-SR=1` 高G/ノルム異常は separate issue として切り分ける
	- [ ] `max|rhoG|`・異常インデックス・近傍成分の追跡を別チケットで継続する

- 実務ルール（当面）:
	1. 本件の結論形成に `HSE-SR=1` の異常値を混ぜない
	2. 本件の再現・比較 run は `HSE-SR=0` 固定で統一する
	3. `HSE-SR=1` 経路の修正が入った場合のみ、回帰確認として再度 A/B を実施する

13.5 比較・監査の継続更新（HSE-SR=0 固定トラック, 2026-04-24）:

- 実施内容:
	- [x] 強化トレース版バイナリで `HSE-SR=0/1` を同条件再実行
	- [x] ログ: `run_vh_iso_contract_hsesr_0_diag.log`, `run_vh_iso_contract_hsesr_1_diag.log`
	- [x] 本件の比較判断は `HSE-SR=0` を主系列として継続する方針を維持

- first-pass 指標（初回出現、A/B比較）:
	- `DC-HARTREE-rho-in mean`
		- `=0`: `1.054489256019869E-01`
		- `=1`: `1.054489256019869E-01`
	- `DC-HARTREE-vh-out mean`
		- `=0`: `-6.705669969914791E-16`
		- `=1`: `-6.397917175472980E-16`
	- `DG-HARTREE-rho-in mean`
		- `=0`: `1.054489256019868E-01`
		- `=1`: `1.054489256019868E-01`
	- `DG-HARTREE-vh-out-prewall mean`
		- `=0`: `-9.621932880084689E-16`
		- `=1`: `-9.868649107779170E-16`
	- `POISSON-FT-KERNEL`
		- `max|rhoG|`: `=0/1` とも `1.054489256019869E-01`
		- `bad_rhoG=0`, `bad_kernel_out=0`（両ケース）

- 監査上の扱い:
	- [x] A/B 差は軽微で、主系列（DC vs DG 不一致）の主因を HSE-SR 分岐へ帰属しない
	- [x] 以降の比較表・監査結論は `HSE-SR=0` 固定結果を正規採用
	- [ ] `HSE-SR=1` 経路は separate bug track として、必要時に独立再現ログを更新

13.6 A/B 差分量の明示（first-pass, 2026-04-24 更新）:

- 対象ログ:
	- `run_vh_iso_contract_hsesr_0_diag.log`
	- `run_vh_iso_contract_hsesr_1_diag.log`

- first-pass 指標の差分（`delta = value(flag=1) - value(flag=0)`）:
	- `DC-HARTREE-rho-in mean`
		- `=0`: `1.054489256019869E-01`
		- `=1`: `1.054489256019869E-01`
		- `delta`: `+0.000E+00`
	- `DC-HARTREE-vh-out mean`
		- `=0`: `-6.705669969914791E-16`
		- `=1`: `-6.397917175472980E-16`
		- `delta`: `+3.078E-17`
	- `DG-HARTREE-rho-in mean`
		- `=0`: `1.054489256019868E-01`
		- `=1`: `1.054489256019868E-01`
		- `delta`: `+0.000E+00`
	- `DG-HARTREE-vh-out-prewall mean`
		- `=0`: `-9.621932880084689E-16`
		- `=1`: `-9.868649107779170E-16`
		- `delta`: `-2.467E-17`
	- `POISSON-FT-KERNEL max|rhoG|`
		- `=0`: `1.054489256019869E-01`
		- `=1`: `1.054489256019869E-01`
		- `delta`: `+0.000E+00`

- 更新判断:
	- [x] A/B 差分は first-pass で数値的に微小
	- [x] 本件主系列は `HSE-SR=0` 固定で継続して妥当
	- [ ] 次段では `HSE-SR=0` 固定のまま、DC/DG の `vh-out` 分布差（offset/shape）を case 横断で表形式化する

13.7 本体監査の進行方針（確定, 2026-04-24）:

- 方針:
	- [x] 以降は迷わず `HSE-SR=0` を基準に、Poisson/境界条件/定義差の本体監査を進める
	- [x] `HSE-SR=1` は separate bug track として独立管理し、本体の結論形成に混在させない

- 本体監査の優先順（`HSE-SR=0` 固定）:
	1. [x] Poisson 出力 `vh-out` の差を `offset`（平均差）と `shape`（centered 相関/centered relL2）へ分離して表化
	2. [x] 境界条件契約（wall/no-wall, periodic wrap, g0拘束）の DC/DG 対応を 1 表に整理
	3. [x] 定義差契約（box index, local/global 座標、注入 `rho` の対応）を再確認し、指標差に対する寄与仮説を更新

- 運用ルール:
	- [x] 比較 run は `HSE-SR=0` を既定値として統一
	- [ ] `HSE-SR=1` の再測定は、独立バグ側で修正が入った時点の回帰確認に限定

13.8 本体監査 3項目の実測表（`HSE-SR=0`, 2026-04-24）:

- 対象データ/ログ:
	- `dc_rs_component_vloc_if1_sp1_fragroot.dat`
	- `dc_rs_component_rho_if1_sp1_fragroot.dat`
	- `dg_vh_from_dc_rho_itt1_if1_rank0.dat`
	- `run_vh_iso_contract_hsesr_0_diag.log`

1) Poisson `vh-out` 差の `offset/shape` 分離（DC `Vh` vs DG `vh_from_dc_rho`）:

| 項目 | 値 |
|---|---:|
| N (比較点数) | 8000 |
| relL2 | 1.364188535 |
| centered_relL2 | 1.226919882 |
| centered corr | 0.037509153 |
| mean offset (DG-DC) | +0.716705408 |
| max abs(DG-DC) | 4.709643676 |

- 判定:
	- `offset`（平均差）は大きい（`+0.7167`）
	- `shape` も不一致が支配的（`centered_relL2=1.2269`, `corr=0.0375`）
	- よって「単純な定数シフトだけでは説明不可」。

2) 境界条件契約（wall/no-wall, periodic wrap, g0拘束）DC/DG 対応:

| 契約項目 | DC | DG | 判定 |
|---|---|---|---|
| wall 境界 (`yn_put_wall_z_boundary`) | `n`（variables.log既報） | `n`（同条件） | 一致 |
| periodic Poisson 経路 | `POISSON-FT-KERNEL` 呼出あり | `POISSON-FT-KERNEL` 呼出あり | 一致 |
| `g0_coef` | `0.0` | `0.0` | 一致 |
| `g0_kernel_out(re,im)` | `(0,0)` | `(0,0)` | 一致 |
| `DG prewall/postwall` 差 | - | `vh-out-prewall == postwall`（mean/l2/min/max 同値） | wall非寄与を確認 |

- first-pass 根拠行（`run_vh_iso_contract_hsesr_0_diag.log`）:
	- `DC-HARTREE-rho-in`: line 61
	- `POISSON-FT-KERNEL`: line 63
	- `DC-HARTREE-vh-out`: line 67
	- `DG-HARTREE-rho-in`: line 152
	- `DG-HARTREE-vh-out-prewall/postwall`: lines 159/160

3) 定義差契約（box index / local-global / 注入 `rho`）再確認と寄与仮説更新:

| チェック項目 | 実測 | 判定 |
|---|---|---|
| index 範囲 (DC rho) | `ix,iy,iz = 1..20` | DG と一致 |
| index 範囲 (DG vh isolation) | `ix,iy,iz = 1..20` | DC と一致 |
| key 交差数 | `8000 / 8000` | 1:1 対応 |
| 注入 `rho` の比較元 | `dc_rs_component_rho_if1_sp1_fragroot.dat` | 固定済み |

- 仮説更新（寄与順）:
	1. 高優先: Poisson solver 呼出前後の `rho` 分布差（DC-HARTREE-rho-in と DG-HARTREE-rho-in の `l2/min/max` 差）が shape 不一致へ寄与。
	2. 中優先: 離散化/分解（DG 経路の再配置・分割）由来の周波数成分差（`||rhoG||_2`, `||kernel_out||_2` 差）。
	3. 低優先: wall/g0 拘束の扱い差（今回データでは一致し、主因候補から後退）。

13.9 Poisson 段別分解（1本目, `HSE-SR=0` 固定, 2026-04-24）:

- 目的:
	- `rho_in -> rhoG -> vhG -> Vh_raw -> Vh_final` の 5 段で、差が最初に立つ段を特定する

- 実装追加（今回）:
	- `src/poisson/poisson_periodic.f90`
		- `[POISSON-CONTRACT] POISSON-FT-RHOG`
		- `[POISSON-CONTRACT] POISSON-FT-VHG`
	- `src/poisson/hartree.f90`
		- `[POISSON-CONTRACT] DC-HARTREE-vh-out-postwall`

- 実行条件:
	- 入力: `.namelist_dg_rt_c64_nt1_rsdump.tmp`
	- 注入 rho: `dc_rs_component_rho_if1_sp1_fragroot.dat`（既存ファイル再利用）
	- env: `SALMON_HSE_SR_HARTREE=0`, `SALMON_POISSON_CONTRACT_TRACE=1`, `SALMON_DG_VH_ISOLATION_FROM_DC_RHO=1`
	- ログ: `run_vh_iso_contract_hsesr_0_stage5.log`

- first-pass 段別比較（DC vs DG）:

| 段 | 指標 | DC | DG | DG-DC |
|---|---|---:|---:|---:|
| `rho_in` | mean | 1.0544892560E-01 | 1.0544892560E-01 | -9.714E-17 |
| `rho_in` | l2 | 1.4414986521E+01 | 1.4694742588E+01 | +2.798E-01 |
| `rho_in` | max | 3.6063094903E-01 | 4.9121277234E-01 | +1.306E-01 |
| `rhoG` | mean_abs | 1.5636481718E-04 | 2.3097865225E-04 | +7.461E-05 |
| `rhoG` | l2 | 1.1089428281E-01 | 1.1208890380E-01 | +1.195E-03 |
| `rhoG` | max_abs | 1.0544892560E-01 | 1.0544892560E-01 | -9.714E-17 |
| `vhG` | mean_abs | 3.6623993688E-04 | 6.4568172894E-04 | +2.794E-04 |
| `vhG` | l2 | 3.9784635441E-01 | 4.8903721078E-01 | +9.119E-02 |
| `vhG` | max_abs | 1.8884892079E-01 | 2.2125337094E-01 | +3.240E-02 |
| `Vh_raw` | mean | -6.7056699699E-16 | -9.6219328801E-16 | -2.916E-16 |
| `Vh_raw` | l2 | 6.0843536823E+01 | 7.2823500867E+01 | +1.198E+01 |
| `Vh_raw` | max | 1.9560891714E+00 | 2.6998830005E+00 | +7.438E-01 |
| `Vh_final` | mean | -6.7056699699E-16 | -9.6219328801E-16 | -2.916E-16 |
| `Vh_final` | l2 | 6.0843536823E+01 | 7.2823500867E+01 | +1.198E+01 |
| `Vh_final` | max | 1.9560891714E+00 | 2.6998830005E+00 | +7.438E-01 |

- 判定（成功条件）:
	- [x] 差が最初に立つ段は `rho_in`（`l2`, `max` で既に差）
	- [x] `rhoG` 以降でも差は継続・増幅
	- [x] `Vh_raw` と `Vh_final` は同値（wall off 条件のため）

- 補足:
	- stage run の最終終了コードは `138` だが、first-pass の段別指標はログに採取済みで本判定には十分。
	- `relL2/corr` は実空間 `Vh`（`dc_rs_component_vloc` vs `dg_vh_from_dc_rho`）では既報、`rhoG/vhG` は現状 aggregate 指標比較（full field relL2/corr は別途 dump 拡張で対応）。

13.10 rho 注入経路の段別監査（Poisson停止, `HSE-SR=0` 固定, 2026-04-24）:

- 目的:
	- `dc_rs_component_rho_if1_sp1_fragroot.dat -> DG rho_in` のみを追跡し、差が
		1) ファイル読込直後からあるか
		2) DG 再配置で初めて入るか
	  を特定する

- 実装追加（今回）:
	- `src/rt/dg/rt_dg_density_hamiltonian_update.f90`
	- `src/rt/dg/rt_dg_integrator_stage_update.f90`
	- 追加ログタグ: `[RHO-IN-AUDIT]`
		- `stage=dc-file-raw`
		- `stage=dg-injected`
		- `stage=dg-rho-in-pre-hartree`
		- `topdiff rank=1..10`

- 実行条件:
	- 入力: `.namelist_dg_rt_c64_nt1_rsdump.tmp`
	- env: `SALMON_HSE_SR_HARTREE=0`, `SALMON_DG_VH_ISOLATION_FROM_DC_RHO=1`, `SALMON_DG_VH_ISOLATION_IFRAG=1`, `SALMON_DG_VH_ISOLATION_ITT=1`
	- 注入元: `dc_rs_component_rho_if1_sp1_fragroot.dat`
	- ログ: `run_rho_in_audit.log`

- 4段比較（first-pass 実測）:
	1. DC ファイル生値 (`stage=dc-file-raw`)
		- `mean=1.309675655976680E-01`
		- `min=2.794763473796665E-02`
		- `max=9.734089396234934E-01`
		- `l2=1.319802033017150E+01`
	2. DG 注入直後 (`stage=dg-injected`)
		- 統計は DC と一致
		- `diff_max=0.0`, `diff_l2=0.0`
	3. DG box/global 再配置後 (`stage=dg-rho-in-pre-hartree`)
		- `mean=7.759275410629161E-02`
		- `min=0.0`
		- `max=7.021602554546760E-01`
		- `l2=1.007009388252944E+01`
		- 対 DC 差: `diff_max=9.734089396234934E-01`, `diff_l2=9.298198449005293E+00`
	4. 代表点（最大差上位）
		- 例: `ix=1,iy=1,iz=20` で `DC=9.734089396234934E-01`, `injected=9.734089396234934E-01`, `rho_in=0.0`
		- 上位10点の多くで `injected=DC` だが `rho_in=0`（または大幅低下）

- index/順序契約の観測:
	- `index-contract core_s=1 1 1, lg_num=24 24 24`
	- `mean_delta_pre_minus_dc=-5.337481149137639E-02`
	- `20^3` raw index をそのまま `rho_in(ix,iy,iz)` と比較すると大差が出る

- 判定（成功条件に対する結論）:
	- [x] ファイル読込直後（`dg-injected`）では差はない
	- [x] 差は DG 再配置後（`dg-rho-in-pre-hartree`）で初発
	- [x] よって本件の次フォーカスは Poisson ではなく、`20^3` raw index と `rho_iso` の global 配置対応（wrap/order 契約）

13.11 raw→global 逆写像監査（本命, `HSE-SR=0`, 2026-04-24）:

- 追加実装（13.10の拡張）:
	- `src/rt/dg/rt_dg_density_hamiltonian_update.f90`
	- `src/rt/dg/rt_dg_integrator_stage_update.f90`
	- 追加ログ/出力:
		- `[RHO-IN-AUDIT] preprocess mean_subtraction/hvol_scaling/normalization`
		- `stage=dg-rho-in-pre-hartree(actual-map)`
		- `stage=dg-rho-in-pre-hartree(expected-map)`
		- 逆写像ファイル: `dg_rho_map_audit_itt1_if1_rank0.dat`

- 実行条件:
	- 入力: `.namelist_dg_rt_c64_nt1_rsdump.tmp`
	- env: `SALMON_HSE_SR_HARTREE=0`, `SALMON_DG_VH_ISOLATION_FROM_DC_RHO=1`, `SALMON_DG_VH_ISOLATION_IFRAG=1`, `SALMON_DG_VH_ISOLATION_ITT=1`
	- ログ: `run_rho_in_audit.log`

- 実測（first-pass）:
	- `preprocess`: `mean_subtraction=off`, `hvol_scaling=off`, `normalization=off`
	- `actual-map`: `diff_max=9.734089396234934E-01`, `diff_l2=9.298198449005293E+00`
	- `expected-map`: `diff_max=9.734089396234934E-01`, `diff_l2=1.116505426384643E+01`
	- `index-contract`: `core_s=1 1 1`, `box_s=-3 -3 -3`, `lg_num=24 24 24`

- 逆写像ファイル観測例（`dg_rho_map_audit_itt1_if1_rank0.dat`）:
	- ヘッダ列:
		- `ix iy iz`
		- `ux/gx actual`
		- `ux/gx expected`
		- `dc_raw`, `rho_at_actual`, `rho_at_expected`
		- `diff_actual`, `diff_expected`
	- 代表点:
		- `ix=1,iy=1,iz=20`: `dc_raw=9.734e-01`, `rho_at_actual=0`, `rho_at_expected=0`
		- `ix=13,iy=1,iz=1`: `dc_raw=7.077e-01`, `rho_at_actual=1.794e-01`, `rho_at_expected=0`

- 判定:
	- [x] 単なる mean/hvol スケーリング要因ではない（3項目とも `off`）
	- [x] `raw 20^3 -> global` の配置契約（特に wrapped order / raw order, 軸対応）を主因候補として確定
	- [ ] 次は `dg_rho_map_audit_*.dat` を使って、どの raw 点がどの global 点へ落ちる/消えるかを規約表（式）として固定する

13.12 配置規約統一の試行（identity 20^3 配置, 2026-04-24）:

- 実装変更:
	- `src/rt/dg/rt_dg_density_hamiltonian_update.f90`
	- `src/rt/dg/rt_dg_integrator_stage_update.f90`
	- 変更内容:
		- raw→global 変換を一旦停止し、`rho_iso(ix,iy,iz)=rho_box(ix,iy,iz)` の identity 配置経路を導入
		- `Vh` 出力も同じ `ix,iy,iz` インデックスで取得

- 実行条件:
	- ログ: `run_rho_in_audit_identity.log`
	- env: `SALMON_HSE_SR_HARTREE=0`, `SALMON_DG_VH_ISOLATION_FROM_DC_RHO=1`, `SALMON_DG_VH_ISOLATION_IFRAG=1`, `SALMON_DG_VH_ISOLATION_ITT=1`

- first-pass 実測:
	- `stage=dg-rho-in-pre-hartree(identity-20x20x20)`
		- `diff_max=9.734089396234934E-01`
		- `diff_l2=9.298198449005293E+00`
	- `index-contract`:
		- `core_s=1 1 1`, `box_s=-3 -3 -3`, `lg_num=24 24 24`, `mg_ie=12 12 12`

- 重要観測:
	- `mg_ie=12 12 12` であり、rank0 の局所所有は `12^3`。
	- したがって rank0 単独の `rho_in` 参照では、`20^3` 全点比較にゼロ化が混ざる（所有外点が未保持）
	- 次段は rank 所有を跨ぐ形で `rho_in` を集約した上での差分評価が必要。

- 現時点の判断:
	- [x] 配置規約統一（identity 経路）の実装自体は反映済み
	- [ ] 「機械誤差まで低下」の判定は、rank 所有集約版の監査で再評価する

13.13 full-box 集約後の再判定（`HSE-SR=0`, 2026-04-24）:

- 実装変更（rank 所有バイアス除去）:
	- `src/rt/dg/rt_dg_density_hamiltonian_update.f90`
	- `src/rt/dg/rt_dg_integrator_stage_update.f90`
	- 変更内容:
		- `rho_in` 監査用に `20^3` の local buffer を各 rank で構築
		- `comm_summation(..., nx*ny*nz, dg_frag%icomm)` で full-box を集約
		- rank0 比較は集約後配列を使用

- 実行条件:
	- 入力: `.namelist_dg_rt_c64_nt1_rsdump.tmp`
	- env: `SALMON_HSE_SR_HARTREE=0`, `SALMON_DG_VH_ISOLATION_FROM_DC_RHO=1`, `SALMON_DG_VH_ISOLATION_IFRAG=1`, `SALMON_DG_VH_ISOLATION_ITT=1`
	- ログ: `run_rho_in_audit_aggregated.log`

- 実測（集約後, first-pass）:
	- `stage=dg-rho-in-pre-hartree(identity-20x20x20,aggregated-fullbox)`
		- `diff_max=0.0`
		- `diff_l2=0.0`
	- `stage=dg-rho-in-pre-hartree(legacy-core-order,aggregated-fullbox)`
		- `diff_max=0.0`
		- `diff_l2=0.0`
	- `mapping`: `mean_delta_pre_minus_dc=0.0`

- 判断:
	- [x] 旧 `diff_max=9.734e-01` は rank0 単独観測による所有分割バイアスが支配
	- [x] full-box 集約後は `DC raw` と `rho_in` が一致し、少なくとも本条件（`ifrag=1`, `core_s=1 1 1`）では配置規約起因の差は非主因
	- [ ] 次段は、`ifrag` を変えて `core_s != 1` 条件でも同じ full-box 監査を実施し、mapping バグの有無を一般条件で再判定する

13.14 Poisson 段別比較（`HSE-SR=0`, 2026-04-24, post-rho-aggregation）:

- 方針（13.13 の結論を反映）:
	- [x] `rho_in` の所有分割バイアスは除去済み（full-box 監査で `diff_max=diff_l2=0`）
	- [x] 本命を `Poisson solver / kernel / gauge / 境界条件` 差へ移す

- 実行条件:
	- 入力: `.namelist_dg_rt_c64_nt1_rsdump.tmp`
	- env:
		- `SALMON_HSE_SR_HARTREE=0`
		- `SALMON_POISSON_CONTRACT_TRACE=1`
		- `SALMON_DG_VH_ISOLATION_FROM_DC_RHO=1`
		- `SALMON_DG_VH_ISOLATION_IFRAG=1`
		- `SALMON_DG_VH_ISOLATION_ITT=1`
		- `SALMON_DG_VH_ISOLATION_DC_RHO_FILE=dc_rs_component_rho_if1_sp1_fragroot.dat`
	- ログ: `run_vh_iso_contract_hsesr0_post_rhoagg.log`

- first-pass 段別比較（DC block line 61 vs DG block line 155）:

| 段 | 指標 | DC | DG | DG-DC |
|---|---|---:|---:|---:|
| `rho_in` | mean | 1.054489256019869E-01 | 1.054489256019868E-01 | -9.714E-17 |
| `rho_in` | l2 | 1.441498652059725E+01 | 1.469474258847251E+01 | +2.798E-01 |
| `rho_in` | max | 3.606309490325457E-01 | 4.912127723415918E-01 | +1.306E-01 |
| `rhoG` | mean_abs | 1.563648171755354E-04 | 2.309786522459622E-04 | +7.461E-05 |
| `rhoG` | l2 | 1.108942828097870E-01 | 1.120889037979707E-01 | +1.195E-03 |
| `rhoG` | max_abs | 1.054489256019869E-01 | 1.054489256019868E-01 | -9.714E-17 |
| `vhG` | mean_abs | 3.662399368832726E-04 | 6.456817289369468E-04 | +2.794E-04 |
| `vhG` | l2 | 3.978463544118600E-01 | 4.890372107809040E-01 | +9.119E-02 |
| `vhG` | max_abs | 1.888489207909100E-01 | 2.212533709382038E-01 | +3.240E-02 |
| `Vh_raw` | mean | -6.705669969914791E-16 | -9.621932880084689E-16 | -2.916E-16 |
| `Vh_raw` | l2 | 6.084353682327360E+01 | 7.282350086653017E+01 | +1.198E+01 |
| `Vh_raw` | max | 1.956089171432791E+00 | 2.699883000453326E+00 | +7.438E-01 |
| `Vh_final` | mean | -6.705669969914791E-16 | -9.621932880084689E-16 | -2.916E-16 |
| `Vh_final` | l2 | 6.084353682327360E+01 | 7.282350086653017E+01 | +1.198E+01 |
| `Vh_final` | max | 1.956089171432791E+00 | 2.699883000453326E+00 | +7.438E-01 |

- 解釈:
	- [x] 本比較では `rho_in` の mean は一致だが `l2/max` には差が残る
	- [x] 差は `rhoG`、`vhG`、`Vh_raw` へと連続して残存・増幅
	- [x] `Vh_raw == Vh_final`（wall off 条件）
	- [x] 本命は Poisson 本体（solver/kernel/gauge/境界）で継続

- 次段（本命）:
	1. `rho_in` の同一性を mean だけでなく full-field 指標（`relL2`, corr, top-k diff）で DC/DG 同時に監査する
	2. `rhoG` / `vhG` について full-field dump を追加し、スペクトル空間で `relL2/corr` を直接比較する
	3. `yn_ffte=n/y` 固定 A/B を同じ表で再測定し、solver 実装差と kernel/contract 差を分離する

13.15 スペクトル段 full-field 比較（`HSE-SR=0`, C64, DC rho 注入, 2026-04-24）:

- 実装追加（worktree側）:
	- `src/poisson/poisson_periodic.f90`
		- `SALMON_POISSON_SPECTRAL_DUMP=1` で `rhoG` / `vhG` の full-field dump を出力
		- 出力形式: `# ix iy iz re im`
		- 出力ファイル:
			- `poisson_spectral_DC_FT_rhoG.dat`
			- `poisson_spectral_DC_FT_vhG.dat`
			- `poisson_spectral_DG_FT_rhoG.dat`
			- `poisson_spectral_DG_FT_vhG.dat`
			- `poisson_spectral_DC_FFTE_rhoG.dat`
			- `poisson_spectral_DC_FFTE_vhG.dat`
			- `poisson_spectral_DG_FFTE_rhoG.dat`
			- `poisson_spectral_DG_FFTE_vhG.dat`
	- `src/poisson/hartree.f90`, `src/poisson/poisson_dg_distributed.f90`
		- DC/DG 文脈ラベルを Poisson 側へ渡し、同一フォーマットで分離出力

- 実行条件:
	- 共通 env:
		- `SALMON_HSE_SR_HARTREE=0`
		- `SALMON_POISSON_CONTRACT_TRACE=1`
		- `SALMON_POISSON_SPECTRAL_DUMP=1`
		- `SALMON_DG_VH_ISOLATION_FROM_DC_RHO=1`
		- `SALMON_DG_VH_ISOLATION_IFRAG=1`
		- `SALMON_DG_VH_ISOLATION_ITT=1`
		- `SALMON_DG_VH_ISOLATION_DC_RHO_FILE=dc_rs_component_rho_if1_sp1_fragroot.dat`
	- `yn_ffte=n`:
		- 入力: `.namelist_dg_rt_c64_nt1_rsdump.tmp`
		- ログ: `run_poisson_spectral_ft.log`
	- `yn_ffte=y`:
		- 入力: `.namelist_dg_rt_c64_nt1_rsdump_ffte_y.tmp`（`&parallel` に `yn_ffte='y'` を追加）
		- ログ: `run_poisson_spectral_ffte.log`

- 評価指標（DC vs DG, complex full-field）:
	- `relL2 = ||DG-DC||_2 / ||DC||_2`
	- `corr_abs = corr(|DC|, |DG|)`
	- `maxabs = max |DG-DC|`

| mode | stage | N | relL2 | corr_abs | maxabs |
|---|---|---:|---:|---:|---:|
| FT | rhoG | 1728 | 7.8501982634E-02 | 9.9788337288E-01 | 2.3141304479E-03 |
| FT | vhG  | 1728 | 4.3771321288E-01 | 9.4110806368E-01 | 1.1689651690E-01 |
| FFTE | rhoG | 1728 | 7.8501982634E-02 | 9.9788337288E-01 | 2.3141304479E-03 |
| FFTE | vhG  | 1728 | 4.3249655790E-01 | 9.2360238846E-01 | 6.5773647432E+03 |

- 判定:
	- [x] 最初に差が立つスペクトル段は `rhoG`（`relL2≈7.85E-02` で非ゼロ）
	- [x] 差は `vhG` で増幅（`relL2≈4.3E-01`）
	- [x] `yn_ffte` A/B で `rhoG` 指標は同一、差の増幅は `vhG` 段で顕著
	- [ ] FFTE の `vhG maxabs` 異常大（`~6.58E+03`）は separate bug track で追加切り分け

13.16 `rhoG` 差の根因切り分け（`HSE-SR=0`, FT 主系列, 2026-04-24）:

- 目的:
	- [x] `rhoG` 差が FFT 前処理の単純契約差
		- 配列原点
		- index ordering
		- periodic wrap
		- mean handling
		- normalization
	  のいずれで説明できるかを確定する

- 追加実装（worktree 側）:
	- `src/poisson/poisson_periodic.f90`
		- `SALMON_POISSON_SPECTRAL_DUMP=1` 時に、FFT 直前の実空間入力 `rho` も one-shot dump
		- 出力形式: `# ix iy iz rho`
		- 出力ファイル:
			- `poisson_input_DC_FT_rho.dat`
			- `poisson_input_DG_FT_rho.dat`

- 実行条件:
	- DC:
		- 入力: `.namelist_dc_gs_64_for_dg_lcfo120_rs.tmp`
		- env: `SALMON_HSE_SR_HARTREE=0`, `SALMON_POISSON_CONTRACT_TRACE=1`, `SALMON_POISSON_SPECTRAL_DUMP=1`
		- ログ: `run_poisson_input_dc_ft.log`
	- DG:
		- 入力: `.namelist_dg_rt_c64_nt1_rsdump.tmp`
		- env: `SALMON_HSE_SR_HARTREE=0`, `SALMON_POISSON_CONTRACT_TRACE=1`, `SALMON_POISSON_SPECTRAL_DUMP=1`, `SALMON_DG_VH_ISOLATION_FROM_DC_RHO=1`, `SALMON_DG_VH_ISOLATION_IFRAG=1`, `SALMON_DG_VH_ISOLATION_ITT=1`
		- ログ: `run_poisson_input_dg_ft.log`

- first-pass 指標:
	- `rho_input(DC vs DG)`
		- `relL2 = 1.4182021094E-01`
		- `mean_delta = -2.7755575616E-17`
	- `rhoG(DC vs DG)`
		- `relL2 = 7.8501982634E-02`
		- 最良複素スケール `alpha = 1.0077493551E+00 + i*6.7504555778E-03`
		- `alpha` 適用後 `relL2 = 7.7226147110E-02`

- `rhoG` 差が最大の `G` 点 top-10（`gx gy gz rhoG_dc(re,im) rhoG_dg(re,im) absdiff`）:
	1. `0 2 0   ( 1.2730384765E-02, 3.2467870832E-03) -> ( 1.4764906153E-02, 4.3494758204E-03), 2.3141304479E-03`
	2. `2 0 0   ( 1.2409063265E-02, 3.3560107161E-03) -> ( 1.4400621948E-02, 4.4558880491E-03), 2.2750903577E-03`
	3. `0 0 4   ( 5.5847555380E-03, 4.9106473376E-03) -> ( 5.2103745955E-03, 6.9812360233E-03), 2.1041622075E-03`
	4. `0 0 1   (-5.1459552288E-12, 4.7101852071E-12) -> (-1.5048362472E-03,-1.3670597746E-03), 2.0330726882E-03`
	5. `0 0 2   ( 1.2409047931E-02, 3.3560159772E-03) -> ( 1.3686459621E-02, 4.8920854429E-03), 1.9978213209E-03`
	6. `0 4 0   ( 5.9034733110E-03, 4.8292282763E-03) -> ( 6.4169902302E-03, 6.3764362887E-03), 1.6302000675E-03`
	7. `2 0 1   (-1.0164362153E-12, 1.1150364169E-12) -> (-1.4571202094E-03, 5.7716897024E-04), 1.5672661933E-03`
	8. `4 0 0   ( 5.5847708133E-03, 4.9106434290E-03) -> ( 6.0012827123E-03, 6.4205604613E-03), 1.5663114654E-03`
	9. `0 2 1   (-1.1401075852E-12, 9.3813480105E-13) -> (-1.4231773902E-03, 4.1344598143E-04), 1.4820160119E-03`
	10. `0 0 3  ( 1.8230193589E-12,-4.4164842596E-12) -> (-9.2018411268E-04, 1.0260844687E-03), 1.3782554736E-03`

- FFT 前処理契約の確認結果:
	1. mean handling:
		- `rho_input_mean_delta ≈ 0`
		- `g=0` 成分も既報どおり一致
		- [x] 平均値処理差を主因から外す
	2. normalization:
		- 最良複素スケール `alpha` を当てても `rhoG relL2` は `7.85E-02 -> 7.72E-02` と軽微改善のみ
		- [x] 単純正規化差を主因から外す
	3. 配列原点 / periodic wrap:
		- スペクトル側で整数循環シフト + 最良 `alpha` を総当たりしても最良は `(sx,sy,sz)=(6,6,6)` で `relL2=7.7226147105E-02`
		- shift なし `alpha` と実質同等で、説明力は増えない
		- [x] 単純原点ずれ / wrap ではない
	4. index ordering:
		- 実空間 `rho_input` に対し、軸 permutation 6 通り + 循環シフト総当たりの最良でも
			- `perm=(2,1,0), shift=(12,0,12), relL2=1.4181941160E-01`
		- baseline `1.4182021094E-01` とほぼ同じ
		- [x] 単純 permutation / 並び替えでもない

- 実空間 `rho_input` 差分の支配パターン:
	- 最大差分点は buffer seam の単純シフト形ではなく、`z=1..2` かつ `x/y=1,13,14,23,24` 近傍に集中
	- top-12 例:
		- `(13,24,1): 3.3614744985E-01 -> 1.7189046417E-01, delta=-1.6425698567E-01`
		- `(13,23,1): 3.4493336838E-01 -> 1.8206297981E-01, delta=-1.6287038857E-01`
		- `(24,1,1):  3.3829626680E-01 -> 1.8594578332E-01, delta=-1.5235048348E-01`
	- 軸集約でも `z=1,13,24` と `x/y` 端面側に大きな重みが立ち、全体が一様にずれているわけではない

- `rhoG` top-8 差モードの逆変換による実空間パターン:
	- 強い点は `z=6/7` の 2 面と、`x,y` の 4 隅対称に現れる
	- 最大点例:
		- `(12,12,6)`, `(6,12,6)`, `(6,6,6)`, `(12,6,6)` で `|partial_r| ≈ 1.3324422687E-02`
		- `(1,1,7)`, `(7,1,7)`, `(7,7,7)`, `(1,7,7)` で `|partial_r| ≈ 1.3078788499E-02`
	- [x] 支配モードは「全体の位相傾き 1 本」ではなく、実空間で対称な多点構造を持つ

- 判定:
	- [x] `rhoG` 差は FFT 内で新たに生じたものではなく、FFT 直前の実空間 `rho_input` 差をそのまま反映している
	- [x] 根因は `配列原点`, `index ordering`, `periodic wrap`, `mean handling`, `normalization` の単純契約差では説明できない
	- [x] 次の本命は「Poisson 呼出前に DC と DG が保持している owner-local `rho` 値そのものの差」を upstream で切ること

13.17 DG_STAGE1_RHO 再構築段別監査（`ifrag=1` owner-local slab, 2026-04-24）:

- 目的:
	- [x] `DG_STAGE1_RHO` 生成直前の再構築経路を `basis -> coef -> occ -> assembly` に分解し、
	  `DC_RT_INIT_STAGE1_RHO` との差を段別に定量化する。

- 追加実装（worktree）:
	- `src/rt/dg/rt_dg_density_reconstruct.f90`
		- env: `SALMON_DG_RHO_RECON_STAGE_DUMP=1`
		- env: `SALMON_DG_RHO_RECON_STAGE_IFRAG=<ifrag>`
		- 追加ダンプタグ:
			- `DG_RECON_STAGE_BASIS`
			- `DG_RECON_STAGE_COEF`
			- `DG_RECON_STAGE_OCC`
		- 既存タグ:
			- `DG_STAGE1_RHO`
			- `DC_RT_INIT_STAGE1_RHO`

- 実行条件:
	- 入力: `.namelist_dg_rt_c64_nt1_rsdump.tmp`
	- env:
		- `SALMON_DG_ENFORCE_STATE_TRUNCATION=0`
		- `SALMON_HSE_SR_HARTREE=0`
		- `SALMON_RHO_STAGE_DUMP=1`
		- `SALMON_RHO_STAGE_DUMP_PREFIX=rho_stage_recon_audit`
		- `SALMON_DG_RHO_RECON_STAGE_DUMP=1`
		- `SALMON_DG_RHO_RECON_STAGE_IFRAG=1`
	- ログ: `run_recon_audit.log`

- owner-local slab（rank0, `is=1,1,1`, 1728点）での段別比較:
	- 比較先: `DC_RT_INIT_STAGE1_RHO`

| stage | relL2 | centered_relL2 | corr | maxabs |
|---|---:|---:|---:|---:|
| `DG_RECON_STAGE_BASIS` | `4.007292468965E-01` | `3.845981865602E-01` | `0.9681488020` | `1.254520424336E-01` |
| `DG_RECON_STAGE_COEF`  | `4.975049220604E-01` | `4.903483220016E-01` | `0.9841063805` | `2.001730931885E-01` |
| `DG_RECON_STAGE_OCC`   | `9.737121561785E-02` | `1.908711766599E-01` | `0.9841063805` | `1.326809853219E-01` |
| `DG_STAGE1_RHO`        | `9.737121561785E-02` | `1.908711766599E-01` | `0.9841063805` | `1.326809853219E-01` |

- `top-k diff`（rank0, `DG_RECON_STAGE_OCC` vs `DC_RT_INIT_STAGE1_RHO`, abs差上位5）:
	1. `(12,12,12): DC=3.585317870197E-01, DG=4.912127723416E-01, diff=+1.326809853219E-01`
	2. `(12,11,12): DC=2.712667326213E-01, DG=3.730568686688E-01, diff=+1.017901360475E-01`
	3. `(1,1,1):    DC=2.832299484597E-01, DG=1.857897595499E-01, diff=-9.744018890981E-02`
	4. `(11,12,12): DC=2.856770267324E-01, DG=3.806711401772E-01, diff=+9.499411344480E-02`
	5. `(12,12,11): DC=2.865121558421E-01, DG=3.811106926171E-01, diff=+9.459853677497E-02`

- 追加確認:
	- [x] rank0 owner-local では `DG_RECON_STAGE_OCC == DG_STAGE1_RHO`（`relL2=0`, `maxabs=0`）
	- [x] したがって、`ifrag=1` owner-local 観測では差の初発段は `occ` を掛けた密度組み立て式

- 14.18% 差の式レベル確定:
	- [x] 全 owner-local slab 合成（13824点）では `DG_STAGE1_RHO` 対 `DC_RT_INIT_STAGE1_RHO` が
	  `relL2 = 1.418202109398E-01`（既報値と整合）
	- [x] よって `~14.18%` 差は Poisson 前、`DG_STAGE1_RHO`（密度再構築の `occ` 反映後）で既に立っている
	- [ ] 次段は `DG_RECON_STAGE_COEF` の寄与差を同一 full-field 契約で抽出するため、
	  remote/send 成分を含む段別全体場ダンプへ拡張する

13.18 DG_STAGE1_RHO 最終式の DC 直接突合（`ifrag=1`, rank0 slab top diff 点, 2026-04-24）:

- 目的:
	- [x] `DG_RECON_STAGE_OCC` の最終式
	  `rho_DG(r) = Σ_n occ_n |Σ_i c_i^(DG)(n) φ_i^(DG)(r)|^2 = Σ_ij D_ij^(DG) φ_i^(DG)(r) φ_j^(DG)(r)`
	  を、DC 側の
	  `rho_DC(r) = Σ_n occ_n |Σ_i c_i^(DC)(n) φ_i^(DC)(r)|^2 = Σ_ij D_ij^(DC) φ_i^(DC)(r) φ_j^(DC)(r)`
	  と 1 項ずつ比較し、差分項を確定する。

- 追加実装（worktree）:
	- `src/rt/dg/rt_dg_density_reconstruct.f90`
		- env: `SALMON_DG_RHO_FORMULA_AUDIT=1`
		- env: `SALMON_DG_RHO_FORMULA_AUDIT_IFRAG=<ifrag>`
		- 出力:
			- `dg_rho_formula_audit_coef_if1_rank0000.dat`
			- `dg_rho_formula_audit_basis_if1_rank0000.dat`
			- `dg_rho_formula_audit_state_if1_rank0000.dat`
	- DC 側は既存の `data_dcdft/fragments/000001/{wavefunctions,basis_functions,rgrid_index}.bin` を直接再構成に使用。

- 実行条件:
	- 入力: `.namelist_dg_rt_c64_nt1_rsdump.tmp`
	- env:
		- `SALMON_DG_ENFORCE_STATE_TRUNCATION=0`
		- `SALMON_HSE_SR_HARTREE=0`
		- `SALMON_DG_RHO_FORMULA_AUDIT=1`
		- `SALMON_DG_RHO_FORMULA_AUDIT_IFRAG=1`
	- ログ: `run_recon_formula_audit.log`

- first-pass 判定（指定 3 点）:
	1. `occ weighting`
		- [x] DG audit の occupied state weight は一様に `2.0`。
		- [x] DC 側 `nelec=256, nspin=1` 条件と一致し、`occ` は差分源ではない。
	2. `basis normalization / basis value`
		- [x] fragment 1 の `φ_i(r)` は指定 3 点で DC 再構成値と機械誤差一致。
		- 例: `phi relL2 ~ 1E-16`。
		- [x] 少なくともこの 3 点では `basis value / normalization` は差分源ではない。
	3. `pair sum range`
		- [x] fragment 1 の basis 本数は DC/DG とも `19`。
		- [x] 指定 3 点では DC 側も fragment 1 のみ寄与し、cross-fragment term はゼロ。
		- [x] よってこの 3 点では `pair sum range` は差分源ではない。
	4. `spin/state indexing / coef product`
		- [x] DG の coefficient 行列は DC fragment 1 ファイル値と非一致。
		- `coef maxabs diff = 1.03310727095`, `coef relL2 diff = 8.61235191633E-01`。
		- `D = Σ_n occ_n c_i c_j` で比較すると
		  `relL2(ΔD) = 1.79633721987E-01`, `maxabs(ΔD) = 6.31818459069E-01`。
		- [x] state 列の対応も崩れているが、単なる permutation ではなく `D` 自体が変化している。
		- [x] したがって差分源は `coef product` 項であり、`spin/state indexing` は「単なる列並べ替え」では説明できない。

- 指定点での式レベル結論:
	- `(12,12,12)`
		- `rho_DC = 3.585317870197E-01`
		- `rho_DG = 4.912106016955E-01`
		- `delta = +1.326788146758E-01`
		- 支配 `ΔD_ij φ_i φ_j` 項:
			- `(16,16): +6.193217861266E-02`
			- `(19,16)+(16,19): +4.845700402708E-02`
	- `(12,11,12)`
		- `rho_DC = 2.712667326213E-01`
		- `rho_DG = 3.730594701846E-01`
		- `delta = +1.017927375633E-01`
		- 支配 `ΔD_ij φ_i φ_j` 項:
			- `(16,16): +4.331507171547E-02`
			- `(19,16)+(16,19): +4.165880002561E-02`
	- `(1,1,1)`
		- `rho_DC = 2.832299484597E-01`
		- `rho_DG = 1.858129870924E-01`
		- `delta = -9.741696136727E-02`
		- 支配 `ΔD_ij φ_i φ_j` 項:
			- `(16,19)+(19,16): -6.812520055869E-02`
			- `(16,16): +2.623899733844E-02`
			- `(19,19): -1.700883053398E-02`

- 結論:
	- [x] 指定 3 点では差分は `occ weighting`, `basis normalization`, `pair sum range` ではなく、
	  `D_ij = Σ_n occ_n c_i c_j` の DG/DC 非一致で立っている。
	- [x] よって `DG_STAGE1_RHO` の最終差分項は `ΔD_ij φ_i(r) φ_j(r)` であり、
	  root cause は `coef product` 側にある。
	- [ ] 次段は `dg_frag%coef` が DC file coefficient からどの段（初期化 / overlap / propagation 前処理）でズレるかを upstream 追跡する。

## 14. DG Core Weak-Form State Alignment 仕様（JA, 2026-04-24）

目的:

- RT 側で用いる正準 Hamiltonian と初期状態の整合条件を固定し、buffered DC seed を投入した場合でも「state は core 弱形式 Hamiltonian に整合している」ことを保証する
- startup 時の便宜的な再割当てや relabel に依存せず、DG RT の初期状態を `S_core`, `H_core` に対する一般化固有値問題の解として定義する

正準定義:

- 基底自由度は core 領域に属する DOF を正とする
- overlap は core 弱形式で定義する
	- `S_core(i,j) = ∫_core φ_i(r) φ_j(r) dr`
- Hamiltonian は core 弱形式で定義する
	- `H_core(i,j) = ∫_core φ_i(r) (H φ_j)(r) dr`
- buffer は `Hφ` の評価安定化と stencil/support のためにのみ用い、最終積分領域には含めない
- fragment 内の周期補完は fragment-periodic を正とし、global periodic completion は正準定義に使わない
- nonlocal pseudopotential は self fragment の core 弱形式として評価し、他 fragment 由来の混入を正準定義へ持ち込まない

状態の正準性:

- RT 初期状態は `H_core c = ε S_core c` の解でなければならない
- buffered 領域まで含めた DC 側固有関数は、そのままでは RT 正準状態とみなさない
- DC 由来の buffered state は seed としてのみ扱い、必要なら RT 側で `S_core`, `H_core` に対する一般化 EVP を解いて正準状態へ写像する

startup ポリシー:

- buffered DC dataset に対して、fragment-local generalized EVP の固有ベクトルを global state へ逐次上書きする startup stationary projection は既定で禁止する
- 既定動作は「skip」であり、明示 override がある場合のみ診断目的で有効化してよい
- 許可される startup 補正は以下に限る
	- state label を保存する Lowdin 正規直交化
	- 占有部分空間内の unitary rotation
- 以下は明示的に禁止する
	- `frag_state_begin`, `nocc_frag`, local eigenvector index に基づく sequential overwrite
	- `esp_old -> esp_new` の破壊的 relabel
	- occupied/unoccupied をまたぐ state injection

整合判定:

- 正常系では startup 前後で state label の意味が保存されること
- `C† S_core C - I` の残差が機械精度近傍であること
- Rayleigh 指標が `H_core` に対して自己無撞着であること
- 密度再構成は正準状態集合のみを入力として行い、startup の途中生成物や fragment-local assignment 結果を直接使わないこと

監査要件:

- upstream audit は少なくとも以下を比較対象に含める
	- startup 入力 coef
	- `S_core`, `H_core` で定義された固有ベクトル
	- startup 後 coef
- 指標は少なくとも以下を保持する
	- `relL2`, `maxabs`
	- `D` の `topDdiag`, `topDoff`
	- `C†SC` 残差
	- Rayleigh 残差
	- `esp_old ↔ esp_new` 対応表
	- `topDdiag` を駆動した global state の assignment source
- root cause 切り分けでは rollback の前、すなわち `before-rollback` を第一観測点とする

実装契約:

- `src/rt/dg/rt_dg_fragment_hamiltonian.f90`
	- startup policy の入口
	- `S_core`, `H_core` に整合する初期状態構成の実装責務
- `src/rt/dg/rt_dg_fragment.f90`
	- upstream audit stage と state-label 監査の実装責務
- `src/rt/dg/rt_dg_density_reconstruct.f90`
	- 正準状態のみを用いた密度再構成の実装責務
- `src/gs/dc/lcfo.f90`
	- seed 供給側の実空間/状態情報が必要な場合の補助実装対象

当面の安全運用:

- buffered basis を使う RT startup では stationary projection を既定で無効にし、このルールを維持する
- 次段の本実装は「startup で補修する」のではなく、「RT 正準状態を core 弱形式で再定義する」方向で進める
- 既存の audit/diagnostic は、新実装後も回帰判定に継続使用する
