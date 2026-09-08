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

### 7.4 mixed-Z local production-compatible route

mixed-Z local route は、現時点では guarded production-compatible route として扱う。
通常の本線 flag は `&propagation` namelist に明示する。環境変数は後方互換
または一時 override 用に残すが、本線入力では namelist を優先する。

| 対象 | namelist flag | legacy env override | production-compatible target |
| --- | --- | --- |
| mixed-Z route | `yn_dg_mixed_z='y'` | `SALMON_DG_MIXED_Z=1` | enable mixed-Z production-compatible route |
| propagation | `yn_dg_mixed_z_local_prop_writeback='y'` | `SALMON_DG_MIXED_Z_LOCAL_PROP_WRITEBACK=1` | production-equivalent split backend |
| propagation backend | `dg_mixed_z_local_prop_backend='global_mixed_split_backend'` or `'fragment_local_mixed_split_backend'` | `SALMON_DG_MIXED_Z_LOCAL_PROP_BACKEND` | selected accepted backend |
| field block | `dg_mixed_z_frag_local_field_block='all'` | `SALMON_DG_MIXED_Z_FRAG_LOCAL_FIELD_BLOCK` | fragment-local field block policy |
| density | `yn_dg_mixed_z_local_rho_writeback_wwonly='y'` | `SALMON_DG_MIXED_Z_LOCAL_RHO_WRITEBACK_WWONLY=1` | WW-only production grid target |
| Pz | `yn_dg_mixed_z_local_pz_writeback_total='y'` | `SALMON_DG_MIXED_Z_LOCAL_PZ_WRITEBACK_TOTAL=1` | `WW_full + WP/PW + PP` |
| current | `yn_dg_mixed_z_local_current_writeback_total='y'` | `SALMON_DG_MIXED_Z_LOCAL_CURRENT_WRITEBACK_TOTAL=1` | `para_WW + nonlocal + diamagnetic` |

density は意図的に `WWONLY` production target であり、`TOTAL` target ではない。
Pz と current は production-compatible な `TOTAL` target を使う。

dry-run / series validation には、必要に応じて

```text
SALMON_DG_MIXED_Z_LOCAL_PROP_OBS_SERIES_SOURCE=local_mixedz_writeback
```

を使い、`rho_diff`, `Pz_diff`, `current_diff_norm`, `coef_diff_Snorm`
を確認する。

未採用の探索経路は production 本線から削除する。
削除対象と理由は以下。

| removed route/log | reason |
| --- | --- |
| old Pz WW-only replacement | superseded by Pz TOTAL |
| current source-bridge diagnostic | accepted current TOTAL path fixed |
| current momentum-block diagnostic | exploration complete |
| `DG-CURRENT-BLOCKS` trace | redundant block trace |
| global-slice comparison | rejected convention/source path |

### 7.5 mixed-Z local performance stage

7.4 の acceptance は正しさの acceptance であり、効率化 acceptance ではない。
現行の propagation candidate は

```text
candidate_kind = global_mixed_split_backend
```

であり、production-equivalent な安全経路として扱う。
弱スケーリング向けの本線は、同じ correctness acceptance を保ったまま

```text
candidate_kind = fragment_local_mixed_split_backend
```

へ落とすことである。

Performance stage は正しさ検証と分けて、以下の順に進める。

| stage | 目的 | 主な確認 |
| --- | --- | --- |
| Perf-0 | 計測 | walltime breakdown, heavy kernel call count, MPI reduction count, dense build / zgemm / eigensolve / expdiag count |
| Perf-1 | cache 化 | field-independent `S/Z/H` block, local basis mapping, operator block, observable work arrays |
| Perf-2 | per-step 重複削除 | rho/Pz/current/writeback/series validation が同じ係数・block を二重走査していないか |
| Perf-3 | local 化 | global mixed split backend を fragment-local backend へ置換し、通信を boundary/neighbor に制限 |
| Perf-4 | acceptance | correctness diff は roundoff、かつ walltime, memory, MPI count が改善 |

Performance acceptance でも correctness は必ず維持する。

```text
coef_diff_Snorm = 0
rho_diff = 0
Pz_diff = roundoff
current_diff_norm = 0
```

効率化ログや探索ログは developer/debug diagnostic として隔離し、
default runtime では accepted path のログだけを見せる。

### 7.6 mixed-Z fragment-local field approximation policy

`fragment_local_mixed_split_backend` の field-on 経路は、DG locality を保つため
local block を以下に固定する。

```text
W_owner + P_self + face-neighbor P
```

標準設定は

```text
SALMON_DG_MIXED_Z_FRAG_LOCAL_FIELD_BLOCK=all
```

であり、`all` は self + 6 face-neighbor layout の範囲内での最大 block を意味する。

Scale-3d の metadata-only 診断では、shell=2 は `Z` coverage を改善せず、
`missing_P_coef_count` だけを増やした。そのため、shell=2 P layout 拡張や
global full-Z equivalence を fragment-local backend の必須条件にはしない。

field-on の global full-Z との差は、まず local field truncation error として扱う。
今後の評価は完全一致ではなく、以下の安定性・収束性で行う。

```text
rho/Pz/current finite
no bad/NaN/FATAL
norm/energy drift controlled
field strength dependence
block-size dependence within DG-local layout
```

Run summary は `SALMON_DG_MIXED_Z_PERF_COUNT=1` で
`[DG-MIXEDZ-FRAG-LOCAL-STABILITY-SUMMARY]` として出力する。

### 7.7 mixed-Z response comparison reference policy

Dielectric / conductivity response の最終比較では、Full TDDFT reference と
mixed-Z / split 追加誤差を分けて扱う。

Full TDDFT reference は別条件として固定する。LCFO/DC 結果から Full TDDFT を
行う場合は、DG-Fragment RT ではなく conventional-from-DC 経路を使う。

```text
method     = real-space TDDFT
propagator = Taylor propagation in the conventional RT path
dt         = 0.02 au
route      = yn_conventional_from_dcdft='y'
```

この経路は DC-LCFO fragment data から conventional real-space wavefunction を
再構成し、その後は通常の real-space TDDFT propagation を行う。DG 側の
`yn_dg_fragment_from_dcdft='y'` / `yn_dg_fragment_rt='y'` とは混ぜない。
diamond64 で使う場合は、同じ `data_dcdft` と構造・grid・state 条件を使い、
RT 入力では以下を基本形にする。

```text
&calculation
  theory = 'tddft_pulse'
  yn_conventional_from_dcdft = 'y'
/

&propagation
  propagator = 'middlepoint'
  yn_dg_fragment_rt = 'n'
/
```

通常 RT の `propagator='middlepoint'` 経路では内部の real-space Taylor
propagator が適用される。`time_integrator_dg_fragment='taylor4pc'` は
DG-Fragment RT 専用の指定なので、Full TDDFT reference の propagator 表記には
使わない。

比較では物理時刻をそろえる。`T_max` を固定し、各手法の `dt` に応じて
`nt` を調整する。例えば mixed-Z 側が `dt=0.1 au, nt=100` なら
`T_max=10 au` であり、Full TDDFT 側は `dt=0.02 au, nt=500` とする。

FFT の前処理は物理時間単位で指定する。

```text
skip_time
T_max
damping
window
padding
```

time grid が異なる場合は、各 `dt` のまま FFT して同じ energy grid へ補間して
比較する。field waveform は物理時間で一致させ、impulse/kick では特に
同じ impulse area を使っているかを確認する。

比較表には少なくとも以下を含める。

```text
case | method | propagator | dt_au | block | peak_shift_eV | peak_height_error | rms_error | bad
```

誤差は次の二段に分ける。

```text
Full TDDFT vs global WF+PW all:
  mixed-Z / block basis approximation error

global WF+PW all vs split / fragment-local WF+PW all:
  split / localization additional error
```

このため、fragment-local backend の acceptance では Full TDDFT との差を直接
単独の pass/fail 条件にしない。まず global WF+PW all を中間 reference として、
mixed-Z 近似誤差と fragment-local 化の追加誤差を分離する。

current-derived response の軽量比較には以下を使う。

```text
python3 tools/compare_dg_current_response_fft.py \
  --reference FULL_OR_GLOBAL_REF_rt.data \
  --axis z \
  --skip-time-au SKIP_AU \
  --tmax-au TMAX_AU \
  --damping DAMP_AU \
  --pad-factor PAD \
  --emin EMIN_EV --emax EMAX_EV \
  --table \
  --candidate-meta case=...,method=...,propagator=...,dt_au=...,block=... \
  CANDIDATE_rt.data
```

候補の time grid が reference と異なる場合、各 trajectory は自分の `dt` で FFT
され、candidate spectrum が reference energy grid へ補間される。

### 7.8 Stage 1g/1h final laser-response benchmark design

laser-pulse excitation の最終評価では、無理に誘電関数へ変換するよりも、
polarization-derived HHG / harmonic spectrum を主指標にする。弱場・線形応答を
見る場合だけ、`Pz(t) -> epsilon(omega)` を副指標として使う。
kick/current-derived response は、外場定義や current convention の差を
切り分けるための補助診断として扱う。

```text
Stage 1g:
  laser-field convention bridge
  purpose:
    classify Ac/E waveform, gauge, pulse envelope, polarization direction,
    and, only as needed, current sign, volume normalization, nonlocal current treatment
  target:
    laser waveform and Pz(t) / auxiliary Jz(t) raw waveform
  not final accuracy metric:
    HHG / harmonic spectra are evaluated from Pz(t) in Stage 1h

Stage 1h:
  polarization-derived laser-response benchmark
  final metrics:
    Pz(t) -> HHG / harmonic spectrum
    harmonic peak position
    harmonic peak height
    plateau / cutoff tendency
    spectrum RMS over selected harmonic/energy windows
    Pz(t) -> epsilon(omega) only for weak-field linear-response checks
    Jz(t)-derived spectrum as secondary/current-derived diagnostic
```

論文・報告書用の最終比較は以下の 5 case を基本にする。

| case | basis | fragment | purpose |
| --- | --- | --- | --- |
| Full TDDFT | real-space Taylor | no | reference laser-pulse HHG / polarization response |
| Global WF | WF only | no | basis truncation laser response |
| Global WF+PW | WF+PW all | no | mixed-Z laser response |
| Fragment WF | WF only | yes | fragmentation + WF laser response |
| Fragment WF+PW | WF+PW all | yes | final method laser response |

Current implementation note: field-on `Global WF` with the expdiag formal
length-gauge path is intentionally guarded because the local-only field path
gave nonphysical `P(t)`. Until a supported global WF field path is added,
`Global WF` is a pending case. Do not bypass this guard for benchmark tables.
`Fragment WF` can be evaluated through the `taylor4pc` DG-Fragment route.

誤差は次のように分解する。

```text
Full TDDFT -> Global WF / Global WF+PW:
  basis / mixed-Z approximation error

Global WF -> Fragment WF:
  WF-only fragmentation additional error

Global WF+PW -> Fragment WF+PW:
  final fragmentation additional error
```

Stage 1h の primary acceptance は、Full TDDFT との絶対一致ではなく
`Global WF+PW -> Fragment WF+PW` の追加誤差で判定する。

```text
Stage 1h acceptance:
  accepted target:
    Global WF+PW -> Fragment WF+PW HHG benchmark
  Fragment WF+PW reproduces Global WF+PW harmonic spectrum
  peak_shift_eV is zero or classified as sufficiently small
  cutoff_shift_eV is zero or classified as sufficiently small
  harmonic intensity errors of order 10% are treated as negligible for the
  present DG-fragmentation accuracy target
  bad = F

Stage 1h report:
  rel_rms_error
  hhg_window_rms_error
  harmonic-order relative intensity errors
```

Full TDDFT との差は、basis completeness / observable bridge /
gauge-field convention を含む総合誤差として別枠で評価する。
DG fragmentation の主張では、`Global WF+PW -> Fragment WF+PW`
の差を最重要指標にする。

比較表の標準列は以下に固定する。

```text
case
method
basis
fragment
propagator
dt_au
Tmax_au
field_abs
axis
block
observable
peak_shift_eV
peak_height_error
rel_peak_height_error
rms_error
rel_rms_error
low_energy_rel_rms_error
hhg_window_rms_error
cutoff_shift_eV
harmonic_rel_intensity_errors
bad
```

For the final laser-pulse benchmark, polarization-derived HHG / harmonic
spectra are preferred over dielectric-function conversion or
kick/current-derived response. This removes ambiguity from kick-definition
differences and uses the natural WF+PW observable, `P(t)`.

If conventional Full TDDFT does not directly output `P(t)`, use integrated
matter current only as a Stage-1g bridge diagnostic until the sign, gauge, and
normalization convention is fixed. The final Stage-1h table should clearly mark
whether the reference observable is direct `Pz(t)` or current-integrated `Pz(t)`.

### 7.9 Stage 2 Full TDDFT vs WF+PW discrepancy analysis

Stage 1h で

```text
Global WF+PW -> Fragment WF+PW:
  fragmentation additional error accepted
```

と整理できたため、以後の大きな差は DG fragmentation ではなく
Full TDDFT と global WF+PW 表現の差として扱う。Stage 2 では
fragment-local backend は一旦固定し、Global WF+PW を代表として
Full TDDFT reference との差の物理・数値的起源を切り分ける。

誤差源は以下の layer として分離する。

```text
Layer 0:
  Full TDDFT

Layer 1:
  Global WF+PW

Layer 2:
  Fragment WF+PW (DG)
```

`Layer 0 -> Layer 1` は mixed-Z 近似、observable convention、propagator、
basis truncation を含む。`Layer 1 -> Layer 2` は DG fragmentation、
local field truncation、fragment boundary、communication/localization を含む。
この分離により、Stage 1 は `Layer 1 -> Layer 2`、Stage 2 は
`Layer 0 -> Layer 1` を主対象にする。

Stage 2 は次の 3 phase に分ける。

```text
Stage 2A: observable / convention bridge
  confirm that the same physical quantity is compared
  - laser pulse E(t), A(t), envelope, CEP
  - gauge convention
  - P(t) / J(t) relation
  - current sign and electron-charge convention
  - nonlocal current contribution
  - volume normalization

Stage 2B: propagator convergence
  classify numerical time-integration error
  - Full TDDFT dt dependence, e.g. dt=0.02 -> 0.01 au
  - Full TDDFT reference acceptance is based primarily on Ne conservation
  - Taylor order dependence when available
  - WF+PW expdiag / split dt dependence
  - HHG peak / cutoff / intensity changes vs dt

Stage 2C: basis completeness convergence
  classify mixed-Z / WF+PW representation error
  - WF count
  - PW count / cutoff
  - PW shell radius
  - projector count
  - HHG spectrum as a function of basis size
```

For Full TDDFT, Stage 2B does not require the dt=0.02 and dt=0.01 response
spectra to agree. If dt=0.02 shows a large drift of `Ne` away from 256, a large
change in HHG / P(t) / current-integrated P(t) at dt=0.01 is expected and should
be treated as removal of a dt=0.02 artifact, not as evidence that the dt=0.01
reference is unusable.

Full TDDFT reference acceptance is therefore:

```text
Primary:
  Ne stays sufficiently close to 256 over the selected time window
  max|Ne - 256|, 256 - Ne_min, and 256 - Ne_final are small enough for the run
  no bad/NaN/FATAL

Secondary:
  dt=0.02 vs dt=0.01 response difference documents the time-step artifact
  HHG / P(t) / J-integrated P(t) from the accepted dt are used as reference
```

Once the dt=0.01 Full TDDFT run satisfies the Ne criterion, the remaining
Full TDDFT vs Global WF+PW discrepancy should be routed back to Stage 2A/2C:
observable convention, gauge/current bridge, and basis / mixed-Z completeness.
Do not mark the dt=0.01 Full TDDFT reference as bad only because it differs
substantially from a dt=0.02 run whose `Ne` drift is already large.

Stage 2C is the highest-priority physics benchmark. A convergence curve such as

```text
Full TDDFT
  <- Global WF+PW(shell=3)
  <- Global WF+PW(shell=2)
  <- Global WF+PW(shell=1)
```

would support the interpretation that the remaining Full TDDFT gap is mainly a
basis / mixed-Z completeness error rather than a fragmentation error.

Stage 2C acceptance は以下に固定する。

```text
Primary:
  HHG spectrum converges with basis size
  plateau intensity converges with basis size
  cutoff converges with basis size
  harmonic-order intensities approach the Full TDDFT reference trend

Secondary:
  polarization RMS
  current RMS
  total energy drift
  norm / charge drift
```

HHG では time-domain RMS よりも cutoff、harmonic intensity、plateau shape を
主評価にする。RMS は補助指標として残す。

Stage 2 comparison tables should keep the Stage 1h HHG columns and add the
basis-control fields needed for convergence plots.

```text
case
method
basis
WF_count
PW_count
PW_cutoff_or_shell
projector_count
dt
propagator_kind
observable_source
gauge
volume_normalization
field_abs
axis
block
peak_shift_eV
rel_peak_height_error
hhg_window_rms_error
HHG_peak_positions
HHG_plateau_metric
HHG_cutoff
cutoff_shift_eV
harmonic_intensities
polarization_RMS
current_RMS
energy_drift
norm_charge_drift
bad
```

Stage 2C sweep comparison can be generated from a TSV/CSV manifest with

```text
python3 tools/compare_global_wfpw_hhg_sweep.py \
  --reference FULL_TDDFT_REFERENCE_rt.data \
  --reference-source rt_current_integral \
  --manifest global_wfpw_basis_sweep.tsv \
  --axis z \
  --tmax-au TMAX_AU \
  --emax EMAX_EV \
  --hhg-emin HHG_EMIN_EV --hhg-emax HHG_EMAX_EV \
  --fundamental-ev LASER_OMEGA_EV \
  --harmonic-orders 3,5,7,9,11,13
```

The manifest must contain `path` and may contain the basis-control columns
listed above. Missing `current_RMS`, `energy_drift`, and `norm_charge_drift`
are left blank unless supplied by a separate run summary.

Stage 2C manifest schema is fixed as follows.

```text
Required:
  path
  basis_id
  WF_count
  PW_count
  PW_cutoff_or_shell

Recommended:
  dt
  propagator_kind
  observable_source
  gauge
  volume_normalization
```

Bad classification is available through optional tolerances in the sweep
script. A row should be classified as bad if the input time series cannot be
read, the peak/cutoff is outside tolerance, or the HHG-window RMS is outside
tolerance. Supplied norm/charge drift is recorded for diagnostics; use it as a
bad criterion only when the comparison target explicitly requires that. For the
Full TDDFT reference, the relevant charge diagnostic is the separate Stage 2B
`Ne` acceptance described above, not an automatic bad flag on every response
comparison row.

For plotting, prefer expanded harmonic-order columns such as
`H03_rel_intensity_error`, `H05_rel_intensity_error`,
`H07_rel_intensity_error`, ... over parsing the compact
`harmonic_intensities` string.

Stage 2C convergence plots are generated from the summary table with

```text
python3 tools/plot_global_wfpw_hhg_sweep.py \
  --summary global_wfpw_basis_sweep_summary.tsv \
  --plot-dir plots/stage2c \
  --x-axis PW_count
```

Use `--x-axis PW_cutoff_or_shell` for shell/cutoff sweeps. The plotter writes
one compact TSV per metric and writes PNG when `matplotlib` is available.
Without `matplotlib`, it writes simple SVG files. Bad rows are excluded by
default; use `--include-bad` only for debugging.

Stage 2C-5 plot acceptance:

```text
cutoff_shift_eV vs basis size is generated
HHG_plateau_metric vs basis size is generated
Hxx_rel_intensity_error vs basis size is generated
bad=T rows are excluded by default or explicitly included
output file names include x-axis and metric
```

Full TDDFT agreement is the final physical validation target, but the first
paper-level claim for DG fragmentation remains:

```text
Global WF+PW and Fragment WF+PW give the same HHG response within the accepted
fragmentation tolerance.
```
