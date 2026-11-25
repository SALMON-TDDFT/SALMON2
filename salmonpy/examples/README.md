# Pythonサンプル集

[English](README.en.md)

[初回テスト](#first-test) · [既存計算への追加](#add-to-existing)

<a id="add-to-existing"></a>
## 今ある計算に機能を追加する

導入と基本計算を済ませた方は、ここから目的を選べます。再インストールや初回テストは不要です。手元のモデル・GS・波形・保存結果を使い、追加する機能に必要な準備だけを行います。モデルや状態の対応条件は各ガイドで確認してください。

| 追加したいこと | 手元で使うもの | 追加準備・最短の入口 | 得られるもの |
|---|---|---|---|
| 実験の光学定数からLDを使う | 材料ID・利用する周波数帯 | [光学材料](../OPTICAL_DATA.md)：`load_material(...).model()`。保存済み係数の利用は再フィット不要 | LDモデルとフィット誤差・出典 |
| 誘電率の方向依存を入れる | 等方LDの係数、主軸・材料の向き | [テンソルLD](../LD.md)：`TensorLorentzDrude`。非対角成分とFDTDの条件を確認 | 異方性の分極・電流・誘電率 |
| 電流から誘電率・伝導率を求める | 同じ時刻列の電荷電流jと電場E | [スペクトル解析](../OPTICAL_TRAJECTORIES.md)、[計算例](../WORKFLOWS.md)。完全テンソルには独立な駆動方向が必要 | 光学スペクトルと有効な周波数範囲 |
| 同じ物質に別のパルスを与える | 保存GSと新しい波形 | [初期状態](../INITIAL_STATES.md)、[保持GS](../NATIVE_GS.md)。既知の未収束GSは先に残差を確認 | 独立したTDDFT応答 |
| SCFを続ける・初期推定を再利用する | 保持中のGS、SCFチェックポイント、GSアーカイブ | [SCF操作](../NATIVE_GS.md)。同じ計算の再開と、新しい計算の初期推定を区別 | 混合履歴を含む継続、または独立した新SCF |
| 電流・磁化・スピン流をまとめて読む | 出力ファイル、再開前後のファイル一覧 | [出力読込み](../OUTPUTS.md)：`read_table`／`read_spin_dynamics`。ソルバー不要 | 単位と時刻を保持した配列 |
| LD・SBE・TDDFTを混ぜる | 材料定義・配置・使用可能ランク数 | [混在計算](../COUPLED.md)：`mpi_mixed_native`。対応するFDTD・電子計算ライブラリとMPI | 場と各点の物質応答 |
| 実経過時間で保存して終了する | 既存の更新ループと一括保存関数 | [保存と再開](../CHECKPOINTING.md)：`CheckpointSchedule`。保存時間の余裕を確保 | 次ジョブで再開できる状態一式 |
| 更新ごとの検査を省く | 動作確認済みの同じモデル・入力 | [性能モード](../NATIVE_STEPPING.md)、[LD](../LD.md)、[SBE](../SBE.md)。対応部品で`execution_mode="performance"`を明示 | 検査・退避の負担を減らした時間発展 |
| FDTDの反射・透過・吸収を確かめる | E/H面データ、物質のE/j | [エネルギー収支](../API.md)。面の位置・時刻・符号を揃える | 散乱・回折も含めた収支診断 |
| 格子幅・dtを選び直す | 物質・擬ポテンシャル・評価したい電流 | [収束診断](../CONVERGENCE.md)。計画作成と計算実行を分ける | 電子数保存と電流差の比較 |
| 自作の物質応答へ差し替える | 自作のA→j／E応答と単位・時刻の定義 | [自作コードの接続](../MATERIALS.md)。まず単独応答と適合検査を行う | 既存の場との結合 |

### 保存済みの光学材料をLDとして使う

```python
from salmonpy.optical_materials import load_material

data = load_material('au-johnson-1972')
model, fit = data.model()  # 保存済みの係数。ダウンロード・再フィットなし
print(fit['relative_rms_error'], data.metadata['conditions'])
```

必要なのは導入済みwheelです。フィットの周波数帯と実験条件を確認してから使います。別の実験データへの再フィットには[LDフィッティング](../LD.md)を使い、選ぶ方法によって任意依存のSciPyが必要です。

### 再開で分かれた出力を読む

```python
from salmonpy.analysis import read_table
from salmonpy.outputs import read_spin_dynamics

current = read_table(['job1/rt_rt.data', 'job2/rt_rt.data'])
spin = read_spin_dynamics(['job1/rt_rt_spin.data', 'job2/rt_rt_spin.data'])
```

パスを手元の既存ファイルへ置き換えます。これは後処理だけで、計算を再実行しません。ファイルは計算順で指定し、絶対時刻を保持します。完全一致する重複だけを除き、不一致は通知します。欠測区間は補間しません。電流とスピンの時刻列は別々に保持します。単位・電流の符号は自動変換しません。SALMONの`Jm`から電荷電流`j`を作る場合は`j=-Jm`です。

### 経過時間で保存して終了する

既存のループと保存関数がある場合の差し込み位置です。

```python
from salmonpy.checkpointing import CheckpointSchedule

schedule = CheckpointSchedule(stop_after_seconds=3300, comm=comm)
for step in range(restored_step + 1, final_step + 1):
    advance_one_step()
    if schedule.after_step(step, save_all_states):
        break
```

`comm`は結合計算全体のコミュニケータです。直列では省略します。`advance_one_step()`と`save_all_states(step)`は既存の処理へ結びます。場・LD・電子状態を同じ完了ステップでまとめて保存し、新しい保存先を使ってください。ジョブ制限時間から1ステップ・保存・終了処理の時間を引いて設定します。単一状態のスナップショットだけでは結合計算の再開データになりません。

### AIに追加機能を頼む

「既存の計算に機能を追加」と検索し、必要な機能のタスクIDを選んで詳細だけを読みます。たとえば「このGSを再利用して円偏光の電流を計算。使えるのは4コア・4GiB」のように、既存データと変更点を伝えます。既に分かっている環境・メモリ・CPU条件を再質問する必要はありません。既存データが対応するか不明な場合だけ、必要な項目を確認します。

<a id="first-test"></a>
## 最初のテストを終える

```sh
python -m salmonpy.quickstart
```

対応するwheel導入後、この一コマンドでLD・合成SBEの短い計算を各1スレッドで完了します。
追加入力は不要です。`SALMONpy first calculations: PASS` と表示された保存先を確認してください。
ログと `summary.json` が残ります。電子計算・FDTDの追加プログラムは、初回テストとは別に準備します。

対応するSALMONpy wheelを導入したPython環境で実行します。ソースツリーでも実行できますが、
[Si固定入力例](../SI_INTRO.md)は現行ランタイムで再検証済みです。対応する更新wheelと新しいApple Silicon用ZIPを配布者から受け取った場合に進みます（公開Releaseは未作成）。

LD/SBEにはビルド済みの計算用共有ライブラリが必要です。import時にコンパイルは行いません。
下記の `python` はSALMONpyをインストールした環境のPythonを指します。

まず、使ってよい総メモリ・CPUコア数・目的を決めてください。以下の格子や時間刻みは
APIの使い方を示す例です。メモリ上限や物理量の収束を保証するものではありません。
出力先は作業用ディレクトリを指定してください。生成したデータをコードツリーへ保存しません。

LDと内蔵SBEはwheel同梱のFortran/Cの共有ライブラリを使います。波形の作成・解析はソルバー不要です。
`FieldSamples.write_salmon` の出力をSALMONで使う場合は `ae_shape1='input'`、`unit_system='au'` を指定します。

## まずこの3本

1. **LD単独応答**：`python -m salmonpy.examples.ld_response`。wheelだけで材料応答を試します。
2. **SBE単独応答**：`python -m salmonpy.examples.sbe_response new-sbe-output`。合成模型でA→jと独立状態を試します。実材料の物性例ではありません。
3. **SBEと場の結合**：`python -m salmonpy.examples.field_coupling new-field-output`。追加の計算プログラムなしで、合成模型を場と結合します。

初めて使う場合はこの順で十分です。以下の一覧は目的別の追加例、MPIで計算を結合する例は応用向けです。出力先にはコードツリー外の新しい場所を使ってください。

| 目的 | サンプル | 追加で必要なもの |
|---|---|---|
| 自作の計算コードの保存対応・互換性検査と独立再開 | [custom_checkpoint.py](custom_checkpoint.py) | 新規出力先。外部ソルバー不要 |
| LDの時間応答・誘電率 | [ld_response.py](ld_response.py) | wheel内のLDライブラリ |
| パルス作成・スペクトル | [pulse_spectrum.py](pulse_spectrum.py) | なし |
| 形状・材料ID・PMLとの接触確認 | [geometry_preview.py](geometry_preview.py) | 新規.npz保存先 |
| FDTDの誘電体光伝播 | [fdtd_dielectric.py](fdtd_dielectric.py) | 対応SALMON実行ファイル |
| SBEのA→J・独立状態 | [sbe_response.py](sbe_response.py) | wheel内のSBEライブラリ。合成模型 |
| 場の部品（tr／lo／1次元伝播）とSBEの結合 | [field_coupling.py](field_coupling.py) | wheel内のSBEライブラリ。合成模型、新規出力先 |
| SBEの初期状態読み込み | [sbe_from_salmon.py](sbe_from_salmon.py) | SALMONテキスト出力・格子ベクトル |
| DFTのGS・エネルギー・力 | [dft_si.py](dft_si.py) | 対応SALMON・Siの4価電子PZ/KY擬ポテンシャル |
| GS再利用・TDDFTのA→J | [tddft_response.py](tddft_response.py) | 対応SALMON・保存済みGS |
| GSからバンド計算・描画 | [bands_from_gs.py](bands_from_gs.py) | 対応SALMON・立方Si GS・Matplotlib |
| 既存band.dataの描画 | [plot_bands.py](plot_bands.py) | 指定形式の出力・Matplotlib |
| 配列の保存・読み込み | [outputs_roundtrip.py](outputs_roundtrip.py) | wheel内のLDライブラリ |
| ASE経路からバンド図 | [ase_bands.py](ase_bands.py) | 対応SALMON・立方Si GS・ASE・f90nml・Matplotlib |
| ASE Calculatorのenergy/forces | [ase_energy_forces.py](ase_energy_forces.py) | 対応SALMON・Si PZ/KY擬ポテンシャル・ASE |

## 追加データなしで試す

```sh
python -m salmonpy.examples.ld_response
python -m salmonpy.examples.pulse_spectrum
python -m salmonpy.examples.geometry_preview material-grid.npz
python -m salmonpy.examples.sbe_response new-sbe-output
python -m salmonpy.examples.field_coupling new-field-coupling-output
python -m salmonpy.examples.outputs_roundtrip new-ld-output
```

SBEの合成模型はAPI学習用で、実材料のスペクトルを表しません。結晶への入力はA/c、
出力は電荷電流密度 j です（SALMONの電子物質電流は −j）。LDはE駆動です。
格子ベクトルはDFTでは行、SBEでは列です。各コード内の単位指定を確認してください。

## 外部SALMONを使う

`EXECUTABLE`、`SHA256`、`SI_RPS`、`NEW_GS`などを自分の環境の値に置き換えます。
SHA256にはソース・オブジェクトを監査した実行ファイルの記録値を指定します。
準備方法と対応範囲は [DFT.md](../DFT.md)、[WORKFLOWS.md](../WORKFLOWS.md)、
[PACKAGING.md](../PACKAGING.md) を参照してください。wheelには完全なSALMON本体は含まれません。

```sh
python -m salmonpy.examples.dft_si EXECUTABLE SHA256 SI_RPS NEW_GS --threads 1
python -m salmonpy.examples.tddft_response EXECUTABLE SHA256 NEW_GS NEW_RT --threads 1
python -m salmonpy.examples.bands_from_gs EXECUTABLE SHA256 NEW_GS NEW_BANDS --threads 1
python -m salmonpy.examples.ase_bands EXECUTABLE SHA256 NEW_GS NEW_ASE_BANDS --threads 1
python -m salmonpy.examples.ase_energy_forces EXECUTABLE SHA256 SI_RPS NEW_CALCULATIONS --threads 1
python -m salmonpy.examples.fdtd_dielectric RUNTIME NEW_FDTD
```

各引数付きサンプルの `--help` で必要な入力を確認できます。空白を含むパスは引用符で囲みます。
全コアを自動使用しません。BLASなどのスレッド数も利用環境で調整し、許可された総並列数に
収めてください。DFT/バンドとプロセス分離型 `salmonpy.tddft` のソルバー起動は単一MPIランクです。
内部MPIのTDDFTは別の `mpi_coupled_tddft` 入口を使います。
TDDFT例は同じGSから2回を順次実行し、励起履歴を共有しません。

実材料SBEローダーの例：

```sh
python -m salmonpy.examples.sbe_from_salmon BUNDLE SYSNAME NEW_SBE --bands 32 --kpoints 8 --electrons 32 --lattice lattice-columns-bohr.txt
```

この数値は例です。実際の全バンド数・k点数・電子数・格子に合わせてください。
既存band.data描画例は17物理k点・24バンド・32電子の8原子立方Si専用です。
他の出力では点数やエネルギー基準をコード内で変更します。
力の利用ではegg-box効果と格子・k点・SCF収束を確認してください。Calculator例は
自動的に最適化を開始しません。詳細は [ASE_CALCULATOR.md](../ASE_CALCULATOR.md)。

公式入力例の参照先はAIマニュアルの `sources` に登録したSALMON-inputsです。
必要ファイルとバージョンを確認し、サンプルの物理条件と一致させて使用してください。

SALMON VR用書き出し: [salmon_vr_export.py](salmon_vr_export.py)。詳細は [SALMON_VR.md](../SALMON_VR.md)。

## 多材料・MPIサンプル

大文字の引数は利用環境のパスへ置き換えます。`RUNTIME` は管理済み**直列**計算プログラム、
`MPI_BUILD_JSON` はMPIビルド完了時の `build.json` です。相互に代用できません。
`GS_BUNDLE` はRT入力・擬ポテンシャル・再開データを含む自己完結したディレクトリです。
`RT_INPUT` はその中の入力ファイル名です。保存GSと汎関数・格子などを一致させてください。

| サンプル | 配置・資源 | 出力 |
|---|---|---|
| [fdtd_layers.py](fdtd_layers.py) | 直列FDTD、人工的なLD二層、1スレッド | 観測ログ、入力、実行記録 |
| [fdtd_mpi_layers.py](fdtd_mpi_layers.py) | FDTD 4ランク、各1スレッド、LD三材料・計768物質セル | 観測ログ、入力、MPI分割を含む実行記録 |
| [mpi_response.py](mpi_response.py) | 場役1＋物質群2×2＝5ランク、一群一物質・一点 | 点ID 10/20のA/cと電荷電流密度JのJSON |
| [mpi_tddft.py](mpi_tddft.py) | 場役1＋物質群2×1＝3ランク、直列TDDFTワーカー2個 | 点ID 1/2のA/cと電荷電流密度JのJSON |

```sh
python -m salmonpy.examples.fdtd_layers RUNTIME NEW_LD_DIRECTORY
python -m salmonpy.examples.fdtd_mpi_layers MPI_BUILD_JSON NEW_MPI_DIRECTORY --grid 2 2 1 --threads 1
mpiexec -n 5 python -m salmonpy.examples.mpi_response NEW_CUSTOM_JSON
mpiexec -n 3 python -m salmonpy.examples.mpi_tddft RUNTIME GS_BUNDLE RT_INPUT NEW_TDDFT_JSON
```

FDTD例は**新しい直列Pythonから**起動します。後二つはPython自体をMPI起動します。
`mpi_response` の `MyMaterial.prepare` を自分のローダーに、`MyResponse.create_state` を自分の
応答計算に置き換えられます。群内の通信には渡された `solver_comm` を使います。
`mpi_tddft` の二点は同じGSから作りますが、励起履歴は共有しません。
いずれもメモリ上限の自動検査はありません。全状態・ワーカー分のメモリとBLASスレッドも予算に含めます。

FDTDの `salmonpy-run.json` で終了状態・ランク数・分割・実行ファイルを確認します。
応答JSONの `history` は原子単位の `time`, `A_over_c`, `J_charge` で、行順は `point_ids` です。
これは後処理用データで、波動関数や全内部状態を含む再開ファイルではありません。

LD例は有限サイズの人工構造、独自応答例は人工的な履歴則です。TDDFT例も短いAランプの接続確認で、
物理的に収束した多層膜の計算ではありません。後二つはFDTD伝播を行わず、TDDFT内部MPIも未対応です。
対応範囲は[材料と界面](../FDTD_MATERIALS.md)と[二階層MPI](../DISTRIBUTED.md)を確認してください。

## E/H FDTDと内部MPI TDDFTを結合する

[mpi_coupled_tddft.py](mpi_coupled_tddft.py) は場2ランク＋物質2群×2ランクの6ランク例です。
[呼び出し元プログラム共有ライブラリの準備・時刻と状態の仕様](../COUPLED.md)に従います。

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 mpiexec -n 6 python -m salmonpy.examples.mpi_coupled_tddft BUILD_JSON GS_BUNDLE RT_INPUT NEW_OUTPUT
```

`response.json` は二点のA/E/Jと二つの時計、場ディレクトリは実行記録と観測ログを保持します。
一プロセス一回の計算処理生成に限定され、別の試行は新しいMPIジョブで開始します。

TDDFTでは実用上1次元の伝播・層構造を推奨し、電子応答を結合する3次元はSBEまたはユーザー作成の電流計算コードを推奨します。この40³・二点の例は接続回帰用で、1次元化や3次元TDDFTの標準的な研究条件を示すものではありません。

## ASEからSiの電流を計算

`python -m salmonpy.examples.ase_tddft RUNTIME_DIRECTORY NEW_OUTPUT_DIRECTORY --verify-independence`

ASE・f90nmlと、DFT/TDDFT/PZ対応の監査済み計算プログラムが必要です。[ASEガイド](../ASE.md)を参照。中性電子数を自動決定し、短いSi接続試験を1スレッドで実行します。光学物性の収束例ではありません。

## snapshotの再実行

```sh
python -m salmonpy.examples.snapshot_restore
```

LD・SBEと4種類のPythonの場を検査します。導入済み計算ルーチンが必要で、SALMON本体は不要です。[保存と復元](../SNAPSHOTS.md)を参照してください。

自作A駆動モデルは直列・MPIで同じ仕様を使えます。[作り方と適合チェッカー](../MATERIALS.md)、`python -m salmonpy.examples.custom_material`を参照してください。

物質ID・自作の計算コードの設定検証は [material_config.py](material_config.py)。追加データ不要で、設定だけを検証します。

See [CONFIG.en.md](../CONFIG.en.md) / [CONFIG.md](../CONFIG.md).

入力ファイルを転送用にまとめる例：[prepare_inputs.py](prepare_inputs.py)。[手順](../PREPARATION.md)に従い、ローカルの合成係数JSONと新規出力先を指定します。

構造ファイル・物質ID・幅指定から転送用入力を作る例：[prepare_structure.py](prepare_structure.py)。[手順](../STRUCTURE.md)。ASEが必要です。

物質ID付き設定からLD二物質を実行：[configured_ld.py](configured_ld.py)。新規出力先と同梱計算処理ライブラリが必要です。

物質設定からSiのDFT→TDDFTを実行：[configured_si.py](configured_si.py)。[手順](../DFT_CONFIG.md)。既存計算プログラムとASE/f90nmlが必要です。

## 追加の接続・入力サンプル

| サンプル | 目的 | 必要条件・引数 |
|---|---|---|
| [field_inputs.py](field_inputs.py) | A/E/H入力の仕様と自作応答 | LD kernel |
| [point_ld.py](point_ld.py) | 材料IDと点IDを保持したLD計算 | LD kernel; new output |
| [point_tddft_restart.py](point_tddft_restart.py) | 点ごとのTDDFT再開 | Catalog, runtime, waveform; new output |
| [catalog_tddft.py](catalog_tddft.py) | カタログからTDDFTを開始 | Catalog, runtime, waveform; see --help |
| [mpi_configured_native.py](mpi_configured_native.py) | 設定済みLD／TDDFTのMPI応答 | 3 ranks; catalog, waveform; TDDFT runtime |
| [mpi_recorded_material.py](mpi_recorded_material.py) | 記録付き自作物質／SBEの直列・MPI比較 | 3 ranks; optional SBE catalog |
| [mpi_point_identity.py](mpi_point_identity.py) | MPIで点IDと材料情報を保持 | 3 ranks; new output |
| [mpi_multiscale_tddft.py](mpi_multiscale_tddft.py) | 1次元場と独立TDDFT群を結合 | 3 ranks with --material-ranks 1; 計算処理 build and GS |
| [isolated_response.py](isolated_response.py) | 孤立系のE駆動・双極子出力 | Validated E-worker, receipt, GS bundle, E NPZ |
| [prepare_grid_field.py](prepare_grid_field.py) | FDTDのE/Hを分子入力にする（準備のみ） | E/HのNPZ、パッチの格子情報、分子位置 |
| [fit_near_field.py](fit_near_field.py) | 空間依存E/Hを展開で圧縮する（準備のみ） | E/HのNPZ、中心、次数、任意の減衰率 |
| [dipole_radiation.py](dipole_radiation.py) | 双極子からFDTDへ渡す電流を作成 | No data, or --archive and --cell-volume; no FDTD propagation |

MPIサンプルも通常のPythonから `--help` を表示できます。実計算には mpi4py と互換MPIが必要です。MPIの既定ランク数が予算を超える場合は、そのまま実行せず、対応するランク数指定を確認してください。

[mpi_field_inputs.py](mpi_field_inputs.py): 3ランクで異なるA/E/H仕様を持つ合成物質を計算。入力不一致の更新前拒否も確認します。

[field_sampling.py](field_sampling.py): 合成場の時刻差と明示的な空間重みを使う補間。外部データ・ソルバー不要。

[mpi_packet_coupling.py](mpi_packet_coupling.py): 合成物質と1D／計算処理場のpacket往復結合。3ランク（計算処理は4ランクも可）。

- `molecular_field_selection.py`：分子自身の真空場だけを引く電場選択。人工標本のみ、ソルバー不要。

大きな並列テスト: `mpi_mixed_native` は `--field-shape` に対応します。5〜15ランクの構成と計測の範囲は [COUPLED.md](../COUPLED.md) を参照してください。4GiB・4コアの入門条件は引き続き利用できます。

収束試験の計画のみ：`python -m salmonpy.examples.convergence_plan NEW_JSON`。CONVERGENCE.md参照。例のPP半径は説明用で、物質計算は行いません。

## FDTDはこの順で

1. `fdtd_dielectric RUNTIME NEW_OUTPUT`：検証済みバンドルを`load_runtime`で開き、誘電体の伝播を確認します。
2. `fdtd_layers RUNTIME NEW_OUTPUT`：LD二層と材料ID非依存の界面を試します。
3. `fdtd_mpi_layers MPI_BUILD_JSON NEW_OUTPUT`：MPI対応の監査済みビルド記録で空間並列へ進みます。

`RUNTIME`は実行ファイル単体ではなく`runtime.json`のあるディレクトリです。既存の`load_runtime`を使うアプリもそのまま利用できます。MPI例は別途能力を確認したビルド記録を使う応用例です。いずれも新しい直列Pythonから開始し、既存の出力先を上書きしません。

`fdtd_periodic_slab RUNTIME NEW_OUTPUT [--ld]` は横方向周期境界と格子に沿ったスラブの例です。[材料と界面](../FDTD_MATERIALS.md)に、Cの解析式を使うFresnel比較試験を説明しています。

### TB-mBJの直列GS入口

`DFT(..., xc="tbmbj", cval=None)` または `dft_from_atoms(..., xc="tbmbj", cval=1.2)` を使います。JSONの物質レシピでも `xc` と任意の `cval` を指定できます。係数はGS記録、warm startの整合性検査、TDDFT入力へ引き継ぎます。TB-mBJはモデルポテンシャルのため、全エネルギーと力は未検証の診断値として警告し、計算は続行します。ASE Calculatorによる最適化の対象にはしません。サンプルは `python -m salmonpy.examples.dft_tbmbj --help` で確認できます。今回の確認は合成密度の計算処理XC評価とインターフェース検査であり、実材料の物理精度保証ではありません。

## 誘電関数・光学伝導率（単一偏光）

`salmonpy.analysis.optical_conductivity` は同じ時刻の全電場Eと電荷電流密度jからsigma=j/Eを求めます。結晶A駆動で得た電流にも使えますが、AそのものはEではありません。呼び出し側で時間配置と微分規約を揃えたEを用意してください。SALMON物質電流Jを使う場合はj=-Jへ明示的に変換します。単位はGaussian原子単位のみです。完全なテンソルには独立な駆動方向の計算が必要で、円偏光一回から全テンソルを復元しません。

```python
from salmonpy.analysis import optical_conductivity, dielectric_function
sigma = optical_conductivity(time, charge_current_x, total_E_x,
    current_convention="charge", field_kind="total", units="atomic_gaussian",
    window="hann", drive_rtol=1e-8,
    metadata={"material_id": "Si", "polarization": "x"})
epsilon = dielectric_function(sigma)
# sigma.omega, sigma.values, sigma.valid; epsilonも同じ構造
```

フーリエ変換は既存の `spectrum` と同じexp(+i omega t)、dt倍の和です。時間依存exp(-i omega t)に対応し、epsilon=1+4pi*i*sigma/omegaです。弱い駆動周波数（最大振幅のdrive_rtol以下）と数値変換不能点はvalid=False、値はNaNとなります。DCの誘電関数は未定義です。全駆動ゼロでも解析結果をマスクするだけで、時間発展を停止しません。外場だけの比をバルク誘電関数と呼び替えません。電場と電流には同じ窓を使い、窓・時間長・開始時刻・dt・電場振幅・利用者の来歴をconfigに保持します。窓付き有限時間比は近似で、特に単色駆動以外では窓依存性を確認してください。

`optical_linearity(reference, candidate, rtol=...)` は利用者が渡した異なる振幅の伝導率を共通の有効周波数で比較し、相対RMSを返します。材料・偏光・駆動形状の一致は利用者が確認してください。電流が全てゼロなら相対誤差を定義せず判定も未定義とします。`optical_sum_rule(sigma, expected_weight=..., normalization="...")` は有効な隣接周波数間のRe sigmaを台形積分します。正周波数の全電子系では期待値pi*n/2が基準ですが、擬ポテンシャル・非局所電流・不足バンドに応じた基準は利用者が指定します。欠損帯域を補外せず、物理精度の合格証明にはしません。追加のSCF/TDDFT計算やダウンロードは起動しません。

合成データのサンプル: `python -m salmonpy.examples.optical_response`。実材料の収束済み光学定数の例ではありません。

制約付きLDフィット: `python -m salmonpy.examples.ld_fit --poles 1 2 3`。合成データで極数ごとの残差・精度判定・刻みの目安と来歴を比較します。

実験光学定数の一覧: `python -m salmonpy.examples.optical_materials`。CC0元データ・出典・測定条件と候補係数をオフラインで読みます。

`python -m salmonpy.examples.ld_response --anisotropic`: 合成相反テンソルLDの非対角応答。FDTD伝播例ではありません。

`fdtd_waveplate`: 四分の一／半波長板の合成FDTD例。既定はプレビューのみ、--runで明示開始。FDTD_MATERIALS.md参照。

界面平滑化を試す場合は `fdtd_waveplate` の `--subpixel` を任意指定します。中心周波数に合わせるLD専用近似です。既定は無効で、SBE・TDDFTには適用しません。まず表示のみで確認し、計算には `--run` を明示してください。

保存済み軌跡の時間長・窓依存とk点電流寄与は `optical_postprocess`。表示が既定で、`--run` は後処理だけを実行します。ソルバーを起動しません。OPTICAL_TRAJECTORIES.mdを参照してください。

`mpi_runtime_receipt` はジョブ内の検証I/Oだけを測るサンプルです。表示が既定で、`--run` でもSALMON計算を開始しません。RUNTIME.mdを参照してください。

混在MPIのA／場パケットは担当点だけへ配ります。Jの場への返却順は同じです。通信の情報はcommunication_reportで確認できます。

担当領域へのJ配信に対応する計算プログラムでは、通常のmpi_mixed_native例も場の空間分割に応じてJを配信します。代表ランク以外の結合J配列は担当行だけが有効です。全点の後処理には代表の保存出力を使ってください。詳しくはCOUPLED.mdを参照してください。

準備だけ: `python -m salmonpy.examples.prepare_fdtd NEW_DIRECTORY --mpi-grid 2 1 1` はFDTD入力保存で終了し、計算・投入はしません。PREPARATION.mdを参照してください。

電子計算の準備だけ: `python -m salmonpy.examples.prepare_electrons dft Si_rps.dat NEW_DIRECTORY` または `tddft rt.inp NEW_DIRECTORY --asset pseudo1_rps.dat Si_rps.dat`。計算・GS探索・ジョブ投入はしません。PREPARATION.mdを参照してください。

`fdtd_energy`: 保存済みマクロ点E/j・分離済みE/Hポートから吸収仕事、R/T、収支を診断（ソルバー起動なし）。

`fdtd_periodic_slab --native-planes`: 各MPIプロセスの担当領域から横断面を取得・位置合わせして流束を表示。対応計算プログラムが必要で、未対応なら計算結果を保持して警告します。

混在計算処理サンプルは任意の`--planes-x FRONT BACK --plane-stride 16`で横周期境界の面を確定後だけ取得し、最後にplane-seriesへまとめます。既定は取得しません。

## 分子とFDTDの双方向結合

`molecular_feedback.py` は既存GSを再利用する完全なMPI例です。必要な計算プログラムと実行方法は [分子とFDTD](../DIPOLE_RADIATION.md) を参照してください。SCFやビルドは暗黙に実行しません。

`isolated_response.py`は`--stop-after 7 --direct-checkpoint NEW`で分子を途中保存できます。同じ入力で別途`--resume-direct NEW`を実行すると残りのEだけを計算します。最大ステップ数と実行ファイルは同一にし、revision5の検証記録を使います。

分子と場の一括再開：`molecular_feedback.py --stop-after 30 --checkpoint NEW`で保存し、同じ入力・容量・MPI数の別ジョブから`--resume NEW`で再開します。各回の出力先は別にします。


## 電流・磁化・スピン流を保存済みSOC出力から読む

```sh
python -m salmonpy.examples.outputs_roundtrip NEW_OUTPUT --rt-data material_rt.data --spin-data material_rt_spin.data
```

[詳細と単位](../OUTPUTS.md)。この入口は後処理だけで計算を起動しません。

- `python -m salmonpy.examples.grid_electric_response INPUT.npz NEW_OUTPUT.npz --dt DT`: 保存した電荷密度から電場経路に沿った分極・電流を準備（双方向実行ではない）

- `python -m salmonpy.examples.molecular_charge_density --help`: 対応検証済みの分子TDDFTから密度を明示取得。恒久保存は `--save-density` で選択

- isolated_response の --grid-series は保存したE/Hを逐次渡します。初期入力の幾何条件と直接保存・再開を使います

- 分子の空間依存E/Hと磁気電流の自己真空参照結合は、ISOLATED_TDDFT.md の `test_grid_feedback_mpi.py --orbital-magnetic` を参照してください。固定C₂H₂の接続・再開検証で、物理的収束を保証するものではありません。

- `molecular_grid_response.py`: 保存した空間依存FDTD E/Hから分子だけを励起。直接保存／再開に対応し、電流はFDTDへ返さない。
