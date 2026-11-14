# Pythonサンプル集

[English](README.en.md) · Topic: `examples` · Paired revision: `2026-11-13.5`
対応するSALMONpy wheelを導入したPython環境で実行します。ソースツリーでも実行できますが、
LD/SBEには構築済みのネイティブライブラリが必要です。import時にコンパイルは行いません。
下記の `python` はSALMONpyをインストールした環境のPythonを指します。

まず、使ってよい総メモリ・CPUコア数・目的を決めてください。以下の格子や時間刻みは
APIの使い方を示す例です。メモリ上限や物理量の収束を保証するものではありません。
出力先は作業用ディレクトリを指定してください。生成したデータをコードツリーへ保存しません。

LDと内蔵SBEはwheel同梱のネイティブライブラリを使います。波形の作成・解析はソルバー不要です。
`FieldSamples.write_salmon` の出力をSALMONで使う場合は `ae_shape1='input'`、`unit_system='au'` を指定します。

| 目的 | サンプル | 追加で必要なもの |
|---|---|---|
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
出力は電荷電流密度 j です（SALMONの電子物質電流は −j。2026-11-14 以前は物質電流を返していました）。LDはE駆動です。
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
python -m salmonpy.examples.fdtd_dielectric EXECUTABLE SHA256 NEW_FDTD
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

大文字の引数は利用環境のパスへ置き換えます。`RUNTIME` は管理済み**直列**ランタイム、
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
[ホスト共有ライブラリの準備・時刻と状態の契約](../COUPLED.md)に従います。

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 mpiexec -n 6 python -m salmonpy.examples.mpi_coupled_tddft BUILD_JSON GS_BUNDLE RT_INPUT NEW_OUTPUT
```

`response.json` は二点のA/E/Jと二つの時計、場ディレクトリは実行記録と観測ログを保持します。
一プロセス一回のネイティブ生成に限定され、別の試行は新しいMPIジョブで開始します。

TDDFTでは実用上1次元の伝播・層構造を推奨し、電子応答を結合する3次元はSBEまたはユーザー作成の電流計算コードを推奨します。この40³・二点の例は接続回帰用で、1次元化や3次元TDDFTの標準的な研究条件を示すものではありません。
