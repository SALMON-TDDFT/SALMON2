# Pythonサンプル集

対応するSALMONpy wheelを導入したPython環境で実行します。ソースツリーでも実行できますが、
LD/SBEには構築済みのネイティブライブラリが必要です。import時にコンパイルは行いません。
下記の `python` はSALMONpyをインストールした環境のPythonを指します。

まず、使ってよい総メモリ・CPUコア数・目的を決めてください。以下の格子や時間刻みは
APIの使い方を示す例です。メモリ上限や物理量の収束を保証するものではありません。
出力先は作業用ディレクトリを指定してください。生成したデータをコードツリーへ保存しません。

| 目的 | サンプル | 追加で必要なもの |
|---|---|---|
| LDの時間応答・誘電率 | [ld_response.py](ld_response.py) | wheel内のLDライブラリ |
| パルス作成・スペクトル | [pulse_spectrum.py](pulse_spectrum.py) | なし |
| 形状・材料ID・PMLとの接触確認 | [geometry_preview.py](geometry_preview.py) | 新規.npz保存先 |
| FDTDの誘電体光伝播 | [fdtd_dielectric.py](fdtd_dielectric.py) | 対応SALMON実行ファイル |
| SBEのA→J・独立状態 | [sbe_response.py](sbe_response.py) | wheel内のSBEライブラリ。合成模型 |
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
python -m salmonpy.examples.outputs_roundtrip new-ld-output
```

SBEの合成模型はAPI学習用で、実材料のスペクトルを表しません。結晶への入力はA/c、
出力Jは電子物質電流です（電子の電荷電流は−J）。LDはE駆動です。
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
収めてください。現在のPython DFT/バンド/TDDFT入口は単一MPIランクです。
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
