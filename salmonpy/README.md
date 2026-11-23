# SALMONpy：Pythonで物質応答・電子の動き・光の伝播を計算する

[English](README.en.md)

物質の光応答、電子の動き、光の伝播を単独で計算し、文書化した接続仕様で組み合わせられます。自作の電流計算コードの接続、構造・光パルスの準備、既存結果の解析にも使えます。現在は開発版で、PyPIには未公開です。

| 使い方 | できること | 準備 |
|---|---|---|
| Pythonから材料応答を計算する | LD・SBE、パルス・構造・解析 | 配布用wheel（Pythonパッケージ）を導入 |
| SALMONを別プロセスで実行する | DFT・TDDFT・FDTD・バンド計算 | 対応するSALMON実行ファイルを追加 |
| Pythonから一ステップずつ計算を進める | 電子計算とFDTDの結合、自作コードとの連携 | 対応する共有ライブラリとMPIを追加 |

初めての場合は最初の行から進めます。追加の計算プログラムは[準備手順](RUNTIME.md)で導入します。対応するビルド済み配布物を優先し、ない場合だけ同梱ソースを明示的にビルドします。別のSALMONリポジトリを取得する必要はありません。

## 最初の計算：二つのコマンド

OS・CPU・Python版に合う配布用wheelを用意します（Python 3.10以上）。
現在PyPI未公開なので、下のファイル名は実際に受け取ったwheelへ置き換えます。
ファイルがなければ配布担当者に確認してください。Fortranコンパイラは不要です。

```sh
python -m pip install --only-binary=:all: "PATH_TO_COMPATIBLE_WHEEL.whl"
python -m salmonpy.quickstart
```

`SALMONpy first calculations: PASS` と出れば、二つの計算が完了です。
ローレンツ＝ドルーデ（LD）模型を1000更新、合成2バンド模型の半導体ブロッホ方程式（SBE）を500更新します。
各1スレッドで実行し、追加入力・MPI・SALMON実行ファイル・描画ソフトは不要です。計算中の通信もありません。

表示された保存先には `summary.json`、二つの計算ログ、SBEの結果が入ります。
通常は数秒程度ですが、パッケージの取得・導入時間は環境によります。
保存場所を選ぶ場合は `--output NEW_DIRECTORY` を付けます。既存フォルダを上書きしません。
これは導入確認であり、実材料の物性・収束を確認する計算ではありません。
失敗時はログを読み、別の計算プログラムの構築や導入を自動で始めないでください。

次は下の目的別一覧か[サンプル集](examples/README.md)を使います。
基底状態や電子計算へ進む場合は[DFT](DFT.md)と[計算プログラムの準備](RUNTIME.md)を参照してください。[Si固定入力例](SI_INTRO.md)は現行ランタイムで再検証済みです。新しい対応wheel・ZIPを配布者から受け取って利用します（新ZIPの公開Releaseは未作成）。擬ポテンシャルの個別取得は[こちら](PSEUDOPOTENTIALS.md)。
導入状態だけを調べる場合は `python -m salmonpy.doctor --json` を使います。

## 目的から選ぶ

A/cは光のベクトルポテンシャルを光速で割った量、Eは電場、Hは磁場強度です。GSは、時間発展の初期条件として保存する基底状態を指します。

| やりたいこと | 入口 | 必要な準備 |
|---|---|---|
| 材料の誘電率・Eに対するP/J | `salmonpy.LorentzDrude` | wheel |
| 結晶のA/c→j応答（SBE） | `salmonpy.sbe` | wheel、実材料にはバンド・遷移行列データ |
| 場（j→A）と物質の結合：横波／縦波／薄膜／1次元伝播 | `salmonpy.maxwell` | NumPy、物質側の部品。[FIELDS.md](FIELDS.md) |
| 光パルス・構造・解析・結果読み込み | `fields` / `geometry` / `analysis` / `outputs` | NumPy、描画にはMatplotlib |
| DFTの基底状態・力 | `salmonpy.dft` | 対応する計算プログラム・擬ポテンシャル |
| 結晶のA/c→j応答（TDDFT） | `salmonpy.tddft` | 対応する計算プログラム・整合する基底状態 |
| 光伝播（FDTD） | `salmonpy.fdtd` | 対応する計算プログラム |
| バンド計算・描画 | `salmonpy.bands` | 対応する計算プログラム・基底状態・Matplotlib |
| ASEの構造・経路・energy/forces | `salmonpy.ase` / `salmonpy.ase_calculator` | 上記に加え任意のASE依存 |

[Pythonサンプルと実行方法](examples/README.md)から始めてください。
各機能の対応範囲は [WORKFLOWS.md](WORKFLOWS.md)、[DFT.md](DFT.md)、
[ASE_CALCULATOR.md](ASE_CALCULATOR.md) に記載しています。
wheel内のLD/SBEと、別途準備が必要なSALMON本体は区別してください。

## 並列計算の選び方

以下のランク数・点数はサンプルの構成であり、パッケージ全体の上限ではありません。実際に使えるMPI配置は計算部品ごとの対応条件で確認します。
LD・SBE・TDDFTを混ぜて光伝播を計算したい場合は、まず[通常混在サンプルの手順](COUPLED.md)と`salmonpy.examples.mpi_mixed_native`を参照してください。材料ID・並び・厚み・各電子計算の並列数を指定して組み立てます。入力準備と実行は分け、実計算は計算ノードで起動します。

| 目的 | サンプルの起動方法 | このサンプルの構成 |
|---|---|---|
| 多材料LDの光伝播 | 直列Pythonから `fdtd_mpi_layers` | FDTDが空間を分担し、各MPIプロセスが担当領域のLD状態を管理 |
| 自作A→Jモデルの群内並列 | `mpiexec -n 5` で `mpi_response` | 二つの物質群、各2ランク・一物質・一マクロ点 |
| 保存GSから二点のTDDFT応答 | `mpiexec -n 3` で `mpi_tddft` | 二つの1ランク物質群、各点に独立した直列ワーカー |

[材料とMPI FDTD](FDTD_MATERIALS.md)、[応答サービス](DISTRIBUTED.md)、
[実行コマンドと出力](examples/README.md)を参照してください。後二者の場役は所定のAを配ります。
TDDFT内部MPIとE/H結合には、別入口の[実験的なPythonから直接呼ぶFDTDと電子計算の結合](COUPLED.md)を使います。
標準例は `mpi_coupled_tddft`（場2＋物質2群×2＝6ランク）です。
バッチFDTDの起動をMPIジョブ内で入れ子に呼びません。

## AIに依頼する場合

まず [AI_USAGE.md](AI_USAGE.md) と下記の検索コマンドを使い、既に回答された条件を再利用して、
研究計算の前に不足する「目的・使ってよい総メモリ・CPUコア数」を確認します。最初の合成模型テストには追加の選択は不要です。

```sh
python -m salmonpy.manual
python -m salmonpy.manual レーザーによる電子の動き
python -m salmonpy.doctor --json
```

検索は上位5件の概要と開始案内を返します。タスクIDを選んで`python -m salmonpy.manual --task TASK_ID`で詳細を読みます。一覧は検索結果の`task_catalog`、表示件数の変更は`--limit N`を使います。Pythonでは`search`で概要、`task_details`で一つのタスクを取得できます。
`sample_files` が実行例、`beginner_intake` が資源上限の確認手順です。
環境診断で検出したCPU数を、そのまま利用許可と解釈しません。
対応する計算ではメモリ見積もりと実測を利用できます。見積もりが上限を超えた場合は警告し、自動では停止しません。未測定の条件では上限を保証しません。

外部SALMONを確認する場合は、監査済みの記録値を指定します。

```sh
python -m salmonpy.doctor --executable PATH_TO_SALMON --sha256 VALIDATED_SHA256
```

ハッシュ一致はファイルの同一性だけを示します。対応機能・バイナリ互換性や実行成功は
検証しません。準備・OS別検証状況は [PACKAGING.md](PACKAGING.md) を参照してください。

## LDのPython API

```python
import numpy as np
from salmonpy import LorentzDrude

material = LorentzDrude(
    plasma_frequency=0.5,
    strength=[0.8, 0.2],
    resonance=[0.0, 0.3],
    damping=[0.01, 0.02],
)
E = np.zeros((101, 3))
E[:, 0] = 1e-3
result = material.respond(E, dt=0.1)
# result.time: (101,), result.polarization/current: (101, 3)

with material.new_state(dt=0.1, initial_field=E[0]) as state:
    P, J = state.step(E[1])  # fields at t=0.1; physical P and J
```

値はSALMONの原子単位・Gaussian単位系で渡します。角周波数・減衰率・プラズマ周波数は原子時間の逆数、時間は原子時間、E・P・Jも対応する原子単位です。eVやfsで指定した値は換算が必要です。振動子の強度は無次元で、各振動子は次の式に従います。

```
dP/dt = J
dJ/dt = f * wp**2 * E / (4*pi) - gamma*J - w0**2*P
```

`respond` は最初のEを初期条件とし、P=J=0から計算します。最初のサンプルでは時間を進めません。`step` には次の時刻のEを渡し、内部で保持した前のEと合わせて更新します。dtは各状態で固定です。`respond` を呼ぶたびに独立した状態を作り、別々の `new_state` の計算履歴は混ざりません。

PとJは全振動子の和です。FDTD更新用の補助量Qとは区別してください。

`permittivity(omega)` は時間依存を exp(-i omega t) とする連続モデルの誘電率です。
`epsilon_infinity + sum(f*wp**2 / (w0**2 - omega**2 - i*gamma*omega))`
で評価し、瞬時に応答する背景の `epsilon_infinity` は既定で1です。上の例は等方的な材料です。方向によって異なる応答には `salmonpy.ld.TensorLorentzDrude` を使います。フィッティングと対応する材料条件は[LDガイド](LD.md)を参照してください。有限のdtによる時間発展はこの連続モデルの近似なので、dtを小さくして結果を確認します。

計算中のP・JなどはFortran側で保持され、`step` のたびに計算状態を作り直す必要はありません。使い終わった状態は `close()` または `with` で終了します。同じ状態を複数スレッドから同時に操作しないでください。

SBE、TDDFT、FDTDの実行例は[サンプル集](examples/README.md)から選べます。結晶のTDDFTとSBEはA駆動です。物質ごとに固定する設定と呼び出しごとに渡す値の区別は[入力設定](INPUT_CONTRACT.md)を参照し、確定した設定は `material.config` と `state.config` で確認できます。

## 電子計算の対応範囲

別プロセス型DFTは周期・固定イオン・非スピン偏極のPZ、TB-mBJ、Libxc-PZ/PBEに対応します。
TB-mBJは模型ポテンシャルで、全エネルギーと力は未検証の診断値です。HSEは[別のGS計算入口](NATIVE_GS.md)を使います。
結晶の別プロセス型TDDFTは整合する基底状態と中点伝播を使います。分子のE駆動は[固定イオン・PZの別入口](ISOLATED_TDDFT.md)です。
対応するPython入口があることと、すべての条件・環境で検証済みであることを区別します。

## 実行環境と出力先

以下の `python` はSALMONpyを導入した環境のPythonです。Python 3.10以上を使い、環境によってはコマンド名または実行ファイルの絶対パスを指定します。Windowsで `python3` という名前があることは前提にしません。

サンプルはwheelに含まれ、`python -m salmonpy.examples.NAME` で実行できます。ソースのビルドや `tools/`、開発用テストを使う場合は、リポジトリのルートから実行してください。

例の実行ファイルパス・SHA256は、自分の環境で検証した値へ置き換えます。`executable_path` と `validated_sha256` は利用者が渡す値です。出力先には書き込める新しいディレクトリを使います。相対パスの基準はPythonを起動した作業ディレクトリで、パッケージのインストール先ではありません。

解析にはNumPy、描画にはMatplotlibが必要です。画面を開けない計算ノードでは `import matplotlib; matplotlib.use('Agg')` でファイル保存用の描画方式を選べます。SALMON計算プログラムの対応環境はwheelとは別に確認します。[配布・ビルド](PACKAGING.md)を参照してください。

## ソースから開発する場合

ソース構築のコンパイラ・BLAS要件は [PACKAGING.md](PACKAGING.md) を参照してください。
リポジトリのルートで明示的に構築します。再構築後はPythonを再起動します。

```sh
python -m salmonpy.build_native
python -m unittest discover -s salmonpy/tests -v
```

## AIへの相談と関連資料

`python -m salmonpy.manual --lang ja "TDDFT"` で目的に合う手順を検索できます。
`python -m salmonpy.doctor --lang ja` で実行環境を確認できます。
開発・配布方法は[開発者向け資料](../docs/development/README.md)を参照してください。
公開版SALMON VRへの接続は [SALMON_VR.md](SALMON_VR.md) を参照してください。
