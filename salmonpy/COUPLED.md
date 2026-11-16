# MPI TDDFTとMultiscale1D／E/H FDTDの結合（実験的）

[English](COUPLED.en.md) · Topic: `coupled` · Paired revision: `2026-11-14.1`

PythonのMPIジョブを場群と物質群へ分け、**各物質群に一物質・一マクロ点**を置きます。
物質群の全ランクがSALMONの電子計算に参加し、k点・軌道・実空間の内部分割はその群内に閉じます。
子プロセスや二重のMPI起動は使いません。LDは場群内の既存FDTDが多数の局所状態を管理します。

入口は `NativeTDDFTMaterialProvider`、`NativeTDDFTBackend`、`NativeEHField`、`CoupledStepper`。
後二つは別モジュール `salmonpy.coupled` です。既存のプロセス分離型 `TDDFTResponseBackend` と
バッチFDTDも従来どおり利用でき、それらの制約をこの新しい入口と混同しないでください。

## 次元と応答モデルの選び方

3次元の場・マクロ点配置を扱えますが、TDDFTは各マクロ点に電子状態を持つため、
実用計算では**まず1次元の伝播・層構造を推奨**します。ここでの次元は電磁場・マクロ配置の次元で、
微視的なDFT/TDDFTの電子状態を1次元化する意味ではありません。

**電子応答を結合する3次元計算は、SBEまたはユーザー作成の電流計算コードを使う場合に推奨**します。
TDDFTの3次元計算は禁止しませんが、必要な物理・独立状態数・メモリ・計算時間を確認して選びます。
SBEや独自コードも、モデルの有効範囲と各マクロ点の独立状態、A/c→電荷電流密度の契約を満たす必要があります。
これはモデル選択の推奨であり、未実装のSBE結合アダプターや自動1次元化を対応済みとするものではありません。
LD単独の3次元FDTDは別の利用経路です。

現在の40³格子・TDDFT二点のサンプルは接続回帰用です。一般利用向けの3次元TDDFT推奨例でも、
横方向に一様な1次元多層膜の実装例でもありません。格子の横方向を単に1セルへ減らすだけでは1次元化できません。

## 初回準備と資源

最初に、目的・使える総メモリ・総CPU並列数を確認します。標準例は場2ランク＋物質2群×2ランク＝6ランク、
電子状態は2個です。各ランク1スレッドで始め、OpenMPとBLASの追加スレッドを予算へ含めます。
メモリ使用量の自動見積もり・上限強制はありません。小さなGSで測定してから規模を増やします。

この入口には**実行環境のMPIで明示的に構築した共有ライブラリ**が必要です。
配布用直列ランタイムや通信検査用ブリッジは代用できません。import時にはビルドしません。
現在のビルド手順はmacOS/Linux・GNU互換コンパイラを対象とします。

```sh
python tools/prepare_salmon_runtime.py UPSTREAM_CHECKOUT NEW_SOURCE
python tools/build_embedded_tddft.py NEW_SOURCE NEW_BUILD --jobs 2 \
  --cmake-arg=-DCMAKE_Fortran_COMPILER=mpifort \
  --cmake-arg=-DCMAKE_C_COMPILER=mpicc \
  --cmake-arg=-DBLA_VENDOR=OpenBLAS \
  --cmake-arg=-DOPENMP_FLAGS=-fopenmp \
  --cmake-arg=-DCMAKE_BUILD_TYPE=Release \
  --cmake-arg='-DCMAKE_Fortran_FLAGS_RELEASE=-O2 -ffree-line-length-none -fallow-argument-mismatch'
```

MPIラッパー・BLASの所在は環境に合わせます。Cラッパーも指定したOpenMPオプションに対応する必要があります。
ビルドはMPI・PICを有効にし、ソースとオブジェクトの固有値ドライバー監査を行い、
ライブラリSHA256・ソースハッシュ・構築コマンドを `NEW_BUILD/build.json` に記録します。
wheelにも構築ツールを `salmonpy/_runtime_source/tools` に同梱します。
`importlib.resources.files('salmonpy').joinpath('_runtime_source')` で場所を調べられます。

mpi4pyとライブラリは同じMPIを使用します。共有ライブラリはホスト依存で、他のMPI環境へコピーして動くとは限りません。
Libxcが必要なら対応したビルド設定を明示します。今回の内部分割の数値検証はPZとTB-mBJです。

## 実行例

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 mpiexec -n 6 python -m salmonpy.examples.mpi_coupled_tddft \
  NEW_BUILD/build.json GS_BUNDLE RT_INPUT NEW_OUTPUT
```

`GS_BUNDLE` は保存済みGS、相対パスの擬ポテンシャル、RT入力を含む自己完結した束です。
`RT_INPUT` はその中のファイル名。GS作成条件と汎関数・格子などを一致させます。
このローダーはDFTを自動実行せず、保存形式も変換しません。
`NEW_OUTPUT` は同じ計算群の全ランクから読める共有ファイルシステム上の新規保存先です。

標準例は40³の場格子と二つの隣接する電子マクロセルです。人工的な光パルスによる接続確認で、
収束した多層膜や実験を再現する条件ではありません。`--field-ranks 1 --material-ranks 2` なら計5ランクです。
各物質群の内部分割は `(k, orbital, x, y, z)`、積は群のランク数と一致させます。
標準例は軌道方向を分割します。k点・バンド数・実空間格子に対して有効な分割を選びます。

`ResponseSession` に渡す `MaterialBinding` のprovider/backendは利用者のコードへ交換できます。
初期物質の読み込みとA→Jの物理計算は別の契約です。物理的な**電荷電流密度**を返し、
電子数電流の符号変換を二重に行わないでください。群内通信は渡された `solver_comm` に限定します。

## 1ランクのMultiscale1Dと分散TDDFT

`CoupledStepper.from_multiscale(context, bindings, field=None)` はPythonの
`Multiscale1D` 状態を同じ物質群へ接続します。場群は1ランクに限定し、そのランクで
`model.new_state(dt=dt, steps=steps)` を作り、物質ランクでは `None` を渡します。
全ランクでfactoryと `step()` を呼びます。物質セル1〜Nを
`context.plan.point_ids` の順に対応付けます。同じGSを二つの物質群に渡しても、
時間発展中の状態は独立です。場1ランク＋物質2群×2ランクなら計5ランクです。

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 mpiexec -n 5 python -m salmonpy.examples.mpi_multiscale_tddft \
  NEW_BUILD/build.json GS_BUNDLE RT_INPUT NEW_OUTPUT --material-ranks 2 --steps 400
```

`--material-ranks 1` なら計3ランクで起動します。物質2セル、間隔5原子単位、
真空(300, 50)、dt=0.02原子単位、z偏光の人工パルスを使います。軌道分割は用意したGSに
適合させてください。`response.json` に端点A/j履歴と時刻宣言、`field.npz` にt=0を含む
格子全体のA履歴、`communication.json` に計測を保存します。再開ファイルではありません。
サンプルはビルドやGS準備を自動実行しません。

factoryはネイティブ状態を作る前に物質点のA(0)を検査します。非ゼロなら点と値を示して拒否します。
左側の真空セルを増やすか、入射パルスを遅らせてください。真空中の入射Aは非ゼロでも構いません。
物質点のA(0)を黙ってゼロへ置き換えません。初版は物質点の初期Aがゼロの場合に限定し、
ネイティブ状態は一プロセス一回の生成という制約を保持します。

更新順は「前回j → Multiscale1Dの場更新 → 新しいA → 物質更新」です。
初回フィードバックはゼロです。場と物質の時計はともにn·dt（offset 0）で、
E/H用の半ステップ変換は適用しません。`step()` は場ランクで点順の
`A_over_c` と `J_charge` を返します。元の場状態を保持すれば `grid_A` と `grid_E`
を読めます。後者は中心差分で時刻(n−1)·dtです。
返される各frameの `E_time` にEの時刻を明記します。Multiscale1Dの中心差分Eは
(n−1)·dt、Aはn·dtで、ネイティブE/HのEはn·dtです。配列の添字だけでEとAを同時刻と扱わないでください。
`coupled.response.initial_current` は点IDをキーとする物理的なt=0電流の辞書です。
独自状態が初期電流を公開しない場合は `None` です。この観測量と初回のゼロフィードバックは区別します。
サンプルは `response.json` の `initial_current` に保存します。

終了は `coupled.response`、存在するランクの `coupled.field`、Pythonコンテキストの順です。
時間発展後の失敗では結合を無効にし、群をまたいだ原子的な巻き戻しや再開は保証しません。

`coupled.response.communication_report()` は全ランクで呼び、ランク別の通信計測を取得します。
A配布とj回収を物質計算とは別に測ります。通信呼び出し時間には他ランクの待ち時間も入り得るため、
純粋なネットワーク遅延とは解釈しません。計測の失敗と計算の失敗を分け、取得不能な値はその旨を残します。
単一ホストの計測から複数ノードでの拡張性能を主張しません。

受け入れ検査はSALMONの二点Weylマルチスケール出力に対する格子全体のAと各点のjの比較、
直列ワーカーとの比較、ランク数を変えた一致確認です。合成GSによる検査は接続と更新順の検証で、
物理的な収束検証とは区別します。

## ネイティブE/Hの時間・空間の契約

新しい方式名は `salmonpy.eh-centered-current/1`。従来のWeyl形式のマルチスケールとは
離散化が異なるため、両者のビット一致は主張しません。通常のFDTD更新順と光源の評価時刻は維持します。

| 量 | 場の時計での時刻 |
|---|---|
| n回目更新後のE、従来ログのE/H（Hは時間平均） | n·dt |
| 新しいAとJ、更新後の生のH | (n+1/2)·dt |
| E更新へ渡す前回のJ | (n−1/2)·dt |

電子計算のローカル時刻0は、場の時刻dt/2です。初期E・H・AとフィードバックJは0に固定します。
初期GSの計算電流が非ゼロでも、勝手に初期フィードバックへ代入しません。
更新は「前回JでE/H更新 → EからA更新 → TDDFTでJ更新」の順です。
`A_new = A_old − dt * S(E_n)` で、結晶に渡す量は常にA/cです。

この表は `salmonpy.time/1` の形で宣言しています（`NativeEHField.timing`：入力 j は offset ½・lag 1・t=0 不使用、出力 A は offset ½）。
`CoupledStepper` は全点の物質側の宣言を MPI で集め、電子系の時計（場の時計 + ½、`ELECTRONIC_CLOCK_OFFSET`）に置いて照合し、
食い違いがあれば最初のステップの前に拒否します。照合結果は `stepper.timing` と `response.json` に残ります。[FIELDS.md](FIELDS.md) を参照してください。

Sは既存FDTDの観測位置への補間と同じです。例えばExは `(i,j,k)` と `(i−1,j,k)` の辺の平均、
Ey/Ezはそれぞれy/z方向の隣接辺の平均です。電流はその転置写像で、対応する二辺へ半分ずつ戻します。
各マクロ点が一様格子の一セル体積を代表するため、追加の体積係数はありません。
MPI境界では各場ランクが自分の所有する辺へ加算します。

マクロ点は境界から10セル以上離し、光源面・その周囲と、背景のLD/誘電体の構成則が使われる辺を避けます。
周囲3³セルが真空であることを検査します。LDを別の場所に置くことはできます。

## ネイティブE/Hの状態・出力・制約

- 一プロセスにつきネイティブ状態の生成は一回です。同じプロセスで再初期化・複数状態・スレッド共有はできません。
  同じGSから独立した二群の計算は可能です。別の試行は新しいMPIジョブで始めます。
- 入力・初期状態・dt・内部分割は生成時に固定。各段でAまたは前回Jだけを渡します。
- SALMONの相対I/Oに合わせ、呼び出し中だけ作業ディレクトリを切り替え、戻します。
  メインスレッド限定です。同じプロセスで並行する相対パスI/Oを行わないでください。
- 反復の途中で場または物質を単独に進めると、結合器は同期不一致を拒否します。
  ネイティブ更新後の異常では再試行・巻き戻しをしません。close後、新しいジョブで再開します。
- closeは物質応答 → 場 → Pythonコンテキストの順に全群で実行します。FortranはMPI自体を終了しません。
  ネイティブの停止・通信障害にはジョブ側のタイムアウトが必要で、サンプルは異常時にジョブ全体を停止します。
- 再開、非ゼロ初期場、インパルス光源、回転、異なるマクロ体積、時間サブサイクル、TTMはこの入口では未対応です。

`response.json` は点ID順のA/c、E、前回J、新しい電荷電流と二つの時計を保持します。
場の `salmonpy-coupled-field.json` は写像、MPI分割、ライブラリ、終了状態を記録します。
物質ごとの作業コピーとSALMON出力も残します。大きい密度などはRT入力で選択し、保存間隔を指定します。
これらのJSONは完全な再開ファイルではありません。

ローカル回帰では既存の真空/LD FDTDとの一致、三方向のMPI境界の電流注入、
場1/2ランクの一致、TDDFTのA履歴再実行を確認しています。物理量の格子・時間・GS収束と大規模性能は別途確認が必要です。

## 異なる材料と時間刻み

同じFDTDへ接続するすべての物質状態は、材料・格子・XCが異なってもFDTDと共通のdtを使います。異なるdtの独立状態を作れることと、異なるdtで場と結合できることは別です。現在の結合器は各点の時刻宣言を開始前に照合し、dt不一致を拒否します。整数比（例：物質dtが場dtの半分）も、そのまま接続できません。

Yee結合の半ステップずれと前回電流のlag1は維持します。同じdtでも半ステップの対応が違えば拒否します。サブサイクル、途中のAの補間、返す電流の時間平均は未実装で、自動適用しません。複数のdtが必要な場合は、独立した計算として扱うか、別途検証した時間結合アダプターが必要です。

## LD・SBE・TDDFTの混在に向けた入口

`mpi_mixed_multiscale`は、同じ1次元場にLDと合成SBEを配置する接続サンプルです。
```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 mpiexec -n 3 python -m salmonpy.examples.mpi_mixed_multiscale NEW_OUTPUT
```
場1ランク、物質点2ランクで動かし、非ゼロの応答と独立した直列更新との一致を検査します。点ID10/20と物質キー、A/c・電荷電流j、通信時間を保存します。LDは端点A差分からEを得ます。native E/Hの半陰的LD更新と同じ方式ではありません。

4ランクにして`--runtime RUNTIME --bundle GS_BUNDLE --input rt.inp`を指定すると、管理済み直列TDDFTワーカーを第三点に接続する入口を使えます。Si入門系の実GS（PZ、16³格子、2×2×2 k点）を用いた4ランク・100ステップの三種類混在実行は、直列計算と完全一致しました。実GSと対応ランタイムを明示的に用意してください。これは管理済み直列ワーカーの検証で、TDDFT群内部のMPI並列とは別の経路です。光学物性の収束は保証しません。記録: `docs/mixed-three-material-validation-2026-11-16.md`。

最終構成は、FDTDランク内の多数のLD状態と、別MPI群のSBE/TDDFT状態の混在です。この小規模サンプルの「LD一セルに一ランク」を本番の配置方針にしないでください。実材料の収束や三次元混在対応を保証するサンプルではありません。

## native E/H内のLDと、SBE・内部MPI TDDFT

`mpi_mixed_native`は、場1ランクが二種類のLD領域を持ち、別の1ランクがSBE、別の2ランクが1点のTDDFTを担当する接続試験です。外側4ランクでTDDFTの軌道分割を内部2ランクにします。
```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 mpiexec -n 4 python -m salmonpy.examples.mpi_mixed_native BUILD_JSON SERIAL_RUNTIME GS_BUNDLE NEW_OUTPUT
```
`--tddft-ranks 1`なら外側3ランクです。`--without-ld`はLD領域を真空にする対照条件です。BUILD_JSONは同じMPI環境で明示ビルドした`tools/build_embedded_tddft.py`の記録、SERIAL_RUNTIMEは応答再計算の比較用です。自動ビルド・ダウンロードは行いません。

LDはnative E/Hの既存の係数・補助量・更新順を使います。電子状態はA/c駆動を維持し、受け渡しは電荷電流jです。点ID10はSBE、20はTDDFT、LD領域ID1/2は場の材料格子上で管理します。電子点の周辺はLDと重複させません。半ステップずれの場時計と物質の局所時計を区別します。

この40³格子・電子2点の小さな三次元配置は接続検査用です。実用TDDFTは引き続き一次元伝播を推奨します。巨大な多数LD、分散した複数FDTDランク、全計算再開、物理収束の検証とは区別してください。ネイティブ異常は回復可能ではなく、サンプルはMPIジョブを停止します。

検証済み: 128ステップ、TDDFT内部1/2ランクの場・電流比較、直列応答再計算、LDなし対照、LD単独の従来一致、wheel導入後実行。範囲と誤差は`docs/native-mixed-material-validation-2026-11-16.md`を参照してください。

### FDTD群も並列化する

場をy方向に2分割する例（場2＋SBE1＋TDDFT1＝外側4ランク）:
```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 mpiexec -n 4 python -m salmonpy.examples.mpi_mixed_native BUILD_JSON SERIAL_RUNTIME GS_BUNDLE NEW_OUTPUT --field-grid 1 2 1 --tddft-ranks 1 --compare REFERENCE_OUTPUT
```
REFERENCE_OUTPUTには、同じサンプルを場1ランク・TDDFT1ランク（外側3ランク）で実行した出力を指定します。`--compare`を省略すれば単独実行です。比較は同じステップ数・dt・材料設定・点順を要求し、A/c・E・jの履歴を確認します。必要な外側ランク数は`product(field-grid)+1+tddft-ranks`です。

x/y/z各方向の2分割は、128ステップで場1ランクの結果と完全一致しました。y/z分割では両LD材料が両方の場区画に分かれます。場2ランクとTDDFT内部2ランクを同時に使う5ランク構成は、今回の4コア枠では未検証です。大規模性能・メモリ上限の保証ではありません。記録: `docs/distributed-field-mixed-validation-2026-11-16.md`。

## 物質ID・並び・厚みからMPI構成を準備する

`salmonpy.layers.assemble_native_layers` は `salmonpy.native-layers/1` の配置と、物質IDからLDモデルまたは `MaterialBinding` への辞書を受け取ります。物理パラメータと初期状態は既存の設定ローダーで用意し、配置ではIDを参照します。自作応答コードも明示的な `MaterialBinding` で接続できます。

```python
import json
from salmonpy.layers import assemble_native_layers

layout = json.loads(layout_path.read_text())
assembly = assemble_native_layers(
    layout, materials, grid=grid, dt=dt, steps=steps, source=source,
    field_grid=(1, 1, 1),
    ranks_per_point={"sbe-point": 1, "si-tddft": 2}, rank_budget=4,
)
report = assembly.record()
print(report["required_ranks"], report["layers"])
# parallel_plan を ParallelContext.create、calculation と points を
# NativeEHField、bindings を ResponseSession に渡します。
```

上の変数は事前に用意したモデル、Grid、PlaneWaveと計算条件です。構成の作成はソルバー・MPIの起動やダウンロードをしません。表示された必要ランク数で、環境のMPI起動方法を使います。

`examples/mixed_layers.json` はSBE、Si、真空、Drude、Lorentzの順の例です。長さは明示的にBohr。各層は一意な `id`、`material_id`（真空はnull）、`thickness`、横方向2軸の `transverse_cells` を持ちます。`axis` は光源の伝播軸と一致させ、`origin_cell` は最初の軸方向セルと横方向の基準セルです。厚みは格子幅で切り上げ、要求値・実際の厚み・終端を含まないセル範囲・物質ID・ランクを記録します。横方向の開始位置は `anchor-count//2` です。

電子物質の各格子セルに独立した状態とMPI群を割り当て、LDは場のランク内で管理します。有限の直方体配置であり、無限膜の均質化や代表点への間引きはしません。厚みや断面積を増やすと電子点と必要ランクが増えます。ランク予算超過、未知のID、バックエンド不整合、吸収境界との重なり、LD構成則の近傍に入る電子点を拒否します。真空の隙間も明示指定が必要です。ランク予算はメモリ上限ではなく、`memory_bytes` は未推定です。初期状態準備、MPI起動、物理的な収束確認は別途行います。

上記4ランクの `mpi_mixed_native` コマンドに `--layers 配置JSONのパス --rank-budget 4` を追加できます。このサンプルの既知IDは `sbe-point`、`si-tddft`、`drude`、`lorentz` です。公開の構成関数では任意の明示的な物質IDを使えます。SBE/LDだけの配置ではTDDFTのランタイムとGSは不要で、サンプルの位置引数には未使用の仮パスを指定できます。LDだけならFDTD単独入口を使います。

### 同じ物質から複数の独立点を作る

同梱の配置JSONを読み、Si層の厚みを2から4 Bohrへ変えると、Siが2セルになり、同じ初期状態を使う独立したTDDFT点が2つ作られます。

```python
import json
from pathlib import Path

layout = json.loads(Path("mixed_layers.json").read_text())
for layer in layout["layers"]:
    if layer["material_id"] == "si-tddft":
        layer["thickness"] = 4
Path("two_si_points.json").write_text(json.dumps(layout, indent=2))
```

`mpi_mixed_native` に `--layers two_si_points.json --tddft-ranks 1 --rank-budget 4` を指定して4ランクで実行します（場1＋SBE1＋TDDFT1×2）。物質IDは共通でも、点IDとMPI群と時間発展状態は別です。厚みを増やすと後続の層も移動するため、報告されるセル範囲を確認してください。

SBEを2点にする場合は、上の条件を `sbe-point` に変更します。どちらも128ステップで各点の異なるA/c履歴を独立に再生し、電流一致を確認しています。同一プロセスに複数のnative TDDFT状態を置く検証ではありません。
