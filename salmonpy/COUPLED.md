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
混在サンプルには不確実性付きのメモリ見積もり・計測があります。予算超過は警告のみで上限を強制しません。小さなGSで測定してから規模を増やします。

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

電子物質の各格子セルに独立した状態とMPI群を割り当て、LDは場のランク内で管理します。有限の直方体配置であり、無限膜の均質化や代表点への間引きはしません。厚みや断面積を増やすと電子点と必要ランクが増えます。ランク予算超過は警告のみで、全物質点と必要ランク数を維持します。実行には表示された必要ランク数が必要です。未知のID、バックエンド不整合、吸収境界との重なり、LD構成則の近傍に入る電子点は拒否します。真空の隙間も明示指定が必要です。ランク予算はメモリ上限ではなく、`memory_bytes` は未推定です。初期状態準備、MPI起動、物理的な収束確認は別途行います。

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

### 結合サンプルの再実行による保存／再開

`mpi_mixed_native` の `--save-at 64` は64ステップで停止し、出力先の `checkpoint/` に入力照合情報とA・E・電流履歴を保存します。同じ全ステップ数・配置・ライブラリ・GS・MPI構成を指定し、新しい出力先と `--resume 保存先/checkpoint` で再開します。`--compare 連続実行の出力先` で最終履歴を比較できます。

これは初期状態から保存時刻まで場・LD・電子状態を同じ順序で再計算する再開です。過去のA・E・電流が完全一致することを検査してから続きを計算します。内部状態の直接保存ではなく、過去部分の計算時間が再び必要です。初期状態とライブラリを保持してください。現時点では同じ環境・構成の同梱サンプル専用です。汎用の結合状態アーカイブは未実装です。同梱サンプルの直接復元は次節を参照してください。

### 内部状態を直接保存して再開する

同梱の `mpi_mixed_native` で `--direct-save-at 64` を指定すると、64ステップで停止し、`OUTPUT/direct-checkpoint/` に各ランクの内部状態を保存します。同じ入力・ライブラリ・GS・MPI構成で、新しい出力先と `--direct-resume 保存先/direct-checkpoint` を指定すると、保存時刻までの時間発展を繰り返さずに再開します。全ステップ数は元の128のまま指定します。`--compare 連続実行の出力先` で比較できます。`result.json` の `restored_step` と `native_steps_this_job` が開始時刻と今回の実行ステップ数です。

保存対象はnativeのE/H・PML、LDのP/J/Qと旧場、電子の波動関数バッファ・密度・局所ポテンシャル・Hartree履歴・RT履歴・A・電流、Pythonの積算A・直前電流・時計です。SBEの密度行列も直接保存します。復元時は静的な入力・係数・初期資源を準備してから動的状態を読みます。全ランクのファイル一覧とSHA256を記した完了ファイルがない保存、不一致の入力・ライブラリ・構成は拒否します。

対象は同じ環境・分割を使う周期結晶・中点伝播サンプルです。XCはPZ、TB-mBJ（自動／固定cval）、Libxc-PZ／PBEに対応します。Libxcには対応ライブラリを明示ビルドしたnativeランタイムが必要です。自作バックエンドの汎用保存、別MPI分割、HSEや他のXCは別途対応・検証が必要です。既定の保護付き復元では扱えるステージング失敗時も元の状態を維持します。明示的な無保護経路や未処理のnative異常には再試行保証がありません。大きな状態の保存量とメモリには余裕を持たせてください。

直接保存形式v2はXC名・導関数設定・入力cvalを照合し、Libxcのハンドルは新しい初期化で再作成します。元のランタイムと保存ファイルのハッシュを保持してください。サンプルのRUNTIME引数には管理済みランタイムのほか、明示構築・監査したシリアル参照計算用のシリアル版SALMON実行ファイルも指定できます（MPI版を参照用子プロセスに使わないでください）。別のXCを試す場合は、そのXCで用意したGSとRT入力を使います。

#### TB-mBJ／Libxc-PBEの内部MPI直接再開

同じXCで準備したGS/RTバンドルを使い、外側4ランク（場1、SBE1、TDDFT軌道分割2）で実行します。`--tddft-ranks 2 --rank-budget 4 --layers 配置JSON` を指定し、連続実行、`--direct-save-at 63` で止める保存実行、新しい出力先への `--direct-resume 保存先/direct-checkpoint --compare 連続実行出力先` の順で実行してください。全ステップ数は128のままです。

TB-mBJ（自動係数／cval=1.2）とLibxc-PBEで、この構成の再開後65ステップと連続128ステップのA/c・E・電流がビット単位で一致しました。同じ二段階MPI構成での直接再開の検証であり、MPI分割を変更する再開は対象外です。場2ランクとTDDFT2ランクを同時に使うにはSBEを含め5ランクとなり、今回の4コア枠では未検証です。

直接再開のチェックサム検証はランク間で分担します。全ランクが同じ完了manifestを読み、検証がすべて成功したことを確認してから内部状態を復元します。保存先は全ランクから同じ内容を参照できる共有ディレクトリが必要です。ノードごとの異なるローカルコピーは対象外です。`result.json` の `checkpoint_validation` に各ランクの検証ファイル数・バイト数を記録します。

## 高速PR版HSEを外側MPI／FDTDへ接続する

`FastHSEMaterialProvider` と `FastHSEResponseBackend` は、高速PR版HSEの専用serialワーカーを既存のA/c→電荷電流j契約へ接続します。FDTDは既存のMPI入口、LDは場ランク内、各HSE点は1ランクの物質群が独立した子プロセスを所有します。HSE内部MPIは使用しません。ホストのカレントディレクトリを変更せず、import時にMPI初期化・ビルド・取得はしません。

```python
from salmonpy.fast_hse import FastHSEMaterialProvider, FastHSEResponseBackend
from salmonpy.material_response import MaterialBinding

binding = MaterialBinding(
    "hse-point",
    FastHSEMaterialProvider(bundle_directory, "rt.inp", max_steps=256),
    FastHSEResponseBackend(worker_build_json, existing_work_root, threads=1),
)
# このbindingを既存のassemble_native_layers / ResponseSessionへ渡します。
```

監査済み高速PR版ビルド記録から、ワーカーを新しい外部ディレクトリへ明示構築します。通常版HSEを指定しません。

```sh
python tools/build_fast_hse_worker.py FAST_PR_BUILD_JSON NEW_WORKER_DIRECTORY
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 mpiexec -n 3 python -m salmonpy.examples.mpi_mixed_native FIELD_BUILD_JSON UNUSED_REFERENCE HSE_BUNDLE NEW_OUTPUT --fast-hse-build NEW_WORKER_DIRECTORY/build.json --layers PATH_TO_mixed_hse_layers.json --dt .001 --steps 256 --rank-budget 3
```

`mixed_hse_layers.json` はLD2材料＋SBE＋HSEの配置です。`UNUSED_REFERENCE` はこのHSE専用例では未使用の位置引数で、存在しない仮パスで構いません。HSEの独立参照計算は同じ専用ワーカーで行います。場2ランクでは `mpiexec -n 4 ... --field-grid 1 2 1 --rank-budget 4` とします。各HSE群は1ランク限定です。子ワーカーがその群の割り当てCPUを使い、呼出元ランクは応答を待ちます。OpenMPを増やす場合は別途CPU予算を割り当ててください。

入力は `tddft_response`、`xc='hse06'`、明示的な高速伝播器（`hse_taylor4`／`hse_taylor4_full`／`hse_ptcn`）を含む元GS/RTバンドルです。ワーカーは選択したdtと上限ステップ数だけを設定し、他の物理条件を保持します。内部並列数は1、学習型には対応するdense-k／固定占有／全支持領域の条件が必要です。初期ゼロ場のサンプルは初期HSE電流が1e-12 au以下であることを要求し、基準電流を差し引きません。合成テスト入力は `tools/prepare_fast_hse_coupled_fixture.py` で作れますが、実材料のGSや物理収束の代わりにはなりません。

`--direct-save-at 127` で場/PML/LD/SBEとHSEワーカーの内部状態を同時刻に保存し、新しいジョブで `--direct-resume SAVE_OUTPUT/direct-checkpoint --compare CONTINUOUS_OUTPUT` を指定します。HSEのLOAD応答により受理済みステップ・電流も復元し、保存したA・結合履歴と照合します。過去127ステップを再実行せず、残り129ステップを進めます。共有保存先・同一ワーカー／初期バンドル／入力／MPI構成が必要です。新しいワーカーの初期資源作成と、検証用の独立参照再生は別途行います。

対象はmacOS/Linux向けPOSIXワーカー（`/dev/stdin`使用）です。今回の実測はApple Silicon。このワーカー方式ではHSE内部MPIを使いません。変更した分割への再開、GPU、任意バックエンドの一般保存、native異常後の継続は未対応です。

HSEワーカーへの入力は原子単位（unit_system=au）に限定します。内部パルスを無効にし、ae_shape1はnoneまたは明示的なゼロimpulse、ae_shape2はnone、trans_longiはtrとします。FDTDからのAと内部パルスの二重駆動を拒否します。状態のtimeは受理済みstep×dtです。

## 高速PR版HSEの内部k点MPI

`salmonpy.fast_hse_mpi.NativeFastHSEBackend` は、Python側で用意した物質コミュニケータ内で高速PR版HSEを集合的に実行します。MPIの初期化・終了はホストが管理します。import時にMPI読込やビルドは行いません。1群1点、1プロセス1回のnative寿命、メインスレッド専用です。呼出中は作業ディレクトリを一時変更するため、相対パスを使う別スレッドと併用しないでください。

`tools/prepare_hse_rt_state.py NEW_SOURCE --explicit-routing` で新しいソースを準備し、`tools/build_fast_hse_mpi.py PREPARED_SOURCE NEW_OUTPUT --jobs 2 --cmake-arg=...` で明示ビルドします。ホストと同じMPIのFortran/Cコンパイラ、OpenBLAS/LAPACK、FFTW、必要ならLibxcのパスをRUNTIMEガイドに従って指定します。MPI/HSEを必須とし、GPU・Poisson-FFTW・ScaLAPACK・EigenExa構成は拒否します。記録にはソースとアダプターのハッシュを含めます。実行エンジンは固定した高速PR版のみです。

```python
from salmonpy.fast_hse import FastHSEMaterialProvider
from salmonpy.fast_hse_mpi import NativeFastHSEBackend
from salmonpy.material_response import MaterialBinding

binding = MaterialBinding("hse-point",
    FastHSEMaterialProvider(bundle_directory, "rt.inp", max_steps=256),
    NativeFastHSEBackend(mpi_build_json, existing_shared_work_root))
# assemble_native_layers / ParallelPlanでこの点に2ランクを割り当てる。
```

k点のみの1／2ランク分割に対応し、軌道と実空間は分割しません（`nproc_ob=1,nproc_rgrid=1,1,1`）。バックエンドは`nproc_k`を群サイズに、dtと容量を指定値に設定します。群に十分なk点を用意してください。外場・固定占有・全支持領域の制約は上記と共通です。サンプルでは結合後の独立比較用にserialワーカーのビルド記録も指定します。

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 mpiexec -n 4 python -m salmonpy.examples.mpi_mixed_native FIELD_BUILD_JSON UNUSED_REFERENCE HSE_BUNDLE NEW_OUTPUT --fast-hse-build SERIAL_WORKER_BUILD_JSON --fast-hse-mpi-build MPI_BUILD_JSON --hse-ranks 2 --layers PATH_TO_mixed_hse_layers.json --dt .001 --steps 256 --rank-budget 4
```

場1＋SBE1＋HSE2で4ランクを使い、2種類のLDは場ランク内に保持します。同じ引数で連続実行、`--direct-save-at 127`、新しいジョブで`--direct-resume SAVE_OUTPUT/direct-checkpoint --compare CONTINUOUS_OUTPUT`を順に実行します。HSEは群サイズ・群内ランク・k範囲を含むランク別ファイルを保存し、全ファイルを検査してから数値状態を反映します。初期化は静的な交換／FFT資源を再構築し、過去の時間発展は行いません。同じバイナリ・元GS・入力・分割と共有保存先が必要です。早期の不正・切断ファイルは数値反映前に拒否します。既定の保護付き復元では、扱えるステージング失敗時も元の全状態を維持します。明示した無保護経路や未処理のnative異常には再試行保証がありません。

Apple Siliconで微視的8³格子・4×1×1 k点の合成入力を実測しました。Taylor／学習PT-CNとも場1＋SBE1＋HSE2、dt=.001、256ステップ中127で直接保存し、残り129回の呼出後にA/c・E・jの全履歴がビット単位で一致しました。独立serial高速PR応答との最大HSE電流差は4.07e-16／1.71e-16 auです。接続・再開の回帰検査であり、物理収束・大規模性能の保証ではありません。記録：`docs/fast-hse-kmpi-validation-2026-11-16.md`。

層構成時、バックエンドは任意で `supported_group_sizes=(1,2)` のように対応ランク数を宣言できます。空でない正の整数の重複なしtupleに限定し、未対応ランク数はMPI・初期状態準備より前に拒否します。高速HSEは `(1,2)` を宣言します。宣言しない自作バックエンドは既存のexecution・分割検査を使います。この属性は数値精度・再開機能の認証ではなく、直接ResponseSessionへ接続する場合にも、同じ検査と宣言のランク間照合を物質準備前に行います。

ID付き直接保存の既定manifestは`salmonpy.native-direct/3`で、従来の/2と共通の互換情報として、保存方式、save/restore対応、各MPI群の所有ランク・点ID、場の分割と電子状態のk／軌道／実空間分割を記録します。全ランクで現在の状態と照合してからLOADを呼ぶため、総ランク数だけが同じ分割変更も拒否します。初期資源の作成は検査前に行いますが、検査中の伝播・snapshot/restore/LOAD呼出はありません。旧manifest/1には必要な互換情報がないため新しい経路では復元を拒否します。旧保存データは対応する元の版で扱ってください。nativeペイロードの形式は変更していません。自作direct-state/1コールバックには、明示的なcheckpoint_format、checkpoint_modes（save/restoreのtuple）、checkpoint_decomposition（群サイズを積に持つ正の整数tuple）が必要です。宣言は利用者の契約であり、その実装・失敗時の原子的復元を認証するものではありません。単にsnapshotメソッドがあるだけでは直接再開と判定せず、同梱LD/SBE以外の再生方式を混入させません。

保存対応の宣言とMPI互換性の照合は共通API `salmonpy.checkpointing` を使います。自作バックエンドや直列ホストからの利用は[CHECKPOINTING.md](CHECKPOINTING.md)を参照してください。既存のnative-direct/2の記録形式と再開経路は維持します。

native混在サンプルの場・電子・LD/SBE保存処理も `CheckpointTarget`／`execute_checkpoint` を通します。保存形式とnative-direct/2の互換情報は維持し、コールバック失敗を共有した後は次の処理へ進みません。完了情報の公開失敗も全ランクへ通知します。

通常サンプルの`--direct-resume`は保護付き復元が既定です。場内LD・SBE・native TDDFT／高速PR HSEとPython側の時計・電流を、全群の候補成功後に反映します。古いABIと即時LOADのみの対象は読込前に拒否し、自動で従来方式へ切り替えません。明示的な従来方式は`--unprotected-direct-restore`です。途中失敗からの継続・再試行と保証範囲は[CHECKPOINTING.md](CHECKPOINTING.md)を参照してください。

`--memory-budget-gib 4`は全MPIランク合計4GiBの助言予算です。通常状態と復元候補の同時保持分を報告し、超過時も停止しません。[MEMORY.md](MEMORY.md)を参照してください。

## 混在サンプルのIDと再開

`python -m salmonpy.examples.mpi_mixed_native` は `--run-id`、論理点ID、物質ID、物質定義・初期状態のSHA256、MPI所有ランク、LDの材料番号対応を保存します。`mixed-example` は試験用の既定ラベルです。実験ごとに明示した `--run-id` を使用してください。定義の識別には出力先や一時ディレクトリを含めません。

`identity.json`、通常実行の `result.json` と `point-response` アーカイブに同じ識別情報を引き継ぎます。直接保存の外側の形式は `salmonpy.native-direct/3`、再計算方式は `salmonpy.native-replay/2` です。数値状態の内部形式と更新順は変更しません。直接復元はIDとMPI分割を数値状態の読み込み前に照合し、保存時は全ランクでのID一致を検査します。

IDのない旧保存を使う場合だけ `--legacy-checkpoint-identity` を明示します。これは識別情報を検証できない旧経路であり、新形式のID不一致を回避する指定ではありません。保護付き復元は引き続き既定です。自作ホストでは `ResponseSession` の `point_record` / `point_map` に識別情報を渡してください。IDは呼び出し側の宣言であり、物理精度の認証ではありません。

## native TDDFTの更新前エラー拒否（限定対応）

新しく準備・ビルドした標準embedded TDDFTには `sp_tddft_error_protocol() == 1` があり、以下を結合先の物質MPI群だけで照合します。10: 不正・非有限入力／古いA／入力ファイル不在、11: ランク間不一致、12: 状態未作成・使用済み初期化、13: ステップ上限、14: 未確定の復元候補あり。これらは数値更新前に全ランク同じコードで返り、既存の受理済み状態を変更しません。初期化前の10/11は生存期間を消費せず、Pythonからも入力を修正して初期化を再試行できます。失敗時の作業ディレクトリは診断用に残ります。

Pythonは版付きの契約と全ランク一致を確認できた拒否だけを `InputError` とします。旧ライブラリ・未知のコード・非有限な戻り値は `RuntimeFailure` のままで、状態保持を推定しません。MPI群の全員が同じ入口へ参加する必要があります。無効なポインタ、MPI・プロセス異常、確保失敗は対象外です。本項の更新前検査ではソルバー深部を扱いません（次節の限定した深部復帰を参照）。初期化に入った後は引き続き1プロセス1生存期間です。高速HSEとnative FDTDのエラー契約はこの変更では拡張していません。

この保証は個別TDDFT状態の更新前拒否です。結合ドライバで他の部品が既に進んだ場合、全体の巻き戻しを保証せず、従来どおりセッションを失敗扱いにします。通常計算中の一般的な原子性が完成したとは判断しないでください。

## 大きなMacでの並列・格子テスト

4GiB・4コアの入門条件は、その条件での確認用として残します。大きなマシンでも利用者の実メモリ・コア数を確認し、OSや他の作業に余裕を残してください。64GB・18コアのMacでは、今回は48GiBを警告用予算として、最大15 MPIランク、各ランク1スレッドで測定します。これは推奨値の固定やハード制限ではありません。

`mpi_mixed_native` の `--field-shape 96 96 96` はFDTDの格子点数を変更します。各軸は40以上の偶数で、MPI分割に整合させます。既定値は40³のままで、格子間隔も2 bohrのままです。外側の真空領域が広がるため、異なる格子サイズを同一物理条件の比較に混ぜないでください。既定の物質点の物理座標は保ちます。電子状態のGS格子・k点は変更しないため、96³のTDDFT検証を意味しません。3Dは接続検証で、TDDFTの実用途は引き続き1D推奨です。

FDTD2＋SBE1＋TDDFT2は `mpiexec -n 5 ... --field-grid 1 2 1 --tddft-ranks 2`。FDTD12＋SBE1＋TDDFT2は `mpiexec -n 15 ... --field-grid 3 2 2 --field-shape 96 96 96 --tddft-ranks 2`。`...` には既存のPythonモジュール呼び出しとBUILD・RUNTIME・GS・新規出力先を入れます。`OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1` を指定してください。

ソース配布の `tools/test_mixed_scaling.py BUILD RUNTIME GS NEW_OUTPUT --rank-budget 18 --memory-budget-gib 48` は導入済みwheelを用い、同じ96³の入力で4/5/7/11/15ランクを順番に2回ずつ比較します。各計算のA/E/jを最初の計算と照合します。格子サイズが違う結果は比較拒否します。MPI実装固有の引数は `--launcher-arg` で指定し、計測中は他の重い処理を避けてください。

結果の `coupled_loop_seconds_max` は初期化・保存・直列参照計算を除いた時間発展区間です。通信計測には待ち時間が含まれ、純粋な通信帯域ではありません。FDTD内部のハロー通信は個別に分解していません。メモリは各ランクのピークRSSの和であり、同時刻の総常駐量や全システム使用量ではありません。時間・資源の計測失敗を計算失敗と混同しないでください。短い反復測定だけで一般的な並列効率や物理精度を保証しません。

## TDDFT深部の数値異常と直接巻き戻し

新しく準備・ビルドした標準native TDDFTは `sp_tddft_recovery_protocol() == 1` を持ち、各ステップの数値状態をメモリへ退避します。伝播後の波動関数（20）、密度（21）、Hartree場（22）、XCポテンシャル／エネルギー（23）、電流／全エネルギー（24）を検査し、保護経路の伝播器エラー（25）も返します。検査は物質MPI群全体で合意してから戻ります。エラー時は波動関数・密度・場・伝播履歴・XC状態・カウンタ・管理フラグを受理済み状態へ直接戻し、伝播の再実行をしません。失敗したステップの通常RT出力も追加しません。

Pythonの `NativeTDDFTBackend` が返す状態の `advance(A_old,A_new)` は、版付き契約と全ランク一致を確認できた復帰だけを `NumericalError` で通知します。`A` と `steps` はそのままで、原因を修正して同じ受理時刻から再試行できます。不正な波形を自動的に弱めたり、固定dtを変更したりしません。旧ライブラリや未知のエラーは `RuntimeFailure` の保守的な扱いを維持し、復帰を推定しません。

対象は周期・非スピン分極・固定イオン・中点伝播のPZ、TB-mBJ（自動／固定係数）、Libxc-PZ/PBEです。SOC、+U、jellium、イオン運動／力出力、Gram–Schmidt、projection/unfold、microscopic Maxwell結合などの保護対象外設定は更新前に15を返します。標準バッチ計算は従来の演算経路を使います。既存の外部checkpoint形式も変更しません。メモリ退避・有限性検査の費用が毎ステップ増えるため、資源見積もりには余裕を残してください。予算は警告のみです。

これはTDDFT群の状態保証です。FDTD・他の物質がすでに進んだ混在ステップ全体を自動巻き戻しする契約ではありません。MPI／プロセス停止、確保失敗、未制御の外部ライブラリSTOP、初期化深部のSTOPや出力装置障害は未対応です。試験専用の一部ランク障害注入は `SALMONPY_RECOVERY_TESTING` を明示して別ビルドした場合のみで、本番ビルドには注入用C入口がありません。

## 高速PR版HSEの深部数値異常からの復帰

新しい `NativeFastHSEBackend` 用ライブラリは `sp_fast_hse_recovery_protocol() == 1` を宣言します。対応範囲は周期・非スピン分極・固定イオンの高速PR版HSE06、Taylor4／学習PT-CN、物質群内k点分割のMPI1／2です。通常版HSEへの切替はしません。

各ステップ前に電子状態とACE・因子／PT-CN履歴・交換キャッシュ・カウンタを退避します。FFTプランやMPIコミュニケータなどの所有資源はコピーしません。波動関数・密度・Hartree・XC・端点の有限性検査に加え、対応したACE演算エラー、PT-CNの収束／Gram／射影検査、履歴更新エラーを物質群全体で合意して返します。受理済み状態を直接戻し、過去ステップや交換履歴を再計算しません。失敗試行の通常RT出力を追加せず、A・時計・電流・エネルギーを維持します。診断ログは残ることがあります。

Pythonは宣言された版と全ランクの既知エラーを確認できる場合だけ `NumericalError` を返し、同じ受理時刻から修正した入力で再試行できます。旧ライブラリ・未知エラーは `RuntimeFailure` の保守的扱いを維持します。dtや入力振幅を自動変更しません。退避と検査のメモリ・時間費用は増え、予算超過は警告のみです。

これはHSE物質群の保証で、FDTDや他の物質を含む混在ステップ全体の自動巻戻しではありません。初期化深部、確保失敗、プロセス／MPI停止、未制御の外部STOP、出力装置障害は対象外です。異常を捕捉できない場合に継続を推定しないでください。

検証用 `tools/test_fast_hse_deep_recovery_mpi.py BUILD GS NEW_OUTPUT --steps 32 --fault-step 13` は2ランク、OMP／BLAS1、外側のtimeoutで実行します。`--split` は4ランクの独立2群です。`--overflow-A 1e100` はPT-CNの実数値発散試験用で、科学計算の入力ではありません。`--faults` は `SALMONPY_RECOVERY_TESTING` を明示した別ビルドだけに使い、本番には注入用入口を含めません。検証記録は `docs/fast-hse-deep-recovery-validation-2026-11-17.md` です。合成小規模系の実装検証であり、物理精度の収束保証ではありません。

## native FDTDの数値異常復帰

新しい標準ホストMPIライブラリは `sp_eh_recovery_protocol() == 1` を宣言します。`NativeEHField.advance` は対応する検査付き入口を選び、場・PML・場内の全LD・電流格子・時計をメモリ退避します。E更新後（20）、LD更新後（21）、H更新後（22）、点サンプル／積算A／通常観測出力の演算（23）の異常を場MPI群で合意し、退避した状態へ直接戻します。波の再伝播はしません。既知の復帰では `NumericalError` を返し、PythonのA・E・時計・failedも維持するため、入力を修正して再試行できます。失敗試行の通常観測ファイルは追加しませんが、診断・進行表示は残る場合があります。

対象は現在のnative E/H平面波経路と場内LDです。TTM、線形応答専用経路、ASE／ART、周波数別観測は更新前に拒否します。バッチFDTDは従来経路を維持します。新しい入口を持たない旧バイナリ・未知エラーでは復帰を推定せず、従来の失敗扱いです。初期化・確保・MPI／プロセス停止・外部STOP・出力装置障害は保証外です。場群単体の保証であり、`CoupledStepper`／`ResponseSession`全体を自動巻戻し・自動再試行する機能ではありません。退避と検査の費用を見込み、メモリ予算は引き続き警告に使います。

`tools/test_eh_deep_recovery_mpi.py BUILD NEW_OUTPUT` を1／2ランク、OMP／BLAS1、外側timeoutで実行します。`--split` は4ランクの独立2群です。`--faults` は明示した別の `SALMONPY_RECOVERY_TESTING` ビルドだけに使います。本番には注入入口を含めません。数値オーバーフローの試験値は科学入力ではありません。

## 明示した契約での混在ステップ保護

`CoupledStepper.step_protected(targets)` は自作コードと完全な退避・復元コールバックを接続する任意の入口です。[STEPPING.md](STEPPING.md) と6ランクの `examples/mpi_protected_step.py` を参照。native全体の保護は未対応で、個別ソルバーの復帰から推定しません。

監査済みnative全体の保護は [NATIVE_STEPPING.md](NATIVE_STEPPING.md) の別プロトコルと `--protect-steps` で利用できます。別のHSE GSは `--hse-bundle`、同時TDDFT／HSEの例は `examples/mixed_tddft_hse_layers.json` です。

## 高速PR版HSEの対称性opt-in

NativeFastHSEBackend(...,allow_symmetry=True)は内部kランク1／2・Taylor／PT-CNに対応し、学習には新しい履歴capabilityが必要。
通常混在例は--hse-symmetry --hse-reference-bundle FULL_K_BUNDLEを明示する。
Aが対称性を破ったときのInputErrorから、全体を戻すには--protect-steps／step_protectedを使う。
直接／保護付き再開と2種類の場内LD・SBE・PBE-TDDFT・HSE混在を検証した。
初期状態・場の局所方向・XCの適合性は別途確認する。極小勾配の安定評価と元の失敗例での検証も含め、
[対称性ガイド](FAST_HSE_SYMMETRY.md)を参照。

縮約k点の内部2ランクと同一分割での直接再開に対応します。追加のnative k-MPI capabilityが必要です。分割変更による再開は未対応。詳しくは対称性ガイドを参照してください。

### 警告と実行エラーの区別

メモリ・ランク数の予算、推奨規模や物理精度の注意は警告として扱い、実行可能な計算を止めません。入力や刻み幅を自動変更しません。読み出せない・破損したファイル、配列寸法や実際のMPI分割の不整合、nativeが停止する既知の未対応条件はエラーにします。未検証であることだけを理由に拒否しません。

高速PR入力検査では、有限な占有数の範囲逸脱や電子数との和の不一致は警告のみで、値を切り詰めたり再規格化したりしません。起動時の波動関数検査はバイト数の整合確認にとどめ、巨大ファイルの有限性全走査を毎回追加しません。そのためNaNを含む軌道を起動前に検出する保証はありません。旧GSヘッダーから物理的なGSの同一性を証明することもできません。

SCF未収束は警告のみです。有限な保存状態は後続計算に使え、`GroundStateResult.converged=False` と保存された残差・閾値・反復数で識別できます。

高速PRのnative起動前に、自動生成k格子と保持対称性の適合・代表点数・重み、物質種ごとの擬ポテンシャル参照、GSのヘッダー・占有数・軌道のバイト数を検査します。KYの検査はnativeと同様に継続レコード・余分な末尾列・Eを省略した指数を扱います。他の認識済みPP形式は参照と読み出し可否までで、内容検査はnativeに委ねます。対応GS形式はnative endianの4バイトrecord marker・整数／論理値、binary64占有数、単一global real64／complex128軌道ファイルです。任意の `file_kw`、unfolding、偶数Gamma端点格子は現在の高速PR経路の制約です。GSを作ったMPIランク数と今回のランク数の一致は要求しません。直接RT checkpointの分割互換性検査とは別です。

高速PR版HSEの外場はCartesianの3成分A/cで、x・y・z・円偏光を固定z条件で制限しません。縮約k点では実際のAを各ステップで保持対称性と照合し、不適合なら状態変更前に拒否します。適合群がなければ全k点の初期状態を使います。新コアは一般セルの対称性付き初期化にも対応します。[対称性ガイド](FAST_HSE_SYMMETRY.md)を参照してください。

### 内部RT状態の所有分離（開発用）

保護付きステップの数値退避・出力メタデータと直接復元候補は、プロセス共通ではなく
各 `tddft_context` が所有します。同じ不変の物質入力とコミュニケータを使う、対称性なしの
内部高速PRコンテキストでは、一方の巻き戻し・候補の確定／取消・終了が他方の候補を
消費しません。個別に選んだdt・容量と3成分の波形で交互に進められます。
出力の再接続は空いている番号を取得し、巻き戻し中の状態を終了しても、他の状態が
再利用した出力番号を閉じません。

公開Python／Cバックエンドは引き続き1プロセス1native寿命です。安定した複数状態のCハンドル、
異なる共有入力・任意のコミュニケータ・スレッドからの呼び出しの対応を意味しません。
対称性・初期化・I/O・タイマーの共有状態と、同期的な直列化処理の配列上限は今後の分離対象です。
コンテキスト全体や所有権を持つFFT／MPI／Libxcハンドルをコピーしないでください。
テキスト出力の巻き戻しは、ホストで保存済みのファイル内容を戻してから再接続します。
数値チェックポイントだけではログを切り詰めません。

開発者向け検証は `tools/test_fast_hse_context_owners_mpi.py BUILD FULL_K_FIXTURE NEW_OUTPUT`
を明示的に実行します。内部プローブをコンパイルし、Taylor・学習Taylor・学習PT-CNの独立実行と
交互実行を比較します。同時保護、復元候補の確定／取消、片方の終了、出力番号の再利用を含みます。
既定は独立した2群、各群MPI1／2・OpenMP1／2です。合成入力による状態分離検証であり、
実材料の収束検証ではありません。

再接続に失敗しても退避と出力の所有情報を保持するため、再試行・終了が相手の再利用済み番号に作用しません。

同一モデルの世代番号付きCコンテキストは、別の[nativeハンドル入口](NATIVE_HANDLES.md)を参照してください。Pythonバックエンドの1寿命の契約は維持します。
