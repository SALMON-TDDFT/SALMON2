# A・E・Hの共通入力

ローカル場のデータ、モデルの要求ポート、既存APIへの明示アダプターを分けます。
結晶はA/c、LDはEのままです。長さゲージのモデルはEを宣言し、自作モデルはA・E・Hの必要な
組合せを宣言できます。Hは磁場強度でBではありません。磁気構成則やH→B変換は推定しません。
この型だけでnative孤立系TDDFTや磁気相互作用が有効になるわけではなく、
MPIでは後述の明示的なfield_contractsとstep_fieldsを使います。

`FieldValue(value, time=..., position=..., unit=...)`は有限実数3ベクトルと位置の不変コピーを持ちます。
時刻は常に原子単位です。position_unitは既定bohrで、別単位も明示宣言できます。
samplingは既定pointで、既存TimeSpecのラベルを使います。
`FieldPacket(A=None, E=None, H=None)`は未提供とゼロを区別し、成分ごとに異なる時刻を保持できます。

```python
from salmonpy.field_input import FieldValue, FieldPacket, FieldContract, FieldAdapter

# 既存の結晶/SBE状態。A/cの端点時刻を宣言しているもの。
contract = FieldContract(state.timing, units={'A': 'A/c:au'},
                         position=(0., 0., 0.), gauge='velocity')
adapter = FieldAdapter(state, contract, legacy_quantity='A')
packet = FieldPacket(A=FieldValue([0., 0., 1.e-5], time=state.time+state.dt,
                                 position=(0., 0., 0.), unit='A/c:au'))
j = adapter.step_fields(packet)
```

E駆動ではInputPort('E')のTimeSpec、単位E:au、legacy_quantity='E'を指定します。
長さゲージならgauge='length'、LDなら'phenomenological'を宣言してください。
backendのゲージを自動選択しません。velocityはA/c、lengthはEを要求しますが、宣言は物理実装の検証ではありません。
Aの微分、H/B変換、単位変換は行いません。戻り値も既存のままです（LDはP,J、SBE/TDDFTは電流）。
既存stepへのアダプターは標準の原子単位、端点サンプリング、clock offset/lagなしに限定します。
backendがtiming/dt/timeを公開していれば一致を確認します。元のstep(A)/step(E)も引き続き使えます。

## 自作の複数ポートモデル

TimeSpecの入力をA/E/Hから宣言し、必要な成分と単位の厳密な対応、位置、gauge='custom'
（または適切なゲージ）でFieldContractを作ります。step_fields(packet)を実装した状態を
legacy_quantityなしでFieldAdapterに渡します。単位ラベルは完全一致で、変換機能ではありません。
HのSI:A/m等は呼出側の宣言です。物理的な解釈と異なる物理量間の単位整合性は自作backendの責務です。
サンプルはE/Hを別々に記録し、異なる次元の量を加算しません。

受理時刻index nでの各入力時刻は`t0 + (n + clock_offset + port.offset - port.lag)*dt`です。
既存の半ステップoffsetでHの時間差を宣言できます。時刻比較には丸め誤差の許容幅を使います
（最大dt*1e-3）。位置・位置単位・samplingは完全一致です。
扱うのは一つのモデル位置へ評価済みの場であり、生のYee格子の各位置ではありません。
必要な補間は明示的に先に行い、値を暗黙に別位置へ移しません。
n=0ではuses_t0=Falseのポートは要求しません。

検査は更新前に完了します。必須場の欠落、形状不正、複素数・bool・非有限値、単位・時刻・位置の
不一致を拒否します。余分に渡された未宣言成分は、複数ポートbackendへ渡すpacketから除きます。
既存stepへは書換可能なコピーを渡し、元の場データは不変です。

## 状態と失敗の扱い

FieldAdapterは状態を借用し、closeしません。step_indexの既定は0で、進んだ状態へ接続する場合は
受理済みindexを明示します。接続中にbackendを別途advance/restoreせず、必要なら正しいindexで
アダプターを作り直してください。既存・自作複数ポートのどちらも、公開されたtiming/dt/timeを更新前後に照合します。
timeはt0+(受理index+clock_offset)*dtです。属性がないモデルでは呼出側の宣言に依存します。
物理契約自体はbackendの責務です。入力検査の失敗ではbackendを呼ばず、アダプターを再使用できます。
backendが例外を投げた場合はアダプターを無効化し、backendが巻き戻ったとは扱いません。
MPIや結合全体の原子的更新ではありません。

```sh
python -m salmonpy.examples.field_inputs
```

wheelサンプルはnative LDの直接呼出しとの一致と、自作の時間差付きE/H受渡しを確認します。
native SBEの直接計算との一致もテスト対象です。孤立系native接続とMPI通信は、後述の別APIと検証範囲を参照してください。

## 実験的な孤立系native E接続（ソース作業コピー）

同梱standardランタイムに検証用の別プロトコルを追加しました。
`SALMONPY_TDDFT_WORKER=E1`で`INIT_E dt capacity E0`、`STEP_E Ex Ey Ez`を受け、
`SALMONPY_TDDFT_E_V1 step delta_dipole[3] total_dipole[3]`を返します。
単位は原子単位、双極子はe·bohrです。電子双極子の変化と全電荷双極子はnativeの
rt.dataと同じ定義です。電流密度ではなく、暗黙の時間微分や真空セル体積での除算はしません。

検証範囲は孤立系・固定イオン・非偏極PZ・スカラー擬ポテンシャル・1ランク・中点伝播、
初期E=0です。step nに渡すEはnativeのcalcVboxがindex nで使う値です。
バッチ参照ではこの値を-(A[n+1]-A[n])/dtで準備しますが、worker内ではAを再構成しません。
AETRS、外部双極子源、磁場との結合は未対応です。周期系は従来のA駆動のままです。

開発用の`tools/test_isolated_reference.py`でC₂H₂参照計算を作り、
`tools/test_isolated_worker.py`で専用実行ファイルとの比較、FieldAdapter経由で2つの
子プロセスを交互に更新する独立性、E履歴の再実行を検証します。
`tools/isolated_worker_driver.py`は検証用で、インストール済み公開ローダーや原子的な
snapshot/restore APIではありません。backend異常時は子プロセスを閉じます。
光学物性の収束や同一プロセス内の複数solverを保証するテストではありません。

ソースを検査して明示的にビルドした実行ファイルを使用します。既存の配布済みランタイムが
自動的にE対応するわけではありません。公開ローダーとメモリ上の復元APIはISOLATED_TDDFT.mdを参照してください。
ローカルのバイナリ別検証記録で対応を判定します。ディスク再開はISOLATED_TDDFT.md、
MPI packet通信は次節を参照。ビルド済みcatalog配布は残作業です。ソース内の
`docs/isolated-worker-validation-2026-11-15.md`に検証記録があります。

## MPIでの明示的なA/E/H転送

`ResponseSession(..., field_contracts={point_id: contract, ...})` と
`session.step_fields({point_id: FieldPacket(...), ...})` を使います。
全ランクで同一の契約を指定し、点ごとの必要成分・単位・時刻・位置・ゲージを状態作成前に照合します。
入力mappingを渡すのは場のleaderだけで、他ランクはNoneです。点IDはplanの整数transport IDです。
全点の入力を更新前に検査し、受信側でも不変packetを再構築・検査します。
物質stateはstep_fieldsを実装し、t=0始まりの整数時刻n*dtにおける電荷電流密度jの3成分を返します。
出力契約はpoint標本・offset=0のjを1つ指定し、backendにも宣言があれば一致を検査します。
既存FieldAdapterでstateが宣言する時刻との整合も検査します。余分な成分は転送しません。
従来のstep(A)は維持し、同じsessionで両方式を混在させません。

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 mpiexec -n 3 python -m salmonpy.examples.mpi_field_inputs
```

これは合成電流模型による通信例です。孤立系TDDFTの双極子を直接jとして返すことはできません。
A/E変換、空間・時間補間、双極子から電流への規格化は自動化しません。
FieldContractに合う標本を呼び出し側で明示して用意します。
backend更新の失敗はsessionを無効化し、全物質の巻き戻しは保証しません。
通信計測のfield value bytesは数値3成分のみで、packetのメタデータ通信量を含みません。
