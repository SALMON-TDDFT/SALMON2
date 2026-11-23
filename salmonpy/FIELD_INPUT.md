# A・E・Hの共通入力

各点のA・E・Hと、モデルが必要とする成分を指定して既存の計算へ渡します。
結晶はA/c、LDはEのままです。長さゲージのモデルはEを宣言し、自作モデルはA・E・Hの必要な
組合せを宣言できます。Hは磁場強度でBではありません。磁気構成則やH→B変換は推定しません。
この型だけで孤立系TDDFTや磁気相互作用が有効になるわけではなく、
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
物理仕様自体はbackendの責務です。入力検査の失敗ではbackendを呼ばず、アダプターを再使用できます。
backendが例外を投げた場合はアダプターを無効化し、backendが巻き戻ったとは扱いません。
MPIや結合全体の原子的更新ではありません。

```sh
python -m salmonpy.examples.field_inputs
```

wheelサンプルは計算処理 LDの直接呼出しとの一致と、自作の時間差付きE/H受渡しを確認します。
計算処理 SBEの直接計算との一致もテスト対象です。孤立系計算処理接続とMPI通信は、後述の別APIと検証範囲を参照してください。

## 実験的な孤立系計算処理 E接続（ソース作業コピー）

同梱standard計算プログラムに検証用の別プロトコルを追加しました。
`SALMONPY_TDDFT_WORKER=E3`で`INIT_E dt capacity E0`、`STEP_E Ex Ey Ez`を受け、
`SALMONPY_TDDFT_E_V3 step delta_dipole[3] total_dipole[3]`を返します。
単位は原子単位、双極子はe·bohrです。電子双極子の変化と全電荷双極子は計算処理の
rt.dataと同じ定義です。電流密度ではなく、暗黙の時間微分や真空セル体積での除算はしません。

検証範囲は孤立系・固定イオン・非偏極PZ・スカラー擬ポテンシャル・1ランク・中点伝播、
初期E=0です。step nに渡すEは計算処理のcalcVboxがindex nで使う値です。
バッチ参照ではこの値を-(A[n+1]-A[n])/dtで準備しますが、worker内ではAを再構成しません。
AETRS、外部双極子源、磁場との結合は未対応です。周期系は従来のA駆動のままです。

開発用の`tools/test_isolated_reference.py`でC₂H₂参照計算を作り、
`tools/test_isolated_worker.py`で専用実行ファイルとの比較、FieldAdapter経由で2つの
子プロセスを交互に更新する独立性、E履歴の再実行を検証します。
`tools/isolated_worker_driver.py`は検証用で、インストール済み公開ローダーや原子的な
snapshot/restore APIではありません。backend異常時は子プロセスを閉じます。
光学物性の収束や同一プロセス内の複数solverを保証するテストではありません。

ソースを検査して明示的にビルドした実行ファイルを使用します。既存の配布済み計算プログラムが
自動的にE対応するわけではありません。公開ローダーとメモリ上の復元APIはISOLATED_TDDFT.mdを参照してください。
ローカルのバイナリ別検証記録で対応を判定します。ディスク再開はISOLATED_TDDFT.md、
MPI packet通信は次節を参照。ビルド済みcatalog配布は残作業です。ソース内の
`docs/isolated-worker-validation-2026-11-15.md`に検証記録があります。

## MPIでの明示的なA/E/H転送

`ResponseSession(..., field_contracts={point_id: contract, ...})` と
`session.step_fields({point_id: FieldPacket(...), ...})` を使います。
全ランクで同一の仕様を指定し、点ごとの必要成分・単位・時刻・位置・ゲージを状態作成前に照合します。
入力mappingを渡すのは場のleaderだけで、他ランクはNoneです。点IDはplanの整数transport IDです。
全点の入力を更新前に検査し、受信側でも不変packetを再構築・検査します。
物質stateはstep_fieldsを実装し、ゼロ起点で宣言した整数または半ステップ時刻の電荷電流密度jの3成分を返します。
出力仕様はpoint標本・offset=0または1/2のjを1つ指定し、backendにも宣言があれば一致を検査します。
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

## 明示的な時空間補間

`salmonpy.field_sampling.SampledField(quantity, times, positions, values, unit)`
は不変の点標本を保持します。valuesは[時刻,位置,3]です。
`sample(time=..., position=..., weights=..., temporal='exact'|'linear')` で
FieldValueを得ます。空間の重みは利用者が指定し、非負・合計1・指定位置の再現を検査します。
これは選択したステンシルの補間であり、適切な物理ステンシルの自動判定ではありません。
時間はexactでは保存時刻と厳密一致、linearでは履歴の範囲内だけ許可します。
範囲外への外挿、最近傍への丸め、単位変換はしません。

`sample_packet(contract, step_index=..., histories=..., weights=..., temporal=...)`
は成分ごとの時計を使いFieldPacketを作ります。historiesとweightsは必要成分名のmappingです。
点標本のみ対応し、区間平均への読み替えは拒否します。Hが必要ならH履歴が必要です。

```sh
python -m salmonpy.examples.field_sampling
```

このサンプルは合成場です。A/E標本の取得は次節を参照してください。電流の戻し時刻を含む
結合driverは後述の対応時刻配置に限定します。補間だけで任意の結合が成立するとは扱いません。

## 結合場のA/E標本を取得

`NativeEHField.sample_fields()` と `MultiscaleField.sample_fields()` は、
少なくとも1ステップ進めた使用可能な場から、A/EのSampledFieldを返します。
これはローカルの読取操作で、通信や更新をしません。返り値は元の場から独立した不変コピーです。
行順はfield.point_ids、位置はbohrです。Hのキーはありません。

- FDTDのE/H更新：Aは(n+1/2)dt、Eはn*dt。位置は登録したマクロセルの中心。
- Multiscale1D：Aはn*dt、Eは(n-1)dt。位置は物質セルのx座標、y=z=0。

各時刻のsamplesを保持し、同じquantity・unit・positionsを持つ標本のtimesとvaluesを
連結してSampledFieldを作れば、前節の明示補間へ渡せます。場を進める前の標本や閉じた場は拒否します。
`mpi_coupled_tddft` はfield_samples.npzにA/E、A_time/E_time、positions、point_idsを保存します。
このファイルは原子単位の点標本で、Hや波動関数restartではありません。
Hの取得は次節を参照してください。packet方式の結合は後述のPacketCoupledStepperを参照してください。

## 計算処理 Hの明示的な集合読取

H対応ABIを明示的にビルドした場合、全fieldランクで
`magnetic = field.sample_magnetic_field()` を呼びます。
1ステップ後から、時刻(n+1/2)dtのraw Hをマクロ点へ空間平均したSampledFieldを返します。
これはn*dtへ時間平均された観測ファイルのHとは異なる標本です。
ステップは進めません。観測出力の間隔にも依存しません。
必要な`sp_eh_sample_h`がない旧ライブラリは全fieldランクで明示的に拒否します。
ローカル読取のsample_fields()とは異なり、leaderだけで呼んではいけません。
結合サンプルは`--sample-h`でH/H_timeをfield_samples.npzへ追加します。
磁気応答模型・H→B変換・双方向packet結合を自動的に有効にする機能ではありません。

## 明示packet結合ドライバ

`PacketCoupledStepper(response, field)` を全contextランクで構築し、step()を全員で呼びます。
responseはfield_contracts付きResponseSession、fieldは場ランクのみで指定します。
場のマクロ点と物質点は1対1、位置・単位は一致させます。対応時刻は次のとおりです。

| 場 | 入力port | 戻りj |
|---|---|---|
| NativeEHField | A/H: offset=1/2、E: offset=0、すべてlag=0 | offset=1/2 |
| MultiscaleField | A: offset=0 lag=0、E: offset=0 lag=1、H不可 | offset=0 |

TimeSpecのt0とclock_offsetは0、標本はpoint、dtは共通です。状態のローカル時計はn*dt、
各portとjの物理時刻はoffset/lagを加えて宣言します。backend出力宣言があれば一致が必要です。
Hを使う場合はH対応runtimeを開始前に確認します。未来標本・位置不一致は開始前に拒否します。
初期feedbackは明示的にゼロ。各stepは場(previous j)→標本→物質→次のfeedbackの順です。
計算中の異常はresponseを無効化し、全体の巻き戻しは保証しません。
戻り辞書にfield_time、current_time、成分ごとのSampledFieldとJ_chargeを保持します。

```sh
mpiexec -n 3 python -m salmonpy.examples.mpi_packet_coupling NEW_OUTPUT
mpiexec -n 3 python -m salmonpy.examples.mpi_packet_coupling NEW_NATIVE --field native --build BUILD_JSON --ports eh
```

例の物質は合成電流模型です。孤立系の双極子を直接jとして返さず、磁気応答や自己場を自動推定しません。
既存A駆動用CoupledStepperは維持します。新driverは明示packetを返すbackend向けで、
既存TDDFT/SBEプロバイダーを黙ってpacket形式へ変換しません。一般の時空間補間やサブサイクルは後続です。

### 既存のA駆動binding

`salmonpy.material_response.packetize_A_binding(binding, field_phase='endpoint')`
で既存のA→j bindingを明示的に包み、Multiscale1Dへ接続します。FDTDのE/H更新結合では
`field_phase='yee_half_step'`を選びます。ResponseSessionには上記と同様、対応する
点ごとの`field_contracts`を渡します。provider・物質ID・実行方式・電流符号を維持し、
各状態が別々の受理済みA履歴を持ちます。元の状態は`step(A)`または
`advance(A_old, A_new)`と、端点A/jの時刻宣言を実装します。

物理ポートの時刻対応を明示する処理で、workerのローカル時計や数値更新順は変えません。
初期Aはゼロです。Aを取得できないbackendも、この初期化仕様を守る必要があります。
絶対時刻に依存する自作モデルは、この対応の意味を利用者が確認してください。
`response.timing_records()`に元と変換後の時刻を記録します。E/H変換・補間・双極子の
規格化は行いません。

既存サンプルに`--adapt-a`を付けると明示アダプター、`--legacy`で比較経路を選べます。
`--runtime RUNTIME --bundle GS_BUNDLE`を追加すると二つの独立した管理TDDFT workerを
使います。bundleには`rt.inp`とrestartデータが必要です。3ランク、worker各1スレッドで
実行してください。短い接続回帰であり、物理的な収束検証ではありません。計算処理経路は
さらに`--field native --build BUILD_JSON`が必要です。取得・ビルドは明示操作のままです。

このアダプターはleader_only backend・物質群あたり1ランクに限定します。
群内collective backendは起動前に拒否します。FDTD側のMPI並列とは別の制限です。
