# A・E・Hの共通入力

ローカル場のデータ、モデルの要求ポート、既存APIへの明示アダプターを分けます。
結晶はA/c、LDはEのままです。長さゲージのモデルはEを宣言し、自作モデルはA・E・Hの必要な
組合せを宣言できます。Hは磁場強度でBではありません。磁気構成則やH→B変換は推定しません。
この型だけでnative孤立系TDDFTや磁気相互作用が有効になるわけではなく、
MPI ResponseSessionの通信は現時点ではAだけです。

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
native SBEの直接計算との一致もテスト対象です。孤立系native TDDFT長さゲージとMPI E/H通信は
後続の接続作業であり、この検証で実装済みとは扱いません。
