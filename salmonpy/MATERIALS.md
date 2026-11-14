# 自作物質モデルを直列とMPIで使う

[English](MATERIALS.en.md) · Topic: `materials` · Paired revision: `2026-11-14.2`

物質モデルは一つの直列契約で書きます。`crystal.respond`で波形への応答を求め、
`materials.bind_material`で既存のMPI応答サービスへ接続できます。
MPI側では計算点を分散し、各点に独立した状態を作ります。
単一の直列モデルを複数ランクで並列化する機能ではありません。

## モデル側に必要なもの

初期モデルは読み取り専用で、`backend`に空でないアルゴリズム名・版を持ちます。
`new_state(*, dt, time0=0., initial_A=(0,0,0))`は毎回独立した可変状態を返します。
初期配列を状態同士で共有・変更せず、受け取ったA配列も変更しません。

| 状態の入口 | 契約 |
|---|---|
| `time` | 受理済み時刻。初期値time0、nステップ後はtime0+n·dt |
| `current` | 現時刻の有限な3成分電流密度 |
| `timing` | 明示的な`TimeSpec`。固定dt・開始時刻、端点A入力と端点電流出力 |
| `step(A_new)` | 次時刻のAを受け、同時刻の電流を返す。currentも同じ値に更新 |
| `close()` | 必要な資源を解放。資源がなければ省略可能 |

AはSALMONのAc=A/c、電流は実験室座標・原子単位の**電荷電流密度j**です。
初期電流は物理的な値を返し、単なるプレースホルダにしません。
`TimeSpec(dt, (InputPort('A'),), (OutputPort('j'),), t0=time0)`で宣言します。
内部がSALMONの物質電流Jなら`OutputPort('J_matter')`を明示します。
直列respondと新アダプターは、その宣言に基づいて一度だけj=−Jへ変換します。
単位や符号の物理的意味を数値から推測する機能はありません。

`step`の失敗時は受理済み状態を変更しないよう、計算と検査を済ませてから状態へ反映します。
モデルが内部でMPIを呼ぶ場合、この直列アダプターを使わず後述の集団バックエンドにします。
同じモデル・状態を複数スレッドから同時に呼びません。

## 同じモデルを直列で使う

実行可能な実装は[custom_material.py](examples/custom_material.py)です。
モデルと状態の二つのクラスを、自作の構成則や外部ライブラリ呼出しへ置き換えて使います。
例の履歴則は人工的な接続試験で、実材料のモデルではありません。

```python
from salmonpy.examples.custom_material import ExampleMaterial, waveform
from salmonpy.materials import check_material, bind_material
from salmonpy.crystal import respond

model = ExampleMaterial(gain=2., decay=0.8)
drive = waveform()
report = check_material(model, drive)
response = respond(model, drive)
binding = bind_material(model, key='sample', identity=model.identity)
```

`check_material(initial, drive, rtol=1e-12, atol=1e-14)`は利用者が指定した代表波形だけを使い、
二つの状態を作って、有限値・時刻・入力の所有権・状態の独立性・再現性を確認します。
E波形や明確な不等間隔時刻は作成前に拒否し、累積する時刻ずれも時間発展中に検査します。成功時は不変の`MaterialCheck`を返し、
作成した状態は成功・失敗とも解放します。元のモデルのcloseは呼びません。
二つの状態とその履歴に必要な資源を確保できる小さい試験波形から始めます。

これは物理法則、単位の正しさ、任意の入力、モデルの完全な不変性、失敗の原子性、
MPI集団処理を証明するものではありません。返り値の`not_tested`にも示します。
rtol/atolは同じ手法の再実行の比較用であり、物理精度の目標値ではありません。
bind_materialや本番セッションがチェッカーを暗黙に実行することはありません。

## MPIへの接続

`binding`は既存の`MaterialBinding`です。`ResponseSession`へそのまま渡します。
一つの材料キーには同じbindingインスタンスを再利用してください。
`identity`はアルゴリズム版・全パラメータ・初期データを識別する文字列を、作者が責任を持って指定します。
Pythonのオブジェクトアドレスやランク番号を使いません。異なる設定に同じidentityを付けても
内容の差を自動発見できるとは限りません。実装やパラメータはセッション中に変更しません。

```sh
python -m salmonpy.examples.custom_material
mpiexec -n 2 python -m salmonpy.examples.custom_material --mpi NEW_JSON
mpiexec -n 3 python -m salmonpy.examples.custom_material --mpi ANOTHER_JSON
```

MPI例は場役1ランク、残りを各1ランクの物質群とし、4点を割り当てます。
物質群が複数点を持つことはできます。初期Aは既存サービスと同じ0です。
電流は各点で完成した値を一度だけ回収し、点やランクの数を掛けません。
通常のimportではmpi4pyの読込み・MPI初期化・ランク起動を行いません。MPI実行時だけ
同じMPIに対応した任意依存`salmonpy[mpi]`を使います。手順は[DISTRIBUTED.md](DISTRIBUTED.md)。

新アダプターは明示的な時刻宣言を要求し、古いモデルの宣言省略を推測しません。
`check_material`はtime0/初期Aを試験波形から取り、MPIサービスはtime0=0/初期A=0です。
モデルは両方で使う条件を実装してください。

状態作成・各ステップで契約違反を検出すると例外を返し、失敗したアダプターは再利用できません。
MPIでは全群でセッションが無効になります。既に進んだ他の点を巻き戻す保証はありません。
必ず全ランクからセッションを閉じ、その後にコンテキストを閉じます。

## 内部MPI・ネイティブ計算との境界

内部MPIが必要なモデルは既存のprovider/backendの集団契約とsolver_commを使います。
このアダプターは複数ランクの物質群を拒否します。直列コードを各ランクで重複実行して
内部並列と扱うことはしません。[DISTRIBUTED.md](DISTRIBUTED.md)を参照してください。

ネイティブTDDFTは検証済みランタイムを扱う既存の`TDDFTMaterialProvider`と
`TDDFTResponseBackend`を使います。LDはE駆動のため既存の`LDMaterialProvider`と
`LDResponseBackend`を使います。場と物質のMPI全体の巻き戻し、checkpointの一括適用、
MPI全体の公開C ABIは今回の入口に含みません。直列のネイティブ物質プラグインは[C_ABI.md](C_ABI.md)を参照してください。
