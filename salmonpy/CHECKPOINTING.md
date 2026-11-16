# 保存／再開の対応可否と互換性

[English](CHECKPOINTING.en.md) · Topic: `checkpoint_compatibility` · Paired revision: `2026-11-17.1`

`salmonpy.checkpointing` は、内蔵・自作バックエンドに共通する保存／復元前の検査を提供します。
自分のSAVE／LOAD関数を呼ぶ前に使い、全参加ランクの検査後に独立したJSONメタデータを返します。
直列ホストは`comm=None`、MPIホストは自分のコミュニケータを渡します。ソルバーの内部型や継承は不要です。
importや検査時にMPI初期化、ファイル操作、コンパイル、ダウンロード、時間発展、保存／復元のコールバックは行いません。

## 対応能力の宣言

```python
from salmonpy.checkpointing import CheckpointCapability

capability = CheckpointCapability(
    "my-material-json/1", ("save", "restore"), (1,),
    identity="implementation-and-fixed-model-input-digest",
)
record = capability.to_dict()
assert CheckpointCapability.from_dict(record) == capability
```

宣言は不変です。`format`は版付きの保存形式、`modes`はsave／restoreの一方または両方を含むtupleです。
`decomposition`は正の整数tupleで、積を実際の群サイズに合わせます。軸の順序も契約の一部で、
`(2,1)`と`(1,2)`は異なります。実装・固定モデル入力・ランタイムの識別には`identity`へ検証した識別子や
ハッシュを渡します。空のidentityは既存native-direct/2との互換性のため記録から省略します。
検査対象は明示された宣言です。モデルの識別と保存ファイルの完全性はホスト側で管理します。

## 直列の自作バックエンド

同梱例は、自作モデル、JSON保存、SHA256検証、独立した状態、連続実行との比較を一通り実行します。

```sh
python -m salmonpy.examples.custom_checkpoint NEW_OUTPUT
```

17ステップで保存し、残り15ステップだけ進めて32ステップの連続実行と電流が完全一致します。
モデルを変更した場合は、ローダーを呼ぶ前に拒否します。漸化式は接続確認用の合成データであり、物質の物理予測ではありません。

この例のうち、対応可否・互換性の検査部分は次のように書けます。

```python
from salmonpy.checkpointing import check_checkpoint_compatibility
from salmonpy.examples.custom_checkpoint import UserCurrent

groups = [{"name": "user", "role": "material", "ranks": [0], "point_ids": ["film"]}]
state = UserCurrent()
saved = check_checkpoint_compatibility(
    None, groups, states={"film": state.checkpoint_capability}, operation="save")
# ホスト側で本体と返されたメタデータを保存し、ハッシュを検証する。
candidate = UserCurrent()
check_checkpoint_compatibility(
    None, groups, states={"film": candidate.checkpoint_capability},
    operation="restore", saved=saved,
    required_files=("state.json",), available_files=("state.json",))
# 検査成功後、ホスト側でcandidate.load(検証済みファイルのパス)を呼べる。
```

`available_files`はファイル名の一覧です。APIは安全な相対POSIXパスと必須ファイルの存在を一覧上で確認します。
ファイルの読込やハッシュ照合は行いません。ホスト側で保存内容を認証し、候補が未進行状態であることを確認してから
ローダーを呼びます。上の断片はメタデータの検査例で、実際のI/Oとハッシュ検査は実行可能な同梱例に含まれます。

## MPIホストと失敗の共有

渡したコミュニケータの全ランクが、同じ順序でAPIを呼びます。`groups`は通常のmappingのlist／tupleで、
nameとranks、必要ならroleとpoint_idsを指定します。ランク番号は渡したコミュニケータ内の番号で、全ランクを
重複なく群へ割り当てます。物質群の点IDは整数または文字列で、文字列表現も含めて所有先を重複させません。
場の群はrole=field、点IDなしとします。物質単独利用には場の群を用意する必要はありません。

各ランクは自分の群の`states`と、場のランクでは`field`の宣言だけを渡します。宣言にはオブジェクトと明示的な
mappingのどちらも使えます。群サイズ、点の網羅、群内宣言、全体の所有ランクと保存時メタデータを照合します。
操作と保存情報は全ランクで一致させ、群内ランクの並びも変更しません。

ランクごとに未検査の設定がある場合は、mappingを渡すとAPI内でエラーを集合的に扱えます。
APIに入る前の宣言構築・読込で例外が出た場合は捕捉し、その文字列を`local_error`として渡して全ランクが検査へ入ります。
1ランクのエラーでも、全ランクがLOADへ進む前に拒否します。コミュニケータは有効なものをホストが用意し、
初期化・分割・終了を管理します。APIはrank／size／allgatherを使い、mpi4pyを自分でimportしません。

```python
from mpi4py import MPI
from salmonpy.checkpointing import check_checkpoint_compatibility

comm = MPI.COMM_WORLD  # ホストが初期化済み。この例は3ランク。
groups = [
    {"name": "field", "role": "field", "ranks": [0], "point_ids": []},
    {"name": "material", "ranks": [1, 2], "point_ids": ["film"]},
]
field = {"format": "user-field/1", "modes": ["save", "restore"], "decomposition": [1]}
material = {"format": "user-current/1", "modes": ["save", "restore"],
            "decomposition": [2, 1], "identity": "fixed-model-and-runtime-digest"}
saved = check_checkpoint_compatibility(
    comm, groups, operation="save", field=field if comm.rank == 0 else None,
    states={} if comm.rank == 0 else {"film": material})
```

復元時はoperation=restoreとし、全ランクで同じ保存メタデータを渡します。必須ファイルの一覧はランクごとに用意できます。
エラーは`salmonpy.errors.InputError`（ValueErrorの派生）で通知し、成功時は`salmonpy.direct-compatibility/1`のJSONを返します。

## 既存の混在native計算

混在サンプルも共通APIを使い、ソルバー固有の未進行状態・native ABIの確認を維持します。
既存native-direct/2の互換情報、固定したnativeペイロード、数値の更新順は変えません。
自作のdirect-state/1サンプルアダプターは、`checkpoint_capability`へオブジェクト／mappingを公開できます。
従来のcheckpoint_format／checkpoint_modes／checkpoint_decompositionも使えます。
コールバック、時計、A・電流の契約は[COUPLED.md](COUPLED.md)を参照してください。

現在は同じ分割と宣言された形式・識別子を照合します。保存内容の変換、MPI分割変更時の再配置、
汎用アーカイブのI/O、失敗時の原子的な巻き戻しはバックエンド／ホスト側の担当です。
APIが任意の保存実装を実行・認証することはありません。互換情報のない旧native-direct/1は、元の対応版で扱います。
