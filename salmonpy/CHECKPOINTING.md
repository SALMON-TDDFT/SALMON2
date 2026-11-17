# 保存／再開の対応可否と互換性

[English](CHECKPOINTING.en.md) · Topic: `checkpoint_compatibility` · Paired revision: `2026-11-17.1`

`salmonpy.checkpointing` は、内蔵・自作バックエンドに共通する保存／復元の検査と実行を提供します。
`execute_checkpoint` は全参加ランクの検査後に、明示的に渡したSAVE／LOAD関数を呼びます。
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

## バックエンド独自の保存／復元を実行する

`CheckpointTarget`へ引数なしのコールバックを登録し、`execute_checkpoint`のoperationにsave／restoreを指定します。
パスと保存形式はバックエンドが決めます。引数を渡す場合はクロージャーや`functools.partial`を使えます。
ファイル形式の完全な例は`python -m salmonpy.examples.custom_checkpoint NEW_OUTPUT`で実行できます。

同じ入口から、メモリ上のバックエンド独自トークンも受け渡せます。

```python
from salmonpy.checkpointing import CheckpointCapability, CheckpointTarget, execute_checkpoint

groups = [{"name": "user", "ranks": [0], "point_ids": ["film"]}]
cap = CheckpointCapability("user-memory/1", ("save", "restore"), (1,), identity="model-v1")
state = {"value": 17}
saved_run = execute_checkpoint(
    None, groups, operation="save",
    states={"film": CheckpointTarget(cap, save=lambda: dict(state))})
payload = saved_run["states"]["film"]
candidate = {}
execute_checkpoint(
    None, groups, operation="restore", saved=saved_run["compatibility"],
    states={"film": CheckpointTarget(cap, restore=lambda: candidate.update(payload))})
assert candidate == state and candidate is not state
```

戻り値には`compatibility`のJSON、プロセス内の`field`の戻り値、文字列の点IDをキーにした`states`の戻り値が含まれます。
不透明なトークン、配列、保存情報はコピーやMPI転送を行わず、そのプロセス内に保持します。所有権はバックエンドが管理します。
Noneを含む任意の戻り値を受け付け、失敗は例外で通知します。Falseを返すだけでは失敗扱いになりません。

対応能力とコールバックの存在を全ランクで検査してから処理を始めます。各群の点順を同期したスロットで進めるため、
群ごとの点数が違っても使えます。通常のコールバック例外は全ランクでRuntimeFailureとなり、後続スロットを停止します。
実行済みの状態変更は残ります。MPIを使うバックエンドは、自分の群の通信手順を完了してから例外を返してください。
バックエンド内で待ち続けるランクの解放や、Fortran STOP・プロセス異常からの回復は共通入口では行いません。
影響した候補を破棄するか、nativeジョブを中止するかはホストが判断します。MPI初期化・分割・終了もホストが管理します。
メタデータ・ファイルのハッシュ照合と完了情報の公開はホスト／バックエンドの方針に従い、native例では既存の一覧とmarkerを維持します。

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
ローダーを実行入口へ登録します。上の断片はメタデータの検査例で、実際のI/Oとハッシュ検査は実行可能な同梱例に含まれます。

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

## 復元失敗時に元の状態を維持する

保護付き復元には`restore_checkpoint`を使います。同時に切り替える必要があるローカルな可変状態を、
一つの`CheckpointState`へ保持します。呼出中は計算を停止し、成功後はこのホルダーから状態を参照してください。
`stage_restore()`は独立した候補へ読込・検証を行って返します。元の状態を書き換えたり、候補と可変の数値資源を
共有したりしないことが契約です。元と同じオブジェクト、候補の使い回し、Noneの返却は拒否します。

```python
from salmonpy.checkpointing import (CheckpointCapability, CheckpointTarget,
                                   CheckpointState, restore_checkpoint)

groups = [{"name": "user", "ranks": [0], "point_ids": ["film"]}]
cap = CheckpointCapability("user-memory/1", ("save", "restore"), (1,), identity="model-v1")
from salmonpy.checkpointing import check_checkpoint_compatibility
saved = check_checkpoint_compatibility(None, groups, states={"film": cap})
original = {"step": 5}
live = CheckpointState(states={"film": original})
run = restore_checkpoint(
    None, groups, live=live, saved=saved,
    states={"film": CheckpointTarget(cap, stage_restore=lambda: {"step": 17})})
assert live.states["film"]["step"] == 17
assert run["previous"]["states"]["film"] is original
assert original == {"step": 5}
```

共通入口は全候補を準備し、MPIで失敗を共有し、成功の合意後に各ランクのローカル状態一式を一回の参照代入で切り替えます。
切替中に利用者のコードは呼びません。以前の状態は`run['previous']`へ保持するため、成功後の資源破棄は明示的に行います。
既存の外部変数は旧オブジェクトを指し続けます。バインディングを作り直すか、ホルダーを参照するクロージャーを使ってください。
準備中は元と候補が共存する分のメモリ・資源が必要です。

通常の準備失敗では`AtomicRestoreError`（RuntimeFailureの派生）が`original_preserved=True`を返し、
全ホルダーは元の状態一式を維持します。元の状態から計算を継続できます。任意の`discard_restore(candidate)`で、
準備を試みた候補を逆順の同期スロットで破棄できます。準備中の例外や危険な参照を検出した場合は、
該当する全ランクで`discard_restore(None)`も呼びます。部分的な準備資源とNoneを扱い、例外前に自群の通信手順を
完了してください。候補ファクトリー自身が失敗時の割当資源を片付けてから例外を返す方法も使えます。
破棄時のエラーは`error.cleanup_errors`へ分けて報告します。破棄処理でも元の状態を変更してはいけません。
バックエンドの独立性契約の下で、元のオブジェクトを破棄コールバックへ渡すことはありません。

`execute_checkpoint(..., operation='restore', atomic=True, live=live)`も同じ保護経路です。
従来のexecute_checkpointは互換性のため非原子的な実行を既定とします。通常のLOADだけを持つ対象は、保護経路では
LOADを呼ぶ前に拒否し、自動で従来方式へ切り替えません。現行のプロセス内native FDTD／TDDFT／高速HSE ABIには
独立候補の生成口がないため、同梱native例は`atomic=False`を明示し、従来の失敗時制限を維持します。
LD/SBEはnew_stateとsnapshot復元で候補を作れます。自作JSONファクトリーは実行可能なサンプルを参照してください。
native候補生成への対応は別の残作業です。

保証は、参加MPIランクが正常で、正しい独立候補ファクトリーが扱える例外の範囲です。
プロセス異常、通信断、Fortran STOP、準備中に元を変更・共有するバックエンドからの回復は対象外です。

## 既存の混在native計算

混在サンプルも共通の検査・実行APIを使い、ソルバー固有の未進行状態・native ABIの確認を維持します。
既存native-direct/2の互換情報、固定したnativeペイロード、数値の更新順は変えません。
自作のdirect-state/1サンプルアダプターは、`checkpoint_capability`へオブジェクト／mappingを公開できます。
従来のcheckpoint_format／checkpoint_modes／checkpoint_decompositionも使えます。
コールバック、時計、A・電流の契約は[COUPLED.md](COUPLED.md)を参照してください。

現在は同じ分割と宣言された形式・識別子を照合します。保存内容の変換、MPI分割変更時の再配置、
汎用アーカイブのI/O、失敗時の原子的な巻き戻しはバックエンド／ホスト側の担当です。
共通入口は渡された関数を実行しますが、保存実装の正しさや原子的復元を認証するものではありません。互換情報のない旧native-direct/1は、元の対応版で扱います。
