# 保存／再開の対応可否と互換性

## 現行の混在計算の保存形式

直接保存は`salmonpy.native-direct/4`、履歴再実行型は`salmonpy.native-replay/2`です。自作呼び出し元の識別情報なし状態は現行形式で扱えますが、保存側と復元側で識別情報の有無を一致させる必要があります。

[English](CHECKPOINTING.en.md)

`salmonpy.checkpointing` は、内蔵・自作の計算コードに共通する保存／復元の検査と実行を提供します。
`execute_checkpoint` は全参加ランクの検査後に、明示的に渡したSAVE／LOAD関数を呼びます。
直列呼び出し元プログラムは`comm=None`、MPI呼び出し元プログラムは自分のコミュニケータを渡します。ソルバーの内部型や継承は不要です。
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
`decomposition`は正の整数tupleで、積を実際の群サイズに合わせます。軸の順序も仕様の一部で、
`(2,1)`と`(1,2)`は異なります。実装・固定モデル入力・計算プログラムの識別には`identity`へ検証した識別子や
ハッシュを渡します。空のidentityは記録から省略します。保存側と復元側で宣言を揃えてください。
検査対象は明示された宣言です。モデルの識別と保存ファイルの完全性は呼び出し元プログラム側で管理します。

## 計算処理独自の保存／復元を実行する

`CheckpointTarget`へ引数なしのコールバックを登録し、`execute_checkpoint`のoperationにsave／restoreを指定します。
パスと保存形式は計算処理が決めます。引数を渡す場合はクロージャーや`functools.partial`を使えます。
ファイル形式の完全な例は`python -m salmonpy.examples.custom_checkpoint NEW_OUTPUT`で実行できます。

同じ入口から、メモリ上の計算処理独自トークンも受け渡せます。

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
不透明なトークン、配列、保存情報はコピーやMPI転送を行わず、そのプロセス内に保持します。所有権は計算処理が管理します。
Noneを含む任意の戻り値を受け付け、失敗は例外で通知します。Falseを返すだけでは失敗扱いになりません。

対応能力とコールバックの存在を全ランクで検査してから処理を始めます。各群の点順を同期したスロットで進めるため、
群ごとの点数が違っても使えます。通常のコールバック例外は全ランクでRuntimeFailureとなり、後続スロットを停止します。
実行済みの状態変更は残ります。MPIを使う計算処理は、自分の群の通信手順を完了してから例外を返してください。
計算処理内で待ち続けるランクの解放や、Fortran STOP・プロセス異常からの回復は共通入口では行いません。
影響した候補を破棄するか、計算処理ジョブを中止するかは呼び出し元プログラムが判断します。MPI初期化・分割・終了も呼び出し元プログラムが管理します。
メタデータ・ファイルのハッシュ照合と完了情報の公開は呼び出し元プログラム／計算処理の方針に従い、計算処理例では既存の一覧とmarkerを維持します。

## 直列の自作の計算コード

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
ファイルの読込やハッシュ照合は行いません。呼び出し元プログラム側で保存内容を認証し、候補が未進行状態であることを確認してから
ローダーを実行入口へ登録します。上の断片はメタデータの検査例で、実際のI/Oとハッシュ検査は実行可能な同梱例に含まれます。

## MPI呼び出し元プログラムと失敗の共有

渡したコミュニケータの全ランクが、同じ順序でAPIを呼びます。`groups`は通常のmappingのlist／tupleで、
nameとranks、必要ならroleとpoint_idsを指定します。ランク番号は渡したコミュニケータ内の番号で、全ランクを
重複なく群へ割り当てます。物質群の点IDは整数または文字列で、文字列表現も含めて所有先を重複させません。
場の群はrole=field、点IDなしとします。物質単独利用には場の群を用意する必要はありません。

各ランクは自分の群の`states`と、場のランクでは`field`の宣言だけを渡します。宣言にはオブジェクトと明示的な
mappingのどちらも使えます。群サイズ、点の網羅、群内宣言、全体の所有ランクと保存時メタデータを照合します。
操作と保存情報は全ランクで一致させ、群内ランクの並びも変更しません。

ランクごとに未検査の設定がある場合は、mappingを渡すとAPI内でエラーを集合的に扱えます。
APIに入る前の宣言構築・読込で例外が出た場合は捕捉し、その文字列を`local_error`として渡して全ランクが検査へ入ります。
1ランクのエラーでも、全ランクがLOADへ進む前に拒否します。コミュニケータは有効なものを呼び出し元プログラムが用意し、
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

失敗時に元の状態を保つ復元には`restore_checkpoint`を使います。同時に切り替える必要があるローカルな可変状態を、
一つの`CheckpointState`へ保持します。呼出中は計算を停止し、成功後はこのホルダーから状態を参照してください。
`stage_restore()`は独立した候補へ読込・検証を行って返します。元の状態を書き換えたり、候補と可変の数値資源を
共有したりしないことが仕様です。元と同じオブジェクト、候補の使い回し、Noneの返却は拒否します。

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
計算処理の独立性仕様の下で、元のオブジェクトを破棄コールバックへ渡すことはありません。

`execute_checkpoint(..., operation='restore', atomic=True, live=live)`も同じ保護経路です。
`execute_checkpoint`の既定は`atomic=False`で、元の状態を保つ保証はありません。保護付き復元には
`restore_checkpoint`を使います。即時LOADだけを持つ対象は、保護経路ではLOADを呼ぶ前に拒否します。
LD/SBEはnew_stateとsnapshot復元で候補を作れます。自作JSONファクトリーは実行可能なサンプルを参照してください。

### 計算処理数値候補（実験的）

明示的に再ビルドしたPythonから直接呼ぶFDTD・標準TDDFT・HSEには、checkpoint_protocol=1の候補ABIがあります。
checkpointのmodeは0=保存、1=即時復元、2=独立した数値候補の読込、3=準備済み候補の反映、4=破棄です。
候補は1 計算処理 lifetimeにつき1個です。元の計算を止め、全群で2を完了・合意してから3を呼びます。
失敗時は全群で4を呼び、元の状態から継続できます。HSE履歴も候補に保持し、3ではファイルを再読込しません。
古いライブラリには対応を推測せず、版確認で拒否します。

内部ヘルパー`salmonpy._native_checkpoint.NativeRestoreCandidate`と
`tools/test_native_restore_candidates_mpi.py`が共通APIとの接続を示します。
stage_restore／discard_restoreには候補のstage／discardを結び、共通APIで全候補の成功を確認した後、
各群でpublishを呼び、保存時のPython側A・電流・時計も準備済みデータで合わせます。
一般のホルダー切替だけでは計算処理反映は完了しません。この内部仕様は安定した公開C ABIではありません。

独立しているのは可変の数値候補です。MPI・Libxc・FFT等の初期化済み資源は既存のコンテキストに残り、
複数計算の同時実行を可能にするものではありません。反映後に旧参照で旧計算を継続することもできません。
反映は準備済みメモリを使いますが、一部の代入は割当を伴い、メモリ枯渇・プロセス異常は回復保証外です。
通常のmpi_mixed_nativeへの接続は次節を参照してください。互換なMPI分割・入力・ハッシュの検査も引き続き必要です。

保証は、参加MPIランクが正常で、正しい独立候補ファクトリーが扱える例外の範囲です。
プロセス異常、通信断、Fortran STOP、準備中に元を変更・共有する計算処理からの回復は対象外です。

## 既存の混在計算

混在サンプルも共通の検査・実行APIを使い、ソルバー固有の未進行状態・C/Fortranの呼び出し規約の確認を維持します。
保存データは形式、入力の識別情報、MPI分割、ファイルのハッシュを照合してから読み込みます。
自作のdirect-state/1サンプルアダプターは、`checkpoint_capability`へオブジェクト／mappingを公開できます。
各部品が宣言するcheckpoint_format／checkpoint_modes／checkpoint_decompositionも参照します。
コールバック、時計、A・電流の仕様は[COUPLED.md](COUPLED.md)を参照してください。

現在は同じ分割と宣言された形式・識別子を照合します。保存内容の変換、MPI分割変更時の再配置、
汎用アーカイブのI/O、失敗時の原子的な巻き戻しは計算処理／呼び出し元プログラム側の担当です。
共通入口は渡された関数を実行します。保存内容の正しさは保存側の責任で、即時LOADだけでは失敗時の元状態保持を保証しません。

## 通常の混在サンプルの失敗時に元の状態を保つ復元

`mpi_mixed_native --direct-resume PATH`は失敗時に元の状態を保つ復元を既定にします。初期化済みの場・標準TDDFT・
HSEは候補対応計算プログラムが必要です。場の候補にはPMLと全LD材料のP/J・更新履歴も含まれます。
SBEは新しい独立ハンドルへsnapshotを復元します。Python側の時計・電流・A・応答アダプター・所有リストと
表示用のframes.jsonも反映前に準備し、全MPIランクの成功確認後に反映します。追加のSCFや時間発展は行いません。

同じセッションを停止した状態なら、内部の`examples._direct_checkpoint.checkpoint(..., restore=True)`は
進行済み状態にも使えます。準備失敗は全ランクへAtomicRestoreErrorとして返り、元の場・LD・電子状態・
時計・参照・failedフラグを維持します。全ランクで同じ判断を行い、そのまま計算を続けるか、原因を除いて
同じ保存データで再試行できます。入力・履歴の検査失敗もLOAD前に拒否します。

成功後はsession.statesから現在の物質状態を取得してください。古いSBEハンドルは閉じます。
旧ハンドルの解放エラーはcoupled.checkpoint_cleanup_errorsへ分けて報告し、成功した復元を失敗として扱いません。
初期化済み計算処理資源は引き続き単一lifetimeで、反映時のメモリ枯渇・プロセス異常・通信断は保証外です。

古いABI、プロセス分離版HSE、自作の即時LOADのみの対象は、保護経路で読込前に拒否します。
通常サンプルは保護付き復元だけを使います。自作バックエンドで即時LOADを明示的に実行する場合は、共通APIの`execute_checkpoint(..., operation="restore", atomic=False)`を使います。
この指定には途中失敗時の元状態保持の保証がありません。自動で切り替えません。

検査用の`tools/test_mixed_protected_restore_mpi.py`は通常サンプルの引数を受け取り、後の電子状態の候補読込後に
1ランクを意図的に失敗させます。場内2種類LD・SBE・内部MPI TDDFT／HSEの元状態保持、
元からの継続、再試行後の連続実行との一致を確認します。リポジトリ内で使うときだけ`--kernel LIB`で
明示的にビルド済みLD/SBEライブラリを選べます。wheel利用時は不要です。

進行済みセッションを巻き戻しても、既に書いたSALMONの計算出力ファイルは切り詰めません。保存済みframesと再開後の
step戻り値で時系列を管理してください。ファイル出力を分ける場合は、新しいMPIジョブ・出力先で再開します。

復元用の別状態の同時保持を含むメモリ助言は[MEMORY.md](MEMORY.md)を参照してください。通常混在サンプルの予算超過は警告のみで、失敗時に元の状態を保つ復元を拒否しません。

## 混在サンプルのIDと再開

`python -m salmonpy.examples.mpi_mixed_native` は `--run-id`、論理点ID、物質ID、物質定義・初期状態のSHA256、MPI所有ランク、LDの材料番号対応を保存します。`mixed-example` は試験用の既定ラベルです。実験ごとに明示した `--run-id` を使用してください。定義の識別には出力先や一時ディレクトリを含めません。

`identity.json`、通常実行の `result.json` と `point-response` アーカイブに同じ識別情報を引き継ぎます。直接保存の外側の形式は `salmonpy.native-direct/4`、再計算方式は `salmonpy.native-replay/2` です。数値状態の内部形式と更新順は変更しません。直接復元はIDとMPI分割を数値状態の読み込み前に照合し、保存時は全ランクでのID一致を検査します。

現在のAPIと保存形式を使用してください。混在計算サンプルは常に物質・初期状態の識別情報を保存します。

## 点ID付きのメモリsnapshot

`PointState.snapshot()` は内部計算処理のsnapshotを `salmonpy.point-snapshot/1` の不変な点情報で包みます。`state.restore(saved)` は実験・点・物質定義・初期状態のIDとパッケージ版を照合してから計算処理を呼びます。別の点や実験、裸の計算処理tokenは拒否されます。明示した同一run/point IDであっても、計算処理のモデル・時間刻み等の整合性検査は引き続き必要です。

これはプロセス内の巻き戻し用であり、分岐・複製やディスク保存、MPI全体のトランザクションではありません。snapshot未対応の自作の計算コードでは `SnapshotError` となります。計算処理の復元失敗時の保証はその実装に従います。点ID不一致の拒否では内部のrestoreを呼びません。個別状態のsnapshot／restoreは、その状態だけを扱います。

`python -m salmonpy.examples.point_ld /new/output/directory` は別の点への復元拒否と、同一点での復元後の電流一致も確認します。

## Pythonから直接呼ぶTDDFTの更新前エラー拒否（限定対応）

新しく準備・ビルドした標準embedded TDDFTには `sp_tddft_error_protocol() == 1` があり、以下を結合先の物質MPI群だけで照合します。10: 不正・非有限入力／古いA／入力ファイル不在、11: ランク間不一致、12: 状態未作成・使用済み初期化、13: ステップ上限、14: 未確定の復元用の別状態あり。これらは数値更新前に全ランク同じコードで返り、既存の最後に成功した計算状態を変更しません。初期化前の10/11は初期化から終了までの期間を消費せず、Pythonからも入力を修正して初期化を再試行できます。失敗時の作業ディレクトリは診断用に残ります。

Pythonは版付きの仕様と全ランク一致を確認できた拒否だけを `InputError` とします。旧ライブラリ・未知のコード・非有限な戻り値は `RuntimeFailure` のままで、状態保持を推定しません。MPI群の全員が同じ入口へ参加する必要があります。無効なポインタ、MPI・プロセス異常、確保失敗は対象外です。本項の更新前検査ではソルバー深部を扱いません（次節の限定した深部復帰を参照）。初期化に入った後は引き続き1プロセス1初期化から終了までの期間です。HSEとPythonから直接呼ぶFDTDのエラー仕様はこの変更では拡張していません。

この保証は個別TDDFT状態の更新前拒否です。結合ドライバで他の部品が既に進んだ場合、全体の巻き戻しを保証せず、従来どおりセッションを失敗扱いにします。通常計算中の一般的な原子性が完成したとは判断しないでください。

## TDDFT深部の数値異常と直接巻き戻し

新しく準備・ビルドした標準Pythonから直接呼ぶTDDFTは `sp_tddft_recovery_protocol() == 1` を持ち、各ステップの数値状態をメモリへ退避します。伝播後の波動関数（20）、密度（21）、Hartree場（22）、XCポテンシャル／エネルギー（23）、電流／全エネルギー（24）を検査し、保護経路の伝播器エラー（25）も返します。検査は物質MPI群全体で合意してから戻ります。エラー時は波動関数・密度・場・伝播履歴・XC状態・カウンタ・管理フラグを最後に成功した計算状態へ直接戻し、伝播の再実行をしません。失敗したステップの通常RT出力も追加しません。

Pythonの `NativeTDDFTBackend` が返す状態の `advance(A_old,A_new)` は、版付き仕様と全ランク一致を確認できた復帰だけを `NumericalError` で通知します。`A` と `steps` はそのままで、原因を修正して同じ受理時刻から再試行できます。不正な波形を自動的に弱めたり、固定dtを変更したりしません。旧ライブラリや未知のエラーは `RuntimeFailure` の保守的な扱いを維持し、復帰を推定しません。

監査付きの`execution_mode="checked"`で、対応する計算ライブラリを使う場合の範囲です。

| 電子状態 | 対応条件 |
|---|---|
| 非スピン分極 | 周期・固定イオン・中点伝播のPZ、TB-mBJ（自動／固定係数）、Libxc-PZ/PBE |
| SOC・非共線スピン | 周期・固定イオン・CPUのPZ／Libxc-PZ/PBE。保持GSからの`gs.tddft_model(...)`など、SOC対応の入口を使う |

SOCのSCFで補助磁場を使った場合は、TDDFTを始める前に`state.set_scf_spin_field([0., 0., 0.])`でゼロへ戻し、無磁場でSCFを続けて収束と残差・磁化の変化を確認してください。補助磁場が非ゼロのままでは`gs.tddft_model(...)`は開始を拒否します。磁場の解除だけでは無磁場の収束を確認したことにはなりません。手順は[NATIVE_GS.md](NATIVE_GS.md#socのscfを補助磁場で整える)を参照してください。

SOCでも、対応する数値異常からの復帰は波動関数・密度・磁化／スピン流・電流出力を保持し、同じ完了ステップから再試行できます。CPU空間分割SOCには直交セルと`sp_tddft_soc_spatial_protocol()==1`が必要です。試験条件は[NICKEL_SOC.md](NICKEL_SOC.md)を参照してください。
`execution_mode="performance"`では、この自動検査・自動復帰を保証しません。必要なら明示的な保護付きステップか保存状態へ戻します。
+U、jellium、イオン運動／力出力、Gram–Schmidt、projection/unfold、微視的Maxwell結合などの未対応設定には、この保証を適用しません。退避・検査のメモリと時間の費用を見込み、予算超過は警告として扱います。

これはTDDFT群の状態保証です。FDTD・他の物質がすでに進んだ混在ステップ全体を自動巻き戻しする仕様ではありません。MPI／プロセス停止、確保失敗、未制御の外部ライブラリSTOP、初期化深部のSTOPや出力装置障害は未対応です。試験専用の一部ランク障害注入は `SALMONPY_RECOVERY_TESTING` を明示して別ビルドした場合のみで、本番ビルドには注入用C入口がありません。

## HSEの深部数値異常からの復帰

新しい `NativeFastHSEBackend` 用ライブラリは `sp_fast_hse_recovery_protocol() == 1` を宣言します。対応範囲は周期・非スピン分極・固定イオンのHSE06、Taylor4／学習PT-CN、物質群内k点分割のMPI1／2です。通常版HSEへの切替はしません。

各ステップ前に電子状態とACE・因子／PT-CN履歴・交換キャッシュ・カウンタを退避します。FFTプランやMPIコミュニケータなどの所有資源はコピーしません。波動関数・密度・Hartree・XC・端点の有限性検査に加え、対応したACE演算エラー、PT-CNの収束／Gram／射影検査、履歴更新エラーを物質群全体で合意して返します。最後に成功した計算状態を直接戻し、過去ステップや交換履歴を再計算しません。失敗試行の通常RT出力を追加せず、A・時計・電流・エネルギーを維持します。診断ログは残ることがあります。

Pythonは宣言された版と全ランクの既知エラーを確認できる場合だけ `NumericalError` を返し、同じ受理時刻から修正した入力で再試行できます。旧ライブラリ・未知エラーは `RuntimeFailure` の保守的扱いを維持します。dtや入力振幅を自動変更しません。退避と検査のメモリ・時間費用は増え、予算超過は警告のみです。

これはHSE物質群の保証で、FDTDや他の物質を含む混在ステップ全体の自動巻戻しではありません。初期化深部、確保失敗、プロセス／MPI停止、未制御の外部STOP、出力装置障害は対象外です。異常を捕捉できない場合に継続を推定しないでください。

検証用 `tools/test_fast_hse_deep_recovery_mpi.py BUILD GS NEW_OUTPUT --steps 32 --fault-step 13` は2ランク、OMP／BLAS1、外側のtimeoutで実行します。`--split` は4ランクの独立2群です。`--overflow-A 1e100` はPT-CNの実数値発散試験用で、科学計算の入力ではありません。`--faults` は `SALMONPY_RECOVERY_TESTING` を明示した別ビルドだけに使い、本番には注入用入口を含めません。検証記録は `docs/fast-hse-deep-recovery-validation-2026-11-17.md` です。合成小規模系の実装検証であり、物理精度の収束保証ではありません。

## Pythonから直接呼ぶFDTDの数値異常復帰

新しい標準呼び出し元プログラムが管理するMPIライブラリは `sp_eh_recovery_protocol() == 1` を宣言します。`NativeEHField.advance` は対応する検査付き入口を選び、場・PML・場内の全LD・電流格子・時計をメモリ退避します。E更新後（20）、LD更新後（21）、H更新後（22）、点サンプル／積算A／通常観測出力の演算（23）の異常を場MPI群で合意し、退避した状態へ直接戻します。波の再伝播はしません。既知の復帰では `NumericalError` を返し、PythonのA・E・時計・failedも維持するため、入力を修正して再試行できます。失敗試行の通常観測ファイルは追加しませんが、診断・進行表示は残る場合があります。

対象は現在のFDTDのE/H更新平面波経路と場内LDです。線形応答専用経路、ASE／ART、周波数別観測は更新前に拒否します。バッチFDTDは従来経路を維持します。新しい入口を持たない旧バイナリ・未知エラーでは復帰を推定せず、従来の失敗扱いです。初期化・確保・MPI／プロセス停止・外部STOP・出力装置障害は保証外です。場群単体の保証であり、`CoupledStepper`／`ResponseSession`全体を自動巻戻し・自動再試行する機能ではありません。退避と検査の費用を見込み、メモリ予算は引き続き警告に使います。

`tools/test_eh_deep_recovery_mpi.py BUILD NEW_OUTPUT` を1／2ランク、OMP／BLAS1、外側timeoutで実行します。`--split` は4ランクの独立2群です。`--faults` は明示した別の `SALMONPY_RECOVERY_TESTING` ビルドだけに使います。本番には注入入口を含めません。数値オーバーフローの試験値は科学入力ではありません。

ステップ失敗時の計算処理数値／通常出力の巻戻しは [NATIVE_STEPPING.md](NATIVE_STEPPING.md) を参照。checkpointのLOADをステップ退避の代用にせず、両者の同時利用は拒否します。反映でownerが変わったらstep targetを作り直します。

## SCFメモリ保存を共通APIから使う

計算済みの生存中の `NativeGSState` は `checkpoint_capability` と
`checkpoint_target(snapshot=None)` を公開します。内部MPI群の全ランクが同順序で
呼びます。接続時に計算プログラムのプロトコルと生存SCFのメタデータを確認しますが、保存や数値反復は行いません。
形式は `salmonpy.native-scf-memory/1` です。識別子は物質・モデル・計算プログラムと
生存計算処理所有者のトークンに結び付き、同じ物質入力の別状態とも区別します。
MPI分割は実際の `(k, orbital, x, y, z)` です。

```python
from salmonpy.checkpointing import execute_checkpoint
comm = state.session.comm
plan = [{"name": "scf", "ranks": list(range(comm.size)), "point_ids": ["film"]}]
saved_run = execute_checkpoint(comm, plan, states={"film": state.checkpoint_target()})
with saved_run["states"]["film"] as payload:
    state.continue_scf(steps=5)
    execute_checkpoint(comm, plan, operation="restore", saved=saved_run["compatibility"],
                       states={"film": state.checkpoint_target(payload)})
    state.continue_scf(steps=5)
```

通信するのは互換性メタデータだけです。`NativeSCFSnapshot` は各ランクのFortran
メモリに残り、明示的に終了します。JSONだけでは復元できず、ディスク保存／プロセス
再開でもありません。通常PZ／Libxc-PZ/PBEと全k点CPUHSE06の対応範囲・
数値実装は従来どおりです。

`execute_checkpoint(operation="restore")` は引き続き**直接復元**です。
SCF candidate protocol1を持たない旧計算プログラムでは `stage_restore` がないため、
`restore_checkpoint` は数値コールバック前に拒否します。対応計算プログラムでは次節の
独立した数値候補による準備を使えます。`execute_checkpoint` では後の部品が実行途中で失敗すると、先の部品だけ
復元済みになり得ます。全状態のロールバックを保証するものではありません。
異なる所有者・終了済みトークン・互換性不一致は数値コールバック前に拒否します。
候補作成の代わりに破壊的な復元へ落とさないでください。
`python -m salmonpy.examples.scf_resume BUILD_JSON NEW_OUTPUT` は共通入口を使います。
HSE用ビルドでは `--xc hse06` を付けます。

### SCF候補の保護付き準備とFortran/Cで保持する計算状態への反映

`sp_gs_scf_candidate_protocol()==1` の対応計算プログラムでは、
`state.checkpoint_target(payload)` が独立した**数値候補**を作れます。
別の計算可能なGSソルバーを作る方式ではなく、FFT・XC・MPI資源は元の生存状態に
保持します。Fortranでスナップショットを切り離した領域へコピーし、SCF初期化・軌道
計算・交換演算を行いません。候補がある間は元の状態を予約し、メタデータの参照は
できますが、継続・直接復元・保存・終了は反映または破棄まで拒否します。
呼び出しは同じ内部MPI群の全ランクがメインスレッドで同順序に行います。
メモリ内・同じ生存状態・対応XC／MPI分割の制約は変わりません。

```python
from salmonpy.checkpointing import CheckpointState, restore_checkpoint
holder = CheckpointState(states={"film": state})
prepared = restore_checkpoint(
    comm, plan, live=holder, saved=saved_run["compatibility"],
    states={"film": state.checkpoint_target(payload)})
# 全参加者の準備が成功。この時点まで元の計算処理数値状態は変わりません。
for candidate in prepared["states"].values():
    candidate.publish()  # 内部MPI群のcollective。戻り値は元のstate。
state.continue_scf(steps=5)
```

holderに入るのは準備済み候補で、独立に計算を走らせる状態ではありません。
全計算処理候補を群内の同じ順序で反映してから計算を再開してください。
準備途中で処理可能な例外が出ると、途中の候補も含めて破棄して予約を解除し、
元の数値状態を保持します。元の軌道を継続でき、新しいtargetで復元を再試行できます。
手動準備を中止する場合は全体の準備成功確認前に `candidate.discard()` を呼びます。
一部反映後のエラーを無視して計算を続けないでください。計算処理反映では配列代入による
メモリ確保があり得るため、致命的な確保失敗・プロセス／通信故障・反映途中の任意の
異常を元に戻す保証ではありません。反映／破棄にファイル操作や時間発展はありません。
元スナップショットは準備・反映が終わるまで保持し、その後明示的に終了します。

`python -m salmonpy.examples.scf_resume BUILD_JSON NEW_OUTPUT --protected` で試せます。
全k点のHSEビルドなら `--xc hse06` も付けます。
必要な候補作成機能を持たない計算プログラムは、保護付き復元では読み込み前に拒否します。即時復元への自動切り替えはありません。

計算処理 SCFの同一ビルド・入力・MPI分割ディスク形式は、
`state.checkpoint_target(directory=...)` でこの共通APIに接続できます。
`scf_disk_capability` はファイルを読まない宣言で、stageがSHA256と入力を検査します。
全候補の準備後に各candidate.publish()を実行してください。失敗後の再試行には新しい
targetを作ります。詳細と2状態の例は [NATIVE_GS.md](NATIVE_GS.md) を参照してください。

計算処理 SCFの準備済み候補はclose／with文で明示的に破棄できます。GS状態のcloseは
その状態のPython所有候補を、snapshotのcloseはそのsnapshot由来の候補を先に破棄します。
C ABIで直接準備した外部所有候補は自動破棄しません。すべて同じMPI呼出し順が必要です。

## 定期保存と、実経過時間を指定した保存・終了

チェックポイントの自動出力は既定では無効です。電流・磁化などの観測出力の間隔とは別に設定します。
`CheckpointSchedule` は、完了した時間ステップの番号を受け取り、指定がある場合だけ保存関数を呼びます。

```python
from salmonpy.checkpointing import CheckpointSchedule

# comm は結合計算全体のMPIコミュニケータ。直列なら省略できます。
# 初期化も制限時間に含める場合、この時点で時計を開始します。
schedule = CheckpointSchedule(
    interval_steps=1000,        # 任意：1000ステップ進むごとに保存して続行
    interval_seconds=600,       # 任意：実経過600秒ごとに保存して続行
    stop_after_seconds=3300,    # 任意：実経過3300秒を超えたら保存して終了
    initial_step=0, comm=comm)

# 初期化／復元後、各ランクが同じ完了ステップで呼びます。
for step in range(restored_step + 1, final_step + 1):
    advance_one_step()
    if schedule.after_step(step, save_all_states):
        break
```

`advance_one_step` と `save_all_states(step)` は呼び出し側で用意する関数です。保存先はステップ番号ごとの新規パスにしてください。結合計算では場・LD・電子状態をまとめて保存する既存の入口を使います。単一の電子状態だけを保存して結合計算を終了してはいけません。保存形式は各計算部品が担当します。

- すべての引数を省略すると、保存も追加のループ内MPI通信も行いません。
- 複数の間隔を指定した場合は、どれか一つに達した時点で一度保存します。間隔は最後に成功した保存から数えます。
- 時間は物理的な計算時間ではなく、作成時からの実経過秒です。MPIではランク0の時計で決め、同じ完了ステップで保存します。時間指定がある場合も、ステップの途中では中断しません。
- `after_step` の戻り値が `True` なら、保存が成功しています。呼び出し側はループを抜け、状態を閉じて終了してください。ライブラリがプロセスやMPIを強制終了することはありません。
- 保存失敗は全ランクへ通知され、正常終了にはなりません。同じステップで再試行できます。途中のファイルの扱いは保存部品の仕様に従い、必要なら別の保存先を使います。
- ジョブ制限時間より、最大1ステップと保存処理・終了処理の所要時間だけ早めに設定してください。タイマー作成前の準備時間は含みません。

`mpi_mixed_native`、`mixed_hse_models`、`molecular_feedback` のサンプルでは、`--checkpoint-interval-steps N`、`--checkpoint-interval-seconds SECONDS`、`--stop-after-seconds SECONDS` を指定できます。新規の `checkpoint-00001000` のようなディレクトリへ内部状態を直接保存します。既存のステップ番号指定による保存・終了も利用できます。

出力連結は、別プロセスでのディスク再開が生成したNi SOCの1・2ランク出力で確認しました。電流は全8ステップ、磁化・スピン流は保存されている2時刻で連続実行と完全一致しました。これは保存／再開の接続検証であり、物理的収束の保証ではありません。

### 時間指定で保存し、次のジョブで再開する

準備済みの計算プログラム・GSを使う通常混在サンプルの例です。4ランクは場1・SBE1・TDDFT2の既定配置です。下記のパスは手元の準備済みデータへ置き換え、計算ノードで実行してください。

```sh
mpiexec -n 4 python -m salmonpy.examples.mpi_mixed_native \
  BUILD_JSON RUNTIME GS_BUNDLE first-job --steps 128 \
  --checkpoint-interval-steps 50 --stop-after-seconds 3300
```

保存先は`first-job/checkpoint-00000050`のように完了ステップ番号を含みます。時間指定終了時には実際の保存先が表示されます。そのパスを次のコマンドの`SAVED_CHECKPOINT_DIRECTORY`に指定します。時間制限で止まるステップ番号を事前に仮定しないでください。

```sh
mpiexec -n 4 python -m salmonpy.examples.mpi_mixed_native \
  BUILD_JSON RUNTIME GS_BUNDLE second-job --steps 128 \
  --direct-resume SAVED_CHECKPOINT_DIRECTORY --stop-after-seconds 3300
```

入力、GS、計算プログラム、総ステップ数、MPI配置を同じにし、新しい出力先を指定します。`--resume`は履歴を再実行する別方式なので、時間指定の直接保存には`--direct-resume`を使います。分子結合と`mixed_hse_models`の直接保存には`--resume`を使います。それぞれのサンプルの引数を確認してください。

全ステップが指定秒より早く終わった場合、時間指定による最終保存は発生しません。定期保存は指定間隔に達した分だけです。この例は結合検査用サンプルで、最終ステップまで達した通常実行では独立計算との比較も行います。実行時間を見積もる際に含めてください。

`after_step`は、保存して終了する場合だけ`True`です。保存しない場合と、定期保存して続行する場合は`False`です。保存失敗は例外で通知します。
