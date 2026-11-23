# 同じ物質から独立したHSE計算を複数作る

[English](SHARED_HSE.en.md)

同じ物質に異なる光波形を与えたい場合や、複数のマクロ点で同じ物質を使いたい場合の方法です。原子配置・対称性などの入力は共通にし、波動関数と励起・電流の履歴は計算ごとに分けます。一つの計算を進めても、ほかの計算の状態は変わりません。

`NativeFastHSEBackend.open_model(initial, dt, context)` で共通のモデルを用意し、`model.new_state(...)` で元の基底状態から独立した計算を作ります。物質入力の検証・コピーとMPIへの接続は一度だけ行い、実際の電子状態の初期化は最初の `new_state` で行います。

```python
from salmonpy.fast_hse import FastHSEMaterialProvider
from salmonpy.fast_hse_mpi import NativeFastHSEBackend

initial = FastHSEMaterialProvider(gs_directory, max_steps=100).prepare(context)
model = NativeFastHSEBackend(build_json, work_root).open_model(initial, .001, context)
states = []
try:
    with model:
        a = model.new_state('point-a')
        states.append(a)
        b = model.new_state('point-b', dt=.002, capacity=100)
        states.append(b)
        ja = a.step([0., 0., 1e-4])
        jb = b.step([1e-4, 0., 0.])
        # MPIランクごとに異なる共有ファイル名を使う。
        a.checkpoint(rank_checkpoint_path, False)
        restored = model.new_state('point-a')
        states.append(restored)
        restored.checkpoint(rank_checkpoint_path, True)
        j_next = restored.step([0., 0., 2e-4])
finally:
    for state in reversed(states):
        state.close()
    model.close()
```

`context` は `solver_comm` と `point_ids` を持ちます。呼び出し元のParallelContextも利用できます。
全ランクが同じ順序で、メインスレッドから呼び出してください。全ての状態・モデルのleaseを閉じるまで
コミュニケータを維持します。入力・ライブラリの検査は計算処理寿命消費前です。
importだけではMPI・計算処理ライブラリの読み込み、コンパイル、ダウンロードを行いません。

`new_state(point_id=None, *, dt=None, capacity=None)` は、既定でモデルのdtと準備済み最大ステップ数を使います。
point_idを省くと連番整数、指定した文字列・整数IDは `state.point_id` に保持します。
計算処理ハンドルはこれと別です。各状態はA・電流・時計・容量・HSE履歴・出力basenameを所有し、
`state.model_identity` が共通物質スナップショットを識別します。dtやcapacityの変更は兄弟状態や
モデルの既定値を変更しません。元の入力は再利用でき、検証済み作業コピーは変更しないでください。

状態の `step(A)`、`advance(Aold,Anew)`、`checkpoint(path,False)`、
`checkpoint(path,True)`、`close()` は従来と同じ意味です。
直接読込は、dt・容量・バイナリ・GS／入力・MPI分割が対応する新規状態へ行います。
過去ステップの再計算なしに内部履歴を復元します。保存ファイルと同一性の照合は呼び出し元プログラム側で行います。
ログは新しい状態のラベルで開始し、数値復元だけでは過去ログを統合しません。
通常の失敗時に元の状態を保つ時間更新は、その状態の2本の出力だけを退避するため、共通ディレクトリでも
兄弟状態のログを巻き戻しません。

`model.close()` は新規生成を止め、既存状態は継続できます。最後の状態を閉じた時にMPI接続を解放します。
モデルを開いたまま全状態を閉じても、隠れた完全電子状態を保持せずに再生成できます。
モデルのコンテキストマネージャーは生存する状態を閉じません。各状態を `finally` で閉じてください。
計算状態の初期化を一度も試していないモデルの終了は、未使用のPython予約も解放します。
PythonのGCからMPI終了処理は呼びません。

従来の `open_model` の対象は、プロセスごとに一つの不変な**全k点または対称性による縮約k点**モデル、既存backendのk点MPI1／2、
CPU・メインスレッドの範囲です。従来の `create_state` は1寿命のままで、同じプロセスでは
こちらか `open_model` のどちらかを選びます。handle protocol1とmodel protocol1を備えた
計算プログラムを明示的にビルドしてください。異なる入力の同一プロセス実行は後述の明示的レジストリで扱います。
自動ランク再配置とスレッド実行は別作業です。context_mpi_protocol1では初期化時の子コミュニケータを状態終了時に解放します。
所有者のない共通資源はdetachまで保持し、無制限の生成・終了は保証しません。資源予算と物理収束の助言は警告扱いで、入力を自動変更しません。
AはCartesian A/cの原子単位、jは電荷電流密度（SALMON物質電流の負）、dtは原子時間、energyはHartreeです。
x／y／z・円偏光を指定できます。

対応MPIで [examples/shared_fast_hse.py](examples/shared_fast_hse.py) を明示実行できます。
全k点・内部パルスなしの準備済みGSを指定します。

```sh
mpiexec -n 2 python shared_fast_hse.py BUILD_JSON GS NEW_OUTPUT
```

2種類の偏光を進め、新規状態への直接復元を確認する例です。実材料の物理的収束を保証する例ではありません。
下位C ABIと残る共有資源については[計算処理ハンドル](NATIVE_HANDLES.md)を参照してください。

## 通常の混在MPIセッション

`NativeFastHSEBackend(build_json, work_root, shared_model=True)` で明示的に有効化し、
`assemble_native_layers` に `points_per_group={material_id: 2}` を渡します。
`ranks_per_point` は各群のMPIランク数です。全マクロ点の状態とIDを維持し、同じ物質の点を
登場順にまとめます。離れた層の点も同じ群にできます。省略時は従来の一群一点です。
予算超過は警告で、点数や分割を自動変更しません。既知の単一状態TDDFTは群を共有できませんが、
自作の計算コードは独自の生成器を提供できます。共有Fortran製HSEは全k点、または明示的に指定した同じ対称性のモデルが対象です。

`ResponseSession` は任意の `backend.open_group(initial, dt, context)` をローカルbindingごとに
一度呼びます。戻り値の生成器は `new_state(point_id)` と `close()` を持ち、`None` なら従来の
`create_state` を使います。群の全ランクで同じ順序で呼びます。生成器の寿命はセッションが管理し、
状態を閉じてから生成器を閉じます。一つの終了処理が失敗しても他の状態の終了を試みます。
計算処理の再利用だけでは生存セッションを暗黙に共有しません。自作コードに継承は不要です。

通常サンプルは `--hse-points-per-group 2` を受け付けます。

```sh
mpiexec -n 4 python -m salmonpy.examples.mpi_mixed_native \
  STANDARD_BUILD RUNTIME HSE_GS NEW_OUTPUT \
  --layers mixed_shared_hse_layers.json \
  --fast-hse-build SERIAL_HSE_BUILD --fast-hse-mpi-build MPI_HSE_BUILD \
  --hse-ranks 2 --hse-points-per-group 2 --steps 256 --dt .003 --protect-steps
```

JSONは `importlib.resources.files('salmonpy.examples')` から取得し、その中の
`mixed_shared_hse_layers.json` のパスを渡します。LD・SBE・HSE二点を含み、共有時4ランク、
別群では6ランクです。`mixed_shared_tddft_hse_layers.json` は内部2ランクのPBEも含み、
共有時6、別群8ランクです。場を2ランクにする場合は合計を一つ増やし、`--field-grid 1 2 1` を指定します。
計算プログラムは明示的にビルドし、検査済みGSを用意します。配置JSONだけではGS生成・取得をしません。

`--direct-save-at 127` で保存し、新規ジョブに
`--direct-resume SAVED_OUTPUT/direct-checkpoint` を指定して同じ配置・dt・stepsで再開できます。
一ランクに複数バイナリ状態がある場合は `salmonpy.native-direct/4` とランク・点ごとのファイルを使います。
一点の場合の従来名・形式を維持し、版1〜3も読めます。同じバイナリ・入力・MPI分割が必要です。
失敗時に元の状態を保つ復元は独立候補を準備し、後の点で失敗しても場・LD・SBE・全電子状態を維持し、
過去ステップの再計算なしで再試行できます。プロセスクラッシュの復旧や汎用交換形式の保証ではありません。

Taylor／PT-CNで、共有群と別群、独立した直列応答、連続実行と直接再開の電流・場・エネルギーを
比較しました。測定範囲は検証記録を参照してください。物理的収束の保証ではありません。

## 縮約k点・対称性付き共有モデル

通常セッションは `NativeFastHSEBackend(..., allow_symmetry=True, shared_model=True)`、
明示的な `open_model` は `allow_symmetry=True` を使います。対応計算プログラムには
`sp_fast_hse_shared_symmetry_protocol()==1` が必要です。全k点は従来のmodel protocol1とも互換です。
各状態のセル・原子・kメッシュ・保持対称操作は共通不変です。密度対称化テーブルは生成時に同期して
再構築し、時間発展する波動関数・交換の写像／計算ルーチン・学習履歴は別々に所有します。
メインスレッド・群内同順序の呼び出しが必要です。

通常サンプルは `--hse-points-per-group 2` に
`--hse-symmetry --hse-reference-bundle FULL_GS` を加え、HSE入力に縮約GSを渡します。
直列参照は全k点です。明示モデルの例は `--symmetry` を加えると二つの外場がz方向になります
（zを保つ対称群が必要です）。他方向・円偏光も実際の群が保持すれば利用できます。
不適合なAは最後に成功した計算状態を変えずに拒否し、外場の投影や対称群の自動変換はしません。
Taylor strict／learnedとPT-CN off／learnedを内部MPI1／2で検証し、独立実行、全／縮約k点、
新規状態への直接再開を比較しました。異なる対称群の同一モデル混在とMPI分割変更による再開は対象外です。

付属の局在点を持つ3次元混在配置は128ステップ（`--steps 128`、保存は
`--direct-save-at 63`）で確認しました。256ステップでは145ステップ後に局所場が保持群に
不適合となり、失敗時に元の状態を保つ時間更新が元の状態を保って拒否しました。光源がz偏光でも全点の
局所FDTD場がzのままとは限りません。実際の場・配置に適合する対称性か全k点入力を使い、
微小成分の切捨て・外場の回転や投影・利用者の要求時間の自動短縮はしません。

非対称PBE点を混ぜると、z偏光入射でもHSEの保持群に不適合な場を作ることがあります。
今回のPBE混在入力はその拒否例であり、この入力には全k点HSEを使います。
全k点のPBE／HSE混在は新しい計算プログラムでも回帰確認しています。

## 状態ごとのMPI初期化資源

`sp_fast_hse_context_mpi_protocol()==1` を持つ新規ビルドでは、同期的な状態初期化で作る
コミュニケータに固有の所有者IDを付けます。`state.close()` は通信要求キャッシュ、
XC・交換・FFT資源、その状態のコミュニケータの順に解放します。兄弟状態と呼び出し元プログラムの
コミュニケータは維持します。管理表の空き枠を再利用し、接続中は所有者IDを再使用しません。
明示的な計算処理診断 `sp_fast_hse_mpi_resources(int values[2])` は追跡中の子の数と管理表の
枠数を返します。全ランク・メインスレッドから呼び、一般のMPI資源数やメモリ見積もりとは区別します。

k点並列では実空間群が1ランクなので、共通の密度対称化写像はMPI_COMM_SELFを使い、
最後に初期化した状態の終了が兄弟の通信先を無効化しないようにします。他の格子分割は
渡されたコミュニケータを維持します。所有者のない共通資源は接続終了時に解放します。
初期化中のクラッシュ・確保失敗・MPI異常、無制限の入力設定や並行呼出しの保証ではありません。
後述の明示的レジストリはモデルごとの密度対称化マップも所有・切替します。
従来計算プログラムは従来範囲で動きますが、子はdetachまで残ります。新しい解放方式には
明示的に再ビルドし、固定バイナリを書き換えないでください。

## 異なる入力のHSEモデルを同一プロセスで扱う

`backend.open_registry(context)` が一つの既存MPI群に接続した
`NativeFastHSERegistry` を返します。`registry.open_model(initial, dt)` はモデルごとに
入力を一度検証・コピーし、最初の `new_state` で計算処理側に登録します。その後も元のGSから
独立した状態を作れます。従来の `open_model`／`create_state` の経路は維持します。
接続の初期化から終了までの期間ではレジストリ経路か従来経路を選びます。import時の接続・ビルド・取得はありません。

```python
backend = NativeFastHSEBackend(build_json, work_root, allow_symmetry=True)
with backend.open_registry(context) as registry:
    first = registry.open_model(initial_a, .001)
    second = registry.open_model(initial_b, .002)
    states = []
    try:
        a = first.new_state("material-a")
        states.append(a)
        b = second.new_state("material-b")
        states.append(b)
        ja = a.step([0., 0., 1e-4])
        jb = b.step([0., 0., -7e-5])
        a.checkpoint(rank_checkpoint_path, False)
        restored = first.new_state("restored-a")
        states.append(restored)
        restored.checkpoint(rank_checkpoint_path, True)
    finally:
        for state in reversed(states):
            state.close()
```

セル・k点表現・時間発展／学習設定・対称性表の異なる入力を扱えます。
既存Python providerは原子単位（`unit_system=au`）を要求し、この制約は維持します。単位表の所有分離は単位指定の拡張を意味しません。
入力値・単位表・対称性作業領域と密度対称化マップを各モデルが所有し、step／checkpoint／
transaction／closeの前に対象モデルへ切り替えます。電子状態・XC・交換・FFT・MPIの所有者を
丸ごとコピーする方式ではありません。全ランクで同じモデルを選び、メインスレッドから
同じ順序で呼びます。作業ディレクトリを同期的に切り替えて元へ戻すため、スレッド同時実行には未対応です。

モデルを閉じると新規状態の作成を停止し、既存状態は継続できます。レジストリを閉じても
最後の状態を閉じるまで計算処理接続を保ちます。状態を明示的に閉じてから呼び出し元プログラムが管理するMPIを終了してください。
未初期化モデルの終了で接続全体の予約を解除しません。モデル選択不一致・ディレクトリ変更失敗は
計算処理呼出し前に集合検査します。一般のSTOP・確保失敗・MPI／プロセス異常の復帰保証は増やしていません。

`sp_fast_hse_registry_protocol()==1` を持つ、明示ビルド・監査した計算プログラムが必要です。
C仕様はnative_adapters/fast_hse_handles.hにあります。C呼び出し元プログラムは対象モデルの作業ディレクトリを
指定し、正の状態IDと負のモデルIDを不透明な世代付き識別子として扱います。
同じMPI群・k点並列1／2ランクの範囲です。異なるMPI分割への復元や、通常の物質ID／ランク計画からの
異種HSEレジストリ自動組立ては別作業です。従来の同一モデルの組立ては維持します。

既存サンプルに `--second-initial OTHER_GS --symmetry` を付けると、二つの準備済みモデルを使います。
`--symmetry` は入力の実際の対称性設定を許可する指定で、各保持群がx／y／円偏光を含む実際のAを
引き続き検査します。投影・切捨ては行いません。検証は既存・合成GSによる実装整合性の確認で、
物性の収束保証ではありません。資源予算は従来どおり警告扱いです。

## 通常の異種モデル層組立て

各HSEのMaterialBindingで `NativeFastHSEBackend(..., registry_model=True)` を指定し、
`assemble_native_layers(..., model_groups={"hse-models": ["hse-a", "hse-b"]})` へ渡します。
物質ごとのprovider／初期GSと、マクロ点ごとの状態を保持します。同じHSE群では内部ランク数、
ライブラリ、作業ルート、対称性の許可設定を揃えます。標準TDDFTとHSEは別のFortran/Cで計算するMPI群を使います。
既定のランク計画と `shared_model=True` の経路は維持します。

`salmonpy.examples.mixed_hse_models` は通常の層組立て、FDTD内のLD、SBE、二つのHSE物質IDと
保護付き直接保存／再開を接続する例です。materials JSONは既存の
`salmonpy.native-tddft-materials/1` 形式を使い、この例がHSE providerを明示選択します。
ファイル名からbackendを推測しません。物質ごとに `[2,1,1,1,1]` などのk点並列を指定します。
同梱 `mixed_hse_models_layers.json` にID・並び・厚みを定義し、`--layers` で置き換えられます。
`--separate` は比較用の独立MPI群を作ります。40³場・8³電子格子の小さい3D接続試験で、物性精度の例ではありません。

```sh
mpiexec -n 4 python -m salmonpy.examples.mixed_hse_models FIELD_BUILD HSE_BUILD MATERIALS_JSON NEW_OUTPUT
mpiexec -n 6 python -m salmonpy.examples.mixed_hse_models FIELD_BUILD HSE_BUILD MATERIALS_JSON NEW_OUTPUT --separate --compare GROUPED_OUTPUT
```

全体を `--steps 32` として `--save-at 8` で保存し、新しい出力先に
`--resume SAVED_OUTPUT/checkpoint --compare CONTINUOUS_OUTPUT` で再開します。
物質ごとの内部ランク数と入力は同じに保ちます。モデルごとに一意な出力名を使い、対称性入力用の
ファイル番号が他状態と衝突しない最新計算プログラムを指定してください。

入力雛形は `examples/mixed_hse_models_inputs.json` です。bundleを手元の既存GSへ書き換えて使います。

## 呼出しごとのファイル基準

directory scope2/3対応runtimeでは、HSEの通常・共有・登録モデルの呼出しで
呼び出し元プログラムCWDを変更しません。各workフォルダを一時的なファイル解決の基準とし、
失敗時も解除します。scope1／旧runtimeは従来経路です。GSと同じ共通呼出し処理を
使いますが、同時スレッド実行の対応ではありません。直接再開の数値形式は維持します。

HSEの共有／登録状態は、新しいtimer protocol1で `state.timing_seconds`（total、initialization、density、hartreeの秒）を物質MPI群から照会できます。各ハンドルで計測を保持し、兄弟・呼び出し元プログラムの計算時間とpark中の時間を除きます。完了した区間のみで、totalは0の場合があります。計測値を数値checkpointに含めず、復元で使用済み時間を巻き戻しません。新しい照会だけが対応計算プログラムを必要とし、旧計算プログラムの通常計算は維持します。従来singletonには別のlegacy timer protocol1があります。同時スレッド実行は未対応です。

新規計算プログラムでは、初期化の固定設定を確認する記録もモデルが所有します。登録モデルは各モデルに、同一物質の共有生成器はその生成器に保持し、別モデルを選ぶたびに共通記録を上書きしません。成功した初期化だけが記録を更新し、失敗した生成では維持します。検査項目・C入口・数値保存形式は従来どおりです。従来の単独経路には互換用の記録を残し、全共有変数の分離完了とはしません。

モデルを閉じて登録枠を別物質へ再利用するときは、その枠の初期化記録だけを新規状態へ戻します。生存中の別モデルの記録は維持します。

保存／復元時の配列サイズ検査に使う作業値も、呼出しごとのローカル値になりました。状態から計算した上限を各読み書き処理に明示的に渡します。従来の上限式・不正サイズ検査・保存形式は同じです。これは同時スレッド実行への対応を意味しません。

HSEの単一状態計算処理経路も、新しいlegacy timer protocol1では初期化（対称性あり／なし）・更新・終了の計測を呼び出し元プログラムから分離します。`state.timing_seconds` は共有状態と同じ4項目を返します。古い計算プログラムでも通常計算は維持し、新しい照会だけ対応版が必要です。計測値は数値checkpointに含めず、MPI・更新順・初期化から終了までの期間の制約は変えません。

HSE保存処理の学習履歴g・b・coeffも、スカラーIEEE検査で全要素を調べます。従来の履歴条件・拒否判定・保存形式を維持します。今回の有限値／NaN／正負Inf試験は合成履歴で行い、実電子状態のディスク再開やFujitsuコンパイルを新たに検証したものではありません。

`model_settings`は必須です。単一状態・別プロセス計算も含め、呼び出し元が設定を所有し、暗黙の共有設定へ落とす処理はありません。プロセス内の選択変数は残り、スレッド並行実行への対応を意味しません。

`model_settings`は必須です。単一状態・別プロセス計算も含め、呼び出し元が設定を所有し、暗黙の共有設定へ落とす処理はありません。プロセス内の選択変数は残り、スレッド並行実行への対応を意味しません。

## HSEの性能モード

`NativeFastHSEBackend(..., execution_mode="performance")` は呼び出し元が
MPIを管理する更新を切り替えます。既定は `"checked"` です。
選択は状態ごとに固定されます。単一状態・共有モデル・異なるモデルの登録に
対応します。別プロセス型には、下記の `FastHSEResponseBackend` の同名引数を使います。

新旧Aの合意、Pythonの有限値走査、電子状態・交換状態の自動退避、数値異常
復帰用の診断合意・走査を省きます。初期化・ライブラリ照合・ローカルな有効期間／
容量／形状・明示的な保存／復元・観測の検査は残ります。ソルバーの破綻判定、
FFT・Hartree・交換の計算に必要な通信、学習準備状態の通信も維持します。
縮約に残した対称性とAの整合条件は計算式の前提なので維持し、status34なら
更新前に拒否して有効なAで再試行できます。それ以外の更新失敗はPythonの
状態を失敗扱いとし、自動復元を保証しません。Fortranの異常でジョブが終了する
場合もあります。

学習用データを取り込む処理は、両モードでエネルギー・出力処理の前に保ちます。
状態が所有する一時フラグで処理順と異常復帰診断を分け、更新後に戻します。
プロセス共有の実行モード変数は増やしません。Taylor／PT-CN、全／縮約k点、
内部1／2ランクの合成入力で、電流・エネルギー・学習履歴とメモリ内復元後の
一致を確認しています。物理的収束や速度の評価ではありません。

`mpi_mixed_native --hse-execution-mode performance` は
`--fast-hse-mpi-build` を必要とし、HSE更新だけを
切り替えます。場と外側の結合処理は、それぞれの設定の検査を維持します。

## 共有・登録HSEの状態ごとのモード

バックエンドの `execution_mode` がモデルの既定値になります。
`model.new_state(point_id, execution_mode="checked")` または `"performance"`
で状態ごとに上書きし、省略時はモデルの既定値を使います。
`registry.open_model(..., execution_mode=...)` で登録モデルごとの既定値も
指定できます。セッションが登録管理器を再利用しても各物質設定の選択を維持します。
通常・性能の状態を混在でき、電子状態・交換資源・学習履歴・A／電流・出力履歴を
それぞれ保持します。プロセス共有のモード変数は追加しません。

性能更新でもハンドルの世代・有効期間の境界と、縮約対称性に関するAの数学的な
条件は維持します。診断用のトークン合意と自動復帰用コピーを省きます。
モデル入力・作業領域の切替、物理通信、学習処理順は維持します。通常モードの
状態は従来の検査・異常復帰経路を使います。明示的なメモリ／ディスク保存と
復元の検査・負荷は残り、復元先の状態のモードは変えません。性能更新の異常時に
自動復元を保証するものではありません。

`shared_fast_hse --execution-mode performance` はpoint-aを指定モード、
point-bを通常モードにし、point-aの再開も例示します。`--second-initial` を
指定すると、異なる入力モデルの登録を使います。通常の `mpi_mixed_native` の
`--hse-execution-mode performance` も共有HSE点に対応し、場・外側の結合検査は
それぞれの設定を維持します。

## 別プロセス型HSEの使い分け

`FastHSEResponseBackend(..., execution_mode="performance")` で子プロセスの性能モードを選びます。既定は `"checked"`、選択は状態ごとに固定です。Python本体のMPI・作業ディレクトリ・スレッド設定を計算から分離できるので、この経路は維持します。1物質点に1つの直列CPUワーカーを使い、内部MPIを使う大規模計算や縮約k点には `NativeFastHSEBackend` を選びます。別プロセス型は従来から全k点の範囲で、縮約入力は初期化できません。性能モードによってこの物理対応範囲は変えません。

新しいワーカーは、初期化後・更新前に `MODE_PERFORMANCE` を受理し、`SALMONPY_HSE_PERFORMANCE 1` を返す必要があります。未対応の実行ファイルには黙って通常モードへ下げず、状態を閉じて未対応を知らせます。ビルドは明示的に `tools/build_fast_hse_worker.py` で行い、importではビルドしません。

更新時のA・電流・エネルギーの有限値走査と旧Aの照合、Fortran側の追加診断・自動退避を省きます。初期化、実行ファイルのSHA、通信形式・ステップ番号、形状・容量・有効期間、明示的な保存／復元の検査は維持します。数値演算と学習／出力の処理順は変更せず、既存の `rt_step(..., checked=...)` を使います。ワーカーの異常後に自動再起動・履歴再計算はしません。他の独立状態と呼び出し元を終了する処理は追加しませんが、その失敗した状態の継続は保証しません。

直接チェックポイントを復元する場合、復元先が選んだモードを維持します。保存形式はバックエンド側に委ね、旧Aの復元やファイルの認証は既存のチェックポイント契約に従います。`mpi_mixed_native --hse-execution-mode performance` は、MPI版と別プロセス版のどちらのHSEバックエンドにも選択を渡します。外側の結合モードは別指定です。
