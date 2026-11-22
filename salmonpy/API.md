# 公開Python APIと互換性

## C／Fortranから呼ぶ前の共通照会

`sp_host_abi(capacity, values)`は、計算プログラムの版と基本型を調べる入口です。
`native_context.h`を使い、`int64_t`を10個以上用意します。成功時は0を返し、
先頭10個だけを書きます。容量不足は13を返し、配列を変更しません。
MPIの初期化・接続、ソルバーの初期化、配列の確保や通信は行いません。

| C配列の添字 | 内容 | 現行値 |
|---|---|---|
|0|共通呼び出し規約の版|1|
|1～5|C int、C double、int64、Fortran INTEGER、real(kind=8)のビット数|32、64、64、32、64|
|6|コミュニケータの渡し方|1：MPI_Comm_c2fの結果をint64へ拡張。MPI_Commそのものは渡さない|
|7|MPIの初期化・終了を行う側|1：呼び出し元|
|8|状態を操作するスレッドの条件|1：メインスレッドから順に呼ぶ|
|9|I/O操作に必要なディレクトリ処理の版|3。I/O機能を含むことの証明ではない|

Fortranからの配列添字は1始まりです。基本型と版が一致しても、個別操作・MPI実装・
物理モデル・保存形式が対応する保証にはなりません。それぞれの照会が必要です。
組み込みライブラリのPython接続は、SHA256照合後・MPI接続前にこの内容を検査し、
一部ランクだけ不一致でも担当群全員で拒否します。旧版への切替や自動ビルドはしません。
自作バックエンドの独自契約を、この照会へ強制的に合わせるものではありません。

明示的に読み込んだ信頼できるライブラリを調べるPython関数は
`salmonpy.api.query_host_abi(library)`です。不一致は`RuntimeFailure`、
組み込みMPI接続での群内の不一致は`InputError`として通知します。
計算せずに別プロセスで照会する例：

```sh
python -m salmonpy.doctor --native-build BUILD_JSON --component fdtd --route single
```

`single`は一つの状態、`shared`は同じモデルの複数状態、`registry`は異なるモデルの
登録を表します。旧`legacy`名の別名受入れはありません。


[English](API.en.md)

利用できる関数・クラスと、エラーや保存データの扱いを説明します。目的から使い方を選ぶ場合は、先に[はじめに](README.md)または[サンプル集](examples/README.md)を参照してください。

導入した版は `salmonpy.__version__` で確認できます。公開前の開発版のため、APIは変更されることがあります。計算プログラムと同梱マニュアルを組にして使ってください。C/Fortranから呼ぶ場合は[C_ABI.md](C_ABI.md)を参照します。

## 公開範囲

`salmonpy.api.PUBLIC_API` がモジュール別の公開名一覧です。指定モジュールからimportします。
偶然importされている名前や、先頭がアンダースコアの実装名は公開APIではありません。
戻り値の状態・結果オブジェクトは各ガイドに記載したメソッド・プロパティを公開します。
内部の場状態クラスを直接構築しません。`StateSnapshot` は型を示す名前で、利用者が構築するものではありません。

| モジュール | 主な用途 |
|---|---|
| `ld`, `sbe`, `crystal`, `materials` | 物質・初期状態・A駆動応答 |
| `maxwell`, `timing` | 場の生成・結合・時刻宣言 |
| `fields`, `geometry`, `analysis`, `outputs` | パルス・構造・解析・出力 |
| `dft`, `bands`, `tddft`, `fdtd`, `runtime` | 外部ソルバーの明示的実行・準備 |
| `ase`, `ase_calculator`, `salmon_vr` | 任意の連携機能 |
| `distributed`, `material_response`, `embedded_tddft`, `coupled` | Pythonと同じプロセスで計算し、MPIで分担するAPI（実験的） |
| `errors`, `snapshots`, `api` | 例外・変更できない保存状態オブジェクト・公開一覧 |

トップレベルの `__all__` は版、`LorentzDrude`、例外型に限定します。
任意依存を一括importしません。import時のコンパイル・取得・Fortran/C共有ライブラリの読み込み・MPI初期化はありません。
snapshot／restoreは標準出力へ書きません。必須Python依存は引き続きNumPyです。
LD/SBE計算ルーチンや外部計算プログラムの明示的準備、任意依存の要件は従来どおりです。

## 現行APIを使う方針

SALMONpyは未公開のため、後方互換分岐・旧引数の別名・旧形式の変換・廃止予告期間は設けません。現行の文書・API・計算プログラム・保存側を組み合わせて使ってください。単位・電流符号・時刻配置は明示します。再現のため、検証したパッケージ版を固定してください。

保存・時刻・計算プログラム・snapshotの形式はそれぞれ版を持ちます。パッケージ版だけから再開可否を推測せず、形式とモデル・実行条件を照合してください。

結晶は **A/c → 電荷電流密度j** を受け渡します。内部の電子物質電流は J = -j です。基本はtransverseで、Longitudinalは明示的に選択します。メモリ上のsnapshotは同じパッケージ版・同じ元モデルを要求します。

## 例外の仕様

| 型（`salmonpy.errors`） | 意味 | 従来の捕捉型 |
|---|---|---|
| `SalmonpyError` | ライブラリが明示的に送出する例外の親 | `Exception` |
| `InputError` | 値・互換性・状態遷移の不正 | `ValueError` |
| `InputTypeError` | 型の不正。`InputError`の子でもある | `TypeError` |
| `SnapshotError` | 不正・不適合な復元。状態は保持 | `ValueError` |
| `NumericalError` | 検出した非有限値などの数値異常 | `FloatingPointError` |
| `RuntimeFailure` | 計算プログラムの欠如・失敗・結果の未確認 | `RuntimeError` |

既存の `RuntimeUnavailable` と `OutputSaveError` は名前を維持し、`RuntimeFailure`を継承します。
従来の組込み例外型での捕捉も利用できます。OS・ファイルの例外（`FileExistsError`を含む）、
プロセスの時間切れ、NumPy変換、任意依存の不足、利用者コールバックの例外は元の型を保つ場合があります。
外部例外を一律に包みません。収束報告がないことは計算プログラムの失敗であり、数値発散の証明ではありません。

例外の型だけから計算全体の巻き戻しを推測しません。内蔵LD・SBE・Pythonの場のstep／restoreは
失敗時に最後に成功した計算状態を保持します。MPI応答・結合の失敗はセッションを無効化し得ます。
FortranのSTOP／abortを回復可能なPython例外にする変更ではありません。

プロセス分離`TDDFTState`は`save_checkpoint(directory)`と
`restore_checkpoint(directory, max_bytes=None)`も公開します。
同じ内容の初期モデルを明示的に用意し、Zarr履歴から別プロセスで再開できます。
初期GS本体は保存しません。詳細は[状態保存](SNAPSHOTS.md)を参照してください。

直列モデルの共通仕様、適合チェッカー、MPI点分散への接続は[MATERIALS.md](MATERIALS.md)を参照してください。

`fast_hse` の `FastHSEMaterialProvider`／`FastHSEResponseBackend` は明示ビルドした高速PR版serialワーカーを外側MPIの物質応答へ接続するプレビュー入口です。COUPLED.mdを参照してください。HSE内部MPIを初期化しません。

`fast_hse_mpi.NativeFastHSEBackend` は高速PR版HSEの実験的な呼び出し元プログラムが管理するMPI入口です。k点のみ1／2ランクの物質群と、同じ分割での直接保存／復元に対応します。旧create_stateはメインスレッド・1プロセスにつき初期化1回です。open_modelはNativeFastHSEModelを返し、同一の全k点モデルの状態を共有接続で独立に扱います。SHARED_HSE.mdとCOUPLED.mdを参照してください。

`checkpointing.CheckpointCapability` と `check_checkpoint_compatibility` は、自作の計算コードにも使える保存／復元前の共通宣言・検査です。[CHECKPOINTING.md](CHECKPOINTING.md)を参照してください。

`checkpointing.CheckpointTarget` と `execute_checkpoint` は検査後に計算処理独自の保存／復元関数を呼び、不透明な戻り値をプロセス内で返します。

`CheckpointState`／`restore_checkpoint` は、独立候補を全ランクで準備した後にまとめて公開する失敗時に元の状態を保つ復元です。`AtomicRestoreError`で元のホルダー維持と候補の破棄エラーを区別します。

メモリ助言の公開入口はsalmonpy.memoryのMemoryEstimate／MemoryReport／memory_advice／process_memoryです。予算超過は警告のみで、実行可否や上限の保証ではありません。詳細は[MEMORY.md](MEMORY.md)。

`stepping.StepTarget`・`guarded_step`・`StepRecoveryError` は明示した完全退避仕様によるステップ保護のプレビュー入口です。STEPPING日英ガイドを参照してください。

`native_steps.native_step_targets` は、監査済みFortran/Cで保持する計算状態と通常出力を共通の失敗時に元の状態を保つ時間更新へつなぐプレビューAPIです。NATIVE_STEPPING日英ガイドを参照してください。

### ソルバー初期化時の原子入力

標準ソルバーと高速PR版HSEは、原子数・元素数・座標（Cartesian／分数）・元素IDを、内部の
値オブジェクト `s_ion_initial_config` に取り込みます。`initialize_ions` は選択した状態の
格子から座標を解決し、独立した所有配列を作ります。形状不一致・非有限座標・範囲外元素ID・
下限が1でない入力配列は、対象状態を変えずに拒否します。明示的な処理は共有入力を読み書きしません。
既存互換経路では下流処理のための分数座標→Cartesian共有入力正規化を維持し、正常な計算結果を保ちます。

これは内部の所有境界であり、新しいPython APIや安定したC ABIではありません。セル・k点・
対称性・コミュニケータなどの依存は残り、異なる全入力のGS／TDDFTを同一プロセスへ自由に
配置できることは意味しません。メモリ確保失敗や初期化中のSTOPからの復帰も未対応です。

### セル・実空間格子の明示的な解決

標準ソルバーにも高速PR版HSEと同じ値型・処理 `s_cell_config`、`s_cell_geometry`、
`s_rgrid_config`、`resolve_cell`、`init_grid_whole_config` を揃えました。両方とも
入力の取り込みだけを行う `cell_config_from_input` を使い、明示的な解決処理は共有入力を
読みません。正常入力のセル指定の優先順位・演算順・周期／孤立系の添字・epsilonを含む丸め・
Maxwell用整合長を維持します。拒否時は対象の形状・格子添字配列を変えません。
これはセル値と全体格子添字の解決で、全座標・並列分割・対称性・FFT資源の初期化ではありません。

既存互換経路には共有入力の取り込みと、下流処理用の `al` 正規化が残ります。異なるセルの
GS／TDDFTを同一プロセスで自由に初期化できる段階ではなく、新しいPython／C公開入口もありません。
実験用RT設定版にも同じ値型を持たせ、準備時に更新した置換前後のハッシュを検査します。

### 初期軌道構成の明示的な選択

標準版・高速PR版の初期化で、計算モード・周期性・スピン／SOC・状態数／電子数・温度・
k点判定を内部の `s_orbital_initial_config` に取り込みます。`resolve_orbital_layout` は
共有入力を読まずにスピン数・軌道数・実数／複素数格納を選びます。セルの直交性を明示的に渡し、
高速PR版は所有済みの占有設定・実際のk格子と、従来と同じ取り込み時のhybrid判定を使います。
標準版のファイル指定k点も、この判定には従来の `-1` 格子印を維持します。

RTの `temperature < -1e-12` 分岐など、正常入力の既存GS／RT選択規則は維持します。
不明なモード・非有限温度・使えない構成は、結果オブジェクトを変えずに拒否します。
これは初期配列構成の選択であり、占有容量・物理的妥当性・SCF入力全体の所有や、計算状態の初期化からの
異常復帰を保証するものではありません。新しいPython／C公開入口もなく、k点・対称性・
コミュニケータなどの共有依存は引き続き残ります。

HSEは高速PR版のみです。`NativeFastHSEBackend(build_record, work_root, allow_symmetry=True)` で
対称性を明示的に許可できます。内部kランク1／2・周期・非スピン分極・CPU・full-qのTaylor／PT-CNに対応し、共有の対称化設定は不変、1プロセスにつき初期化1回です。既定値Falseと既存C入口は
対称性を拒否します。内部2ランクはsp_fast_hse_symmetry_kmpi_protocol()==1が必要。Taylor strict／learned・PT-CNにはsp_fast_hse_symmetry_history_protocol()==1が必要です。学習PT-CNは既存方式どおりinterval=1・predictor_corrector=nを指定します。

新C入口 `sp_fast_hse_init_symmetry` は既存initと同じ引数で、
`sp_fast_hse_symmetry_protocol()==1` を宣言します。全ランクの有限なstatus34とこの宣言を
確認した場合だけ、Pythonは状態を変えない `InputError` を返し、Aを直した再試行を許します。
未知・旧ABIの失敗は継続しません。FDTD結合では全体の再試行に `step_protected` が必要です。
直接／候補復元は共通checkpoint経路を使います。

[対称性ガイド](FAST_HSE_SYMMETRY.md)には、混在サンプル、対応範囲と極小勾配XCの安定評価と回帰検証を
記載しています。接続・保存復元のテストを物理収束の証明と解釈しないでください。

縮約k点の内部2ランクと同一分割での直接再開に対応します。追加の計算処理 k-MPI capabilityが必要です。分割変更による再開は未対応。詳しくは対称性ガイドを参照してください。

NativeFastHSEBackendは計算状態の初期化前に周期セル・格子・原子入力を読み取り検査し、群内でInputErrorを合意します。入力を直してproviderを作り直せます。再開バイナリの互換性や全計算状態の初期化の安全性は別仕様です。詳しくは対称性ガイドを参照してください。

SCFの閾値未達は `RuntimeWarning` として通知し、有限で読み出せる結果を後続へ渡します。
`GroundStateResult.converged` と保存JSONの `converged` が収束状態を表します。
保存JSONの `status='complete'` は実行完了を意味し、収束済みを意味しません。

異なる入力を同一プロセスで順次切り替える標準TDDFTは [TDDFT_REGISTRY.md](TDDFT_REGISTRY.md) を参照してください。

[保持GSハンドルとPython接続](NATIVE_GS.md)。

`salmonpy.gs_input.GSInputMetadata` / `describe_gs_input` と保持状態の
`state.input_metadata()` で不変の宣言入力を参照できます。任意のmetadata extraは
明示呼び出し時だけ読み込み、source_idは物理的GS同一性や復元互換性を証明しません。
詳細は保持GSガイドを参照してください。

`salmonpy.preparation.preflight_inputs(directory, manifest_sha256=None, registry=None)`は、
計算を開始せずJSON化可能な入力受付レポートを返します。信頼する検証関数の登録は任意で、
省略時はファイル整合性だけが対象です。ジョブの終了コードとランク別記録はPREPARATION.mdを参照してください。

## 計算処理実行経路の仕様情報

```python
from salmonpy.api import native_contract
contract = native_contract('fast_hse', route='shared')
assert contract['model_settings_owner'] == 'explicit model'
assert contract['runtime_verified'] is False
```

`native_contract`は組み込み`tddft`／`fast_hse`の`single`・`shared`・`registry`経路、
および`fdtd`の`single`経路を説明します。ライブラリ読込・計算・新しい制限は行いません。
単一状態の明示設定と明示的モデル所有、MPIの呼出し元所有、メインスレッドの逐次操作、collectiveの順序、電流符号、失敗保証範囲を区別します。
標準TDDFTのlegacy C入口は物質電流J=-jを返し、Python側は電荷電流jへ変換します。
高速HSEと標準handleのC入口は電荷電流jです。照会がローカルかcollectiveかは個別ヘッダーに従い、診断を一律collective扱いしません。

これは宣言された仕様の情報で、実計算プログラムの能力確認や安定C ABIの宣言ではありません。`runtime_verified`は常にfalseです。
実際の計算プログラム・保存再開の対応可否は別途照会します。未知の自作コンポーネントや経路名には`known=False`・`status='unverified'`を返し、
自作の計算コードを拒否しません。その仕様は計算処理自身が定めます。結果は独立した辞書で、呼出し側で注記しても後の照会に混入しません。
既存の計算呼出し、数値状態、legacy経路の省略時処理、計算可能な条件は変えません。C／Fortran呼び出し元プログラムの安定ABI適合検証は別の残課題です。

### 明示指定した実計算プログラムの版照会

```python
from salmonpy.doctor import probe_native_protocols
report = probe_native_protocols('audited-build.json', component='tddft', route='shared')
# errors・照会結果・protocol_matchを確認。runtime_readyは未確認のままです。
```

CLIは`python -m salmonpy.doctor --native-build audited-build.json --component tddft --route shared --json`です。
`native_contract`や独立静的チェッカーとは異なり、明示指定された計算処理ライブラリを、記録のSHA256と照合して別プロセスで読み込みます。
子側でも読込前にhashを再確認します。hashは同一性であり、信頼性やビルド監査の証明ではありません。
読込時に初期化コードが動く可能性があるため、信頼する監査済みライブラリだけを対象にしてください。
既知の引数なし版照会関数だけを呼び、ソルバー状態の作成、MPI初期化・接続、ビルド、インストール、計算は要求しません。
未知の自作経路は記録ファイルも開かず、コードを読み込みません。

`protocol_match`は選択した版照会だけの一致です。関数の欠如や新しい版は記録するだけで、自動的な計算処理拒否や別機能への置換はしません。
`runtime_ready`は未確認のままで、操作関数のシグネチャ、物理的対応範囲、保存形式・MPI分割の互換性の保証ではありません。
子プロセスの異常は呼出し元を落とさず記録し、`library_loaded=None`は子側が読込完了の有無を報告できなかった意味です。
渡されたコードを読み込まない`check_job`とは別の、任意の明示的診断として扱います。

CとISO_C_BINDING Fortranの版照会例は`tests/native/test_native_protocol.c`と、
`tests/native/test_io_ownership.f90`の`SALMONPY_ABI_PROTOCOL_PROBE`モードです。
開発用の実行入口は`tools/test_io_ownership.py BUILD_RECORD NEW_OUTPUT --abi-contract`で、
高速HSEの場合は`--fast-hse`を追加します。計算プログラムhashと全リンク対象のアプリケーションオブジェクトを監査し、
返却値の一致とMPIが未初期化のままであることを確認します。これは版照会の相互呼出し検証で、モデル・状態・step全体のABI適合確認ではありません。
C/Fortranの呼び出し規約は引き続き実験段階です。

### 状態の作成・更新・終了を確認する呼び出し元プログラム試験

`tools/test_io_ownership.py BUILD_RECORD NEW_BUILD_DIR --abi-state --fixture EXISTING_PZ_FIXTURE --state-output NEW_RESULTS_DIR --ranks 2` は、監査済みの通常TDDFTライブラリをCとISO_C_BINDING Fortranから呼びます。高速PR版HSEのビルド記録と既存の全k点HSE入力には`--fast-hse`を付けて同じ検査を実行できます。HSE試験では各実行のk点並列数だけを設定し、渡されたGSは変更しません。縮約k点入力にはさらに`--state-symmetry`を付けます。渡された`sym.dat`を参照し、対称性対応のHSE初期化を使い、この入力の保持対称操作に合うz方向を駆動します。代表k点が実際に出力されたことも確認します。この試験入力の方向指定であり、全k点でのx・y・z駆動を制限しません。MPIの初期化・終了は呼び出し側が所有します。2状態と個別実行の一致、直前Aが不一致のときの状態保持、スロット再利用後の古いハンドルの拒否、新しい状態の独立した初期化を確認します。既存GSはその場で読み、変更されていないことを確認します。GSのコピーや電子状態チェックポイントの保存はしません。結果の保存先は明示してください。

ディレクトリスコープの共通宣言は、両ハンドルヘッダーから取り込む`native_context.h`にあります。macOSではMPI呼び出し元プログラムのリンクにも既存の`tools/build_embedded_tddft.py::mpi_link_command`を使ってください。OpenMPIのFortran COMMON特殊値が実行ファイル内に別配置されると、共有ライブラリとの間でMPI_IN_PLACE処理の判定が変わる場合があります。試験ランナーはこの共通処理を使い、リンク引数と実行ファイルのハッシュを記録します。版の問い合わせだけではこの問題は検出できません。

確認範囲は通常PZと高速PR版HSEの全k点／縮約k点入力、C／Fortran呼び出し元プログラム、MPI 1／2ランク、複数状態と個別実行の電流一致です。HSE固有の引数と拒否コード（直前Aの不一致はHSEで1、通常TDDFTで10）を共通試験内で扱い、公開ABIや数値ソルバーは変更しません。縮約入力は対称性対応PT-CNと学習型交換を使い、全8点を4代表点で扱います。他XC・形状・他環境の検証は別途必要で、公開C/Fortranの呼び出し規約は引き続き実験的です。この小規模入力は接続と状態所有を検査するもので、実材料の収束を保証しません。

### 呼び出し元プログラム試験のPZ入力を準備する

`--abi-state`の標準PZ入力は`rt.inp`、`data_for_restart/`、`pseudo1_rps.dat`を必要とします。`theory='tddft_pulse'`、周期・固定イオン・非スピン分極、`middlepoint`、`trans_longi='tr'`を使い、試験用RT入力から`&parallel`を除いてください。同じランナーが1／2ランクを試すため、固定した並列条件を埋め込まないひな形を提供します。GSのセル・原子・格子・k点・電子数・擬ポテンシャル・XCは入力と一致させてください。

```sh
python tools/prepare_pz_state_fixture.py GS_INPUT EXISTING_RESTART PSEUDO NEW_FIXTURE
python tools/test_io_ownership.py BUILD_JSON NEW_BUILD --abi-state --fixture NEW_FIXTURE --state-output NEW_RESULTS --ranks 2
```

SALMONpy導入済みのPython環境で実行します。`GS_INPUT`には元のGS入力または同じGSから検証したRT入力を指定します。この準備ツールはPZ・周期・単一元素・インライン原子座標の試験用であり、GSを計算・コピー・ダウンロードしません。再開ディレクトリと擬ポテンシャルは既存ファイルへの参照です。小さなGS入力の例は`examples/pz_host_gs.inp`（Si1原子、8³格子、2×1×1 k点）です。必要なら利用者が計算ノードで明示的にGSを準備し、次にこのツールへ渡します。収束した物性の例ではありません。

`BUILD_JSON`にはビルド記録のファイル、または`build.json`を含むディレクトリを指定できます。`--parser-host`はGSを使わず、`status`単独／`native_host=.false.`で通常検査を維持し、明示trueだけが場用検査へ切り替えることをMPIで検査します。

計算プログラムバンドルからは`SalmonExecutable.from_runtime(directory)`を使えます。これは`runtime.load_runtime`へ委譲し、全ファイルと能力記録を検証します。バイナリ単体とハッシュを手作業で組み合わせる必要はありません。

## 誘電関数・光学伝導率（単一偏光）

`salmonpy.analysis.optical_conductivity` は同じ時刻の全電場Eと電荷電流密度jからsigma=j/Eを求めます。結晶A駆動で得た電流にも使えますが、AそのものはEではありません。呼び出し側で時間配置と微分規約を揃えたEを用意してください。SALMON物質電流Jを使う場合はj=-Jへ明示的に変換します。単位はGaussian原子単位のみです。完全なテンソルには独立な駆動方向の計算が必要で、円偏光一回から全テンソルを復元しません。

```python
from salmonpy.analysis import optical_conductivity, dielectric_function
sigma = optical_conductivity(time, charge_current_x, total_E_x,
    current_convention="charge", field_kind="total", units="atomic_gaussian",
    window="hann", drive_rtol=1e-8,
    metadata={"material_id": "Si", "polarization": "x"})
epsilon = dielectric_function(sigma)
# sigma.omega, sigma.values, sigma.valid; epsilonも同じ構造
```

フーリエ変換は既存の `spectrum` と同じexp(+i omega t)、dt倍の和です。時間依存exp(-i omega t)に対応し、epsilon=1+4pi*i*sigma/omegaです。弱い駆動周波数（最大振幅のdrive_rtol以下）と数値変換不能点はvalid=False、値はNaNとなります。DCの誘電関数は未定義です。全駆動ゼロでも解析結果をマスクするだけで、時間発展を停止しません。外場だけの比をバルク誘電関数と呼び替えません。電場と電流には同じ窓を使い、窓・時間長・開始時刻・dt・電場振幅・利用者の来歴をconfigに保持します。窓付き有限時間比は近似で、特に単色駆動以外では窓依存性を確認してください。

`optical_linearity(reference, candidate, rtol=...)` は利用者が渡した異なる振幅の伝導率を共通の有効周波数で比較し、相対RMSを返します。材料・偏光・駆動形状の一致は利用者が確認してください。電流が全てゼロなら相対誤差を定義せず判定も未定義とします。`optical_sum_rule(sigma, expected_weight=..., normalization="...")` は有効な隣接周波数間のRe sigmaを台形積分します。正周波数の全電子系では期待値pi*n/2が基準ですが、擬ポテンシャル・非局所電流・不足バンドに応じた基準は利用者が指定します。欠損帯域を補外せず、物理精度の合格証明にはしません。追加のSCF/TDDFT計算やダウンロードは起動しません。

合成データのサンプル: `python -m salmonpy.examples.optical_response`。実材料の収束済み光学定数の例ではありません。


追加入口: salmonpy.preparation.prepare_fdtd(calculation, directory, mpi=None)は新規のPythonから直接呼ぶFDTD入力群を保存し、Pathを返すだけです。FDTDMPI.rank_argumentsとcommand(executable)はargvの構成・プレビュー用です。実行は自動で始まりません。詳細はPREPARATION.mdを参照してください。

入力準備専用: prepare_dft(model, directory, scf=None, save_restart=True)とprepare_tddft(input_file, directory, files=None)はPathだけを返し、SCF／TDDFT worker／MPI／ジョブ投入を起動しません。PREPARATION.mdを参照してください。


## SCF未収束時の確認待ち（2026-11-21）

`DFT.run()` と `NativeGSSession.run()` は、反復上限で未収束なら警告し、
`SCFReviewRequired` を送出して後続計算への自動継続を止めます。
失敗した実行としてデータを消すことはしません。バッチの保存済み結果／計算処理の生存状態は
例外の `result` に保持し、`report_directory` に残差CSV・任意のPNG・report.jsonを保存します。
MatplotlibがなくてもCSVは利用できます。図の生成失敗と計算失敗を混同しません。
レポートの保存に失敗しても確認待ちは維持し、元の出力を確認するよう警告します。
通常のGS保存設定を維持し、`save_restart=False` なら新たなGS再開データは保存しません。

```python
from salmonpy.errors import SCFReviewRequired
try:
    gs = model.run(engine, 'new-gs-output')
except SCFReviewRequired as review:
    gs = review.result
    print('Review SCF residuals:', review.report_directory)
    # ここで終了し、残差の推移・閾値・入力条件を利用者が確認する。
    # このexcept節からTDDFTを自動で開始しない。
```

残差が下がっているか、停滞／振動しているかを確認し、反復上限・混合条件などを見直します。
確認後に未収束状態を使うと利用者が明示した場合だけ `allow_unconverged=True` を指定できます。
AIが自動的にこの指定で再実行してはいけません。既存結果の `GroundStateResult.open()` や
残差の読取・保存・warm startは引き続き可能です。既知の未収束GSからの
`band_initial_state()`／`tddft_initial_state()` にも同じ確認待ちを適用します。

計算処理の `run(..., steps=K)` と `continue_scf(steps=M)` で意図的に反復境界へ
止めた `paused=True` の状態は正常な途中状態です。全反復上限に達した未収束終了だけが
確認待ちになります。確認待ち後もFortran/Cで保持する計算状態は閉じず、保存や診断に利用できます。
ランク間で `allow_unconverged` の指定を揃えます。計算の更新順・収束閾値は変更しません。
この方針は以前の「未収束でも警告のみで自動継続」を置き換えます。


## FDTD／マルチスケールのエネルギー収支

`analysis.absorbed_energy_density(time, E, j, units='atomic_gaussian',
current_convention='charge', field_kind='total')` は各マクロ点の `∫j·E dt` を
Fortranで台形積分します。E/jは同時刻・同位置の `(時刻数, マクロ点数, 3)`、
返り値は点ごとのHartree/bohr³です。SALMONの物質電流Jは明示的にj=−Jへ変換します。
この値は可逆な物質の蓄積を含む正味の仕事で、散逸熱だけではありません。
TDDFTの電子エネルギー変化やLDの分極蓄積をもう一度加算しません。

`analysis.reflection_transmission` は分離済みの入射・反射・透過E/H
（それぞれ `(時刻数, 3)`）からR(ω)/T(ω)を返します。法線は入射側から透過側へ。
複素Poynting流束の比なので偏光を1方向に固定せず、Eだけの振幅比に置き換えません。
area_weights省略時は同じ面積の一様な平面ポートが対象です。全電場をそのまま反射波にせず、参照計算などで
**場を分離してから**渡します。Yee格子のE/Hは時刻・空間位置を揃えてください。
弱い入射スペクトルの点は `valid=False`・NaNとし、クリップや偽の0を返しません。

```python
from salmonpy.analysis import absorbed_energy_density, reflection_transmission
from salmonpy.fdtd import energy_balance

work = absorbed_energy_density(time, E_points, charge_j_points,
    units='atomic_gaussian', current_convention='charge', field_kind='total')
ports = reflection_transmission(time, E_inc, H_inc, E_ref, H_ref, E_tr, H_tr,
    normal=(0, 0, 1), units='atomic_gaussian', window='none')
check = energy_balance(work, point_volumes, ports,
    material_ids=material_ids, transverse_area=area,
    field_energy_change=delta_field_energy, other_outgoing_energy=other_outflow,
    units='atomic_gaussian', rtol=.01)
print(check.residual, check.relative_residual, check.within_tolerance)
# ports.omega, ports.R, ports.T, ports.valid; workは各マクロ点の順序を保持。
```

収支は `入射−反射−透過−Σ(work×担当体積)−場のエネルギー変化−その他の流出`。
体積はbohr³、面積はbohr²、総エネルギーはHartreeです。物質IDは計算由来のものを維持し、
同じIDの複数点も個別に扱います。必須の場のエネルギー変化には瞬時背景ε∞の蓄積を含め、
物質への仕事を二重計上しません。ポートが拾わない側面流出・PML吸収等は別途計上します。
PMLの外へ置いたポートとPML吸収を同時に加えないなど、領域境界を揃えます。

R/Tは選んだ窓のスペクトルですが、総収支には窓を掛けないE/Hから積分した流束を使います。
R/Tの単純平均や、駆動帯域の `1−R−T` だけで全吸収を推定しません。非線形な周波数変換も
時間領域流束には含まれます。有限時間の残留場・遅い放出・観測不足は診断に現れます。
既定1%は検査の目安で物理精度の認定ではありません。超過は警告し、計算を止めたり
ソルバー状態を変更したりしません。ゼロ入射・非有限入力など計算不能な診断入力は拒否します。
再開データでは重複する境界サンプルを除いて時系列を連結し、同じ始終状態で比較します。
追加のMPI通信・自動ダウンロード・計算開始はなく、担当ランクで解析できます。
非一様な周期面は下記area_weights入口で扱います。閉曲面の自動配置や回折次数別の分解は行いません。

保存済みNPZを使う例は `python -m salmonpy.examples.fdtd_energy INPUT.npz NEW_OUTPUT`。
入力キーと単位はそのサンプルの説明を参照します。ソルバーは起動しません。


### 実際の結合計算での確認（2026-11-21）

Apple Silicon・監査済み計算処理実行ファイルで、実際のLD＋E/H FDTDと
2点TDDFT＋Weyl FDTDを実行しました。入射比の収支残差の絶対値は、LDで
0.0108%、既に進行していた半分の格子幅・dtで0.00270%、TDDFTで0.000711%でした。
これは結合・観測処理の確認で、実材料の光学物性の収束認定ではありません。
利用者の判断に従い、これ以上の精度探索は行っていません。

LDは横方向周期の一様な膜で、実測Eから同じ計算処理 LD計算ルーチンでjを再構成しました。
再構成の最終電流は計算処理のLD保存データとrtol=1e-10・atol=1e-18で一致しています。
生のYee Hを空間方向に位置合わせし、時間中心化済みHを再度時間平均してはいません。
TDDFTは既存GSをその場で読み、2マクロ点・8 MPIランク・dt=0.003 auで実行し、
E/H/jと場エネルギーをSALMONの計算出力から直接取得しました。新しいGSや波動関数は保存していません。

入射を遅らせて初期場を空にした主検証に加え、初期充填波を持つ条件でも
初期場のエネルギーを差し引いた収支を確認しました。有限時間のR/Tは入射の切り詰めで
見掛け上1を超え得るため、初期場と観測時間を必ず確認します。

ソース取得者向けの再解析入口は `tools/analyze_fdtd_energy.py` です。
計算済み検証データだけを読み、ソルバーやMPIは起動しません。検証用C補助ライブラリを
明示的に構築・監査します。通常のwheel利用では保存済みNPZを `examples.fdtd_energy` に渡せます。
実物質計算は必須CIに追加していません。


### 横方向周期境界：横断面全体の流束と回折

非一様な面では `reflection_transmission(..., area_weights=dA)` を指定します。
E/Hは `(時刻数, Ny, Nz, 3)`、`dA`は `(Ny, Nz)`（単位bohr²）です。
`(時刻数, 面要素数, 3)` と1次元の重みでも使えます。全横断面・周期セル1個を
覆う面積重みを渡し、周期端点を二重に数えません。**各点のE×Hを先に計算してから
面積和を取ります。EとHを平均してから外積を作ってはいけません。**

```python
from salmonpy.analysis import surface_flux, reflection_transmission
from salmonpy.fdtd import energy_balance

# 各E/H: (N, Ny, Nz, 3)。同位置・同時刻へ補間済み。
# dA: (Ny, Nz)。一様格子なら各要素にdy*dzを設定。
ports = reflection_transmission(time, E_inc, H_inc, E_ref, H_ref, E_tr, H_tr,
    normal=(1, 0, 0), area_weights=dA, units='atomic_gaussian')
area = ports.config['surface_area']
check = energy_balance(point_work, point_volumes, ports,
    material_ids=material_ids, transverse_area=area,
    field_energy_change=delta_field_energy, units='atomic_gaussian')
net = surface_flux(time, E_total, H_total,
    normal=(1, 0, 0), area_weights=dA, units='atomic_gaussian')
# net.power: 面全体の符号付き電力[Hartree/au時間]
# net.energy: その時間積分[Hartree]。net.fluenceは面積平均[Hartree/bohr²]。
```

R/Tも各点の複素Poynting流束を面積積分してから比を取ります。面内平均で消える
位相を持つ回折次数、斜めに進む散乱光の法線方向成分も、入力に解像されていれば
含みます。真空参照との差を取る場合は **E/Hそれぞれを引いて反射場を作り**、
強度を引いたものを反射場として渡しません。斜めの回折波にH=法線×Eを仮定せず、
実際のHを使います。観測面は構造・PMLから離し、前後で同じ周期セル面積を使います。
有限時間の残留場・反応性の近接場・弱い駆動帯域では符号や比が不安定になり得るため、
validマスク、警告、時間窓と観測面位置を確認します。負の値を勝手にクリップしません。

横方向が同じ周期条件なら、対応する側面の符号付き流出は周期セルで相殺します。
前後の面積積分、物質への仕事、場の始終エネルギーで収支を作れます。外部からの
横方向光源・非周期側面・未観測の追加経路がある場合は別途計上が必要です。
`surface_flux` 自身はMPI通信を行いません。分割面では担当群でpower/energyを加算でき、
R/Tは分割ごとの比を平均せず、全横断面で評価します。閉曲面の自動配置や回折次数ごとの
分解、空間補間は自動実行しません。計算処理の通常point観測とplane出力は別であり、
バッチFDTDの自動取得は下記PlaneMonitorで指定できます。任意の外部配列も渡せます。旧point観測には200点上限が残ります。

Fortranでの面積・時間積分を使い、既存の一様ポート入口は変更しません。
エネルギーの二重計上を防ぐため `transverse_area` は記録された面積と一致させます。
NPZサンプル `examples.fdtd_energy` も任意の `area_weights` を読み取れます。

実Pythonから直接呼ぶFDTDの横方向周期・損失なし誘電体格子（z方向一様）の確認では、
32点で横断面全体を積分し、収支残差は入射の約0.0107%でした。場を先に平均すると
入射の約10.1%に相当する透過分を取りこぼしました。最終場の蓄積も計算処理の
電磁場保存データから評価しています。これは局所的な結合・観測の実装検証であり、
一般の3D材料の光学精度や多ノード性能の保証ではありません。


### 計算処理 FDTDの面取得と自動位置合わせ

バッチの `FDTD(..., planes=(PlaneMonitor(axis, position), ...))` で面を指定できます。
軸は0/1/2、位置はbohrです。最も近い格子面を選び、要求値と実位置を記録します。
横方向は周期境界、面と隣接面は通常の計算領域内でPMLから離して指定します。
旧来の点観測の200点上限とは別の経路なので、大きい面も点観測へ展開しません。

```python
from dataclasses import replace
from salmonpy.fdtd import PlaneMonitor
from salmonpy.geometry import Geometry

# calculationは横方向周期・通常のFDTD設定。管理ノードでは構築/prepareだけで終了。
monitored = replace(calculation,
    planes=(PlaneMonitor(0, x_front), PlaneMonitor(0, x_back)))
# 以下のrunは計算ノード/許可された計算環境で明示的に呼ぶ。
reference = engine.run(replace(monitored,
    geometry=Geometry(monitored.geometry.grid), materials=(), smoothing=None),
    'new-vacuum-reference')
result = engine.run(monitored, 'new-material-run')
plane = result.planes[0]
print(plane.time, plane.E.shape, plane.position)
net = plane.flux()                         # 取得・位置合わせ済みの面全体の流束
ports = result.reflection_transmission(reference)
check = result.energy_balance(reference, point_work, point_volumes,
    material_ids=material_ids, field_energy_change=delta_field_energy)
```

`PlaneData` のE/Hは `(時刻数, N1, N2, 3)`、area_weightsは `(N1,N2)` です。
Eは隣接する2本の辺、Hは4面をFortranで補間し、周期端も折り返して同じ空間位置へ揃えます。
計算処理が保存用に中心化済みのHをそのまま使い、時間平均は重ねません。
面と直前の隣接面を持つランクだけが、観測面ごとにストリームを追記します。
所有タイルを読み込んで全面を組み直す処理に追加のMPI集団通信はありません。
バッチの新規FDTDと、下記の明示的な再開区間の読み取り・連結に対応します。
embedded NativeEHFieldは下記の明示的なsample_planesスナップショットで取得します。従来の保存/再開の数値仕様は変更しません。

`plane.time` は0と観測間隔ごとの時刻です。stepsがobservation_strideで割り切れない場合、
最終観測は最後のstride境界です。物質への仕事と場の始終差もこの時間範囲へ合わせます。
面データだけから材料の蓄積や仕事をゼロと推定しません。これらは利用者が観測から明示します。
参照計算は同じ格子・入射・時刻設定と面配置を使います。E/Hは返却後読み取り専用です。
保存済み出力は `read_native_planes(directory, calculation, ranks=...)` でも再取得できます。

対応計算プログラムには `SALMONPY_NATIVE_PLANES_V1` が必要です。旧計算プログラム、面出力の失敗や
不完全なデータは警告して診断不可と記録し、正常終了したFDTDの点観測と計算結果は保持します。
`result.energy_balance` は不完全な面から収支を捏造しません。importや準備時に自動ビルドせず、
必要な場合だけ明示的に対応計算プログラムを取得/構築します。観測のrawデータ量の目安は
面ごとに `8×6×2×N1×N2×(1+steps//stride)` bytesです。時刻タグとタイルヘッダーが加わります。

実回折格子では観測なし/ありの既存点出力128個と、全電磁場保存データが従来実行ファイルと
バイト一致しました。横方向2ランク分割でも自動取得から収支まで一致し、残差は約0.0107%。
面と隣接面が別ランクになる法線方向分割でも、非ゼロの位置合わせ済みE/Hが逐次版と
ビット一致しました。これらは局所実装検証で、多ノード性能や光学物性の収束保証ではありません。

### バッチ再開の面観測データ

既存計算処理の再開処理・更新式は変更しません。再開側には、保存ステップより後の
観測間隔に沿ったサンプルだけが入ります。時刻ゼロの補充やHの再平均はしません。
計算処理再開入力は元の`nt_em`と観測間隔を維持し、利用者が計算ノードで明示的に実行します。
次の読み取り・連結処理は計算を起動しません。

```python
from salmonpy.fdtd import read_native_planes, join_plane_segments, FDTDResult
from pathlib import Path
# calculationは元の総ステップ数・形状を保持。
# early_dirは保存ステップで終わる完結した区間の出力。
early = read_native_planes(early_dir, calculation, ranks=2, stop_step=checkpoint_step)
late = read_native_planes(restart_dir, calculation, ranks=2, start_step=checkpoint_step)
planes = join_plane_segments(early, late)
result = FDTDResult(Path(restart_dir), (), planes)
budget = result.energy_balance(reference, work, volumes,
    material_ids=ids, field_energy_change=stored_delta)
```

`stop_step`の既定は元の最終ステップです。その時刻までの全サンプルが必要で、余分な
データも拒否します。書き込み中のファイルは渡さず、完結した区間を保存してください。
連結時は面の形状・全体時刻・観測間隔・総予定ステップ・保存ステップの接続を検査し、
欠落・重複・設定不一致を診断のエラーとして通知します。再計算・補間・重複の黙示削除は
しません。元配列は変えず、連結後も読み取り専用です。参照計算・仕事・場の蓄積差は
連結した全期間に合わせます。同じプロセス内での計算実行では下記の確定後だけ取得する入口を使います。

### 同じプロセス内での計算FDTD：確定後だけ面を取得する

同じ`FDTD.planes`を指定し、場の群の全ランクから`field.sample_planes()`を明示的に呼びます。
取得できるのはステップ0か`observation_stride`の観測時刻です。Fortranの`sp_eh_plane_read`は
担当する保存済みE/Hを読むだけで、MPI・ファイル出力・状態変更を行いません。Python側は面と
エラー情報を場の群の一度のgatherで集め、場のrank0だけに読み取り専用PlaneDataを返します。
他ランクはNoneです。全体通信・毎ステップの自動取得・再計算は追加しません。Hの再平均も不要です。
古いライブラリや診断不能時は場のrank0だけが警告し、Noneを返して計算を継続できます。

保護付きstep/復元の確定・候補解放後に呼んでください。未確定の候補がある間は警告してNoneを返し、
途中の状態を公開せずファイルも作りません。任意の`output_directory=NEW_DIR`で単一面サンプルを
既存の上書き禁止アーカイブへ保存できます。保存失敗も警告にとどめ、取得済みサンプルは返します。
計算状態と既存出力は保ちます。これは独立した観測記録であり、後日の巻き戻しで既存記録を削除する
機能ではありません。別の計算分岐には新しい保存先を使います。クラッシュ/OOM/MPI障害までの保証はしません。

```python
from salmonpy.fdtd import combine_plane_samples, FDTDResult
snapshots = []
for step in range(steps):
    coupled.step_protected(targets)  # 全参加群から呼ぶ
    if field is not None and (step + 1) % snapshot_interval == 0:
        sample = field.sample_planes()  # 場の全ランクから呼ぶ
        if sample is not None:
            snapshots.append(sample)   # 場のrank0だけ
if snapshots:
    planes = combine_plane_samples(*snapshots)
    result = FDTDResult(output_directory, (), planes)
    budget = result.energy_balance(reference, work, volumes,
        material_ids=ids, field_energy_change=stored_delta)
```

スナップショットは設定が同じで時刻が等間隔・昇順である必要があります。補間や重複削除はしません。
流束・スペクトルには最低2サンプルが必要で、場の変化を解像する取得間隔と同じ時刻の参照計算を使います。
仕事・場の蓄積差は取得した期間に合わせて明示的に渡します。面積×サンプル数に比例したメモリを使うので、
間隔を選んでください。混在サンプルの`--planes-x FRONT BACK --plane-stride 16`は横方向を周期境界にし、
確定したステップだけで取得して最後に一つの`plane-series`アーカイブへ保存します。既定の実行は変えません。

LDの実計算で直列/MPI2の直接・失敗時に元の状態を保つ復元と再試行を検証し、場の履歴と点出力は一致しました。
従来バッチ版とのE/H差は最大8e-19未満（比較条件rtol1e-12/atol1e-14）です。保存された数値状態は一致します。
以前からある未使用の座標配列領域は独立プロセス間で値が不定なので、全MPIランクの保存ファイル全体の
バイト一致とは区別しています。

座標保存の追補（2026-11-21）：新しい計算プログラムでは未使用の座標配列領域をゼロ初期化します。
有効な座標・halo・更新式は変えず、直列/MPI2で従来の場・電流・点出力に一致し、
保存/復元後も全ランクのチェックポイントがバイト一致しました。上記の未使用領域の差は修正前の記録です。


## 共通操作：引数・寿命・失敗時の扱い

`native_contract(component, route=...)['host_operations']` で、次の6操作の
契約を参照できます。同じ宣言からPythonのC引数型・戻り値型を設定します。
必要な関数をすべて解決してから操作を呼び、一つのランクで関数が欠けても
接続前に群全体で拒否します。既存の接続前合意を使い、時間更新中の通信は
増やしません。自作バックエンドは独自の契約を使えます。型宣言だけでは
バイナリの実装を証明できないため、監査済みビルドと呼び出し試験も必要です。

| C関数 | 引数 → 戻り値 | 通信・前提 | 効果と失敗 |
|---|---|---|---|
| `sp_attach` | `int64_t MPI_Fint` → `int` | 指定した有効な群内コミュニケータで集団呼び出し、未接続 | 通信群を複製。成功0、ハンドル幅・範囲の不一致−1、それ以外はMPI状態値 |
| `sp_detach` | なし → `int` | 接続群内で集団呼び出し、所有する状態はすべて解放済み | 複製した通信群を解放。呼び出し元のMPIは終了しない。成功0またはMPI状態値 |
| `sp_info` | `int* rank, int* size` → void | 接続中、ランク内だけ | 呼び出し元のバッファへ出力。エラー戻り値なし |
| `sp_barrier` | なし → void | 接続中、群内で同じ順番の集団呼び出し | 群内同期。エラー戻り値なし |
| `sp_sum` | `double` → `double` | 接続中、群内で同じ順番の集団呼び出し | 全ランクへ合計値を返す。エラー戻り値なし |
| `sp_library` | `char* buffer, int capacity, int* active` → void | ランク内だけ、接続前も使用可能 | capacity ≥1。末尾NULのMPI情報（容量不足なら切り詰め）と、MPI初期化済み・未終了ならactive=1 |

すべてメインスレッドで逐次呼び出します。ポインタは有効な呼び出し元所有の
領域を指し、戻った後にバッファを保持しません。元の通信群は接続中保持します。
info/sum/barrierは場・電子状態を変更しません。任意のMPI異常、プロセス終了、
無効ポインタは巻き戻し保証の範囲外です。voidや数値だけを返す関数から回復可能な
エラーコードは取得できません。Pythonは終了後の操作を拒否し、二度目のcloseは
何もしません。C関数のdetachには有効な接続が必要です。

ソルバー固有の生成・登録・複製、step/info、checkpoint/transaction、closeは
`standard_tddft_handles.h`、`standard_tddft_models.h`、`fast_hse_handles.h`
の契約に従います。各操作の集団呼び出し順、トークン、電流符号、対応条件を
共通表だけで置き換えることはできません。全操作・他環境での保証は残作業です。


## TDDFT・高速HSEの状態操作

組み込みインターフェースの
`native_contract(component, route=...)['state_operations']` は、必須C関数の
引数順・戻り値・通信範囲・所有権・モードを返します。同じ宣言から
single/shared/registryのPython型指定を設定します。必須関数をすべて解決してから
型指定・状態操作を行い、プロトコルと計算条件の検査も別途維持します。
任意のタイマー・電子数診断・メモリ内GS取り込みは、それぞれの機能問い合わせが
必要で、この一覧だけでは対応済みと判定できません。自作バックエンドは独立です。

ポインタ表記以外はCの値渡しです。`int64`は符号付き64bitの不透明トークンで、
メモリアドレスではありません。`string`はNUL終端したファイルシステムのバイト列、
`vector3`は連続したfloat64の直交座標3成分、`decomposition5`は連続したint32の5成分、
`history8`はint32の8成分です。`int*`・`int64*`・`double*`はスカラー出力ポインタ、
`buffer`は書き込み可能な呼び出し元の領域です。バッファは呼び出し後に保持しません。
プロトコル問い合わせはランク内のみ、状態操作はメインスレッドで、所有する群の
全ランクが同じ順に、モデルのI/Oディレクトリを指定して逐次呼び出します。
成功時の出力だけを使ってください。非ゼロの戻り値が常に巻き戻し成功を意味する
わけではなく、操作ごとの復帰契約に従います。

| 操作 | 寿命・振る舞い |
|---|---|
| handle_init / model_register | 世代付きの正の所有状態トークンを返す。登録は負のモデルトークンも返す。入力・GS・擬ポテンシャルは不変。各状態を一度だけ解放する。 |
| model_new / model_spawn | 元のGSから独立した状態を作る。他の状態の励起は引き継がない。TDDFTはdt・容量・並列分割・初期Aを継承し、HSEはdt・容量を呼び出しごとに渡す。 |
| handle_step | 原子単位の旧・新A/c。旧Aは受理済み端点と厳密に一致させる。状態ハンドルとHSEは電荷電流j、通常TDDFTのsingle C関数は物質電流J=−j。HSEは全エネルギーHartreeも出力する。 |
| handle_info / handle_label | 受理済みのstep/A/j、または状態専用の出力名を取得。名前バッファの容量はNULを含み、不足時13。 |
| handle_checkpoint | 0保存、1新しい状態への読込、2準備、3確定、4取消。パスはランク別でもよい。実行ライブラリ・入力・並列分割の整合性は呼び出し元が検査。直接読込だけでは計算全体の原子的復元にはならない。 |
| handle_transaction | 通常TDDFT：0捕捉、1巻戻し・ストリーム終了、2解放、3再開。HSE：0捕捉、1巻戻し・ストリーム終了、2再開・flush・解放（失敗時保持）、3再開。呼び出し元が巻戻しと再開の間に当該状態所有の出力だけを復元する。 |
| handle_close / model_close | 状態トークンを消費／モデルの新規生成を終了。既存の子状態は最後の解放までモデル資源を保持。再利用されたスロットでも古いトークンは無効。detach前に所有状態を解放する。 |

shared/registryには、11ランク間引数不一致、12古いトークン・接続寿命不一致、
13バッファ・容量不足などの拒否があります。他の番号は操作・計算条件ごとに
異なります。正確なC宣言は同梱ヘッダーを参照してください。singleの必須一覧は
現行の初期化・更新・終了の入口を対象とし、任意の全操作を保証するものでは
ありません。任意のMPI異常・STOP・無効ポインタからの復帰は保証しません。
他の計算条件・プラットフォームでの保証には別の検証が必要です。

各項目の`parameters`には引数名・型・in/outを順番どおり記録します。
`decomposition_order`はk点・軌道・格子x/y/zの順です。通常TDDFTのsingle C出力は
`J_matter`と明記し、電荷電流jと区別します。
