# 公開Python APIと互換性

[English](API.en.md) · Topic: `public_api` · Paired revision: `2026-11-14.3`

版は **0.1.0a5** です。`salmonpy.__version__` と配布メタデータに同じ版を示します。
プレビュー版であり、安定した1.0 APIを宣言するものではありません。物質C ABIのプレビュー1は[C_ABI.md](C_ABI.md)を参照してください。

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
| `distributed`, `material_response`, `embedded_tddft`, `coupled` | 既存の寿命制限がある実験的MPI・埋め込みAPI |
| `errors`, `snapshots`, `api` | 例外・不変の状態トークン・公開一覧 |

トップレベルの `__all__` は版、`LorentzDrude`、例外型に限定します。
任意依存を一括importしません。import時のコンパイル・取得・ネイティブ読込み・MPI初期化はありません。
snapshot／restoreは標準出力へ書きません。必須Python依存は引き続きNumPyです。
LD/SBEカーネルや外部ランタイムの明示的準備、任意依存の要件は従来どおりです。

## 互換方針

配布するAPI・挙動の変更には新しいパッケージ版を付けます。同じリリース内で単位・電流符号・
時刻配置・配列名を黙って変更しません。プレビュー版間では明記した破壊的変更を行うことがあります。
ホスト側は検証した版を固定し、移行説明を確認してください。将来の安定版では廃止予告期間を
設けますが、現在の実験的MPI API全体に安定性を保証するものではありません。
公開一覧への掲載と、機能の成熟度は別です。

アーカイブ・時刻・ランタイム・snapshotのスキーマ版は独立です。パッケージ版だけから
保存形式の互換性を推測しません。メモリ上のsnapshotは同じパッケージ版・同じ元モデルを要求します。

a5への移行で更新式・電流規約は変更しません。結晶は **A/c → 電荷電流密度j**、
内部の電子物質電流は J = -j です。過去のa4開発中に規約が変更されたため、古いデータは
記録されたメタデータで解釈し、一律に符号反転しません。基本はtransverseで、Longitudinalは明示的選択です。

## 例外の契約

| 型（`salmonpy.errors`） | 意味 | 従来の捕捉型 |
|---|---|---|
| `SalmonpyError` | ライブラリが明示的に送出する例外の親 | `Exception` |
| `InputError` | 値・互換性・状態遷移の不正 | `ValueError` |
| `InputTypeError` | 型の不正。`InputError`の子でもある | `TypeError` |
| `SnapshotError` | 不正・不適合な復元。状態は保持 | `ValueError` |
| `NumericalError` | 検出した非有限値などの数値異常 | `FloatingPointError` |
| `RuntimeFailure` | ランタイムの欠如・失敗・結果の未確認 | `RuntimeError` |

既存の `RuntimeUnavailable` と `OutputSaveError` は名前を維持し、`RuntimeFailure`を継承します。
従来の組込み例外型での捕捉も利用できます。OS・ファイルの例外（`FileExistsError`を含む）、
プロセスの時間切れ、NumPy変換、任意依存の不足、利用者コールバックの例外は元の型を保つ場合があります。
外部例外を一律に包みません。収束報告がないことはランタイムの失敗であり、数値発散の証明ではありません。

例外の型だけから計算全体の巻き戻しを推測しません。内蔵LD・SBE・Pythonの場のstep／restoreは
失敗時に受理済み状態を保持します。MPI応答・結合の失敗はセッションを無効化し得ます。
FortranのSTOP／abortを回復可能なPython例外にする変更ではありません。

プロセス分離`TDDFTState`は`save_checkpoint(directory)`と
`restore_checkpoint(directory, max_bytes=None)`も公開します。
同じ内容の初期モデルを明示的に用意し、Zarr履歴から別プロセスで再開できます。
初期GS本体は保存しません。詳細は[状態保存](SNAPSHOTS.md)を参照してください。

直列モデルの共通契約、適合チェッカー、MPI点分散への接続は[MATERIALS.md](MATERIALS.md)を参照してください。

`fast_hse` の `FastHSEMaterialProvider`／`FastHSEResponseBackend` は明示ビルドした高速PR版serialワーカーを外側MPIの物質応答へ接続するプレビュー入口です。COUPLED.mdを参照してください。HSE内部MPIを初期化しません。

`fast_hse_mpi.NativeFastHSEBackend` は高速PR版HSEの実験的なホストMPI入口です。k点のみ1／2ランクの物質群と、同じ分割での直接保存／復元に対応します。旧create_stateはメインスレッド・1プロセス1回のnative寿命です。open_modelはNativeFastHSEModelを返し、同一の全k点モデルの状態を共有接続で独立に扱います。SHARED_HSE.mdとCOUPLED.mdを参照してください。

`checkpointing.CheckpointCapability` と `check_checkpoint_compatibility` は、自作バックエンドにも使える保存／復元前の共通宣言・検査です。[CHECKPOINTING.md](CHECKPOINTING.md)を参照してください。

`checkpointing.CheckpointTarget` と `execute_checkpoint` は検査後にバックエンド独自の保存／復元関数を呼び、不透明な戻り値をプロセス内で返します。

`CheckpointState`／`restore_checkpoint` は、独立候補を全ランクで準備した後にまとめて公開する保護付き復元です。`AtomicRestoreError`で元のホルダー維持と候補の破棄エラーを区別します。

メモリ助言の公開入口はsalmonpy.memoryのMemoryEstimate／MemoryReport／memory_advice／process_memoryです。予算超過は警告のみで、実行可否や上限の保証ではありません。詳細は[MEMORY.md](MEMORY.md)。

`stepping.StepTarget`・`guarded_step`・`StepRecoveryError` は明示した完全退避契約によるステップ保護のプレビュー入口です。STEPPING日英ガイドを参照してください。

`native_steps.native_step_targets` は、監査済みnative状態と通常出力を共通の保護付きステップへつなぐプレビューAPIです。NATIVE_STEPPING日英ガイドを参照してください。

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
これは初期配列構成の選択であり、占有容量・物理的妥当性・SCF入力全体の所有や、native初期化からの
異常復帰を保証するものではありません。新しいPython／C公開入口もなく、k点・対称性・
コミュニケータなどの共有依存は引き続き残ります。

HSEは高速PR版のみです。`NativeFastHSEBackend(build_record, work_root, allow_symmetry=True)` で
対称性を明示的に許可できます。内部kランク1／2・周期・非スピン分極・CPU・full-qのTaylor／PT-CNに対応し、共有の対称化設定は不変、1プロセス1寿命です。既定値Falseと既存C入口は
対称性を拒否します。内部2ランクはsp_fast_hse_symmetry_kmpi_protocol()==1が必要。Taylor strict／learned・PT-CNにはsp_fast_hse_symmetry_history_protocol()==1が必要です。学習PT-CNは既存方式どおりinterval=1・predictor_corrector=nを指定します。

新C入口 `sp_fast_hse_init_symmetry` は既存initと同じ引数で、
`sp_fast_hse_symmetry_protocol()==1` を宣言します。全ランクの有限なstatus34とこの宣言を
確認した場合だけ、Pythonは状態を変えない `InputError` を返し、Aを直した再試行を許します。
未知・旧ABIの失敗は継続しません。FDTD結合では全体の再試行に `step_protected` が必要です。
直接／候補復元は共通checkpoint経路を使います。

[対称性ガイド](FAST_HSE_SYMMETRY.md)には、混在サンプル、対応範囲と極小勾配XCの安定評価と回帰検証を
記載しています。接続・保存復元のテストを物理収束の証明と解釈しないでください。

縮約k点の内部2ランクと同一分割での直接再開に対応します。追加のnative k-MPI capabilityが必要です。分割変更による再開は未対応。詳しくは対称性ガイドを参照してください。

NativeFastHSEBackendはnative寿命前に周期セル・格子・原子入力を読み取り検査し、群内でInputErrorを合意します。入力を直してproviderを作り直せます。再開バイナリの互換性や全native初期化の安全性は別契約です。詳しくは対称性ガイドを参照してください。

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

## native実行経路の契約情報

```python
from salmonpy.api import native_contract
contract = native_contract('fast_hse', route='shared')
assert contract['model_settings_owner'] == 'explicit model'
assert contract['runtime_verified'] is False
```

`native_contract`は組み込み`tddft`／`fast_hse`の`legacy`・`shared`・`registry`経路、
および`fdtd`の`legacy`経路を説明します。ライブラリ読込・計算・新しい制限は行いません。
プロセス互換設定と明示的モデル所有、MPIの呼出し元所有、メインスレッドの逐次操作、collectiveの順序、電流符号、失敗保証範囲を区別します。
標準TDDFTのlegacy C入口は物質電流J=-jを返し、Python側は電荷電流jへ変換します。
高速HSEと標準handleのC入口は電荷電流jです。照会がローカルかcollectiveかは個別ヘッダーに従い、診断を一律collective扱いしません。

これは宣言された仕様の情報で、実ランタイムの能力確認や安定C ABIの宣言ではありません。`runtime_verified`は常にfalseです。
実際のランタイム・保存再開の対応可否は別途照会します。未知の自作コンポーネントや経路名には`known=False`・`status='unverified'`を返し、
自作バックエンドを拒否しません。その契約はバックエンド自身が定めます。結果は独立した辞書で、呼出し側で注記しても後の照会に混入しません。
既存の計算呼出し、数値状態、legacy経路の省略時処理、計算可能な条件は変えません。C／Fortranホストの安定ABI適合検証は別の残課題です。

### 明示指定した実ランタイムの版照会

```python
from salmonpy.doctor import probe_native_protocols
report = probe_native_protocols('audited-build.json', component='tddft', route='shared')
# errors・照会結果・protocol_matchを確認。runtime_readyは未確認のままです。
```

CLIは`python -m salmonpy.doctor --native-build audited-build.json --component tddft --route shared --json`です。
`native_contract`や独立静的チェッカーとは異なり、明示指定されたnativeライブラリを、記録のSHA256と照合して別プロセスで読み込みます。
子側でも読込前にhashを再確認します。hashは同一性であり、信頼性やビルド監査の証明ではありません。
読込時に初期化コードが動く可能性があるため、信頼する監査済みライブラリだけを対象にしてください。
既知の引数なし版照会関数だけを呼び、ソルバー状態の作成、MPI初期化・接続、ビルド、インストール、計算は要求しません。
未知の自作経路は記録ファイルも開かず、コードを読み込みません。

`protocol_match`は選択した版照会だけの一致です。関数の欠如や新しい版は記録するだけで、自動的なバックエンド拒否や別機能への置換はしません。
`runtime_ready`は未確認のままで、操作関数のシグネチャ、物理的対応範囲、保存形式・MPI分割の互換性の保証ではありません。
子プロセスの異常は呼出し元を落とさず記録し、`library_loaded=None`は子側が読込完了の有無を報告できなかった意味です。
渡されたコードを読み込まない`check_job`とは別の、任意の明示的診断として扱います。

CとISO_C_BINDING Fortranの版照会例は`tests/native/test_native_protocol.c`と、
`tests/native/test_io_ownership.f90`の`SALMONPY_ABI_PROTOCOL_PROBE`モードです。
開発用の実行入口は`tools/test_io_ownership.py BUILD_RECORD NEW_OUTPUT --abi-contract`で、
高速HSEの場合は`--fast-hse`を追加します。ランタイムhashと全リンク対象のアプリケーションオブジェクトを監査し、
返却値の一致とMPIが未初期化のままであることを確認します。これは版照会の相互呼出し検証で、モデル・状態・step全体のABI適合確認ではありません。
native ABIは引き続き実験段階です。

### 状態の作成・更新・終了を確認するホスト試験

`tools/test_io_ownership.py BUILD_RECORD NEW_BUILD_DIR --abi-state --fixture EXISTING_PZ_FIXTURE --state-output NEW_RESULTS_DIR --ranks 2` は、監査済みの通常TDDFTライブラリをCとISO_C_BINDING Fortranから呼びます。高速PR版HSEのビルド記録と既存の全k点HSE入力には`--fast-hse`を付けて同じ検査を実行できます。HSE試験では各実行のk点並列数だけを設定し、渡されたGSは変更しません。縮約k点入力にはさらに`--state-symmetry`を付けます。渡された`sym.dat`を参照し、対称性対応のHSE初期化を使い、この入力の保持対称操作に合うz方向を駆動します。代表k点が実際に出力されたことも確認します。この試験入力の方向指定であり、全k点でのx・y・z駆動を制限しません。MPIの初期化・終了は呼び出し側が所有します。2状態と個別実行の一致、直前Aが不一致のときの状態保持、スロット再利用後の古いハンドルの拒否、新しい状態の独立した初期化を確認します。既存GSはその場で読み、変更されていないことを確認します。GSのコピーや電子状態チェックポイントの保存はしません。結果の保存先は明示してください。

ディレクトリスコープの共通宣言は、両ハンドルヘッダーから取り込む`native_context.h`にあります。macOSではMPIホストのリンクにも既存の`tools/build_embedded_tddft.py::mpi_link_command`を使ってください。OpenMPIのFortran COMMON特殊値が実行ファイル内に別配置されると、共有ライブラリとの間でMPI_IN_PLACE処理の判定が変わる場合があります。試験ランナーはこの共通処理を使い、リンク引数と実行ファイルのハッシュを記録します。版の問い合わせだけではこの問題は検出できません。

確認範囲は通常PZと高速PR版HSEの全k点／縮約k点入力、C／Fortranホスト、MPI 1／2ランク、複数状態と個別実行の電流一致です。HSE固有の引数と拒否コード（直前Aの不一致はHSEで1、通常TDDFTで10）を共通試験内で扱い、公開ABIや数値ソルバーは変更しません。縮約入力は対称性対応PT-CNと学習型交換を使い、全8点を4代表点で扱います。他XC・形状・他環境の検証は別途必要で、公開native ABIは引き続き実験的です。この小規模入力は接続と状態所有を検査するもので、実材料の収束を保証しません。

### ホスト試験のPZ入力を準備する

`--abi-state`の標準PZ入力は`rt.inp`、`data_for_restart/`、`pseudo1_rps.dat`を必要とします。`theory='tddft_pulse'`、周期・固定イオン・非スピン分極、`middlepoint`、`trans_longi='tr'`を使い、試験用RT入力から`&parallel`を除いてください。同じランナーが1／2ランクを試すため、固定した並列条件を埋め込まないひな形を提供します。GSのセル・原子・格子・k点・電子数・擬ポテンシャル・XCは入力と一致させてください。

```sh
python tools/prepare_pz_state_fixture.py GS_INPUT EXISTING_RESTART PSEUDO NEW_FIXTURE
python tools/test_io_ownership.py BUILD_JSON NEW_BUILD --abi-state --fixture NEW_FIXTURE --state-output NEW_RESULTS --ranks 2
```

SALMONpy導入済みのPython環境で実行します。`GS_INPUT`には元のGS入力または同じGSから検証したRT入力を指定します。この準備ツールはPZ・周期・単一元素・インライン原子座標の試験用であり、GSを計算・コピー・ダウンロードしません。再開ディレクトリと擬ポテンシャルは既存ファイルへの参照です。小さなGS入力の例は`examples/pz_host_gs.inp`（Si1原子、8³格子、2×1×1 k点）です。必要なら利用者が計算ノードで明示的にGSを準備し、次にこのツールへ渡します。収束した物性の例ではありません。

`BUILD_JSON`にはビルド記録のファイル、または`build.json`を含むディレクトリを指定できます。`--parser-host`はGSを使わず、`status`単独／`native_host=.false.`で通常検査を維持し、明示trueだけが場用検査へ切り替えることをMPIで検査します。

ランタイムバンドルからは`SalmonExecutable.from_runtime(directory)`を使えます。これは`runtime.load_runtime`へ委譲し、全ファイルと能力記録を検証します。バイナリ単体とハッシュを手作業で組み合わせる必要はありません。

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
