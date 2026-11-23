# SALMONpy：AI向け利用・実装ガイド

## No backward-compatibility fallback

Use current inputs/APIs only. Runtime preparation/build accepts one output directory. Shared-library directory protocol3 is required; never emulate an older runtime with Python chdir. Require GS restart_saved, SCF rebuild/2 codec/state revision, direct/4 or replay/2 mixed checkpoints. Do not infer missing metadata, use the removed unidentified-checkpoint option or accept old worker replies without measured electron number. For HSE pass an owned model_settings; never use implicit global settings. Fortran inputs use exx_mlwf_* and exx_k_algorithm, without old aliases. RT restarts require rho0_s.bin. Current custom backends and numerical/platform choices remain supported; these are not version fallbacks.

利用者が行いたい計算から入口を選びます。最初に `manual.json` の検索と[サンプル集](examples/README.md)を使い、既に分かっている条件を再質問せず、目的・使用可能な総メモリ・CPU数の不足分だけを確認してください。

| 目的 | 使う入口・ガイド |
|---|---|
| 材料応答・パルス・構造・解析を単独で使う | `ld`・`sbe`・`fields`・`geometry`・`analysis` |
| 保存GSから電子状態を進める | `tddft`、[利用手順](WORKFLOWS.md) |
| 電子計算とFDTDを一ステップずつ結合する | [結合手順](COUPLED.md) |
| 同じ物質の独立した計算を複数作る | [TDDFT](SHARED_TDDFT.md)・[高速HSE](SHARED_HSE.md) |
| 異なる物質入力を同じMPI群で使う | [TDDFTの登録](TDDFT_REGISTRY.md)・[HSEの登録](SHARED_HSE.md) |
| 保存・再開や失敗時の状態保持を使う | [保存・再開](CHECKPOINTING.md)・[時間更新](NATIVE_STEPPING.md) |

各ガイドの実行手順・使用するライブラリの条件を確認してください。開発・検証の補足や `docs/plans/` の過去の「未接続」「今後」の記述を、現在の対応機能一覧として扱いません。公開API一覧に名前があることと、使用する計算ライブラリで動作確認済みであることは別です。存在しないAPIや未測定の物理精度を案内しないでください。

## 物質応答の境界

結晶の公開インターフェースはA駆動のみです。TDDFTとSBEは、内部状態と新旧Aを受けて時間発展しJを返す交換可能な物質応答として設計します。結晶にE駆動の選択肢を追加しないでください。SBEは既存の演算順を保った共有Fortran計算ルーチンに接続しています。

E,Hを扱うFDTDとA駆動の結晶の間には結合処理を置きます。単位、ゲージ、時間位置、符号、座標回転、Jの体積規格化は現行マルチスケール実装に準拠して確定し、アダプターで明示します。Aの積分履歴は結合側、波動関数・密度行列は物質側が保持します。

ドルーデ＝ローレンツはE駆動の構成則として独立させます。結晶用A駆動APIに合わせるための数式変更はしません。

## 今回利用できるFortran境界

`native/lorentz_drude.f90` の `s_ld_state` がP、J、旧J、旧E、係数、次の電場更新用の補助量Qを所有します。公開関数は `ld_allocate`、`ld_coefficients`、`ld_store_old`、`ld_update`、`ld_finalize` です。具体的な配列形状・呼び方は `tests/native/test_ld.f90` と、外部に準備したSALMONソースの `src/maxwell/fdtd_eh.f90` を参照してください。

`ld_coefficients` の累積先は呼び出し前にゼロ初期化します。`ld_store_old` と `ld_update` は既存演算順を保つため分割E成分を受け取ります。`rj*_fdtd_ld` は物理電流の単純和ではなく補助量Qです。

`fdtd_eh_kernel.f90` は差分演算 `eh_fd` を提供します。`fdtd_eh` は引き続き時間ループ、境界、通信、観測、再開ファイルを調整します。実験的な単一ステップ・結合入口は `salmonpy.coupled.NativeEHField`。COUPLED.mdの範囲に限定します。

## 変更してはいけない更新順

旧J/E保存 → E更新（LD係数を含む）→ 光源 → 旧Q注入 → E通信 → LDのJ/P/Q更新 → 線形応答 → 旧H保存 → H更新 → H通信 → 観測・保存。

LDは電流を返すだけでは分離できません。E更新係数への寄与も必要です。係数を固定したまま時間刻みを変更しないでください。独立計算には別々の状態インスタンスを使います。

## 検証と禁止事項

ビルド前に `python tools/check_forbidden_eigensolvers.py --source-root native` を実行します。構成・増分ビルド時の検査を外しません。新規実行ファイルは全アプリケーションオブジェクトを検査し、再利用時は記録されたSHA256と照合します。禁止された固有値ドライバーへのフォールバックはありません。PZHEEVDの検証には固有値、残差、直交性、ワークスペース容量を含めます。

IEEE判定は各要素へのスカラー呼び出しに限定します。GSは利用者が指定した保存先に保存して再利用し、不要な検証用複製は収集しません。計算結果は利用者が指定した書き込み可能な出力先へ保存し、コードツリーに格納しません。

変更後は独立LDの旧演算との比較とFDTD全体の比較を行います。単一スレッドのビット一致と並列の数値許容差検証を区別します。再開互換性や未検証の計算機での動作を推定だけで保証しないでください。

## 設定を決めるタイミング

INPUT_CONTRACT.md に従います。材料係数はモデル作成時、dtと初期条件は状態作成時、Eは各stepで指定します。`material.config` と `state.config` は確定値の読み取り専用レコードです。設定変更は新しいモデル・状態の作成で行い、private属性を書き換えません。`respond` は毎回新規状態です。設定レコードはチェックポイントではありません。

## 配布

本体計算プログラムの準備は RUNTIME.md / RUNTIME.en.md と manual.json の `runtime` タスクに従います。既存の監査済み計算プログラムを再利用し、次にOS・CPU・必要機能に合う配布物を探します。対応品がない場合だけ明示的な `python -m salmonpy.runtime build` を案内します。importやinstall失敗を理由に自動ビルドしません。独立SALMONpyの非公開Releaseと期限付きActions配布物を区別し、公開URLを作り話で補わないでください。Libxcを必要とする計算には `--require libxc` を付け、PZ版に勝手に代替しません。取得・検証した `load_runtime` の戻り値を既存の `executable=` に渡します。

wheelの構築と修復・OS別検証は PACKAGING.md に従います。検証済み環境とCI設定のみの環境を混同しません。wheel利用者にはFortranコンパイラを要求せず、ソース構築時のみ要求します。

構造・波形・解析の単独入口は WORKFLOWS.md に記載しています。FDTDの一ステップAPIはCOUPLED.mdの条件のみ例示できます。周期構造、斜め入射、未対応のTDDFT計算方式を対応済みとして例示しないでください。共通Acos2パルスとFDTD側のCEP規約の差はアダプターで0.25周期補正します。


## SBEの接続

`salmonpy.sbe.SBEInitialState.from_salmon` は原子単位・単一スピンの
`_k.data`, `_eigen.data`, `_tm.data` を読みます。`bands` はファイルの全バンド数、
`active_bands` は任意の伝播バンド数です。占有数は現行SBEと同じ閉殻規約です。
`new_state(dt=..., initial_A=...)` ごとに独立した密度行列を生成します。
結晶はA/c駆動のみ。中点Aで四次Taylor伝播、新しい端点AでJ評価の順を保持します。
補正次数と非局所補正は初期モデルで固定し、dtは状態作成時に固定します。
SBEの長さゲージ・バイナリ再開は未対応です。四次Taylor法なので厳密な
ユニタリー性や自動誤差制御を保証しません。APIは刻み幅を自動変更しません。


## TDDFTの呼び出しと終了

`salmonpy.tddft.TDDFTInitialState.from_salmon` は自己完結した入力ディレクトリを
複製し、ハッシュ付きの初期スナップショットを所有します。初期GSのバイナリ解釈は
SALMONの既存ローダーが行います。保存用バンドルは1組を状態間で共有し、各 `new_state` は専用プロセスと独立した作業コピーを所有し、
波動関数、Vh履歴、予測子修正子、XCや通信の共有変数を他の計算と共有しません。
`with` で解放してください。共通 `respond` は成功・失敗とも自動的に解放します。
PZ／TB-mBJ／Libxc-PZ／Libxc-PBE、周期系、固定イオン、非スピン分極、中点伝播、CPU一MPIランクが対象です。
実行ファイルは検証済み配布物を取得するか、明示的に構築・監査し、SHA256を指定します。wheelだけでTDDFTが動くとは
説明しないでください。dt、最大ステップ数、初期Aは状態作成時に固定します。
共通の出力は電荷電流密度 j です（2026-11-14 から）。SALMONの電子物質電流 J = −j は状態の内部に保持し、`state.matter_current` で参照できます。符号を二重に反転しないでください。
異常終了時は受理済みA履歴を再実行して復元し、復元失敗時は状態を閉じます。
実GSによる検証と合成波動関数による接続回帰を区別してください。

Libxc-PZ/PBEには `USE_LIBXC=ON` の外部SALMONが必要です。Libxc全汎関数の対応ではありません。GSとRTで汎関数・擬ポテンシャルを一致させてください。現ローダーはGS生成条件との整合を自動判定しません。PBEの検証はヘッダーでPBE・NC・スカラー相対論・NLCCありを確認したSi UPFを使用しました。詳細は docs/libxc-validation-2026-11-13.md を参照してください。

## 後処理用アーカイブ

OUTPUTS.md に従います。`crystal.respond(..., output_directory=...)`、
`TDDFTState.finish(directory)`、`outputs.OutputArchive` と `outputs.save_arrays` が
実在する入口です。TDDFTのfinishは状態を閉じ、初期ファイルハッシュと選択した
本体出力を残します。大きい分布は初期入力のanalysisで選択・間引きます。
波動関数チェックポイントとprojection/unfoldingは未対応です。NPZアーカイブを
restartファイルとして渡さないでください。別ノード通信やコア固定は未実装・未検証です。

## バンド図

BANDS.mdを参照してください。`bands.BandStructure.from_salmon` と `.plot` が入口です。
固有値の単位、実際のCartesian k、セルの選択、余剰点、収束対象バンド、エネルギー基準を
明示してください。経路上の価電子帯上端を全域のVBMと呼ばないでください。

バンドの計算入口は `BandInitialState.from_salmon` → `run(BandPath, ...)`。
初期束の外の新規保存先、対応した監査済み実行ファイルを使い、SCFが反復上限で未収束なら警告とSCFReviewRequiredで一旦止め、残差履歴を利用者が確認するまで後続へ進めません。実行完了とSCF収束を区別し、`GroundStateResult.converged` と `ground-state.json` の `converged` を確認して未収束を収束済みと説明しません。残差・閾値・反復数は保存されます。破損出力・非有限値・計算処理実行失敗は引き続きエラーです。

ASE連携は任意の独立モジュール `salmonpy.ase`。ASE.mdを参照。
`run_bands` は構造を照合し、ASEのk点列をそのまま計算します。Calculator・力は独立したASE_CALCULATOR.mdの対応範囲を参照してください。

GS作成は `salmonpy.dft.DFT`、ASE変換は `salmonpy.ase.dft_from_atoms`。DFT.mdに従います。ASE入口の電子数はKY/NC-UPFの価電子数から中性系として決定・照合します。直接DFTでは電子数を明示します。軌道数・格子条件・k点は必須で、収束保証のある自動既定値はありません。
`salmonpy.ase.tddft_from_atoms` はGSを保存して既存TDDFT初期状態を返します。ASE.en.mdとexamples/ase_tddft.pyも参照してください。結晶はA/c駆動、戻り値は電荷電流j。状態は独立で、計算プログラムの自動取得・コンパイルはしません。
収束済み `GroundStateResult` からバンド／TDDFT初期状態を作れます。PZ／Libxc-PZ／Libxc-PBE・KY／NC-UPFのCalculatorと本体の力はASE_CALCULATOR.mdを参照してください。

`salmonpy.ase_calculator.SalmonCalculator` は任意依存の独立モジュール。energy/forcesのみ、PZ／Libxc-PZ／Libxc-PBEとKY／NC-UPF・単一ランクが対象です。Libxc対応計算プログラムと擬ポテンシャルの整合を確認します。
NLCC付きUPFの力は検証記録なしではCalculatorが計算前に拒否します。信頼できるバンドルをload_runtime(..., required=['nlcc-forces'])、または監査済みビルドと検証記録をload_validated_buildで読み込んだ場合のみCPU・非スピン分極NC・PZ/Libxc-PZ/PBEの力を有効にします。SalmonExecutable(path,sha256)だけでは無効です。記録はバイナリ・既知パッチ・微分検査・固定Si8実行検査に結び付け、移動後にも実行検査します。互換パッチのXC重み・全密度勾配・明示的コア力の修正は旧バイナリに自動適用されません。実装回帰と物理的な収束を区別し、メタGGA・OpenACC・スピン分極Libxcの力を検証済みと説明しないでください。
`warm_start=True` は直前の互換なGSの波動関数だけを再利用し、既定Falseです。セル・計算条件が変われば新規SCF、元データが破損すればエラーです。旧原子位置・混合履歴は再利用しません。
結晶TDDFTの基本はtransverse（trans_longi=tr）。Longitudinalは明示的な場の選択です。孤立系・non-collinear・Ehrenfestは初版の目標として残しますが、今回のASE APIの対応済み機能には含めません。
粗い格子のエネルギー有限差分はegg-box効果を含むため正解と仮定せず、剛体並進と格子収束を調べます。

## 初心者への最初の案内

manual.jsonの `beginner_intake` を必ず適用します。新しい計算を選ぶ前に、未回答の
「使ってよい総メモリ」「使ってよいCPU並列数」「やりたいこと」を短くまとめて聞きます。
既に回答された条件は再質問せず、「不明」も受け付けます。装置の全メモリ・全コア数を
自動的に利用許可と解釈しません。目的に応じて構造規模・精度・所要時間の希望を補います。

最初に、今できること／追加準備が必要なこと／未対応と、予算内で始められる小規模例を示します。
計算規模の上限は実測・暫定見積もり・不明を区別します。未測定の原子数上限などを断言しません。
必要なら予算に収まる小規模試行でピークメモリ・時間・保存量を測定し、条件を記録して更新します。
見積もりはメモリ強制制限ではありません。自動容量推定器・ハードメモリ制限は現状未実装です。

プロセスとスレッドを合計したCPU使用量、全子プロセスを含むメモリ、保存容量を扱います。
メモリの余裕は許可された上限の内側に確保します。現在のPython起動が単一MPIランクである
機能を、任意のMPI並列数で実行可能と案内しません。収まる計算規模と物理的な収束精度は別です。

実行可能な用途別Pythonサンプルと必要な入力は [examples/README.md](examples/README.md) を参照してください。AI検索の各taskの `sample_files` からも辿れます。

## 初見の環境診断と検索

`python -m salmonpy.doctor --json` でライブラリの読み込み可否、追加依存の所在、
任意指定した外部実行ファイルのハッシュ照合を確認できます。計算・コンパイル・導入は行いません。
外部SALMONの `hash_matches` を機能対応や実行成功と読み替えないでください。
`manual.search` は常に `getting_started` と `task_catalog` を返します。一致しない検索語では
目的別一覧を示して目的を確認し、対応不能とは即断しません。

本リポジトリはSALMONpy専用です。同梱本体ソースはruntime/sources/solverの一系統です。ビルド用srcは `runtime/README.md` の手順でリポジトリ外に準備します。

## 並列経路を選ぶ

LD中心の光伝播には `fdtd_mpi` タスクを選び、既存MPI FDTDの各ランク内で多数のLD状態を管理します。
起動元は直列Pythonです。MPI版の完成ビルド記録が必要で、配布用直列計算プログラムとは区別します。
所定のAに対する独自応答は `distributed_response`、二点TDDFTは `distributed_tddft` を選びます。
これらはPythonを一度MPI起動し、標準例では一物質群に一物質・一マクロ点を割り当てます。
内部MPIのTDDFTとE/Hフィードバックには `native_coupled` タスクとCOUPLED.mdを使います。

## 場の部品（j→A）

`salmonpy.maxwell` の場の部品は `step(j_n) → A_{n+1}` です。`Transverse`／`Longitudinal`／`ThinFilm`／`Multiscale1D` は SALMON の tr／lo／2d／1次元多スケールと同じ式・順序・初期条件です。部品間のインターフェースでは電荷電流 j を受け渡し、物質の内部は SALMON の物質電流 J を保持します（j = −J）。crystal プロトコルは j を返します。J_matter を宣言した部品だけ `maxwell.couple` が変換します。場の部品に J を直接渡さないでください。t=0 の電流は場に戻しません。1次元では x が伝播方向です。座標回転、背景誘電率のある1次元媒質、2次元・3次元の Weyl 伝播、再開は未対応として扱い、対応済みと例示しないでください。時刻配置は `salmonpy.timing`（`salmonpy.time/1`）で宣言し、`couple()` が結合前に検査します。新しい部品には `timing` を実装し、検査を回避する変換を書かないでください。詳細は FIELDS.md を参照してください。
これは一群一点・一プロセス一回の計算処理寿命に限定した実験的APIです。
既存FDTDバッチと応答サービスを並べるだけで結合したと説明しません。
Pythonから直接呼ぶFDTDと電子計算の結合のE時刻n·dtと電子時刻(n+1/2)·dtを混同しません。再開、回転、サブサイクルは未対応です。
実行前に不足する資源予算を確認し、例のランク・スレッド・状態数と、未測定のメモリ使用量を明示します。
機能ごとの実行例・出力形式は examples/README.md / README.en.md にあります。

## 次元選択の推奨

電子応答を結合する計算では、TDDFTは負荷が大きいため、まず1次元の伝播・層構造を提案します。
3次元を推奨するのはSBEまたはユーザー作成のカレント計算コードを使う場合です。
3次元TDDFTが不可能とは説明せず、必要性と状態数・メモリ・時間の予算を確認します。
ここでの次元は場・マクロ配置で、微視的な電子状態の次元ではありません。
現行40³・二点サンプルは接続回帰用です。自動1次元化、未実装SBE結合アダプター、
周期境界を対応済みとして案内せず、実在する入口の制約を確認してください。LD単独FDTDは別の経路です。

## 資源計測と計算の成否を分ける

導入・計算・数値検証・資源計測の結果を別々に記録します。任意の計測器や
`time -l`、sysctl、監視APIの失敗だけで、成功した計算を失敗と判定しません。
計算プロセス自身の終了コード・ログ・出力検証を根拠に判定し、非ゼロ終了や
数値検証の失敗を計測問題として除外しないでください。ラッパーの終了コードしか
分からない場合は計算成否を不明とし、既存出力を保存して、新しい出力先で再試行します。

取得できないピークメモリ等は unavailable / not_measured と null で記録し、0や推定値で
代用しません。計測対象（親プロセスのみ／子を含むジョブ全体）・単位・取得方法を残します。
計測に失敗しても予算内という保証はできません。ハード上限が必要な実行では、有効な
制限手段を確認するまで実行・拡大しません。初回と環境修正後の結果は分けて残します。

wheelからの利用テストはソース側の `docs/wheel-user-test.md` / `.en.md` と
`tools/test_wheel_user_workflow.py` を使います。指定wheelを新規venvに導入し、
LD・合成SBE・波形・形状・保存の公式例を実行します。ソースビルドへ自動切替せず、
wheelハッシュと実際の読込元を記録します。これはスモークテストであり、先行Notebookの
追加の解析解・保存則検証や、外部SALMON・MPIの検証を含みません。

## 公開仕様とメモリ上の状態保存（0.1.0a5）

`public_api` と `snapshots` タスクを参照。公開名は `salmonpy.api.PUBLIC_API`、
版は `salmonpy.__version__`。例外は `salmonpy.errors` の親と分類を使い、
既存の組込み例外による捕捉も維持します。OS・依存・コールバックの例外がすべて
SalmonpyErrorになるとは説明しません。詳細はAPI.md/API.en.md。

プロセス分離TDDFT・LD・SBE・Pythonの4種類の場が `snapshot()`／`restore(token)` に対応します。
同じモデルオブジェクト、同じdt・開始時刻と場の容量・駆動条件を要求します。
保存はメモリ上の不変コピーで、失敗したrestoreは状態を変えません。保存数と状態サイズの
分だけ追加メモリが必要です。トークンは編集・シリアライズせず、Zarr GSや出力NPZと混同しません。
手動結合では同じ最後に成功した計算時刻で場と物質を両方保存します。MPI全体の原子的復元、可変dt、
Pythonと同じプロセスで実行するTDDFTや計算処理E/Hの復元は未対応です。プロセス分離TDDFTのディスク再開は下記の専用APIで対応し、
出力アーカイブや計算処理波動関数checkpointとは区別してください。
examples/snapshot_restore.pyは実行可能な例です。結晶の基本はtransverse・A/c駆動です。

プロセス分離TDDFTのrestoreは別ワーカーでA履歴を再生し、各端点電流をrtol=1e-12・atol=1e-14で照合後に切替えます。失敗時は元ワーカーを保持します。同一初期モデル・dt・初期A・容量・実行ファイル・スレッド数を要求し、復元コストは履歴長に比例、検証中は2ワーカー分の資源が必要です。SNAPSHOTS.mdを参照。

## TDDFTのディスク再開

TDDFTState.save_checkpoint(directory)とrestore_checkpoint(directory, max_bytes=None)を使います。
新規の保存先のみ、Zarr v3の無圧縮little-endian履歴とSHA256 manifestです。GS・擬ポテンシャル本体は含みません。
利用者が同じ初期バンドル・実行ファイル・dt・初期A・capacity・threadsを用意し、別Pythonプロセスから再開できます。
全履歴再生の費用は残り、物質電流を内部保存し公開値は電荷電流です。既定上限は受信モデルの最大ステップ数×48バイトと形式の管理量から算出し、max_bytesは利用者の追加予算です。全RSS上限ではありません。
破損・不一致・再生失敗では現状態を保持します。場やMPIを含む一括再開は未対応です。SNAPSHOTS.md/enを参照してください。

## 自作物質モデル

materials.check_material(initial, drive)で利用者指定A波形について独立性・再現性・有限値・端点時刻・入力所有権を検査します。
materials.bind_material(initial, key=..., identity=...)は同じnew_state/step(A)モデルをMaterialBindingへ変換します。明示TimeSpecが必要です。
MPIは1ランク/物質群の計算点分散であり、単一の直列モデルの内部並列化ではありません。identityは作者が版・全設定・初期データを識別して指定します。
宣言されたJ_matterのみjへ変換し、物理精度やfailure atomicityを証明しません。本番でチェッカーを暗黙実行しません。
Pythonから直接呼ぶTDDFTは既存のmanaged provider、LDは既存E駆動adapterを使います。MATERIALS.md/enとexamples/custom_material.pyを参照してください。

## 自作計算処理物質のC ABIプレビュー1

c_material.CMaterial(library, sha256, configuration_bytes)は公開ヘッダinclude/salmonpy_material.hのABI1を読みます。
モデル/状態は独立int64 handles、固定dtの端点Ac→電荷電流j、明示TimeSpec、statusコード、失敗時原子性が仕様です。
CMaterialは共通respond/check_material/bind_materialで使え、モデルclose後も既存状態が閉じるまで計算処理モデル解放を延期します。
ライブラリは信頼できる実行コードで、SHA照合は障害隔離ではありません。ロード済みファイルは置換せず版ごとに別ファイル、全計算処理呼出しはprocess-wide直列化、Pythonモデル/状態の同時呼出しは未対応です。
C++/Fortran呼び出し元プログラム例は明示build、ソース/全object監査を保持。Cプラグインの製造解検証であり、LD/SBE永続handle化や複数Pythonから直接呼ぶTDDFT contextを実装済みと説明しません。
C_ABI.md/enを参照。MPI collective C ABI、波動関数checkpoint、GPUは別作業です。

## SBE永続ハンドル

SBEInitialStateのモデル配列を一度計算処理へ準備し、状態は独立したaccepted/trial密度を持ちます。
通常stepはA/J3成分だけ渡し、density/snapshot時だけ全密度をexportします。close/context managerと補助GC解放に対応。
初期モデルはdata-onlyのasdict/pickle/deepcopyを維持し、nativecacheは外部weak-key登録で移送しません。
shared velocity kernels/演算順序/既存OMP外側k分割は不変。常駐メモリはモデルnativecopyと各状態のrho2枚を必要とします。
明示サイズoverflow/allocateエラーは拒否しますが、未変更kernelの自動work領域OOMを回復可能とは説明しません。
既存固定nativebinaryを上書きせず、新規build/installで検証できます。SBE.md/enを参照。LD/TDDFTのhandle移行を実装済みとしません。

## LDの状態保持

LD.md／LD.en.mdに従います。係数と極ごとの履歴は状態ごとに保持し、stepではEと総P/Jだけを受け渡します。close／withで解放し、respondとMPI応答サービスも所有状態を解放します。固定のメモリ上限はありません。snapshot／restoreは同一モデル・固定dtでのみ利用し、状態そのもののpickleは案内しません。共有LD計算ルーチンと物理的P/Jの規約は不変です。

TDDFTの保存用GSは状態の参照で保持し、初期モデルclose後も既存状態の再生に使えます。初期モデルclose後のnew_stateは拒否し、最後の状態解放後に保存用も削除します。N状態でN+1組、restore検証中は追加の作業用コピーが必要です。論理コピー量と実ファイルシステム上の物理占有量は区別します。

## GS/2配列アーカイブ

`salmonpy.gs_archive.save_ground_state`／`GroundStateArchive.open` は明示的に渡した周期・非スピン分極のGS配列のZarr v3入出力です。GS_ARCHIVE.md／GS_ARCHIVE.en.mdに従い、k,band,z,y,xのセル周期u、reduced逆格子k、原点0・端点なしの格子を使います。SALMONpy DFT生成のΓ点GS/1はGroundStateResult.export_archiveでGS/2へ書き出せます。SCF収束の再認定は行いません。GS/2からのTDDFTはΓ点・直列GNU計算プログラムに固定した監査済みGammaRestartConverterを明示指定した場合だけ利用できます。擬ポテンシャルは参照だけで自動取得しません。Fortran stream helperはチャンクの実数I/Oだけで、完全なSALMONアダプターではありません。

Γ点GS変換ではgs_native.GammaRestartConverterを使い、対象実行ファイルSHA256・実コンパイラ/Fortran_FLAGS・擬ポテンシャルhash・GFORTRAN_CONVERT_UNITを照合します。MPI／多k／別ABIの汎用変換とは説明しません。合成軌道での実計算プログラム接続試験と物理的なGS検証を区別します。

GS書き出しは新しいexporter付きGammaRestartConverterを明示指定します。実数/複素数の直列Γ点のみ、原子単位のgs_eigen.dataと計算処理占有数の一致を必須にします。密度は軌道と占有数から再構成し、混合履歴の保存と呼びません。provenanceに由来・元GS等のhash・局所チャネルを記録し、再入力時の局所チャネルも検査します。実GSを収集せず合成データで接続試験することと物理検証を区別します。

複数k点は新しいSerialRestartConverterを明示指定し、座標軸に沿う直交セル・直列・非スピン分極のみです。GammaRestartConverterはΓ限定を維持します。元GS書き出しは検証済み通常kmeshとgs_k.dataを照合し、再入力はhash付きkpoints.datをfile_kwへ渡します。密度・電子数はk点重みを含めます。波動関数/占有数と、k点テキスト・重み再規格化による電流丸め差を区別し、一般メッシュの電流比較はrtol1e-12/atol1e-14を使います。

## 本体ソースの管理（2026-11-15）

本体ソースと変更済みの完成形はruntime/sources/solverに統合しています。外部SALMONのcloneや互換パッチ適用を案内せず、`python -m salmonpy.runtime build NEW_BUILD --jobs N`を明示的ビルド入口にします。HSEは高速PR版のみ。標準版でUSE_HSEを有効にして代替しません。元コミット・ライセンス・全ファイルのハッシュを保持し、import時はビルドしません。旧patch_sha256は由来と過去記録の互換用識別子で、新しいソースの認証ではありません。新規数値変更はsource-tree/replacementsのハッシュと対応する検証を更新します。

## 高速HSEの本体接続（開発用）

RUNTIME.md/enの計算処理 owner経路に従います。`--explicit-routing`で準備した高速版は、
明示的な検証時に`SALMONPY_HSE_CONTEXT=1`でTDDFT本体へ独立HSE ownerを結合します。
FFT・ACE・履歴の交互実行検査と、本体電流/エネルギーの32ケース比較を区別します。
共通かつ不変のglobal設定、直列・非スピン分極・全支持の範囲です。TDDFT全体の複数状態・
異なるglobal XC/dt・MPI/Wannierの完全な独立性や、Pythonからの新規HSE APIを案内しません。
owner/結合済みsystemはコピーせず、unbind後clearで解放します。資源付き状態のbytes=-1・
status=2は計測不可です。GitHub Actionsは利用枠節約のため、新規起動せずローカル検証を優先します。

## TDDFT全体の独立状態の開発用検査

RUNTIME.md/enの「状態ごとの固定dtを持つTDDFT全体の独立状態」に従います。tddft_context.f90は
計算処理の全step引数を所有し、dt以外は同じ設定・同じGSから独立初期化した2状態の交互実行を検査します。
異なるA波形、逆順、片方の破棄/再生成、別プロセス参照との配列比較、本体Acos2参照を使います。
rt_stepの返すjは電荷電流、エネルギーは毎ステップ、単位は原子単位です。
内部Fortranの開発用入口であり、既存Python TDDFTの制約を解除したとは説明しません。
XC作業配列・modulo表は各状態へ移しましたが、global入力・初期化設定の共有は残ります。異なるXC/grid/k、MPI、同時呼出しは
未対応です。保存した基本設定の変更や非有限Aは進める前に拒否しますが、全計算処理失敗の原子性や
RSS上限は保証しません。traceは合成検査専用で、checkpointではありません。

内部Fortranのrt_open(..., delta_t=...)では状態ごとの固定dtを指定できます。有限・正の値のみ、
省略時は入力dtを取得します。異なるdtの別プロセス参照と旧global dtを変更した交互実行を比較します。
適応刻み・途中のdt変更・FDTDサブサイクル・公開Python APIの拡張とは説明しません。


XC作業配列はsystem%xc_work、周期境界表・最適化可否はstencilが所有します。
calc_xcのworkspace省略時は呼び出し内で確保します。下位の異形状PZ/配列検査を、異なる格子/XCの
TDDFT全体の対応と読み替えないでください。共通格子/k点/XCの保護は維持しています。
AVX/SVE最適化経路・MPI・GPUはこのローカル検証の対象外です。


OpenMPのステンシル/電流/力の分割フラグは状態ごとに固定します。内部rt_openのopenmp_modeは
省略時に従来入力、指定時はauto/orbital/rgridです。プロセスのスレッド数設定・MPIランク分割・
同時実行対応とは別です。異なる分割方法の交互実行試験を、異なる格子/XCや力の検証と説明しません。


内部Fortranの作成時設定はs_rt_run_configです。rt_input_configで従来入力から値を取得し、
rt_open(config=...)へ渡します。dt・最大ステップ数・3種類のOpenMP分割指定をコピーして固定します。
configと旧delta_t/openmp_modeの混在は拒否します。容量の整数検査をメモリ保証と扱わず、
格子/XCなど入力全体が独立化されたとは説明しません。上限到達後のstepは更新前に拒否します。


s_rgrid_configとinit_grid_whole_configは、実空間格子の境界・halo・添字配列のみを初期化する内部入口です。
座標配列やMPI分割表は後続処理で扱います。明示的な設定での異形状・孤立境界の検査を、
異なる格子/XCや孤立系TDDFT全体の対応と読み替えないでください。不正値は既存格子を変更せず返します。


resolve_cellはs_cell_configからセル形状と格子間隔を決める内部ルーチンです。結果は値として独立し、
不正値・体積ゼロでは変更しません。計算処理互換層のglobal al正規化とk/pseudo/XCの依存は残るため、
セル設定全体の独立化や異なるセルのTDDFT対応済みとは案内しません。


内部Fortranのk点・radial擬ポテンシャル・XC入力はs_kpoint_config / s_pp_input_config / s_xc_configへ分離しました。
RUNTIME.md/enの「入力の所有」に従い、値設定はコピー、初期化済みLibxc/HSEハンドルはコピー禁止です。
calc_xcのpre_scf_active省略は本計算を意味します。互換ラッパーは既存SCF段階を明示的に渡します。
下位部品の異設定・共有入力変更テストを、異物質TDDFT全体・MPI・同時呼出し対応と説明しないでください。
SOC/+U専用処理・SCF/アルゴリズム・MPI/ファイル単位の共有は残り、共通grid/k/XCガードを維持します。
検査範囲はdocs/shared-input-validation-2026-11-15.mdで確認し、未検証PP形式や物理収束を保証しません。


スカラー擬ポテンシャルの格子配置はs_pp_grid_configをppgが保持し、原子数・原子種はsystemから取得します。
init_ps(config=...)の独立経路はscalar/+Uなし。作成後の設定・原子数・種対応変更は拒否し、位置変更だけ従来順序で更新します。
SOC/+U専用状態、MPI資源、タイマー、異物質TDDFT全体の共有解消済みとは案内しません。
検証はdocs/pp-grid-validation-2026-11-15.mdとtools/test_pp_grid_config.pyを参照してください。

PP格子を更新する間、読み込んだradial数値表と格子形状は不変です。表内容は毎stepで再hashせず、設定・構成のみ検査します。


内部HSEのexx_set_phase(system, pre_scf)とexx_pre_scf(system)で評価フェーズを明示所有できます。
RT初期化は本計算フェーズを所有し、共有pre-SCFフラグ変更の影響を受けません。
明示的にpre-SCFから戻す際は現軌道でexx_refresh後に交換作用を使います。
GS/DCのpre-SCF閾値・連続回数・必要回数・開始反復も各systemが所有し、集団通信と切替・再構築順は維持します。
DC通信・物質情報交換、適応交換制御、従来交換キャッシュは共有のため、全SCF制御の独立化・並行実行対応とは案内しません。
小さいGSの4条件比較とDC/MPIの8条件比較は切替・数値の回帰検査で、物理収束や複数ソルバーの同時独立実行の保証ではありません。
DC比較はMPI有効ビルド記録を渡し、tools/test_pre_scf_runtime.pyに--dc-ranks 2または4を指定します。
空間分割HSEではwrite_gs_restart_data='no'が必要です。MacではMPI CラッパーがApple Clangを選ぶ場合があり、
GNUの-fopenmpを使う構成ではCコンパイラもGNUに合わせます。MPIソケット起動の制限は計算失敗と区別します。
検査はtools/test_pre_scf_state.pyとtools/test_pre_scf_runtime.py、記録はdocs/pre-scf-validation-2026-11-15.mdを参照します。

評価フェーズの変更はbind済み独立HSE ownerの資源・履歴を無効化します。共有交換キャッシュの転用はしません。

GSの適応交換準備・交換固定・支持領域変更フラグと局在化検査はsystemに結び付いたHSE ownerへ振り分けます。
未結合なら従来の共有ownerを使います。完全なGS独立化・残る設定の分離・同時呼び出し対応とは案内しません。
内部APIはexx_set_scf_flagsと3つの問い合わせ関数、exx_check_localization(system=system)です。

通常GSの限定構成ではHSE資源をGS内部で生成・結合し、結果/checkpoint出力後に解除・解放します。
判定はexx_bind_gs_systemで、status=0は結合成功、2は対象外で従来owner維持、1は結合衝突です。
対象は構造最適化なしの通常DFT、周期・スカラーHSE06、全支持領域CPU交換、Wannier/対称性縮約なし、
整数占有（temperature<0）、各電子分割1ランクです。DC/有限温度/空間分割は対象外です。
ログEXX_GS_OWNER explicitで選択を確認できます。GS全体の独立化や並行実行対応とは説明しません。

GSの同一入力・同一プロセス2回連続実行と元コードの別プロセス実行は、
tools/test_gs_lifecycle.pyでスカラー結果と全生成restartファイルを比較します。
4条件で一致していますが、異なる入力の併存、同時実行、TDDFTへのrestart読込の検証とは区別してください。
比較後にrestart実体は削除し、ハッシュを保持します。

内部GSコンテキストsalmonpy_gs_contextは、実引数をTARGET付きで宣言して安定したアドレスを保ちます。
生存中のs_gs_contextをコピーせず、破棄前にgs_closeを呼びます。gs_runの再実行はstatus=1で拒否、正常終了は0です。
入力・並列・タイマーなどの共有は残るため、明示入力のGS API完成や並行実行対応とは案内しません。
複数状態保持の検証は同一入力のスカラーHSEに限定します。

内部GSの収束指標を含む5つのSCF設定はgs_input_configでs_gs_run_configに取り込み、gs_run(config=value)へ渡します。
max_iterationsは正かつ整数最大値未満、threshold/pre_scf_thresholdは有限・非負、pre_scf_stepsは正です。
不正値は初期化前にstatus=1、正常設定はコピーして使います。省略時は共有入力を一度取り込みます。
convergence_metricはrho_dne/norm_rho/norm_rho_dng/norm_pot/norm_pot_dng。指標と閾値の内部単位を一致させます。
GS混合はs_gs_run_config%mixingにmethod/mixrate/alpha_mb/beta_p/nmemory_p/auto_mixing/update_ratio/initial_no_mix_stepsを保持します。
初期化前にコピーし、適応係数と履歴は各状態で独立です。係数は有限、nmemory_pは1〜21、無混合回数は0以上・整数最大値未満。
methodはsimple/broyden/pulay/simple_potential/simple_dm、auto_mixingはy/n。simple_dmにはspin-orbitが必要です。
内部Fortranのみ。占有数用電子数・構造・物質入力の分離完了、DFT-MD設定の独立化や並行GS対応とは説明しません。
同じ物質で異なる混合条件の二状態を元SALMONと比較します。物理収束の保証ではありません。

明示pre-SCF設定も入力読込と同じ適用条件を再検査します。restart/MD/最適化、checkpoint指定や密度以外の収束などを
設定構造体で迂回しないでください。不正な実行設定は内部status=1、従来実行入口では非ゼロ終了になります。

## 物質IDと設定カタログ

`material_config`タスクとCONFIG.md/enに従います。実装済みは`salmonpy.materials/1`の
設定検証・不変化だけ。明示登録したPython resolverを使い、JSONからimportや計算起動をしません。
初期状態は参照IDとsource/prepareの形だけを検証し、読み込み・物理整合の検証済みとは扱いません。
構造と幅の独立アダプターは後述のstructure_gridを参照。全計算の`salmonpy.simulation/1`、MPI配置、
出力/restartへのID継承は後続です。スーパーセルは別の構造ファイルを参照する設計です。
MateriAppsを構造取得APIと呼ばず、実際の提供元と出典を記録します。数値例を物性データと扱いません。

## オフライン計算ノードへの案内

計算ノードのインターネット接続を前提にしません。未確認なら接続可否と準備・転送先を確認し、
既知の施設制約は再質問しません。RUNTIME.md/enの手順で取得と実行を分け、計算先に合う
wheelと依存、ローカルカタログと計算プログラムZIP、構造/スーパーセル・擬ポテンシャル・入力・
初期状態と全参照データを事前準備します。版・出典・hashを記録し、転送後の照合と実計算ノードでの
小規模確認を案内します。実行時のDB検索・URL取得、暗黙ダウンロードや自動ビルドを提案しません。
不足物を列挙して準備段階へ戻します。ログインノードの通信可否も推測しません。
オフライン一式の自動収集・依存解決、ユーザーcallbackの通信遮断は未実装です。

明示的なローカルファイル収集は`offline_inputs`タスクの`salmonpy.preparation.prepare_inputs`、
転送後の照合は`verify_inputs`です。PREPARATION.md/enに従い、登録済みresolverと
コピー先相対名→元ファイルの一覧を指定します。汎用収集器が自動で依存探索・構造解釈・格子幅解決・科学的妥当性検証を行うとは説明しません。検証器はbackendを呼ばず通信しません。
manifestのhashを別経路で保存し、出力は入力一式の外へ置きます。

## 構造ファイルと格子幅

`structure_grid`タスクとSTRUCTURE.md/enを参照。StructureInputはASE Atomsまたは
明示形式cif/extxyz・非負frame番号でローカル構造を読み、独立コピーと出典hashを保持します。
CIF対称性展開はASE読込に含まれますが、追加スーパーセル化・標準化・折り返しはしません。
grid_from_spacingは正の軸方向直交セルに限定し、幅の単位はangstrom/bohrを明示。
ceilと最低点数/明示倍数で切り上げ、実際の幅と点数を保存します。MPI数から変更しません。
ASEのdft_from_atoms/tddft_from_atomsはgridまたはgrid_spacing＋grid_spacing_unitを選択し、
両指定を拒否。直接DFT/Calculatorは従来どおり点数を渡します。格子生成を収束・資源保証と扱いません。
examples.prepare_structureが既存オフライン収集APIへの明示登録例です。自動backend登録、
擬ポテンシャル選択、実DFTの準備完了とは説明しません。混合・部分占有は代表元素で計算せず拒否します。直接ASE幅指定は採用点数のみ保持し、希望値等を保存する場合はヘルパーのレポートを記録します。

## 組み込み物質の設定接続

`builtin_materials`タスクとBUILTIN_MATERIALS.md/enに従います。builtin_material_registryを
明示作成し、resolve/load_materials→material_from_catalog→new_stateへ接続します。
LD/SBEに加えDFT設定接続とGS→TDDFTの橋渡しは下記を参照。LDはE、SBEはA/c駆動で、符号や更新順は既存と同じです。
SBEは既存の3種類のSALMONテキストをhash付きで読み、変更時は再解決を要求します。
モデル/状態にmaterial_idとdefinition_digestを保持しますが、point/run ID、MPI割当、
出力/restartの自動ID継承やinitial_statesレシピの選択は未対応。record()を結果と明示保存します。
状態作成で導入済み計算処理が必要です。登録/解決/モデル作成に計算処理起動・通信・コンパイルはありません。

## DFT物質設定からTDDFTへ

`dft_material_config`タスクとDFT_CONFIG.md/enに従います。salmonpy.dft/dftを明示登録し、
ローカル構造・PP・XC・states・kmesh・点数/幅から既存DFTを作成します。UPF元素は原子番号と照合し、
KY元素やPPのXC由来まで自動認定したとは説明しません。固定入力と実行時SCF/threads/runtimeは分離します。
ConfiguredMaterial.run→ConfiguredGroundState.tddft_initial_state→new_stateでIDと定義hashを継承。
既存GS検証とtransverse・A/c→j、独立workerを維持します。warm startは結果ラッパーを展開して既存検査を通します。
任意GSのJSON選択・point/run ID・ディスクID自動保存・MPI接続は後続です。configured_siは短い接続検査です。

## 保存GSをIDで選ぶ

`saved_gs_selection`とINITIAL_STATES.md/enに従い、tddft_from_catalogを明示呼出します。
GS/1のsource相対パスとground-state.jsonのSHA256、DFT material_idを指定。全計算処理ファイルと
構造/PP/ライセンスを管理ノードで収集・転送し、verify_inputsで照合します。GS生成時と同じ
実行ファイルhashを要求し、機種間restart互換性やSCF再実行は暗黙に提供しません。
material/initial_state IDと各digestはモデル・状態に保持し、record()を結果に保存します。
既存の「initial_states未対応」は汎用factoryと任意レシピに関する制約です。この入口は
完了したSALMONpy GS/1のみ対応し、prepare、GS/2、point/run IDやMPI自動割当は対象外です。
サンプルのAはA/c、返す電流は電荷電流密度j = -J_SALMONです。収束・メモリ保証とは区別します。

## 直列の点・run対応表

`point_runs`とPOINT_RUNS.md/enを参照。PointRunに順序付きpoint_id→準備済みモデルを明示し、
一つのrun/pointに一状態だけ作成します。物質と初期状態のID/hashを保持し、save_arraysで
配列のpoint_axesを明示して保存、open_point_arraysで対応表を検査します。LDはE、結晶はA/c。
自作モデルにもrecord/new_stateで接続できますが、隠れた状態共有・実行能力は認定しません。
既存ガイドのpoint/run未対応は従来API・MPI・restart自動継承の制約です。この入口は直列出力
に限定し、時刻同期、点の自動並べ替え、checkpoint再開やMPI割当を提供しません。

## MPIの利用者ID

`mpi_point_identity`とMPI_IDENTITY.md/enに従い、ResponseSessionまたは
CoupledStepper.from_multiscaleへ同一point_record/point_mapを全ランクから明示指定します。
整数通信番号は既存のまま、利用者の文字列point_idへ明示対応します。run_idをランクごとに生成しません。
prepare/create前に全ランクの一致を検査し、自群のSolverContext.point_recordsでIDを参照できます。
session.save_arraysは集団呼出しで、引数は場leaderだけです。open_point_arraysで対応表も検証します。
metadataは宣言値の整合性検査でありproviderの実体・物理・restart互換性の認定ではありません。
従来ガイドのMPI ID未対応は省略時の自動継承やrestartについての制約で、この入口は明示的な継承です。

## 記録付きモデルのMPI binding

`recorded_material_binding`とRECORDED_MATERIALS.md/enを参照。設定済みSBE/自作端点Aモデルを
bind_recorded_materialで接続し、response_identityでIDを手で写さず点記録を生成します。
物質群は1ランク。全ランクで同じrun_id/point_mapを指定し、モデルをセッションより長く保持します。
recordは準備と状態生成で再照合しますが、虚偽のhash宣言は検出しません。LDのE変換は明示の専用経路、
TDDFTは管理計算プログラム専用providerを使います。設定済みTDDFTの専用factoryは後続です。

## 設定済みTDDFT/LD専用binding

`configured_native_binding`とCONFIGURED_NATIVE.md/enを参照。bind_configured_tddftは
管理runtimeを明示しsnapshotを借用します。binding/prepareではコピーもworker起動もしません。
bind_configured_ld_endpointsは既存のA端点差分→E変換を明示的に選ぶ入口です。
response_identityでIDを継承し、物質群は各1ランク、TDDFTは独立子worker。準備モデルはセッション解放後にclose。
群内計算処理 MPIや半陰的E/H結合・資源上限・物理収束の保証へ拡張解釈しません。

## 点ID付きTDDFT checkpoint

`point_checkpoint`とPOINT_CHECKPOINTS.md/enを参照。PointStateのsave_checkpoint/restore_checkpointは
既存TDDFT履歴に厳密なrun/point/material/initial ID照合を重ねます。再開時は同じrun_idを明示。
別点/runへの分岐やMPI全体の再開は未対応。初期GS/runtimeは別途検証して用意し、archiveから自動取得しません。
復元は履歴再計算で最大2worker、max_bytesは内部archive予算でRSS保証ではありません。
LD/SBEや自作backendのディスクserializerとして案内しません。生TDDFT archiveはPointStateへ直接渡せません。

## 共通A/E/Hローカル入力

`field_input`とFIELD_INPUT.md/enを参照。FieldPacketは任意A/E/HをFieldValueで持ち、
Noneとゼロを区別。TimeSpec・単位・ゲージ・位置をFieldContractで宣言し、FieldAdapter.step_fieldsへ渡します。
AはA/c、HはBではありません。結晶はAのまま、LDはE。長さゲージはE必須ですが、計算処理孤立系の接続完了ではありません。
既存stepへの橋はlegacy_quantityを明示し、A/c:auまたはE:auの端点値だけを渡します。
自作step_fieldsには宣言成分だけ渡し、半ステップや単位・位置を検査します。補間・微分・単位変換はしません。
MPIのA/E/H転送は後述の明示的なpacket経路を使います。磁気物理の対応はbackendごとに確認します。状態は借用、backend例外でadapterは無効です。

## 孤立系E-workerの現在の範囲

FIELD_INPUT.mdの検証用E3プロトコルは、ソースから明示的にビルドしたstandard runtimeで
PZ・固定イオン・中点伝播の双極子を比較する開発用接続です。結晶のTDDFTInitialStateにEを
渡したり、検証用driverを公開初期状態ローダーとして案内したりしないでください。
公開入口はISOLATED_TDDFT.mdのload_isolated_runtimeとIsolatedTDDFTInitialStateです。
revision5のローカル検証記録に結び付いたバイナリに限定します。savepoint系は同じプロセスへの直接退避、snapshot/restoreは受理履歴の保存・再実行です。
既存配布runtimeの自動対応、電流密度への暗黙変換、孤立系ローダーのMPI応答サービスへの直結は未対応です。孤立系のsave_checkpoint/restore_checkpointは専用E/dipoleスキーマで履歴再実行を行います。GSとruntimeは別途指定し、同一容量・dt・位置・版などを照合します。周期系形式と混同しないでください。

分子双極子からFDTDへの一方向接続はDIPOLE_RADIATION.mdを参照してください。電荷双極子、
FDTDセル体積、配向・個数を明示し、最初の区間電流をゼロにする保持時間を隠さないでください。
双方向の自己場帰還や静的双極子の静電場を計算済みとは案内しないでください。

## MPI場packet
ResponseSessionのfield_contractsとstep_fieldsで明示的なA/E/Hを転送できます。
FIELD_INPUT.md/enとmpi_field_inputs.pyを参照。戻り値は電荷電流密度jです。
孤立系の双極子を自動で規格化したり、場の時空間補間を推測しないでください。
既存step(A)との混在・MPI全体のrollback・計算処理結合driverへの自動接続は未対応です。

場の補間にはfield_sampling.SampledField/sample_packetを使えます。
利用者が選ぶ空間の凸重みとexact/linear時間方針を明示し、範囲外へ外挿しません。
成分の時計・単位を保ち、Hを推定で補いません。計算処理結合は後述の明示された時刻配置に限定します。

NativeEHField/MultiscaleFieldのsample_fields()で1ステップ後からA/E標本を取得できます。
計算処理ではA=(n+1/2)dt、E=n*dt、1DではA=n*dt、E=(n-1)dtです。
Hをゼロで補わないこと。sample_fields自体は履歴を蓄積せず、結合更新もしません。

H対応C/Fortranの呼び出し規約ではsample_magnetic_field()を全fieldランクから明示的に呼べます。
raw Hは(n+1/2)dt、保存用の時間平均Hはn*dtです。両者を混同しないでください。
旧ABIのH欠落を0で補わず、対応runtimeを明示ビルドするかHを必要としない経路を選びます。

PacketCoupledStepperは宣言済みの直接標本に限りpacket往復結合を提供します。
NativeEHではA/Hとjがoffset1/2、Eが0。1DではA/jが0、Eがlag1です。
t0/clock_offsetは0、位置・単位は厳密一致。一般補間や既存電子providerの自動変換はしません。
FIELD_INPUTの表とmpi_packet_coupling例を参照。双極子をjとみなさないでください。

既存A駆動bindingのpacket化はpacketize_A_bindingを明示的に呼びます。
field_phaseはendpoint（1D）またはyee_half_step（FDTDのE/H更新）を必須指定。
元の端点A/j宣言、ゼロ初期Aを要求し、電荷電流の符号を再反転しません。
絶対時刻依存の自作モデルへ機械的に適用せず、時刻対応を確認します。
response.timing_recordsに元/変換後の宣言が残ります。FIELD_INPUT.md/enと
mpi_packet_couplingの--adapt-a、--runtime/--bundleを参照してください。

packetize_A_bindingはleader_only・一物質群1ランクに限定し、collective物質backendを拒否します。

分子への電場選択はdipole_source.select_molecular_fieldで方針を明示します。
provided_driveまたはenvironment_scattered_onlyを使い、対象自身だけの真空参照と
配置・観測・源履歴の一致を宣言します。宣言は計算条件の独立検証ではありません。
DIPOLE_RADIATION.md/enとmolecular_field_selection例を参照してください。
純粋な標本選択であり、native/孤立系TDDFTの双方向結合・放射反作用の完成とは説明しません。

分子から電磁場への応答の受け渡しは開発用振動子で順序検証済みですが、保持したEによる
結合仕事の差が1次で減るだけで、結合全体の保存則は未達です。振動子単体の約2次精度と
混同せず、孤立系TDDFTの双方向対応を宣言しません。新しい検証スキルの一般判定も、
SALMON固有の更新順・固定dt・資源予算やユーザー所有の物理収束判断を上書きしません。

内部GSの収束規格化はs_gs_run_config%normalization（nelec/nelec_spin）に取り込み、状態へコピーします。
非分極はnelec、分極は正のnelec_spin和を優先し、DCは状態内elec_num_totを優先します。
非DCのrho_dneの分母0以下は初期化前にstatus=1。通常は取り込んだ値を変更しません。
これは電子数・電荷・占有数を変更するPython APIではありません。物質入力の全面独立化とは案内しません。
二倍の分母で残差だけ半分にする診断と、通常条件で元SALMONとの一致確認を区別します。

### 状態が保持する占有数入力

内部`s_dft_system%occupation_input`はnelec、nelec_spin、nstate、spin、spinorbit、temperature、redistribution_delayを初期化時にコピーします。初期充填・通常ne2mu・通常GS再分配条件・ギャップの占有バンド番号はその値を使います。未初期化systemには従来入力の互換経路があります。温度は内部原子単位です。

生存中の値は変更しません。gs_runの新しい占有数上書きAPIではありません。DC専用ソルバー・並列配置やHSEの他の設定等には共有入力が残り、任意の電子数・温度の同時GS実行を保証しません。

検証: `tools/test_occupation_inputs.py BUILD_RECORD NEW_OUTPUT`、`tools/test_gs_occupation_runtime.py BUILD_RECORD REFERENCE_RECORD NEW_OUTPUT`。前回の未コミットファイルは再起動後に見つからず、会話から再構成しました。復旧後の証拠は`docs/gs-occupation-recovery-validation-2026-11-16.md`です。

### HSEの状態別占有数入力

HSE交換の電子数・温度・状態数（nelec/temperature/nstate）は、同じ`occupation_input`を参照します。初期化前の従来経路は維持し、演算式・更新順は変更しません。

独立性の検証範囲は、直列CPU・周期境界・非スピン分極・整数占有・全支持の高速HSEです。電子数4/2の二状態について、共有入力を書き換えた単独・交互実行、解放・再作成、作用と交換エネルギーを検査します。有限温度GSは従来版との一致を検査しますが、有限温度/Wannierの複数状態独立性は未検証です。nstateのWannier分岐も新たな独立性検証の対象ではありません。

外側のTDDFT設定一致ガードは維持します。HSE以外の共有入力が残るため、電子数や温度が異なる任意のTDDFT状態を同時に作れる新APIとは案内しません。生存中の占有数設定は変更しません。
検査は`tools/test_hse_explicit_state.py BUILD_RECORD NEW_OUTPUT --native-resources`（必要なら`--threads 2`）。検証記録は`docs/hse-occupation-inputs-validation-2026-11-16.md`です。

### TDDFTの電子数診断

`subdip`の電子数チェックも、状態の`occupation_input.nelec/nelec_spin`を使います。2スピンで正のスピン電子数合計があればそれを優先し、その他はnelecを使う従来規則を維持します。上下限の厳密な不等号、出力、停止処理は変更しません。未初期化状態は従来入力を参照します。
これは広い範囲の異常検知（指定電子数の0.1倍から10倍の間）であり、精密な電子数保存や物理的収束を保証する検査ではありません。共有I/O・イオン情報や外側TDDFTガードは残ります。
検査: `tools/test_occupation_inputs.py BUILD_RECORD NEW_OUTPUT --dip`。検証範囲は`docs/dip-occupation-inputs-validation-2026-11-16.md`を参照してください。

## LD・SBE・TDDFTの混在に向けた入口

`mpi_mixed_multiscale`は、同じ1次元場にLDと合成SBEを配置する接続サンプルです。
```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 mpiexec -n 3 python -m salmonpy.examples.mpi_mixed_multiscale NEW_OUTPUT
```
場1ランク、物質点2ランクで動かし、非ゼロの応答と独立した直列更新との一致を検査します。点ID10/20と物質キー、A/c・電荷電流j、通信時間を保存します。LDは端点A差分からEを得ます。FDTDのE/H更新の半陰的LD更新と同じ方式ではありません。

4ランクにして`--runtime RUNTIME --bundle GS_BUNDLE --input rt.inp`を指定すると、管理済み直列TDDFTワーカーを第三点に接続する入口を使えます。Si入門系の実GS（PZ、16³格子、2×2×2 k点）を用いた4ランク・100ステップの三種類混在実行は、直列計算と完全一致しました。実GSと対応計算プログラムを明示的に用意してください。これは管理済み直列ワーカーの検証で、TDDFT群内部のMPI並列とは別の経路です。光学物性の収束は保証しません。記録: `docs/mixed-three-material-validation-2026-11-16.md`。

最終構成は、FDTDランク内の多数のLD状態と、別MPI群のSBE/TDDFT状態の混在です。この小規模サンプルの「LD一セルに一ランク」を本番の配置方針にしないでください。実材料の収束や三次元混在対応を保証するサンプルではありません。

## FDTDのE/H更新内のLDと、SBE・内部MPI TDDFT

`mpi_mixed_native`は、場1ランクが二種類のLD領域を持ち、別の1ランクがSBE、別の2ランクが1点のTDDFTを担当する接続試験です。外側4ランクでTDDFTの軌道分割を内部2ランクにします。
```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 mpiexec -n 4 python -m salmonpy.examples.mpi_mixed_native BUILD_JSON SERIAL_RUNTIME GS_BUNDLE NEW_OUTPUT
```
`--tddft-ranks 1`なら外側3ランクです。`--without-ld`はLD領域を真空にする対照条件です。BUILD_JSONは同じMPI環境で明示ビルドした`tools/build_embedded_tddft.py`の記録、SERIAL_RUNTIMEは応答再計算の比較用です。自動ビルド・ダウンロードは行いません。

LDはFDTDのE/H更新の既存の係数・補助量・更新順を使います。電子状態はA/c駆動を維持し、受け渡しは電荷電流jです。点ID10はSBE、20はTDDFT、LD領域ID1/2は場の材料格子上で管理します。電子点の周辺はLDと重複させません。半ステップずれの場時計と物質の局所時計を区別します。

この40³格子・電子2点の小さな三次元配置は接続検査用です。実用TDDFTは引き続き一次元伝播を推奨します。巨大な多数LD、分散した複数FDTDランク、全計算再開、物理収束の検証とは区別してください。計算処理異常は回復可能ではなく、サンプルはMPIジョブを停止します。

検証済み: 128ステップ、TDDFT内部1/2ランクの場・電流比較、直列応答再計算、LDなし対照、LD単独の従来一致、wheel導入後実行。範囲と誤差は`docs/native-mixed-material-validation-2026-11-16.md`を参照してください。

### FDTD群も並列化する

場をy方向に2分割する例（場2＋SBE1＋TDDFT1＝外側4ランク）:
```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 mpiexec -n 4 python -m salmonpy.examples.mpi_mixed_native BUILD_JSON SERIAL_RUNTIME GS_BUNDLE NEW_OUTPUT --field-grid 1 2 1 --tddft-ranks 1 --compare REFERENCE_OUTPUT
```
REFERENCE_OUTPUTには、同じサンプルを場1ランク・TDDFT1ランク（外側3ランク）で実行した出力を指定します。`--compare`を省略すれば単独実行です。比較は同じステップ数・dt・材料設定・点順を要求し、A/c・E・jの履歴を確認します。必要な外側ランク数は`product(field-grid)+1+tddft-ranks`です。

x/y/z各方向の2分割は、128ステップで場1ランクの結果と完全一致しました。y/z分割では両LD材料が両方の場区画に分かれます。場2ランクとTDDFT内部2ランクを同時に使う5ランク構成は、今回の4コア枠では未検証です。大規模性能・メモリ上限の保証ではありません。記録: `docs/distributed-field-mixed-validation-2026-11-16.md`。

## 混在層の構成入口

配置には `salmonpy.layers.assemble_native_layers` と `examples/mixed_layers.json` を使用する。物理設定を配置JSONへ重複させず、準備済みの物質IDとモデル／MaterialBindingを渡す。Bohr単位、厚み切り上げ、全電子セルの独立点化、必要ランク数と実際の厚みを利用者へ示す。ランク予算はメモリ制限ではない。LDと電子点の禁止された近接を真空で黙って変更しない。構成関数はMPI・GS計算・ダウンロードを開始しない。自作コードも明示的なMaterialBindingで接続する。詳細と実行例はCOUPLED.md／COUPLED.en.md。

## 計算処理混在サンプルの直接保存／復元

`examples/mpi_mixed_native.py` の `--direct-save-at N` と `--direct-resume PATH` を使う。入力・計算処理ライブラリ・GS・分割・全ステップ数は同じに保ち、新しい出力先へ復元する。場/PML/LD/SBE/周期PZ／TB-mBJ／Libxc-PZ・PBE中点TDDFTの内部状態と積算A・直前電流・時計を直接復元する。`restored_step` と `native_steps_this_job` で過去の時間発展を再計算していないことを確認する。参照用の独立シリアル計算は検証のため別途実行される。保存は同じ環境・MPI構成の同梱サンプル用であり、汎用計算処理保存やHSE・未検証XCの直接復元を推定しない。詳細はCOUPLED日英ガイド。

TB-mBJ（自動／固定cval）とLibxc-PBEは、場1＋SBE1＋TDDFT内部軌道MPI2の4ランクでも、直接保存63→残り65ステップの再開が連続128ステップとビット単位で一致した。分割を同じに保つこと。場2＋SBE1＋TDDFT2は5ランク必要で、4コア枠では未検証。

HSEの実行エンジンは高速PR版のみ。通常版HSEや比較用ベースラインを利用者向けエンジンに追加しない。高速HSEの直接保存／復元はACE・因子履歴・専用伝播器・計算処理資源を含む対応が必要で、中点伝播のXC許可リストへ名前だけを追加しない。

高速PR版HSEは `native_adapters/tddft_context.f90` の内部serial Fortran入口で直接保存／復元に対応した。`rt_checkpoint(ctx,path,restore,status)` と `rt_accepted_step(ctx)` を使い、新しいプロセスで同じ元GS・入力から初期化した後に内部数値状態を読み、残りだけ進める。ACE・因子／PT-CN履歴とFFT資源の再構築を含む。`tools/test_hse_direct_checkpoint.py` は1／2スレッド、2／3 k点、4伝播設定、6／7分割の32ケースで電流・エネルギー・内部トレースが完全一致した。これは呼び出し元プログラムが管理するMPI/FDTD結合の高速HSE対応や安定したPython/C APIを意味しない。通常HSEは追加しない。詳細はRUNTIME日英ガイド。

高速PR版HSEの外側MPI／FDTD接続は `fast_hse.FastHSEMaterialProvider` と `FastHSEResponseBackend` を使う。専用serialワーカーを明示ビルドし、各HSE点は1ランク群にする。HSE内部MPIではない。場1／2ランク＋SBE1＋HSE1と場内LDを結合し、共有保存先へ直接保存／復元する。`mpi_mixed_native --fast-hse-build ... --layers mixed_hse_layers.json` の手順はCOUPLED日英ガイド。初期電流は実測値を返し、差し引かない。ゼロ場サンプルでは1e-12 au以下を要求。入力の高速伝播器と学習条件を守り、通常HSE・Gammaの未対応学習モードへフォールバックしない。全ランクの保存検査後にHSE LOADを行い、counter/current/Aを照合する。

内部k点MPIの高速PR版HSEは `fast_hse_mpi.NativeFastHSEBackend` を明示ビルドした呼び出し元プログラムが管理するMPIライブラリへ接続する。従来のcreate_stateは物質群1／2ランク、軌道／実空間分割なし、1群1点、メインスレッド・1プロセスにつき初期化1回。明示的なopen_modelの全k点共有状態はSHARED_HSE.mdに従う。場1＋SBE1＋HSE2と場内LDの4ランク例は `--fast-hse-mpi-build ... --hse-ranks 2`。同じ分割・バイナリ・共有保存先でランク別直接保存／復元し、過去ステップの再計算をしない。静的初期化と独立serial参照比較は別途行う。全shardの検査後に反映するが、後段I/O失敗や群間での原子的巻戻しは保証しない。失敗後はジョブを中止する。入力、制約、計測範囲はCOUPLED日英ガイドに従う。

層構成では任意の `backend.supported_group_sizes` を検査する。高速HSEは(1,2)。対応ランク数の宣言は空でない正の整数の重複なしtuple。未宣言の自作の計算コードへ未確認の能力を推定しない。数値・再開の能力検査は別途必要。

対応ランク数の任意宣言はResponseSessionの直接接続でも物質準備前に集合検査する。宣言を全群で照合し、tupleの順序差は同一視する。Aパケットラッパーも宣言を保持する。

直接保存サンプルは計算処理-direct/2の互換情報を全ランクで照合してからLOADする。保存方式・save/restore宣言・所有ランク／点・場と電子の内部分割を比較する。同じ総ランク数でも分割が違えば拒否する。manifest/1は情報不足のため新経路では復元不可、元の対応版を使う。任意のsnapshotを直接保存とみなさず、再生方式を混入させない。自作direct-state/1にはcheckpoint_format／checkpoint_modes／checkpoint_decompositionを明示する。汎用保存APIの完成や原子的復元を意味しない。

保存対応の共通入口は `checkpointing.CheckpointCapability` と `check_checkpoint_compatibility`。自作コードは通常の群mappingと宣言を渡せる。直列はcomm=None、MPIは既存のcommを渡し、全ランクで操作と保存情報を一致させる。ローカル設定エラーはmappingの検査かlocal_errorで集合的に扱う。identityで実装・固定モデル入力・計算プログラムを識別する。APIは宣言・MPI分割・必要ファイル名の照合のみで、ファイル本体のハッシュ検査、候補の寿命検査、I/O・LOADは呼び出し元プログラム側が担当する。CHECKPOINTING日英ガイドとcustom_checkpoint例を参照。旧計算処理-direct/2の互換情報は変更しない。

保存／復元の実行入口は `checkpointing.CheckpointTarget` と `execute_checkpoint`。引数なしのsave／restore関数を渡し、パスや保存形式は計算処理が決める。戻り値はcompatibilityのJSONとプロセス内のfield／statesで、保存トークンや配列をMPI転送・コピーしない。失敗は例外で通知し、Falseなどの戻り値はエラー扱いしない。MPIでは群の点数差を同期スロットで扱い、コールバック失敗後は次のスロットへ進まない。実行済み変更の巻き戻しや計算処理内部の通信停止からの回復は保証しない。完了markerやハッシュ認証は呼び出し元プログラムの方針で管理する。custom_checkpoint例は実際のJSON保存と復元も共通入口から呼ぶ。

元の状態を維持する復元はCheckpointStateとrestore_checkpointを使う。各Targetはstage_restoreで独立候補を読み、全ランクの成功合意後にホルダー一式を切り替える。ホルダー外の別名は旧状態を指すため、計算はlive.field／live.statesを通して再開する。元へLOADする関数をstage_restoreに登録せず、可変の数値バッファを共有しない。初期化済み不変計算処理資源は専用候補仕様に従う。discard_restoreは失敗した準備についてNoneも受け、自群の通信手順を完了する。候補破棄の失敗はAtomicRestoreError.cleanup_errors、元の維持はoriginal_preservedで確認する。既存execute_checkpointは既定で非原子的。保護経路で未対応を検出したら自動で従来経路へ下げない。再ビルドしたPythonから直接呼ぶFDTD／TDDFT／高速PR HSEには実験的な数値候補ABI（mode2準備／3反映／4破棄）がある。内部ヘルパーとMPI試験はCHECKPOINTINGを参照。全群の候補成功後に別途計算処理反映と準備済みPythonメタデータの更新が必要で、ホルダー切替だけでは完了しない。初期化済み不変資源の共有は維持するが数値候補は分離する。複数計算処理 lifetimeの同時実行ではない。通常mpi_mixed_nativeの--direct-resumeは候補対応ABIでの保護復元が既定。場内全LD・SBEの独立ハンドルとPython側の時計・電流・アダプター・所有リスト・framesを反映前に準備する。進行済みセッションの準備失敗でもfailedを立てず、元を維持して全ランクで継続／再試行できる。成功後はsession.statesを参照し、旧SBE解放エラーはcheckpoint_cleanup_errorsへ分ける。旧ABI・プロセス分離HSE・自作direct-onlyは事前拒否。通常サンプルの旧復元オプションは削除済み。共通APIのprotected=Falseは自作バックエンド用の明示的な非原子操作であり、自動fallbackには使わない。古いABIは版確認で拒否し、反映時のメモリ枯渇やプロセス異常は回復保証外。

通常混在のライブ巻戻しは既存SALMONの計算出力ファイルを切り詰めない。保存framesと以後のstep戻り値を使うか、
新規MPIジョブ・出力先で再開する。数値の一致とappend済みファイル全体のバイト一致を混同しない。

## メモリ助言（停止しない）
既知の4GiB・4コア等は聞き直さず再利用する。通常混在では--memory-budget-gibを設定し、MPIランク×OMPスレッドを割当コア数内にする（BLASも1スレッド）。初期化前概算・初期化後最大RSS・復元前の候補同時保持と一時コピーをmemory-advice.jsonで確認する。超過時は警告と超過量を示し、勝手に停止・格子変更・ランク変更しない。強制上限はユーザーが撤回済み。MEMORY日英ガイドとmemory.MemoryEstimate／memory_adviceを参照。未取得部分をゼロにしない。Noneとunknownは不明。低い見積もりを実行可能保証にせず、倍率は計算処理 heap上限の証明ではない。最大RSSは現在値でも子ワーカー合計でもない。資源計測・警告・助言出力の失敗を計算失敗と区別し、failedを立てない。自作の計算コードは明示したバッファ量と根拠を渡し、見積もり失敗でも集合呼出を省略せず不明を渡す。

## 混在サンプルのIDと再開

`python -m salmonpy.examples.mpi_mixed_native` は `--run-id`、論理点ID、物質ID、物質定義・初期状態のSHA256、MPI所有ランク、LDの材料番号対応を保存します。`mixed-example` は試験用の既定ラベルです。実験ごとに明示した `--run-id` を使用してください。定義の識別には出力先や一時ディレクトリを含めません。

`identity.json`、通常実行の `result.json` と `point-response` アーカイブに同じ識別情報を引き継ぎます。直接保存の外側の形式は `salmonpy.native-direct/3`、再計算方式は `salmonpy.native-replay/2` です。数値状態の内部形式と更新順は変更しません。直接復元はIDとMPI分割を数値状態の読み込み前に照合し、保存時は全ランクでのID一致を検査します。

旧版向けの引数・再開形式・欠落項目の補完は提供しません。現在のAPIと保存形式を使用してください。混在計算サンプルは常に物質・初期状態の識別情報を保存します。

## 点ID付きのメモリsnapshot

`PointState.snapshot()` は内部計算処理のsnapshotを `salmonpy.point-snapshot/1` の不変な点情報で包みます。`state.restore(saved)` は実験・点・物質定義・初期状態のIDとパッケージ版を照合してから計算処理を呼びます。別の点や実験、裸の計算処理tokenは拒否されます。明示した同一run/point IDであっても、計算処理のモデル・時間刻み等の整合性検査は引き続き必要です。

これはプロセス内の巻き戻し用であり、分岐・複製やディスク保存、MPI全体のトランザクションではありません。snapshot未対応の自作の計算コードでは `SnapshotError` となります。計算処理の復元失敗時の保証はその実装に従います。点ID不一致の拒否では内部のrestoreを呼びません。裸の状態に対する従来のsnapshot APIは変わりません。

`python -m salmonpy.examples.point_ld /new/output/directory` は別の点への復元拒否と、同一点での復元後の電流一致も確認します。

## Pythonから直接呼ぶTDDFTの更新前エラー拒否（限定対応）

新しく準備・ビルドした標準embedded TDDFTには `sp_tddft_error_protocol() == 1` があり、以下を結合先の物質MPI群だけで照合します。10: 不正・非有限入力／古いA／入力ファイル不在、11: ランク間不一致、12: 状態未作成・使用済み初期化、13: ステップ上限、14: 未確定の復元用の別状態あり。これらは数値更新前に全ランク同じコードで返り、既存の最後に成功した計算状態を変更しません。初期化前の10/11は初期化から終了までの期間を消費せず、Pythonからも入力を修正して初期化を再試行できます。失敗時の作業ディレクトリは診断用に残ります。

Pythonは版付きの仕様と全ランク一致を確認できた拒否だけを `InputError` とします。旧ライブラリ・未知のコード・非有限な戻り値は `RuntimeFailure` のままで、状態保持を推定しません。MPI群の全員が同じ入口へ参加する必要があります。無効なポインタ、MPI・プロセス異常、確保失敗は対象外です。本項の更新前検査ではソルバー深部を扱いません（次節の限定した深部復帰を参照）。初期化に入った後は引き続き1プロセス1初期化から終了までの期間です。高速HSEとPythonから直接呼ぶFDTDのエラー仕様はこの変更では拡張していません。

この保証は個別TDDFT状態の更新前拒否です。結合ドライバで他の部品が既に進んだ場合、全体の巻き戻しを保証せず、従来どおりセッションを失敗扱いにします。通常計算中の一般的な原子性が完成したとは判断しないでください。

## 大きなMacでの並列・格子テスト

4GiB・4コアの入門条件は、その条件での確認用として残します。大きなマシンでも利用者の実メモリ・コア数を確認し、OSや他の作業に余裕を残してください。64GB・18コアのMacでは、今回は48GiBを警告用予算として、最大15 MPIランク、各ランク1スレッドで測定します。これは推奨値の固定やハード制限ではありません。

`mpi_mixed_native` の `--field-shape 96 96 96` はFDTDの格子点数を変更します。各軸は40以上の偶数で、MPI分割に整合させます。既定値は40³のままで、格子間隔も2 bohrのままです。外側の真空領域が広がるため、異なる格子サイズを同一物理条件の比較に混ぜないでください。既定の物質点の物理座標は保ちます。電子状態のGS格子・k点は変更しないため、96³のTDDFT検証を意味しません。3Dは接続検証で、TDDFTの実用途は引き続き1D推奨です。

FDTD2＋SBE1＋TDDFT2は `mpiexec -n 5 ... --field-grid 1 2 1 --tddft-ranks 2`。FDTD12＋SBE1＋TDDFT2は `mpiexec -n 15 ... --field-grid 3 2 2 --field-shape 96 96 96 --tddft-ranks 2`。`...` には既存のPythonモジュール呼び出しとBUILD・RUNTIME・GS・新規出力先を入れます。`OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1` を指定してください。

ソース配布の `tools/test_mixed_scaling.py BUILD RUNTIME GS NEW_OUTPUT --rank-budget 18 --memory-budget-gib 48` は導入済みwheelを用い、同じ96³の入力で4/5/7/11/15ランクを順番に2回ずつ比較します。各計算のA/E/jを最初の計算と照合します。格子サイズが違う結果は比較拒否します。MPI実装固有の引数は `--launcher-arg` で指定し、計測中は他の重い処理を避けてください。

結果の `coupled_loop_seconds_max` は初期化・保存・直列参照計算を除いた時間発展区間です。通信計測には待ち時間が含まれ、純粋な通信帯域ではありません。FDTD内部のハロー通信は個別に分解していません。メモリは各ランクのピークRSSの和であり、同時刻の総常駐量や全システム使用量ではありません。時間・資源の計測失敗を計算失敗と混同しないでください。短い反復測定だけで一般的な並列効率や物理精度を保証しません。

## TDDFT深部の数値異常と直接巻き戻し

新しく準備・ビルドした標準Pythonから直接呼ぶTDDFTは `sp_tddft_recovery_protocol() == 1` を持ち、各ステップの数値状態をメモリへ退避します。伝播後の波動関数（20）、密度（21）、Hartree場（22）、XCポテンシャル／エネルギー（23）、電流／全エネルギー（24）を検査し、保護経路の伝播器エラー（25）も返します。検査は物質MPI群全体で合意してから戻ります。エラー時は波動関数・密度・場・伝播履歴・XC状態・カウンタ・管理フラグを最後に成功した計算状態へ直接戻し、伝播の再実行をしません。失敗したステップの通常RT出力も追加しません。

Pythonの `NativeTDDFTBackend` が返す状態の `advance(A_old,A_new)` は、版付き仕様と全ランク一致を確認できた復帰だけを `NumericalError` で通知します。`A` と `steps` はそのままで、原因を修正して同じ受理時刻から再試行できます。不正な波形を自動的に弱めたり、固定dtを変更したりしません。旧ライブラリや未知のエラーは `RuntimeFailure` の保守的な扱いを維持し、復帰を推定しません。

対象は周期・非スピン分極・固定イオン・中点伝播のPZ、TB-mBJ（自動／固定係数）、Libxc-PZ/PBEです。SOC、+U、jellium、イオン運動／力出力、Gram–Schmidt、projection/unfold、microscopic Maxwell結合などの保護対象外設定は更新前に15を返します。標準バッチ計算は従来の演算経路を使います。既存の外部checkpoint形式も変更しません。メモリ退避・有限性検査の費用が毎ステップ増えるため、資源見積もりには余裕を残してください。予算は警告のみです。

これはTDDFT群の状態保証です。FDTD・他の物質がすでに進んだ混在ステップ全体を自動巻き戻しする仕様ではありません。MPI／プロセス停止、確保失敗、未制御の外部ライブラリSTOP、初期化深部のSTOPや出力装置障害は未対応です。試験専用の一部ランク障害注入は `SALMONPY_RECOVERY_TESTING` を明示して別ビルドした場合のみで、本番ビルドには注入用C入口がありません。

## 高速PR版HSEの深部数値異常からの復帰

新しい `NativeFastHSEBackend` 用ライブラリは `sp_fast_hse_recovery_protocol() == 1` を宣言します。対応範囲は周期・非スピン分極・固定イオンの高速PR版HSE06、Taylor4／学習PT-CN、物質群内k点分割のMPI1／2です。通常版HSEへの切替はしません。

各ステップ前に電子状態とACE・因子／PT-CN履歴・交換キャッシュ・カウンタを退避します。FFTプランやMPIコミュニケータなどの所有資源はコピーしません。波動関数・密度・Hartree・XC・端点の有限性検査に加え、対応したACE演算エラー、PT-CNの収束／Gram／射影検査、履歴更新エラーを物質群全体で合意して返します。最後に成功した計算状態を直接戻し、過去ステップや交換履歴を再計算しません。失敗試行の通常RT出力を追加せず、A・時計・電流・エネルギーを維持します。診断ログは残ることがあります。

Pythonは宣言された版と全ランクの既知エラーを確認できる場合だけ `NumericalError` を返し、同じ受理時刻から修正した入力で再試行できます。旧ライブラリ・未知エラーは `RuntimeFailure` の保守的扱いを維持します。dtや入力振幅を自動変更しません。退避と検査のメモリ・時間費用は増え、予算超過は警告のみです。

これはHSE物質群の保証で、FDTDや他の物質を含む混在ステップ全体の自動巻戻しではありません。初期化深部、確保失敗、プロセス／MPI停止、未制御の外部STOP、出力装置障害は対象外です。異常を捕捉できない場合に継続を推定しないでください。

検証用 `tools/test_fast_hse_deep_recovery_mpi.py BUILD GS NEW_OUTPUT --steps 32 --fault-step 13` は2ランク、OMP／BLAS1、外側のtimeoutで実行します。`--split` は4ランクの独立2群です。`--overflow-A 1e100` はPT-CNの実数値発散試験用で、科学計算の入力ではありません。`--faults` は `SALMONPY_RECOVERY_TESTING` を明示した別ビルドだけに使い、本番には注入用入口を含めません。検証記録は `docs/fast-hse-deep-recovery-validation-2026-11-17.md` です。合成小規模系の実装検証であり、物理精度の収束保証ではありません。

高速HSEのNumericalErrorは物質群の受理状態維持を意味する。ResponseSession全体を自動再試行しない。呼び出し元プログラムは場と他物質の進行を管理する。未知エラーや旧ABIでは継続しない。

## Pythonから直接呼ぶFDTDの数値異常復帰

新しい標準呼び出し元プログラムが管理するMPIライブラリは `sp_eh_recovery_protocol() == 1` を宣言します。`NativeEHField.advance` は対応する検査付き入口を選び、場・PML・場内の全LD・電流格子・時計をメモリ退避します。E更新後（20）、LD更新後（21）、H更新後（22）、点サンプル／積算A／通常観測出力の演算（23）の異常を場MPI群で合意し、退避した状態へ直接戻します。波の再伝播はしません。既知の復帰では `NumericalError` を返し、PythonのA・E・時計・failedも維持するため、入力を修正して再試行できます。失敗試行の通常観測ファイルは追加しませんが、診断・進行表示は残る場合があります。

対象は現在のFDTDのE/H更新平面波経路と場内LDです。線形応答専用経路、ASE／ART、周波数別観測は更新前に拒否します。バッチFDTDは従来経路を維持します。新しい入口を持たない旧バイナリ・未知エラーでは復帰を推定せず、従来の失敗扱いです。初期化・確保・MPI／プロセス停止・外部STOP・出力装置障害は保証外です。場群単体の保証であり、`CoupledStepper`／`ResponseSession`全体を自動巻戻し・自動再試行する機能ではありません。退避と検査の費用を見込み、メモリ予算は引き続き警告に使います。

`tools/test_eh_deep_recovery_mpi.py BUILD NEW_OUTPUT` を1／2ランク、OMP／BLAS1、外側timeoutで実行します。`--split` は4ランクの独立2群です。`--faults` は明示した別の `SALMONPY_RECOVERY_TESTING` ビルドだけに使います。本番には注入入口を含めません。数値オーバーフローの試験値は科学入力ではありません。

## 共通のステップ保護

STEPPING日英ガイドとmanual.jsonのprotected_stepを参照。StepTargetのcapture／rollback／releaseは完全な状態と副作用を扱う利用者宣言で、トークンをMPI転送しない。CoupledStepper.step_protectedは各local ownerの対象を一つずつ要求し、結合メタデータも戻す。再試行はStepRecoveryError.original_preservedがTrueかつcleanup_errorsが空の場合だけ、原因修正後に明示する。監査済み計算処理全体の退避・通常出力保護はnative_steps.native_step_targetsの別仕様で接続する。個別の数値復帰プロトコルから全体の安全を推定しない。

## 計算処理混在ステップと出力の保護

NATIVE_STEPPING日英ガイドを参照。対応する新規ビルドとnative_step_targetsを明示し、step_protectedを全ランクで呼ぶ。通常例は--protect-steps。受理済み未flush行を先にflushし、失敗試行のRT／点観測追記を巻き戻して同じRTユニットを開き直す。仮の追記は一時的に見え得るため同時外部アクセス・クラッシュ原子性を保証しない。作業フォルダ直下の専用通常.dataのみ。密度などの未監査出力を黙って無効化せず事前拒否。未知計算処理失敗・復元／解放失敗ではfailedを復活させない。targetはlive ownerに結び付くためcheckpoint反映で状態が置換されたら作り直す。Cpp候補とstep退避の同時利用は拒否。HSEは高速PRのみ。別GSは--hse-bundleで明示し、HSE固有の入力ハッシュと初期状態IDを保持する。メモリは退避同時保持の未較正費用と実測ピークを助言し、超過で停止しない。過去時間発展の再生はせず、再試行するのは拒否されたステップだけ。

## HSEの本筋と対称性

HSEは高速PR版のみを使い、対称性・star／orbit再利用・ACEを本筋とする。
FAST_HSE_SYMMETRY日英ガイドとfast_hse_symmetryタスクを参照。通常版へ切り替えたり、標準側の
別の縮約演算に合わせて高速版を変更しない。

NativeFastHSEBackend(...,allow_symmetry=True) と新しい対応計算処理ライブラリを明示する。
既定Falseは対称性を拒否。対称性opt-inは内部k1／2・CPU・周期・非スピン分極・full-q・factor-off
Taylor strict/learned・PT-CN off/learnedには対応するhistory capabilityを確認する。旧経路は1プロセスにつき初期化1回・共有対称化設定不変。open_modelの複数状態は現状では対称性なし・全k点のみ。内部2ランクにはsp_fast_hse_symmetry_kmpi_protocol()==1を確認する。
新C入口sp_fast_hse_init_symmetryとsp_fast_hse_symmetry_protocol()==1が必要。旧バイナリから
対応を推定しない。有限で一致したstatus34だけが状態維持のInputErrorで、Aの修正後に再試行可。
外場の自動投影・微小成分の切捨て・通常HSEへのfallbackは禁止。FDTD全体にはstep_protectedを使う。

通常混在例は--hse-symmetryと--hse-reference-bundleで明示した等価な全k点初期状態を使う。
実際の局所Aが群を破ると保護付きで拒否される。未適合な物質・配置の例を無理に通さない。
合成検証は--orbital-phase .23を明示するが、利用者のGSを黙って変更してはならない。

極小勾配のLibxc HSEの相殺は、現在の完全ソースで同じWPBEh式をu=s²で安定評価する。
全k参照ワーカーとMPIライブラリを両方新ソースから明示的に構築・監査する。protocol1だけで
数値修正を推定しない。test_hse_low_gradient.pyは80桁エネルギー差分、
test_fast_hse_low_gradient_mpi.pyは元のphase0・全／縮約k点の比較に使う。
probe_hse_low_gradient.pyは生のLibxc診断で、wrapper修正後も感度Trueになり得る。
これをwrapperの失敗と取り違えない。密度・勾配・Aをfloor/smoothせず、許容誤差を緩めない。
元の失敗例と通常回帰の結果・ハッシュを保持する。数値一致は物理収束の保証ではない。
内部MPI2で全／縮約k点と直接再開を検証済み。保存時と同じ分割を使い、分割変更はしない。次は初期化範囲の検証。分子は2026-11-22に再開し、双方向FDTDの中点結合検証を先行します。

高速HSEのsym.datは、計算状態の初期化前に形式・有限値・保持された群の構造を事前検査する。
InputErrorなら元入力を直しproviderを作り直して再試行可。既存providerの入力を黙って変更しない。
続くNumPy-onlyの幾何検査で、有効セル・原子配置／種別・格子と対称性整合を検査する。
k-star・擬ポテンシャル内容・再開バイナリとの互換性や全初期化異常復帰は保証しない。
詳しくはFAST_HSE_SYMMETRY日英ガイド。

幾何検査はSALMONと同じセル優先順位・格子解決・単位を使う。検査中にceil・セル正規化・原子移動などの入力補正をしない。
形式不正の拒否と物理収束を区別する。collective create_stateへ全solverランクで参加する。

高速PR版HSEの外場A/cはCartesianの3成分で、z限定ではありません。x/yや円偏光には実際の全時刻のAとセル・原子配置を保つ対称群を選ぶか、対応する全k点初期状態を使います。一般セルとz非保存群には新コアのcapabilityが必要です。波形の投影・微小成分の切捨て・群の自動変更はしません。詳細はFAST_HSE_SYMMETRY.md／FAST_HSE_SYMMETRY.en.mdを参照してください。

対称性付き高速PRのTaylor strict／learnedとPT-CNには `sp_fast_hse_symmetry_history_protocol()==1` を持つ明示ビルドを使います。学習PT-CNは既存方式どおり exact interval=1・predictor_corrector=n。設定を自動変更しません。履歴・ready・予測状態を保存し、最後に成功した計算時刻を再生せずに再開します。異なるモデルを同一プロセスで自由に混在できるとは説明しません。

## Python共有HSEモデル

NativeFastHSEBackend.open_model(initial,dt,context)で入力を一度だけ検証・コピーし、NativeFastHSEModel.new_state(point_id,dt=...,capacity=...)で元GSから独立した状態を作る。最初のnew_stateが計算状態の初期化を行う。状態ごとのstep/checkpoint/closeと出力ラベルを使い、直接復元は対応する新規状態へ行う。model.closeは生成を止めるだけで、生存状態は最後に閉じるまで接続を保持する。finallyで全状態を明示closeする。空のままのモデルや全状態終了後も、開いているモデルから再生成でき、隠れた完全電子状態を保持しない。model protocol1の新しい計算プログラムが必要。既存のMPI点配置や異種物質の同一プロセス対応を自動変更しない。SHARED_HSE日英ガイドとexamples/shared_fast_hse.pyに従う。

同一物質・セル・kメッシュ・対称操作の共有HSEは allow_symmetry=True で縮約k点にも対応する。shared_symmetry_protocol1の明示ビルドが必要。外場が実際の保持群に適合することを検査し、群や外場を自動変更しない。通常混在セッションはshared_model=Trueとpoints_per_groupを明示する。

context_mpi_protocol1の新しい高速PR計算プログラムでは、rt_open初期化で作った子コミュニケータは状態所有とし、state.closeで通信要求・XC／FFT資源の後に解放する。共通資源・呼び出し元プログラムが管理するMPIは維持する。異種入力対応や全MPI資源の独立化、初期化異常復帰まで完了とは説明しない。

標準TDDFTのcontext_mpi_protocol1は初期化〜初期電流の子コミュニケータを状態終了時に回収し、pending復元用の別状態も破棄する。一計算処理寿命／プロセスの制約は維持する。呼び出し元プログラムが管理するMPI継続と多重状態対応を混同しない。MPIラッパーの内部コンパイラ指定もビルド記録に残す。

標準TDDFTのhandle/model protocol1は、同じ全k点物質・dt・容量・MPI分割・初期Aの独立Fortran/Cで保持する計算状態を扱う実験的C入口。返すjは電荷電流で、旧sp_tddft_init/stepの物質Jとは異なる。保護／候補／A／保存上限は状態所有。Python create_stateの一寿命制約は維持し、共有生成器・異種入力対応とは説明しない。STANDARD_HANDLES日英と対応ヘッダーを参照。

## 収束診断と数値処理の配置

`CONVERGENCE.md`／`CONVERGENCE.en.md`を参照してください。Pythonは入力・計画・呼び出し・記録を担当し、数値処理は原則Fortran/Cに置きます。dtの基本指標は実測した電子数の減少率1e-5です。増加側も記録し、目安超過だけで計算を止めません。NaN/Inf等はその試験だけ直ちに終了します。設定nelecを測定値の代わりに使いません。新しく明示ビルドした通常XCのTDDFT状態はelectron_numberを公開し、run_dt_scanは既定で使用します。内部MPIでは全ランクが同順序で取得します。旧本体の測定口を推測せず、未対応なら更新を明示します。格子試験はSCF収束済みの隣接格子が指定精度を満たしたら0.1Å以前でも終了します。高速HSE／孤立系E駆動の測定はまだ対象外です。

格子収束はrun_current_grid_scanを主入口とし、各格子で電子数保存1e-5を満たしたdtの電流を比較します。共通物理時刻・3成分の相対RMS誤差5%以下を既定基準にします。エネルギーの一致で電流収束を代用せず、SCF未収束やdt未解決を確定扱いにしません。基準を変える場合は利用者の指定を記録します。

開発優先順位は4（SCFの完全保存・再開）→5（共有状態・初期化・異常復帰・C/Fortranの呼び出し規約）→3（NLCC・単位換算）です。1（高速HSEノルム）・2（物理格子収束確認）は保留し、分子作業は2026-11-22に再開しました。単一・固定イオン・E駆動を対象にFortran振動子で中点結合を検証し、その後に孤立系TDDFTを接続します。共通ソースのSCF進捗の状態所有は、完全なSCF継続APIの実装とは区別してください。

通常PZ/Libxc-PZ/PBEのSCFはsession.run(..., steps=K)で反復境界に停止し、continue_scf(steps=M)、snapshot()/restore()で混合履歴を保ったメモリ内継続ができます。実測した初期化／実行反復カウンタを確認し、波動関数だけのwarm startを完全なSCF復元と呼ばないでください。スナップショットは同じ生存状態・ビルド・MPI分割に限定し、ディスク／プロセス再開を実装済みと案内しないでください。新APIでも終端のSCF未収束は結果を保持して確認待ちにします。意図的なpaused=Trueの途中停止は対象外です。

高速PR版HSE06も同じSCFメモリ保存APIに対応します。全k点CPU・k点MPIでHSE対応ビルドを使い、低速HSEへ落とさないでください。混合履歴・PBE事前SCFの進捗・交換演算子の数値キャッシュを保持し、FFT/backend資源は複製しません。復元時の初期化／過去反復の再計算は0です。本来の将来のPBE→HSE切替処理と復元時の再計算を区別してください。HSE移行後の未収束は警告して保持します。対称性・ディスク再開の最新範囲はNATIVE_GS各節を参照してください。SOC/+U/DCは未対応です。

SCFのcheckpoint_capability／checkpoint_targetを使うと、共通execute_checkpointから
通常XC／高速PR版HSEのメモリ保存・直接復元を呼べます。互換性JSONだけを保持しても
再開できず、各ランクのNativeSCFSnapshotと生存所有者が必要です。トークンをMPI通信・
pickle・ディスク再開に転用しないでください。対応するcandidate protocol1があれば
restore_checkpointで独立した数値候補を準備でき、全参加者の成功確認後に各candidateを
publishしてから継続します。holderの中身は準備済み候補で、独立したGSソルバーでは
ありません。準備途中の処理可能な失敗では元の状態を保持して候補を破棄できますが、
反映時の致命的なメモリ／プロセス故障の一般的ロールバック保証ではありません。
旧計算プログラムは保護付き準備を拒否します。直接復元は引き続き非原子的です。
SCF未収束は警告して結果を保持します。

SymmetryOperationsをDFTへ渡すとsym.datを含む資産を物質ごとに保持し、縮約k点HSEの
SCF継続・保護付き準備に使えます。CPU・全q点・内部k点MPI1／2の範囲です。
操作の自動探索や一般の結晶群への対応を推測しないでください。後のRTでは外場を保存する
部分群が必要ですが、偏光全体のz限定ではありません。旧版の対称力でCartesianへの
戻し忘れを再現し、Fortranで修正しました。対称化した力だけは意図的な単位補正です。
固定格子Siの全点／縮約点のE・密度・力は記録した許容差内ですが、構造最適化等の精度認定ではありません。
旧計算プログラムにsp_gs_cartesian_force_protocol宣言がなければ、力を警告付きで保持して
基底未確認と表示します。Hartree/bohrを推測せず、暗黙ビルドやPython自動変換をしないでください。

## GS保存・再利用の基本方針（2026-11-19追記）

利用者の最新指示ではGSを保存して再利用する経路を基本とします。以前のGSデータ収集禁止を、保存機能の開発・利用の禁止へ拡大解釈しません。収束済みGSの再利用と、混合履歴を含む反復途中の完全再開は別の機能です。snapshot()は同じ生存ownerのメモリ内保存です。新しい対応計算プログラムではsave_scf()とrun(restart=...)による同一ビルド・同一入力・同一MPI分割のディスク再開を使えます。異なる計算プログラムへの移行は未対応です。MPI再分割は後述の対応protocolの範囲で明示します。保存先は指定された結果領域とし、不要な複製は作りません。通常XCの明示的対称性付きSCFもメモリ内保存・失敗時に元の状態を保つ復元に対応します。

SCFデータ保存時には残差履歴も確認します。新しい計算処理 GSのsave_scf()はdiagnostics/以下にCSVと任意のPNGを出します。未収束は警告してデータを保持し、グラフの有無とSCF成功・失敗を混同しません。Matplotlibがなくても保存済み数値データとCSVは使えます。履歴の計算・保持はFortran側で行います。

複数の計算処理 SCFをディスクからまとめて戻す場合はcheckpoint_target(directory=...)と共通restore_checkpointを使います。全候補をpublishしてから再開し、失敗後はtargetを作り直します。restore_scfの複数回直接呼出しを全体の原子的復元と扱いません。

計算処理 SCFの準備済み復元用の別状態を使わない場合はclose/discardまたはwith文で破棄します。GS状態のcloseと元のsnapshotのcloseも、それぞれが所有するPython候補を先に破棄します。MPI各rankで同順に呼び、C呼び出し元プログラムが直接作った予約やGCによる自動MPI終了まで保証すると説明しません。

通常XCのSCFディスク再開では、対応計算プログラムと新形式のglobal k範囲メタデータがある場合、run(restart=...,repartition=True)でk-only MPI1↔2の変更を検証済みです。同じbinary・物理入力・物質IDを使い、Fortranが配列を再配置します。高速HSEのk／軌道分割変更は後述の対応protocolで扱います。通常XCの空間分割は後述の専用入口を使い、HSE空間GSと通常XCの混合配置は後述の専用入口を使い、HSEのΓ点軌道×空間混合は後述のprotocol2入口を使います。既定と共通失敗時に元の状態を保つ復元の同一分割条件を勝手に緩めません。担当開始k点がGaussian乱数seedに含まれるため、MPI分割を変えた別の初期軌道生成を比較基準にせず、保存元と同じ初期状態からの連続実行と比べます。

## SCFのk点MPI数を変える再開

`NativeGSSession.run(restart=..., repartition=True)` は同じ計算プログラム・物質・
物理入力を保ち、k点だけのMPI分割を変える明示的な入口です。通常XCはprotocol1/2、
密なACEを使う高速PR版HSEはk-repartition protocol2が必要です。確認済みは1↔2ランク、
HSEの全k点／縮約k点とPBE準備段階です。FortranがACE・キャッシュ・履歴を再配置し、
SCFを再実行しません。HSEの軌道分割は後述の専用入口を使い、HSE空間GSは後述の専用入口を使います。packed/distributed ACE、異なる計算プログラムは
未対応です。失敗時に元の状態を保つ復元のMPI分割は従来どおり同一にします。詳細はNATIVE_GS.md／
NATIVE_GS.en.mdを参照し、物理収束の検証と混同しないでください。

通常XCの純粋な軌道並列も同じ明示的再開入口を使います。両側は
`decomposition=(1, comm.size, 1, 1, 1)` とし、orbital-repartition protocol1/2と
新保存のorbital_boundsを要求します。高速HSEの軌道並列は後述のprotocol2を使います。k／軌道間の配置変更は
未対応です。サンプルは両側に `--orbital-parallel` を指定します。

高速PR版HSEのGamma軌道並列SCFは、新計算プログラムの
`sp_gs_hse_orbital_protocol()` の版1または2 がある場合、同じMPI分割への直接保存／再開を
使えます。サンプルは両側に `--xc hse06 --orbital-parallel` を指定し、再開側に
`--restart` を渡します。ランク数変更には後述のprotocol2を使います。ACEとGamma交換作業状態を
各状態が所有し、混合履歴と一緒に直接復元します。詳細と対象条件はNATIVE_GSの最新節を
参照してください。複数k点の軌道HSEや軌道並列GS→RTを対応済みとは案内しません。

高速PR版HSEの軌道並列数を変える場合は、保存側・再開側とも
`NativeGSBackend(..., decomposition=(1, comm.size, 1, 1, 1), hse_orbital=True)`
とし、再開時に `run(restart=..., repartition=True)` を使います。1ランクでも
同じGamma交換経路を選ぶため、hse_orbitalは省略しません。HSE-orbital protocol2と
orbital-repartition protocol2が必要です。密なACE・タイルなしgaugeのCPU全支持
Gamma SCFを対象とし、空間／混合配置や軌道GS→RTを対応済みとは案内しません。
サンプルの `--xc hse06 --orbital-parallel` はこの明示設定を行います。

通常XCの純粋な空間分割変更は、spatial-repartition protocol1の新計算プログラムで
`NativeGSBackend(..., decomposition=(1,1,nx,ny,nz))` と
`run(restart=..., repartition=True)` を使います。新保存のgrid_boundsが必要です。
Fortranが密度・ポテンシャル・軌道／CG内部値と混合履歴を再配置し、SCFは再実行しません。
周期haloは内部値からコピーします。新分割の静的FFT／通信資源と履歴を持たないXC／
FFT作業バッファの準備を、GS再計算と混同しないでください。比較は格子内部で行います。
サンプルは両側に `--spatial-parallel` を指定します。対象はPZ／Libxc-PZ／PBEで、
HSE空間GSと通常XCの混合配置は後述の専用入口を使います。HSE混合は後述のΓ点protocol2入口を使い、別計算プログラム移行は未対応です。

高速PR版HSEの空間SCF再開は、保存側・再開側とも
`NativeGSBackend(..., decomposition=(1,1,nx,ny,nz), hse_spatial=True)` を使い、
再開時は `run(restart=..., repartition=True)` とします。1ランクでも同じGamma交換
経路を使うため指定を省略しません。HSE-spatial compute1とspatial-repartition2が必要です。
Fortranが密なACEと交換作業配列の格子行を直接再配置し、SCF再計算を行いません。
検証済みは1↔2、PBE準備段階、16³2↔4と二方向分割です。CPU・全支持・Gamma・
直交セル・固定イオン・全占有スピン対・タイルなしgaugeの範囲です。hse_orbitalとは
同時に指定しません。サンプルは `--xc hse06 --spatial-parallel` を両側に指定します。
複数k点×空間は後述のhse_kspace経路を参照し、packed/tiled/sparse表現・別計算プログラム移行・空間HSE GS→RTを
対応済みとは案内しません。これは実装検証で、物理収束の認定ではありません。

通常XCの混合並列SCF再開はmixed-repartition protocol1の対応計算プログラムで、
`NativeGSBackend(..., decomposition=(k, orbital, nx, ny, nz))` と
`run(restart=..., repartition=True)` を使います。積はMPI群のサイズに合わせますが、
物理的なk点・格子やSCF条件を変えません。保存rankのgrid/orbital/k範囲を併せた
5次元被覆でFortranが再配置します。密度・ポテンシャル・混合履歴はk／軌道群に
複製済みなので二重に足しません。再初期軌道やSCF反復で状態を作り直しません。
例は `--decomposition K O X Y Z` として公開し、物理k格子は固定です。
この通常XC入口はPZ／Libxc-PZ／PBEです。HSE混合は次のΓ点protocol2入口を使い、異なる計算プログラムへの移行は未対応です。

高速PR版HSEのΓ点SCFで軌道×空間を混合する時は、保存元・復元先とも
`NativeGSBackend(..., decomposition=(1,O,X,Y,Z), hse_spatial=True)` を使う。
計算 `sp_gs_hse_spatial_protocol()==2`、混合再配置
`sp_gs_scf_m_repartition_protocol()==2` が必要。物理kmeshは(1,1,1)、
完全占有スピン対・直交セル・CPU・密なACE・全支持の範囲。SCF混合履歴とHSEの
ACE／交換作用キャッシュを直接移し、`executed_iterations==0` で読み戻した後に
continue_scfする。再初期化した波動関数や交換作用で置き換えない。
サンプルは `--xc hse06 --decomposition 1 O X Y Z` と
`--save-scf`／`--restart --repartition`。複数k点×空間は後述のhse_kspace経路を参照し、圧縮／分割形式・
異なるruntime間の移行まで対応したと解釈しない。詳細はNATIVE_GS日英ガイド。

高速PR版HSEの全k点×空間SCF接続は、両側でhse_kspace=Trueと
decomposition=(K,1,X,Y,Z)を明示する。compute sp_gs_hse_kspace_protocol1/2/3と
mixed-repartition3を要求する。全一様k格子・CPU・完全占有スピン対・直交セル、
対称性縮約／軌道分割なし。hse_spatial／hse_orbitalやΓ点Wannier入力と混ぜない。
最後に成功した計算状態は空間局所だが、交換更新の全格子FFT作業は空間ランクで複製される。
メモリ削減や交換FFTの空間高速化と説明しない。既存の高速k点計算ルーチンを使い、
低速HSEには戻さない。再開はrun(restart=...,repartition=True)、過去SCF／交換の
再計算0を確認する。サンプルは--xc hse06 --hse-kspace --decomposition K 1 X Y Z。

HSE縮約k点×空間SCFはhse_kspace=Trueの計算protocol2以上を要求する。
既存のSymmetryOperationsと固定した元の一様k格子・代表点／重みを使い、
k群1/2、軌道分割は後述のcompute3。再開はmixed-repartition3。局所格子の整合検査を省かず、
格子行／代表点スライスと混合履歴を直接移す。無効化される幾何学的位相キャッシュを、
電子状態の再計算と混同しない。交換FFT作業は空間ランクに複製される。
サンプルは--hse-kspace --symmetry c2zと明示decompositionを併用する。

複数k点の軌道×空間HSEはhse_kspace=Trueのcompute3と
decomposition=(K,O,X,Y,Z)を使う。O<=states、CPU・完全占有スピン対・直交セル。
全／縮約k点のどちらも同じmigration3で格子行／軌道列／代表点スライスと混合履歴を
直接再開する。交換更新の全格子・全軌道作業領域はr/o群で複製される。
軌道／空間MPIを増やすだけで交換FFTやピークメモリが改善すると説明しない。
圧縮／分割／疎形式のSCF状態を生成・再開できると推定しない。サンプルは既存の
--hse-kspace --decomposition K O X Y Zと必要に応じ--symmetry c2zを使う。

SCFを同じソースの再ビルドへ再開する場合だけ、保持GSのrunに
compatible_runtime=Trueを明示できる。新保存のruntime_contractと、記録済みsolver／
adapterソース指紋・Fortran stream表現・HSE/FFTW種別・OS/CPUの一致が必須。
一致しない／旧保存で情報がない場合にハッシュを書き換えて回避しない。既定は同一
binary SHA256のまま。MPI変更は対応範囲のrepartition=Trueを別途明示する。
scf_restart_originとexecuted_iterations==0を確認し、再保存時は元を上書きしない。
サンプル再開側は新BUILD_JSONと--compatible-runtimeを使う。これは同じソース・
今回の同一Mac/toolchainでの実測であり、ソース版・OS跨ぎや汎用checkpointへの
対応を推定しない。generic restore_checkpoint／TDDFTの引数と混同しない。

異なるソース版ではcompatible_runtime="declared"／--declared-runtimeを明示し、
scf_compatibility.jsonの正確な保存元→復元先の登録だけを使う。逆向き・推移・任意の
同形式版を推定しない。Trueは同一ソース限定、既定は同一SHA。codecの型／生成IO
指紋、計算処理表現、OS/CPU、数値状態版も一致が必要。初期登録は状態版照会の追加
だけで演算／保存項目に変更がない通常XC・高速PR HSEの組。数値意味が変わる版には
適用しない。scf_restart_originのcompatibility_declarationと再計算0を確認する。
汎用checkpointやTDDFTの引数ではない。ソース差をレビューし、実測した組だけを
根拠文書付きで登録する。OS跨ぎ・保存形式変換は未対応。

GSの計算時間はstate.timing_secondsで照会する。全ランクが同順序で呼び、結果は
rank-localの完了区間秒数（total/initialization/scf_iterations/write_results）。
MPI合計・最大や物理精度の指標ではない。別GSの作成／継続で混ざらず、数値復元で
時間を巻き戻さない。別プロセス再開では新計測。旧runtimeは従来のGS計算を維持し、
新照会のみtimer protocol1が必要。同時スレッド実行／TDDFTタイマー分離を推定しない。

モデル切替では単位名と出力換算係数を同時復元する。境界は原子単位のまま、既存
SALMONの換算定数・式を保持する。DOS/PDOS/band/AVSは呼出しごとの空きファイル番号
を使い呼び出し元プログラムのファイルを閉じない。各状態の作業フォルダと入力の出力先設定を使う。
ファイル番号を数値checkpointに入れず、復元で出力を巻き戻すと説明しない。
作業ディレクトリのプロセス共有、同時スレッド実行の制約が解消したと推定しない。

保持GSではdirectory scope1/2/3対応runtimeを自動検出し、作成／継続／復元の際に
os.chdirを使わず状態フォルダを計算処理パス解決に渡す。相対ファイルは状態基準、
絶対パスは維持。例外・他ランク準備失敗でも範囲を解除する。論理入力／数値保存に
一時prefixを追加しない。旧runtimeの従来経路は維持する。旧runtime向けwrapperの
CWD切替や同時スレッド実行まで解消したと推定せず、計算処理呼出しはメインスレッドで
順次行う。通常のrun/continue_scfの呼出し方は変えず、import時にビルドしない。

通常TDDFT（通常／共有／異入力登録）と高速PR HSE（通常／共有／登録）は、scope2/3
runtimeで呼び出し元プログラムCWDを変更しない共通計算処理呼出しを使う。scope1はGSだけの対応で、
TDDFTでは旧切替へ戻す。import時にビルドせず、対応runtimeは明示的に準備する。
相対ファイルは各work基準、絶対パスは維持。準備のrank合意と例外時解除を共通化し、
数値の保存形式・時刻・電流符号は維持する。FDTDはscope3で同じ経路を使い、scope1/2では従来の切替を維持する。
主スレッド逐次の制約は維持し、同時計算処理実行まで対応したと説明しない。

標準TDDFTの共有／登録ハンドルは新しいtimer protocol1でstate.timing_secondsを物質MPI群から集団照会できる。total/initialization/density/hartreeはrankローカル完了区間の秒。totalは未完了の区間を含めず、生存中0の場合がある。呼び出し元プログラム・兄弟状態・park中の時間を除き、数値checkpointには入れず復元でも計測値を巻き戻さない。旧runtimeの通常計算を止めず、新照会だけは明示的な対応runtimeを要する。高速HSEやlegacy singletonのtimer分離まで完了と説明しない。

高速PR HSEもtimer protocol1対応の共有／登録ハンドルからstate.timing_secondsを集団照会できる。単位・完了区間・park中除外・復元で計測を巻き戻さない条件は標準TDDFTと同じ。legacy singletonは別のlegacy timer protocol1で扱い、新診断のために旧runtimeの通常計算を拒否しない。

標準Pythonから直接呼ぶTDDFTの終了時、呼び出し元プログラムのFortran unit30を閉じない修正が新規runtimeにある。共有ハンドルとlegacy singletonの両経路が対象。数値状態／更新順／保存形式は変更しない。古い固定runtimeをその場で改変せず、必要時だけ明示ビルド／対応配布を利用する。その他I/O所有権や一般異常復帰の完了とは説明しない。

標準TDDFT／高速PR HSEは新規runtimeで更新中も既存の状態MPI ownerを有効にし、遅延生成された子communicatorを状態終了時に回収する。戻る経路では所有範囲を解除し、他状態・呼び出し元プログラムのMPIを維持する。Pythonからの呼出し方／分割／数値checkpointは同じ。一般のMPI障害復帰・自由な同時計算処理実行と混同しない。

高速PR HSEの初期化整合性ガードは新runtimeでモデル所有。登録モデルごと／同一モデル生成器ごとに保持し、切替で共通記録を上書きしない。検査条件を緩めず、成功した初期化だけ記録を更新。C API／数値checkpointは不変。legacy単独の互換記録やその他共有設定まで分離完了と説明しない。

高速HSEのcheckpoint配列サイズ上限は呼出しローカルで、各serializerへ明示的に渡す。上限式・不正サイズ検査・stream形式は従来どおり。ユーザーのメモリ予算や新しい計算制限を追加する変更ではない。同時計算処理実行まで対応と説明しない。

標準TDDFT singletonも新runtimeでinit/step/closeの計測を状態所有にする。legacy timer protocol1時にstate.timing_secondsを物質群から照会できる。古いruntimeの計算を拒否せず新照会だけ対応版が必要。HSE singletonは対応する別のlegacy timer protocol1を確認し、同時スレッド対応と混同しない。

高速PR HSE singletonも新runtimeのlegacy timer protocol1でinit/init_symmetry・step・closeの計測を呼び出し元プログラムから分離する。state.timing_secondsは4項目のrankローカル完了秒数。共有ハンドル用protocolと別なので能力照会を確認する。旧runtimeの通常計算を止めず、新照会のみ対応版が必要。MPI障害復帰や同時実行まで対応と説明しない。

FDTDとHSE履歴検査のIEEE問い合わせはスカラーに限定し全要素を走査する。従来の閾値・MPI拒否・拒否時状態保持を維持し、新たな停止条件と説明しない。R0検証はGNU/macOSとMPI、合成履歴であり、Fujitsuコンパイル・実電子状態の新たなディスク再開検証まで完了と主張しない。

R1対応Pythonから直接呼ぶFDTD closeはpending復元用の別状態とstep_savedを解放する。途中終了で完走時解析を強制しない。任意のsp_eh_resourcesはfield/step_saved/pendingの3確保フラグのみで総メモリ・MPI解放の保証ではない。R2のMPI所有権対応と混同しない。

R2 Pythonから直接呼ぶFDTDはMPI子groupをfield所有にしcache終了後に回収する。sp_eh_mpi_resourcesはlive/slots/field_ownedの3値で、slots残存をリークと判定しない。親接続を外す前のlive数と別所有者の保持、detach後の0を検証する。既存sp_eh_resourcesの3確保フラグとは別。分割・更新順・数値checkpoint・一場寿命制約は同じで、一般MPI障害復帰や複数場同時実行と説明しない。

相反テンソルLDはTensorLorentzDrudeを使う。強度は(poles,3,3)の実対称正半定値、背景は(3,3)の実対称正定値。方向未知は従来どおり独立1成分の等方モデル。単独応答・MPI材料入口・process-local snapshot/restoreは対応。Pythonから直接呼ぶFDTDは主軸が格子に沿う対角テンソルに対応し、非対角電場更新は未接続。input_textとinput_assetsの組を保持し、scalarへ縮約して通してはいけない。permittivity_tensorで成分を明示し、A端点のJは分散極だけ、背景行列は場側が所有する。非相反応答・自動回転・一般テンソルフィット・新しいディスク再開形式までは実装済みと説明しない。

波長板サンプルfdtd_waveplateは--runを省くとプレビューのみで計算を開始しない。人工的な無損失板の検証であり、実酸化物の材料データや実験的精度保証と説明しない。carrierのStokes診断はFortran、比較は真空中平板のFresnel・多重反射込み。

相反非対角のPythonから直接呼ぶFDTDはFDTD_MATERIALS.md/enの入口を使います。Yee補間と連立E更新、分散J→P→Q、直接チェックポイントに接続済みです。等方／対角と分けて資産を生成し、古い計算プログラムへ渡さないでください。tensor_maxiter到達は警告、非有限試行は状態保持して再試行です。回転波長板はexamples.fdtd_waveplateの--angle-degreesを指定し、既定はプレビューのみです。別のA-only場への異方性背景や反対称磁気光学まで対応したと解釈しないでください。


## 全体通信を最小限にする

ジョブ全体のcollectiveは初期化、必要な異常合意、明示的な最終集約に絞ります。
物質点ごと・計算プログラム読み込みごとの全体通信を追加しません。担当グループと
点所有者への受け渡しを優先し、診断結果は原則として代表ランクだけへ集めます。
既存の通信が生存を確認できる場合、追加のBarrierや生存確認allreduceを入れません。
同じ段階の状態・データ検査を一つにまとめても、有限性検査、物理的に必要な総和、
異常時の一貫した拒否、保存／復元・再試行の保証は維持します。

検証専用mpi_runtime_receiptは初期bcast一回と最後のroot gather一回で終了します。
計算プログラムの各再利用はMPI通信をしません。NativeEHField.advanceのPython側検査は
担当場グループ内で六回から二回へ統合しています。計算処理内部の通信はこの回数に
含めません。ResponseSessionのAは所有点分だけをScattervで配り、場パケットも担当点の
辞書だけをscatterします。担当点が無い場グループと単一ランク群への追加配布は
省きます。Jは所有点分を代表へ回収します。対応するPythonから直接呼ぶFDTD結合では、
Yee辺への加算先を持つ場ランクだけへ担当点Jを配ります。
物質群のprivateなresponse._Aは担当行だけが現在値です。全点の履歴には場の代表が
保存した出力を使い、物質群の非担当行を全点データとして扱わないでください。
計算処理 Jの担当範囲はFortranのsp_eh_current_maskで取得し、境界をまたぐ点は両側へ配信します。
場の代表は全点Jを保持します。他の場ランクの結合Jは担当行だけが有効です。
全点解析には代表の出力を使います。旧ABIとMultiscale1Dは従来配信を維持します。
保存形式・MPI互換条件とNativeEHField.advance(J)の全点入力仕様は変更しません。
詳細はCOUPLED.mdの担当領域に応じたJ配信を参照します。隣接合意を統合した段階では、場の結果と応答入力の事前検査をまとめました。その段階の1ステップで
A結合は4→3、場パケットはHなし7→4／Hあり7→5です。1回は階層合意で、計算処理内部や
保護付き保存・復元の通信は別です。H採取前の合意と物質更新後の合意は削除しません。
独自の公開stepを上書きした応答は、初期化時に全群で従来の経路を選択します。
一ランクだけ異なる応答でも合意回数は揃います。通常のstepに原子的保証を
追加する変更ではないため、全体の状態保持・再試行には既存の保護付き入口を使います。


### A/J輸送に同梱する異常情報
標準の結合経路では、A／場パケットとJのヘッダに異常有無・通知元ランクだけを載せます。
全群が異常を知った時だけ、そのランクをrootに専用領域で説明文を一回Bcastします。
正常ステップでは説明文Bcastなし。専用領域はParallelContextで一度だけ複製・解放し、
呼び出し元のMPIを解放しません。詳細はCOUPLED.mdのA/J同梱輸送を参照してください。
独立した階層合意は直前版のA3／Hなし4／Hあり5から、現在は1／2／3です。
代表間のPython MPI集団通信も5／6／7から4／5／6へ減っています。計算処理内部、初期化、
保護付きトランザクションを含む総通信回数や通信時間の比率ではありません。
J後の結果通知、H前の合意、受信検査と保護付き再試行は残します。
正常時の説明文Bcastがないことと、通常stepの原子性は別です。全体を戻す場合は
明示的な保護付き入口を使います。独自の公開応答入口と単独ResponseSessionは従来経路です。


### HPC用FDTD入力の準備と起動は分ける
管理ノードではprepare_fdtdまたはexamples.prepare_fdtdを使い、入力保存で終了します。
この入口はSALMON、mpiexec/srun、ビルド、ジョブ投入を呼びません。準備完了を計算開始の
許可や指示と解釈しないでください。計算の起動は別途、ユーザーのサイトのジョブ内で行います。
FDTDMPI.commandはプレビューのみ。rank_argumentsの{ranks}置換以外のシェル解釈はありません。
既定のrunとMPIランク内二重起動拒否は変更しません。既にMPIアプリ内ならCOUPLEDの
caller-owned同じプロセス内での計算入口を使います。DFT/TDDFTのバッチ投入共通化はまだ別課題です。


### DFT／TDDFTの管理ノード準備
prepare_dftはfresh DFT入力と不変PP／対称性ファイルのみを保存します。save_restartの
既定Trueは将来の計算でGSを保存する設定です。準備そのものはGSを生成しません。
prepare_tddftは明示した計算処理 RT入力とfilesだけを保持し、GS探索・読込・snapshot化、
dtや波形などの書換をしません。計算ノード側の既存GSを管理ノードへ要求しません。
どちらもprepared-not-runのPathを返して終了し、実行ファイル、MPI初期化、SCF、worker、
ビルド、ジョブ投入へ自動継続しません。典型例はexamples.prepare_electronsです。
計算が依頼されていない入力準備でTDDFTInitialState.from_salmonやnew_stateへ進まないで
ください。実行・転送は指定された計算ノード／ジョブで別操作です。これらは完全な物理・
restart整合性検査ではなく、HSE高速PR限定や既存実行入口の制約を解除しません。


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

## FDTDのエネルギー保存診断

API.mdのエネルギー収支節を参照。マクロ点Eとcharge jの仕事はabsorbed_energy_density、
分離したE/Hポートはreflection_transmission、総収支はfdtd.energy_balanceを使います。
E/Hの位置・時刻、jの符号、担当体積・面積、場の始終エネルギーを利用者の設定から確認します。
足りない観測量をゼロと推測して保存則が成立したと説明しません。E/JやAだけからR/Tを作らず、
物質の可逆蓄積を熱と呼ばず、TDDFT/LD蓄積を仕事と二重計上しません。診断超過は警告のみ。
非一様な周期面はreflection_transmissionのarea_weightsとsurface_fluxで全点の局所流束を積分します。面内平均E/Hの外積で代用しません。閉曲面の自動観測や回折次数別の分解は未実装です。

実結合の収支確認済み範囲はAPIの日英「実際の結合計算での確認」節を参照します。
LDは観測Eから同じ計算処理計算ルーチンでjを再構成し最終計算処理電流と照合、TDDFTは2点・8ランク
のFDTDのE/H更新/j出力を直接使っています。初期充填波がある場合の初期エネルギーを無視せず、
有限時間のR/T>1をクリップしません。実材料の光学収束まで確認済みと説明しないでください。

横方向周期セルの全面積分では、E/Hの時空間位置合わせ、面積重み、周期端点の重複、
参照場の引き算を確認します。surface_flux.power/energyは総量、PortSpectrumのfluenceは
面積平均でありenergy_balanceが面積を一度掛けます。MPI面分割では比を平均しません。
バッチFDTDの自動面取得はPlaneMonitorを明示し、旧point観測の200点制限とは区別します。

PlaneMonitorを指定したバッチFDTDは所有タイルを計算処理から取得し、FortranでYee E/Hを
空間位置合わせしてresult.planesへ返します。時間中心化済みHを二重に平均しません。
result.reflection_transmission(reference)／energy_balanceで診断へ接続し、物質仕事と場の始終差は
明示的な観測値を使います。旧計算プログラムや面診断失敗をFDTD計算失敗と混同しません。
この入口はfresh batch用であり、NativeEHFieldや再開ストリームの自動統合を対応済みと案内しません。
管理ノードではprepare_fdtdまで、runは利用者が計算環境で明示的に呼びます。


## 同じプロセス内での計算FDTDの面取得

NativeEHField.sample_planesは場の全ランクから、保護付きstep/復元が確定し候補を解放した後に明示的に呼びます。場のrank0だけに不変PlaneDataを返し、診断失敗は警告+Noneで計算状態を変えません。Fortranの局所読み出しは状態更新/MPI/IOなし、Pythonは場の群のgather一度だけです。観測時刻で取得し、combine_plane_samplesで等間隔時系列にして既存の流束/R/T/収支へ渡します。output_directoryは単一サンプルの任意保存で上書き禁止、失敗時も既存出力と計算を保ちます。確定前に呼んで投機状態を公開しないこと。後の復元が以前の独立観測アーカイブを自動削除するとは案内しません。混在サンプルの--planes-x/--plane-strideは明示的な横周期境界と確定後の間引き取得を選択します。仕事・場の蓄積差は観測期間に合わせて明示します。

新しいPythonから直接呼ぶFDTDでは未使用座標領域をゼロ初期化し、計算で使う座標・haloは従来式のままです。既存のソース/保存ファイルを一括書換えしません。再現用の小さな明示的テストはtools/test_fdtd_coordinates.pyで、SALMONジョブ・MPI・GSは実行せず、任意指定コンパイラで依存のない座標ルーチンだけを検査します。GNU検証済み、Fujitsu未検証です。


## Pythonから直接呼ぶFDTDの局所格子エラー

新しいE/Hの1106は既存の最小局所格子条件による初期化拒否で、正常な計算の停止条件を増やしません。全領域境界表を使うため新たな合意通信もありません。sp_eh_init_recovery_protocol()==1と全ランク一致1106を確認した場合だけ、NativeEHFieldは失敗接続と初期化予約を解放します。同じ親MPIを保ち、利用者が設定を直して新しい保存先で再試行できます。ジョブ投入後の割当・入力自動修正はしません。旧版・未知コードや任意STOP/OOM/MPI障害まで回復可能と案内しないでください。


明示的なMPI配置では1105（Maxwell分割の既存条件・接続群サイズ不一致）、1107（実際の配置処理の未実装分岐）も、プロトコル1と全ランク一致が確認できた場合に再試行可能です。1106と同様、失敗接続・初期化予約を解放し、親MPIを保ちます。配置名は外側でホワイトリスト化せず、実際のmapperの分岐に任せます。利用者が追加した計算コードや共有ライブラリの方式を先回りで拒否しないこと。自動分割・Fujitsuトポロジー・任意の確保/MPI異常まで対応済みとは案内しません。

### Recoverable material-ID CUBE input errors

FDTD updates of E and H initialization returns1115 for unavailable, malformed, truncated,
nonfinite or grid-count-mismatched SALMON material-ID CUBE input. Reported tiny
reader staging-buffer allocation failure returns1116; this allocation path has
not been fault-injected and does not cover general OOM or solver allocations.
The reader owns its file unit and checks every value used in integer conversion,
including the final partial record. Existing streaming order and material-ID
conversion are preserved; this is not a general electronic-density CUBE reader.

For external shape input only, one field-group initialization agreement precedes
halo exchange, so a read failure on one rank returns on all ranks. No per-step
collective or generated-shape collective is added. Known unanimous statuses and
recovery protocol1 permit cold cleanup and retry with a new work directory on
the same caller communicator; status-absent standalone behavior is retained.

Local GNU MPI1/2 tests cover retry, caller MPI survival, direct resume and exact
agreement with prior generated-shape fields, point output and checkpoint bytes.
A host-owned unit600 survives rejected shape reads. Ordinary point headers,
point appends, final observation information, plane snapshots/integrals and the
material-ID CUBE writer now open a local newunit and close only that owned unit.
Filenames, record formats, ordering and numerical/communication operations are
unchanged. The generated-shape writer no longer closes the caller's requested
unit for an unimplemented output format. Existing helper signatures are retained.

GNU MPI1/2 tests keep a host-owned unit600 open throughout initialization,
32 steps, direct checkpoint restore and finalization. Field arrays, point files,
checkpoint bytes and generated shapes match prior controls; enabled plane output
also matches. This is file-unit ownership isolation, not output-error recovery:
write failures may still use existing Fortran error handling. Linear-response/spectral/ASE/ART units are addressed below; TTM has been removed;
legacy restart/helper ownership is addressed below. No claim covers arbitrary device failure,
external STOP, process failure or all initialization allocations.

### Legacy EM restart/helper file-unit ownership

The common_maxwell text/binary restart readers and writers now open their own
local newunit, including per-rank data and metadata files and the remapping read
path. txtfile_copy_em independently owns its input and output units instead of
borrowing700/701. Existing helper signatures retain ifn for source compatibility;
it no longer selects an actual file unit in these helpers. Paths, formats,
metadata, array order, copy trimming and communication/reconstruction operations
are unchanged. No additional startup or per-step collective is introduced.

The existing I/O fixture has an explicit --fdtd-restart-io mode: synthetic EM
arrays only, no SCF/GS data. GNU MPI1/2 tests preserve calling program units600/601/700/701
while exercising text values, real/complex1–5D binary arrays, single/all modes,
same-layout and cyclically exchanged target-domain reads, inactive payloads and
text copying. Control/occupied outputs and prior-runtime outputs match byte for
byte. This tests the helpers and does not certify every complete legacy batch
restart configuration or arbitrary rank counts. Normal LD+FDTD/direct-resume
field, point, plane, shape and checkpoint data also match prior controls.

Linear-response/spectral/ASE/ART units are addressed below; TTM has been removed.
This change isolates file-unit ownership; missing/corrupt restart data, IO device
failure and general allocation failure retain their existing handling. No general
IO-error recovery, Fujitsu execution or multi-node certification is claimed.

### Spectral output file-unit ownership

E/H finalization owns local newunit handles for linear-response trajectories and
spectra, frequency-resolved E/H/J planes, and ASE/ART spectra. The two simultaneous
windowed/unwindowed ASE/ART files use two local handles, not base-unit arithmetic.
The periodic linear-response response file is explicitly closed after its last
record; its earlier implicit close at process termination is no longer needed.
No module-wide file selector is added. Fourier transforms, normalizations, field
updates, loop/file ordering, names, record formats and MPI operations are retained.

The existing solver I/O tool has --fdtd-spectrum-io, using actual small LD+FDTD
lifecycles with no electronic GS. GNU MPI1/2 controls cover isolated and periodic
impulse response, isolated non-unit-background raw response, Fourier planes, and
separate ASE and ART cases with three frequencies and both window variants. Caller-owned
units600/601/602 remain connected to their original files. Named spectral outputs
are closed before return; expected numeric outputs are nonempty and finite.
Unoccupied, occupied and previous-runtime outputs are compared byte for byte.
These are output/calling program-contract regressions in atomic units, not convergence or
material-property certification. FDTD updates of E and H's existing supported-mode limits are
not lifted: the batch calling program used here enables spectral modes, not the in-process
step interface. No extra runtime collective, automatic build or startup run is
introduced. TTM has been removed; general IO-error recovery remains separate.

TTM is excluded from SALMONpy: unpublished implementation, automatic input-file activation, build targets, examples and utilities are removed from the maintained solver. A leftover `ttm.inp_ttm` file does not activate any calculation. Historical provenance records are retained.

### Recoverable automatic MPI factorization failure

FDTD updates of E and H returns1117 when the existing automatic MPI distributor cannot assign
its rank count. The current automatic domain routine factors2/3/5; this is a limit
of that algorithm, not a blanket rank-count restriction. A supported explicit
spatial partition remains usable (tested with7 ranks). No Python/outer factor
whitelist is added and the distribution algorithm is unchanged. Correct the
partition and create a new work directory on the same caller communicator.

The internal Fortran helper has an optional status (0success,1unassigned prime
factor,2no allowed non-orthogonal k/orbital split). Without status it retains the
legacy diagnostics/STOP behavior. Matching Fortran modules must be used when
building callers; this is source compatibility, not stable Fortran binary ABI.
The C initialization signature and recovery protocol1 are unchanged. Python
releases a failed lifetime only for a recognized unanimous1117 with protocol1.
The existing cold cleanup/status agreement is reused; no new collective or
per-step work is added. This does not handle allocation, Fujitsu topology, MPI
failure or every specialized initialization stage.

Local GNU tests cover unchanged automatic choices on MPI1/2 and helper failure/
continuation on1/2/7. solver MPI7 repeats raw/public automatic-failure rejection,
keeps parent MPI alive, releases owned children and retries explicit decomposition;
32-step fields, point bytes and direct checkpoint/resume match the prior runtime.
This is local contract/regression evidence, not multi-node scaling certification.

### Retrying rejected LD tensor/subpixel input

The existing solver readers already return1110 for missing/invalid diagonal or
reciprocal-tensor material sidecars and1113 for a missing/invalid reference
subpixel configuration. For these cold errors the solver initializer destroys the
candidate and releases its child resources. Python now recognizes these two
unanimous statuses with sp_eh_init_recovery_protocol1, detaches, and releases its
failed lifetime reservation. Correct/stage the material file and retry in a new
work directory on the same caller communicator. Error messages identify the
material input instead of reporting a generic solver rejection.

No Fortran routine, input validation, coefficient arithmetic, file format,
communication or solver binary changed in this correction. The existing audited
runtime is reused by recorded SHA256. Unknown/mixed statuses or an absent/wrong
recovery protocol do not enable retry. Handled coefficient/Bloch-mode returns are addressed below; arbitrary allocation,
IO-device, MPI or process failure remain outside the bounded guarantee.

The existing E/H recovery tool accepts --material-input-mode diagonal|tensor|subpixel
and --material-input-admission. GNU MPI1/2 repeat missing, invalid-header, truncated
and nonfinite-file rejection twice each through the public Python constructor;
candidate/attachment resources are released and the parent communicator remains
usable. Valid material initialization then completes32 steps and direct restore.
Compare against the pre-change Python package with the same solver binary; this
is an input/lifecycle regression, not a new material-convergence certificate.

### Retrying handled coefficient/special-mode initialization errors

Python also recognizes existing solver cold returns1111 (nonfinite diagonal/
tensor LD update coefficients),1112 (nonzero Bloch phase with the unsupported
off-diagonal tensor mode), and1114 (invalid subpixel update coefficients).
Only unanimous codes with initialization recovery protocol1 release the failed
lifetime. The solver initializer already deallocates its candidate, releases its
send/receive cache before its owned MPI children, and returns; the correction is
Python's classification and diagnostic text, not a solver numerical change.

Correct the input/mode and use a new work directory on the same communicator.
This does not implement Bloch/off-diagonal tensor evolution, silently set Bloch
phase to zero, relax finite checks or change coefficient formulas. Files already
written during rejected initialization remain in the rejected work directory;
the rejection is not a promise of filesystem rollback. Unknown/mixed statuses,
missing protocol, trapped floating-point exceptions and allocation/device/MPI/
process failures remain outside this bounded retry guarantee.

--coefficient-admission in the existing E/H tool uses material-input-mode
(diagonal,tensor,subpixel) to exercise1111/1112/1114. Finite but oversized test
rates provoke the existing coefficient guards; a solver periodic test input
provokes the existing unsupported Bloch/tensor return. GNU MPI1/2 repeat each
failure and then run valid32-step/direct-resume controls. Reader rejections are
combined in the same fixture. The old SHA256-verified solver binary is reused;
no solver rebuild, extra numerical calculation in Python or new collective.

### Retrying handled early grid/boundary/time errors

FDTD updates of E and H already returns1101 for inconsistent cell/spacing/grid-count input,
1102 for incompatible periodic/boundary faces,1103 for an invalid duration/
step-count/time-step specification, and1104 when dt is not below its existing
CFL bound. Python now releases the failed initialization reservation for these
recognized unanimous returns only when initialization recovery protocol1 exists.
The solver cold candidate/owned-child cleanup remains the same. Correct the
input and retry with a new work directory on the same caller communicator.

No new input restriction, automatic repair, CFL threshold change, solver binary
change or collective is introduced. This aligns Python's error classification
with existing solver behavior. Unknown/mixed statuses, missing protocol and
other unhandled failures remain conservative. Rejected directories are retained;
this does not promise rollback of staged or diagnostic files. No automatic job
resubmission or solver execution is added to preparation/checking workflows.

The existing --grid-init-retry fixture now includes these four public-constructor
rejections twice each. Combined with --geometry-admission, GNU MPI1/2 complete35
rejection/retry checks per run, verify parent MPI/child cleanup, then32 steps and
direct checkpoint restore. Field/point/shape/checkpoint bytes match earlier valid
controls. This is a local lifecycle regression, not all initialization/IO-error
recovery or cross-platform/physical-convergence certification.

FDTD called directly from Python core-allocation recovery:1118 is a retryable cold-init return only
when all field ranks agree and the runtime declares the initialization recovery
protocol. It covers material-ID/scratch and core split-E/H/coefficient/saved-E/H
3D allocations, not all LD/PML/auxiliary arrays or OS/MPI failures. Budget excess
remains advisory; do not introduce a hard memory-budget stop. Failure injection
belongs only in a separate `SALMONPY_RECOVERY_TESTING` build and the explicit
`--allocation-admission` validation; never in production or automatic import.

solver cold-init1119 extends recoverable allocation failure to the initial
`ld_allocate`'s24 pole/current/polarization/coefficient/history arrays. Require
unanimous1119 and the declared initialization recovery protocol for retry.
The owner cleans partial state before retry; do not reuse a partial LD object.
This does not cover material/coefficient scratch, tensor auxiliary, subpixel
reallocation or PML/observation allocations. Keep memory-budget excess advisory.
`--ld-allocation-admission` uses only an explicit test runtime; production must
not export the fault hook. No additional per-step collective is introduced.

solver cold-init1120 covers three material coefficient vectors and six boundary
E/H coefficient vectors. Permit retry only with unanimous known1120 and the
runtime's init-recovery protocol. Share the material gate with LD input checking:one extra field-group startup
agreement with LD,two without LD,and no extra step collective. PML bulk3D coefficients are already in1118; small
automatic PML workspaces and arbitrary auxiliary/OS/MPI failures are excluded.
Do not claim blanket PML/OOM recovery or change budget warnings into hard stops.
`--aux-allocation-admission` belongs only to an explicit testing runtime.

Optional solver point-output failures (obsN_at_point_rt.data header/append) warn
with SALMONPY_POINT_OUTPUT_WARNING and continue physical evolution. Never turn
these optional I/O errors into an automatic abort. A later observation retries
append after the destination/header is repaired. Missing/partial rows are not
replayed or repaired; check diagnostic completeness before interpretation. No
extra MPI collective or caller-unit ownership change. Other outputs,required
restart writes and filesystem-wide failures are outside this bounded support.

Optional plane-integral header/append I/O errors warn with
SALMONPY_PLANE_INTEGRAL_WARNING and continue physical evolution. Keep all existing
integration/MPI sums even after root output failure;do not skip peer collectives
or add error-only agreements. Future writes retry after path/header repair;
missing/partial rows remain incomplete. Check diagnostic coverage before using
warned integral/energy-balance results. Other output families remain separate
scopes;do not apply this warning-only policy to required checkpoints implicitly.

Optional obs0_info.data finalization errors emit SALMONPY_OBSERVATION_INFO_WARNING
and continue normal state/MPI cleanup. Separate completed physics from missing/
incomplete visualization metadata. Keep existing data/checkpoints;do not label
this optional-write warning as a numerical failure or automatically rerun physics.
Check metadata availability before analysis. No automatic post-finalization retry
is provided. Required restart-save errors remain a separate policy and must not
be silently accepted under this optional-output rule.

Optional time-resolved plane snapshot errors emit SALMONPY_PLANE_SNAPSHOT_WARNING
and continue other components/planes and physical evolution. Keep every existing
plane aggregation/MPI collective even after root output failure. Later snapshots
use their time-stamped files;missing/partial frames are not fabricated or repaired.
After path repair,a user-driven post-restore continuation can write the same
frame again;restoration itself does not replay steps. Check frame coverage before
analysis. Snapshot staging-array allocation,frequency maps and required restart
errors remain separate scopes. No extra MPI communication/global warning state.

Batch FDTD frequency-resolved plane-map I/O errors warn with
SALMONPY_FREQUENCY_MAP_WARNING and continue other windows/planes/energies and
cleanup. Keep every existing Fourier accumulation/MPI aggregation;no extra
collective or recomputation. Missing/partial maps are not fabricated or retried
after finalization. Verify coverage before analysis. This is the batch path:
FDTD updates of E and H's unsupported diagnostic return15 is unchanged;do not claim solver
protected stepping/checkpoint support for these accumulators. Other output and
required restart errors remain separate policies. The explicit test option
--fdtd-spectrum-io --plane-map-faults focuses the existing calling program probe on two
energies,three planes,two windows and real optional-output faults;it is not CI
physics validation or an automatic import-time calculation.

Batch LR/ASE/ART trajectory/response/spectral I/O failures warn with
SALMONPY_OPTIONAL_SPECTRUM_WARNING and continue other outputs and cleanup. Keep
Fourier processing and every integration/MPI reduction even if a file cannot be
written;ASE/ART windows own independent units/status/names. No extra collective
or recomputation. Missing/partial tables are not fabricated or retried after
finalization:check required trajectory/energy/window coverage before analysis.
FDTD updates of E and H unsupported-diagnostic return15 is unchanged;this is batch support.
Required restart-save and allocation failures remain separate policies.
The explicit --fdtd-spectrum-io --spectrum-faults probe reuses LR isolated,
periodic,non-unit-background,ASE and ART cases with actual file faults;do not
turn this local validation into mandatory CI physical-observable reproduction.

Required saved solver state saves are not optional diagnostics. Nonzero from
sp_eh_checkpoint(path,0) means the save failed even if other ranks produced new
files. Existing files are not overwritten and live numerical state/clocks remain
unchanged. Catch/report failure;repair or choose a new empty destination and retry.
Preserve committed data and quarantine failed-attempt files;do not silently reuse
an incomplete set or delete prior data. Per-rank distinct paths and higher-level
all-success/integrity publication are required;raw calls to the Fortran/C library are not atomic
multi-file archives. Existing solver failure/retry is verified,no new format or
numerical implementation was necessary. Do not generalize this result to legacy
batch restart writers,mid-stream disk-full/close-only errors or process/MPI loss.

### Required legacy batch EM saves

`data_for_restart` text/binary helpers require explicit `io_status`; no unpublished-API backward-compatibility fallback is provided. Caller initializes the phase status to zero; first local failure skips subsequent local saves. Adjacent exclusive `.salmonpy-pending` staging checks open/write/close before rename. Failed staging evidence is retained. Batch callers perform one owning-group agreement per metadata/checkpoint phase and exit nonzero on required-save failure; ordinary time steps gain no agreement.

Do not interpret per-file publication as an atomic archive: files/ranks can contain mixed generations after a failed phase, and no fsync guarantee is made. Preserve a separate completed restart directory; quarantine failed staging files before explicit retry. Do not reuse a failed set merely because individual files exist. solver protected restoration remains a separate path. The explicit GS-free `tools/test_em_restart_save.py NEW_OUTPUT` exercises maintained writers under MPI1/2; `--reference-source PREPARED_TREE` additionally compares all supported normal payload bytes.

### Batch restart-directory creation status

Existing restart directories, including valid directory symlinks, return zero status. Creation failure is agreed once within the owning MPI group using a replacement for the existing creation barrier; batch ranks then exit nonzero. Conflicting existing files/targets are preserved. Do not treat directory-existence success as a writability certificate; actual writers handle I/O errors. No new ordinary-step communication, automatic preflight, or backward-compatibility branch is introduced.

### Required batch EM restart readers

Reader `io_status` is mandatory and collective within the owning solver group: initialize it to zero on all participating ranks and call readers in the same order. Checked read-only open/read/close precedes publication; shared-read status gates payload broadcasts, local reads agree failures, and remapped metadata uses one status-bearing header instead of three broadcasts. Workspace status is agreed before payload transfer. This adds restart/setup gates, not ordinary-step communication. No compatibility fallback is supplied.

A batch failed read terminates initialization nonzero before using invalid values. Low-level helpers return an agreed error while preserving caller MPI/file ownership; explicit repair and status reset permit retry. Their failed buffers and a partially initialized batch model must be discarded, not published as valid state. They do not provide whole-model atomic rollback; solver candidate-based protected restoration remains separate. Never rewrite the input archive to suppress errors. In midpoint comparisons, legacy observation-file copying can strip trailing spaces from previously written rows; compare values separately while retaining exact binary payload comparisons.

Explicit local validation: `tools/test_em_restart_read.py BUILD_JSON NEW_OUTPUT --batch`; no GS is required. Keep this physical run out of required CI.

For shared input/output restart paths, new layout metadata is published only after every required old payload is read successfully; a failed read must leave that input set unchanged. This does not make later multi-file checkpoint saves atomic. Keep a separately completed archive.

### Directional LD material allocation errors

The1120 initialization recovery contract also covers partial allocation of
diagonal/full-tensor LD material-input arrays. Agree within the field group
before array transfer, release partial arrays and owned communication resources,
and permit constructor retry only after the existing recovery protocol succeeds.
Invalid material files remain1110. Replacing the existing status broadcast with
a maximum reduction adds no collective; step arithmetic and checkpoint formats
are unchanged. Memory-budget excess only warns. This does not cover arbitrary
allocation failures or OS termination. No legacy API or format fallback is added.

### LD coefficient preparation allocation failures

The1120 constructor-retry contract covers the LD material list, six common
coefficient arrays and three directional coefficient arrays. Skip dependent
coefficient/LD-state work after local allocation failure, agree within the field
group before use, release partial ownership and retry only after the existing
recovery protocol succeeds. Existing status collectives carry the result; no
startup or step collective is added. Nonfinite directional coefficients remain1111.
Budget excess only warns. Do not claim coverage of arbitrary tensor/subpixel or
observation workspace allocation, OS termination or MPI failure.

### Off-diagonal LD grid-coefficient allocation failures

The1120 initialization-retry contract covers per-grid coupling weights,
diagonal coefficients and integer material indices. Allocate before the existing
field-group material-array agreement; initialize values and calculate coefficients
at the original point. No new collective or arithmetic-order change. Release
partial ownership before retry; budget excess only warns. Do not apply this
coverage claim to subpixel allocations, tensor propagation workspaces, optional
observation/image workspaces or OS/MPI failure. No old-version fallback is added.

### Tensor-update workspace failure

Current Python stepwise E/H calls allocate explicit tensor-update workspaces
before physical mutation, using the existing input agreement for field-group
status. Known agreed24 plus the recovery protocol raises RuntimeFailure without
setting field.failed or advancing E/H, LD, time or outputs. Release memory and
retry the same J; do not replay history or change dt automatically. Workspaces
are call-owned, reused during iteration and excluded from snapshots/checkpoints.
No extra collective; budgets only warn. This is not arbitrary OOM/backup-allocation,
OS/MPI or observation-buffer recovery. Standalone processes stop on failure.

### Common host contract before MPI attachment

Current built-in MPI libraries expose sp_host_abi(capacity,int64_values):10 words,
version1, primitive widths32/64/64/32/64, MPI_Fint widened to int64, host-owned MPI,
main-thread sequential state operations, required I/O directory protocol3. Query
locally before attachment; it initializes no MPI/solver and performs no communication.
Use api.query_host_abi on an explicitly loaded trusted library, or doctor
--native-build BUILD_JSON --component fdtd --route single. SHA256/build audit
trust is still required. Built-in NativeParallelContext rejects subgroup-wide
ABI mismatch before attach; no old-version fallback or automatic build. Use
single/shared/registry, not legacy. Unknown user backends remain unverified and
retain their own contract. This descriptor is not operation/profile/convergence
or checkpoint certification; use the corresponding current queries.

## Common operation contracts

Read `native_contract(component, route=...)['host_operations']` before calling
the built-in host bridge. It declares exact C argument/return tokens, collective
scope, lifetime preconditions, buffer ownership and error limits. Python binds
all six required symbols before executing any of them; missing symbols reject
the solver group at admission. This adds no per-step communication. A declaration
is not machine-code validation or full solver ABI certification. Never interpret
void/value-returning operations as recoverable error-code APIs. Keep solver-specific
state/token/checkpoint/profile contracts separate. User backends keep their own
contracts; do not impose this built-in symbol set on them.

## Built-in electronic state-operation contracts

Use `native_contract('tddft' or 'fast_hse', route=...)['state_operations']` for
current argument order/types, ownership and checkpoint/transaction mode meanings.
Never interchange standard/HSE transaction mode2: standard releases, HSE also
reopens/flushes and retains on failure. Protocol/profile queries remain required.
All state calls use the owning group/main thread/matching order/model I/O directory.
Only successful outputs are valid; do not infer arbitrary rollback from nonzero
status. Keep charge j versus standard single-state C matter J signs explicit.
Stage/commit/cancel plus host coordination provide session protection; raw load
alone does not. Optional diagnostics/retained-GS capabilities and custom backends
keep separate contracts.

Use ordered `parameters` names/directions as well as argument types; old_A and new_A share a type. Standard single C current is `J_matter`; decomposition order is k/orbital/grid-x/grid-y/grid-z.

## FDTD operation contracts

Read `native_contract('fdtd',route='single')['field_operations']` and
`native_fdtd.h`. Init grid triples are one-based; Python inputs are zero-based.
The C parameter J is physical charge j, atomic units, no volume multiplier/sign
flip. current_mask/plane_read are local; H sampling and stepping/checkpoint/
transactions are field-group collective. Field transaction mode3 is no-op, not
electronic stream reopening. Keep optional capability checks at use; do not
require unused features or impose this ABI on custom backends.

## GS operation contracts

Use `native_contract('gs',route='single')['gs_operations']`. Single is the entry
family, not a one-state restriction. State info/read are collective (token
agreement), resources/format queries local. Snapshots retain numerical/mixing/
HSE state only within a live owner, not text logs or lifetime execution counts.
Candidate2/3/4 stage/publish/discard; disk0 save,disk2 stage read then candidate3
publish. Authenticate records and coordinate before publication. Keep default
unconverged-SCF review; an explicit developer-fixture acceptance is not physical
convergence. Never require unused optional capabilities or constrain user backends.

## Units / maintained NLCC evidence

Use `unit_contract()` before interpreting SALMON eV inputs versus Python energy
presentation. Preserve existing constants; relative energy conversion difference
is4.3951750e-8, not a reason for silently rewriting inputs. Geometry/ASE length
uses SALMON's retained Bohr value. NLCC CPU unpolarized NC-UPF force evidence
must match current executable/source-tree and completed whole-app reference;
missing/old source identities reject. Check all6 core geometries and both MPI
formats; MPI-linked validated builds need their own2/4-rank evidence. Completed
application records may be reused explicitly for the exact binary/reference,
never inferred for a different build. No material-convergence/GPU guarantee.

NLCC certification scripts/tests are maintainer source-checkout tools requiring audited build objects; wheel users consume the completed validated runtime/receipt through load_validated_build or runtime bundles, without compiling or running those maintainer tests.

## Explicit initial guess in preparation

prepare_dft(...,initial_guess=GroundStateResult) stages only authenticated
info/occupation/wfn files and records required_executable_sha256. It shares
DFT.run checks/input rendering; no execution/automatic GS discovery/submission.
Choose/verify that executable on the compute node in a separate job. Atomic
publication rejects a changed staged seed. This is wavefunction warm start,not
mixing-history SCF resume or general HSE/live/custom-backend seed support.

## 単一分子の双方向FDTD結合

`MolecularFeedback` は呼び出し元が作成した二つの場・一つのE駆動分子状態・MPIを使う公開入口です。`dipole_radiation` を検索し、DIPOLE_RADIATION.mdと完全例 `examples/molecular_feedback.py` を参照します。`start()` の1dt保持と中点E・自己真空場差し引きの条件を省略しません。分子状態はtotal_rootだけ、場は各担当ランクに渡し、初期化・終了は呼び出し元が担当します。取消は公開のsavepoint/restore_savepoint/release_savepointを使う直接復元です。検証記録revision5を必須とし、旧記録や履歴再生へ暗黙に切り替えません。ワーカー喪失時は復元保証なしとして扱います。`StepRecoveryError.original_preserved` を確認せずに再試行せず、暗黙のSCF・ビルド・ダウンロード・自動再試行は追加しません。

孤立系の `savepoint()` は同じワーカー内で一つだけ保持する退避です。明示的に解放します。`snapshot()/restore()` と `save_checkpoint()/restore_checkpoint()` は受理履歴の保存・再実行であり別の用途です。直接退避中に履歴再生は行わず、保持中の履歴アーカイブ作成は拒否します。プロセスを終了した後に使える直接ディスク保存は未実装です。

分子の直接ディスク再開は`save_direct_checkpoint`／`restore_direct_checkpoint`。
`salmonpy.isolated-direct/1`のFortran状態を新しいワーカーへ読み込み、経過済みE履歴を
再計算しない。同じ実行ファイル・GS入力・dt・容量・位置・スレッド数を指定する。
GSは初期化で読み込むがSCFしない。保存は新規ディレクトリのみ。復元候補の失敗時は
元の状態を維持する。既存のsave_checkpointは明示的な履歴再生形式なので混同しない。
分子だけの保存をFDTD結合全体の保存と案内しない。直接形式には現行revision5検証記録が必要。

分子と二つの場の一括保存は`MolecularFeedback.save_checkpoint`／`restore_checkpoint`。
受理後の本計算FDTD n・自己真空参照FDTD n・分子 n−1を明記する。
次のジョブで同じ入力・ライブラリ・GS・dt・容量・MPI分割から作成し、start()を呼ばずに復元する。
両場のFortran内の番号と、分子の更新回数0をverified_stepsで確認する。
形式はsalmonpy.molecular-feedback-direct/1。分子単体のsave_direct_checkpointと区別する。
任意のhistoryはtotal_rootだけが渡す。サンプルは--stop-after／--checkpoint／--resumeを使う。
準備失敗の原子的拒否と、OS/MPI破損・切替時の致命的確保失敗は区別する。
通常ステップに通信を増やさず、検証は明示的な保存・再開時に行う。

六面エネルギー収支は、正の非分散背景ε・μ、等方／対角／非対角LDとLD界面平滑化、非負のOhmic伝導率を確認済み。
実際のYee更新係数からFortranでエネルギーの重みを取得する。場の更新順や値は変えない。
非対角背景とLD界面平滑化では、背景補間の境界交換量をPoynting流束と別に記録する。Ohmic伝導の散逸量はohmic_workとして分離し、LD分極への仕事・流束等と共に残差へ含める。
field_modified_*（離散保存量）とfield_physical_*（Hを中心化した通常のエネルギー）を区別する。
tools/test_diagonal_field_energy_mpi.pyは実計算を行う開発用検査であり、入力検査ではない。

For non-diagonal background field energy observation, use the current rebuilt 19-component observer. Report constitutive_outward_energy separately from Poynting outward_energy; both enter the corrected residual. LD subpixel observation uses the actual smoothed coefficients; check the subpixel summary and reference-frequency approximation separately. Ohmic dissipation is returned separately as ohmic_work and included in residuals; LD polarization work is not solely heat. Rebuild the field and lightweight libraries together for the current contract. The optional query adds owning-field-group halo work; it does not change the physical step or checkpoint.

Field observation profiling: see MEMORY.md / MEMORY.en.md and tools/test_diagonal_field_energy_mpi.py in the source checkout. It executes a calculation; never launch it automatically as management-node input validation. Distinguish completed calculation with partial resource measurements, individual process high-water RSS from simultaneous job peak, sampled-step energy from cumulative balance, and observation maxima from physical checkpoint payload. Grid-count changes here extend vacuum at fixed spacing; this is not the deferred material-grid convergence campaign.

Performance field execution: NativeEHField(..., execution_mode="performance") skips per-step full-state snapshots, finite scans and validation-only group communication for FDTD/embedded LD. Default is checked. Select once and identically within the field group; setup validation, build audits/SHA, local ABI bounds, physical communication and solver convergence/breakdown remain. No automatic preservation/retry guarantee; invalid values may propagate or terminate the job. Explicit protection/diagnostics/checkpoints retain their own costs. This does NOT disable SBE/TDDFT/HSE matter or outer coupling checks. Rebuild the field runtime for sp_eh_step_performance; see COUPLED.md / COUPLED.en.md.

## Standalone performance mode

Pass `execution_mode="performance"` to `new_state(...)` or the batch response
entry (`material.respond(...)` for LD, `crystal.respond(...)` for SBE). The
default is `"checked"`; the choice belongs to each state and cannot be changed
during stepping. LD supports scalar and reciprocal tensor materials; crystal
SBE remains driven by A/c and returns charge current density.

Performance mode uses the same Fortran update arithmetic. LD omits per-step
pole backups and finite-value scans; SBE updates its density matrix directly
without a trial copy or per-step finite-value scans. Automatic rollback and
retry after a failed update are not guaranteed in this mode. Initial input,
library integrity and shape checks remain. Explicit `snapshot()` / `restore()`
and saved output retain their checks. Restoring physical state does not change
the execution mode. This option does not disable checks in TDDFT or external
coupling drivers.

The installed examples accept `--execution-mode performance`. Compare a checked
run before using this mode for a new configuration.

## Performance mode for serial coupling

`Transverse`, `Longitudinal`, `ThinFilm` and `Multiscale1D` accept
`new_state(dt=..., steps=..., execution_mode="performance")`. The same option
on `maxwell.couple(...)` selects performance mode for both the field and every
matter state. Built-in SBE supports this; TDDFT and fast HSE do not yet support
this entry. The default remains `"checked"`.

Setup still validates inputs, capacity, CFL and timing contracts. During stepping
the field omits current/result finite scans, and the serial driver omits current
finite scans and repeated clock comparisons. Shape/capacity checks remain;
update formulas, ordering, sign conversion and physical communication are unchanged.
Explicit snapshots, restoration and saved output retain their checks. Physical
snapshots exclude the mode; restore preserves the receiving state's mode.
Automatic failure rollback is not guaranteed in performance mode. This option
does not select performance mode for MPI sessions or TDDFT/HSE internals.

Custom field/matter models must accept the `execution_mode` keyword in
`new_state` and expose the selected mode on the returned state. A requested
performance mode is rejected if the model cannot honor it; there is no silent
fallback. Ordinary checked custom models need no changes.

Try `python -m salmonpy.examples.field_coupling NEW_DIRECTORY --execution-mode performance`.

## TDDFT performance mode

`NativeTDDFTBackend(..., execution_mode="performance")` selects unchecked
updates for the caller-managed MPI route. The default is
`"checked"`; the backend and returned state's mode are read-only. Supported
verification profiles are periodic fixed-ion scalar PZ, TB-mBJ, Libxc-PZ/PBE
with transverse supplied A/c and middlepoint propagation. Single-state, shared-model and heterogeneous-registry routes are supported.
Process-isolated workers, fast HSE and molecules are not covered by this option.

Initialization, library SHA256 authentication and explicit checkpoints remain
checked. Normal updates omit Python endpoint/finite/clock agreement, directory
admission votes, Fortran argument agreement, full-state backups and recoverable
numerical scans. Physical solver communication and update order remain. Shape,
capacity and lifetime bounds remain local checks. Caller ranks must issue valid
matching operations; nonfinite values can propagate and native failures can end
the MPI job. A returned update error marks the Python state failed; it cannot be
retried under the checked-mode state-preservation promise.

Explicit norm queries, observation, protected operations and output retain their
checks and costs. The outer coupling driver and FDTD have independent settings.
The sample `mpi_coupled_tddft` accepts `--matter-execution-mode performance` to
select this TDDFT route only; its field/outer checks are retained. The C entry
`sp_tddft_step_performance(Aold,Anew,J_matter)` returns SALMON matter current;
Python changes the sign once to charge current. Aold is ABI-only in this route.

Small synthetic GS tests compare checked/performance/previous library for all
four XC profiles on 1/2 internal ranks. Current and measured electron-number
arrays agree exactly, including explicit in-memory restoration. They establish
neither real-material convergence nor general performance gains.

## Per-state modes in shared and registered TDDFT

`NativeTDDFTBackend(..., shared_model=True, execution_mode="performance")`
and `registry_model=True` set the default for states created by normal sessions.
For direct model use, `model.new_state(point_id, execution_mode="checked")` or
`"performance"` overrides that default per state. Omitting the keyword inherits
the model default. `registry.open_model(..., execution_mode=...)` can set each
registered model's default. Each created state's choice is read-only; sibling
states may use different modes. Shared registry reuse preserves each binding's
selected default instead of inheriting a previous binding's mode.

Performance state-token updates still check local handle generation/lifetime,
capacity and ABI bounds. They omit token/endpoint/model-selection diagnostic
votes, numerical scans and automatic backups. Required model/input/workspace
activation and physical solver communication remain. Checked-state operations
keep the previous collective validation. Explicit snapshots/restoration remain
checked and do not change the receiving state's mode. No automatic rollback or
retry guarantee is implied for performance failures.

Synthetic regression covers interleaved checked/performance states, independent
controls, PZ/TB-mBJ/Libxc-PZ/PBE, heterogeneous PZ↔TB-mBJ and Libxc-PZ↔PBE, and
1/2 internal ranks. This is not a real-material convergence or speed benchmark.

## Fast-PR HSE performance mode

`NativeFastHSEBackend(..., execution_mode="performance")` selects the
caller-managed MPI route; the default is `"checked"`. The choice
is fixed on each state. Single-state, shared-model and heterogeneous-registry
routes are supported; serial subprocess workers are not supported by this option.

The performance update omits endpoint agreement, Python finite scans, electronic
and exchange-state backups, and numerical-recovery diagnostic votes/scans.
Initialization, library authentication, local lifetime/capacity/shape bounds,
explicit checkpoints and observations remain checked. Solver breakdown checks,
FFT/Hartree/exchange reductions and learning-readiness communication remain.
The mathematical compatibility of A with the retained symmetry little group
also remains: status34 rejects that input before mutation and permits retry.
Other returned failures mark the Python state failed, without an automatic
rollback guarantee; unhandled native failures can end the job.

Learning endpoint collection stays before energy/output handling in both modes.
An owner-local scheduling flag separates this ordering from recovery diagnostics;
it is reset after each update and is not an execution-mode selector in shared
process globals. Taylor/PT-CN, full/reduced k and internal1/2-rank synthetic
comparisons preserve current, energy and learning history, including explicit
in-memory restoration. These are not physical convergence or speed benchmarks.

`mpi_mixed_native --hse-execution-mode performance` requires an explicit
`--fast-hse-mpi-build`. It selects only HSE updates;
the field and outer coupling retain their independently selected checks.

## Modes of shared and registered fast-HSE states

A backend's `execution_mode` supplies its model default. Use
`model.new_state(point_id, execution_mode="checked")` or `"performance"`
to override per state; omit the keyword to inherit the model default.
`registry.open_model(..., execution_mode=...)` selects a registered model's
default. Session registry reuse preserves each material binding's chosen mode.
Siblings can mix checked/performance updates with separate electronic, exchange,
learning, A/current and output history. There is no process-global mode selector.

The performance handle update keeps local token generation/lifetime bounds and
the symmetry/A mathematical condition, while omitting diagnostic token votes
and automatic recovery copies. Physical activation of each model's input and
workspaces, solver communication and learning update order remain. Checked
siblings retain their full validation and recovery path. Explicit memory or disk
checkpoint operations retain checks/costs; restoring does not change the
receiving state's selected mode. Performance update failures still have no
automatic rollback guarantee.

`shared_fast_hse --execution-mode performance` demonstrates point-a in the
selected mode and point-b in checked mode, including a matching point-a restart.
With `--second-initial`, the same example uses two registered input models.
The usual `mpi_mixed_native --hse-execution-mode performance` also permits
shared HSE points; field/outer-coupling checks have independent settings.

## Periodic TDDFT in a separate process: performance mode

`TDDFTInitialState.new_state(..., execution_mode="performance")` and
`crystal.respond(initial, drive, execution_mode="performance")` select the
periodic A-driven worker mode. `TDDFTResponseBackend(execution_mode=...)` passes
the choice to each worker used by a response session; outer session/field checks
retain their own settings. Default is `"checked"`, fixed per state. Molecules
and fast-HSE subprocess workers are not covered by this option.

The child must accept `INIT_PERFORMANCE` and reply with
`SALMONPY_TDDFT_PERFORMANCE_V1`; an old binary is rejected, without fallback.
Mode is passed only in that child's environment, never installed in the Python
host environment, and does not select embedded-host behavior. Setup, executable
SHA256, protocol shape/sequence, capacity and explicit output checks remain.
The step omits Python/worker endpoint-current finite checks and automatic
reconstruction/replay after a worker failure. A failed performance worker is
closed; do not continue it under the checked-mode state-preservation promise.
Physical propagation, update order, requested electron-number output and solver
breakdown handling remain. No global collective is added.

Accepted history is retained for explicit `snapshot()` / `restore()` and output.
This worker's explicit restore still starts a candidate and **recalculates saved
history**; it is not direct wavefunction restoration. Snapshot validation and
replay-current comparison remain, and the receiving state's mode is preserved
even when restoring a snapshot made in the other mode. Finite/nonfinite output
validation still applies when explicitly snapshotting or saving results. Normal
GS copying/integrity rules are unchanged; this mode does not reduce startup I/O.

Examples `tddft_response` and `mpi_tddft` accept `--execution-mode performance`.
Use a matching explicitly built or distributed worker; import does not build it.

Molecular E-driven TDDFT also accepts `new_state(..., execution_mode="performance")`. Default to checked mode; use performance only when explicitly requested. Explicit saves/restores remain validated and retain the receiving mode. Do not promise automatic failure recovery or unchecked outer FDTD coupling.

Outer MPI: `ResponseSession(..., execution_mode="performance")` and `CoupledStepper.from_multiscale(..., execution_mode="performance")` select outer diagnostics independently of backend modes. Every rank must agree. Per-step finite scans and redundant coupling admission are omitted; callback exception synchronization and physical A/J communication remain. Explicit protected steps still create backups and retain rollback. Do not describe this as zero communication or universal failure recovery.

Outer E/H+LD+two-point TDDFT performance validation covers 3/6 local ranks and 96 steps, including explicit protected rollback/retry after native propagation. Do not generalize this synthetic PZ result to HSE, molecular field packets, other platforms, or physical convergence.

Fast-PR HSE outer performance mode: backend, field and ResponseSession modes are separate explicit choices. Protected steps retain direct internal-state backups including HSE resources. Never remove reduced-k field-symmetry admission, physical exchange communications or necessary callback exception agreement as a performance optimization.
